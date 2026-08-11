args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(prefix, default = NULL) {
  hit <- grep(paste0("^", prefix), args, value = TRUE)
  if (!length(hit)) return(default)
  sub(paste0("^", prefix), "", hit[[1]])
}

max_scrna <- as.integer(get_arg("--max-scrna=", "4"))
if (is.na(max_scrna) || max_scrna < 1L) {
  stop("--max-scrna must be a positive integer")
}

workspace <- Sys.getenv("GSE200996_WORKSPACE")
if (workspace == "") {
  workspace <- normalizePath(getwd(), mustWork = FALSE)
}

manifest_path <- file.path(
  workspace,
  "03_rebuild",
  "manifests",
  "GSE301741_RAW_ROUTE_SAMPLE_QC.csv"
)
out_dir <- file.path(workspace, "03_rebuild", "validation", "GSE301741_pilot")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(manifest_path)) {
  stop("Raw-route QC manifest not found: ", manifest_path)
}
if (!requireNamespace("Seurat", quietly = TRUE)) {
  stop("R package Seurat is required for this pilot.")
}
if (!requireNamespace("Matrix", quietly = TRUE)) {
  stop("R package Matrix is required for this pilot.")
}

qc <- utils::read.csv(manifest_path, stringsAsFactors = FALSE, check.names = FALSE)
scrna <- qc[
  qc$modality == "scRNA" &
    qc$extracted_exists == TRUE &
    qc$h5_status == "ok",
  ,
  drop = FALSE
]

if (!nrow(scrna)) {
  stop("No extracted QC-pass GSE301741 scRNA H5 files found.")
}

patient_timepoints <- aggregate(
  timepoint_harmonized ~ patient_id,
  data = scrna,
  FUN = function(x) paste(sort(unique(x)), collapse = ";")
)
paired_patients <- patient_timepoints$patient_id[
  grepl("pre", patient_timepoints$timepoint_harmonized) &
    grepl("post", patient_timepoints$timepoint_harmonized)
]

add_priority <- function(df) {
  fraction_priority <- match(
    df$sort_or_cell_fraction,
    c("CD3", "CD45pos", "Unsorted", "CD45neg", "CD45ratio")
  )
  fraction_priority[is.na(fraction_priority)] <- 99L
  library_priority <- match(df$assay_or_library, c("5prime", "3prime"))
  library_priority[is.na(library_priority)] <- 99L
  df$.fraction_priority <- fraction_priority
  df$.library_priority <- library_priority
  df
}

select_balanced_pilot <- function(scrna, paired_patients, max_scrna) {
  scrna <- add_priority(scrna)
  selected <- list()
  if (length(paired_patients)) {
    for (patient in sort(paired_patients)) {
      patient_rows <- scrna[scrna$patient_id == patient, , drop = FALSE]
      for (tp in c("pre", "post")) {
        tp_rows <- patient_rows[patient_rows$timepoint_harmonized == tp, , drop = FALSE]
        if (!nrow(tp_rows)) next
        tp_rows <- tp_rows[
          order(tp_rows$.fraction_priority, tp_rows$.library_priority, tp_rows$geo_accession),
          ,
          drop = FALSE
        ]
        selected[[length(selected) + 1L]] <- tp_rows[1L, , drop = FALSE]
        if (length(selected) >= max_scrna) {
          out <- do.call(rbind, selected)
          out$.fraction_priority <- NULL
          out$.library_priority <- NULL
          return(out)
        }
      }
    }
  }

  already <- if (length(selected)) {
    do.call(rbind, selected)$geo_accession
  } else {
    character()
  }
  remaining <- scrna[!scrna$geo_accession %in% already, , drop = FALSE]
  remaining <- remaining[
    order(remaining$patient_id, remaining$timepoint_harmonized, remaining$.fraction_priority, remaining$.library_priority),
    ,
    drop = FALSE
  ]
  needed <- max_scrna - length(selected)
  if (needed > 0L && nrow(remaining)) {
    selected <- c(selected, split(head(remaining, needed), seq_len(min(needed, nrow(remaining)))))
  }
  out <- do.call(rbind, selected)
  out$.fraction_priority <- NULL
  out$.library_priority <- NULL
  out
}

pilot <- select_balanced_pilot(scrna, paired_patients, max_scrna)

read_counts <- function(path) {
  counts <- Seurat::Read10X_h5(path)
  if (is.list(counts)) {
    if ("Gene Expression" %in% names(counts)) {
      counts <- counts[["Gene Expression"]]
    } else {
      counts <- counts[[1]]
    }
  }
  counts
}

summary_rows <- list()
objects <- list()

for (i in seq_len(nrow(pilot))) {
  row <- pilot[i, , drop = FALSE]
  sample_id <- row$geo_accession
  message("Pilot loading ", sample_id, ": ", row$raw_archive_file_name)
  counts <- read_counts(row$extracted_path)
  cells_before <- ncol(counts)
  genes_before <- nrow(counts)
  obj <- Seurat::CreateSeuratObject(
    counts = counts,
    project = "GSE301741_pilot",
    min.cells = 3,
    min.features = 200
  )
  obj$geo_accession <- row$geo_accession
  obj$patient_id <- row$patient_id
  obj$timepoint_harmonized <- row$timepoint_harmonized
  obj$sort_or_cell_fraction <- row$sort_or_cell_fraction
  obj$assay_or_library <- row$assay_or_library
  obj$site <- row$site
  obj$treatment <- row$treatment
  obj$response_harmonized <- row$response_harmonized

  obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(obj, pattern = "^MT-")
  obj[["percent.ribo"]] <- Seurat::PercentageFeatureSet(obj, pattern = "^RP[SL]")

  meta <- obj[[]]
  summary_rows[[i]] <- data.frame(
    geo_accession = row$geo_accession,
    patient_id = row$patient_id,
    timepoint_harmonized = row$timepoint_harmonized,
    sort_or_cell_fraction = row$sort_or_cell_fraction,
    assay_or_library = row$assay_or_library,
    cells_before_create = cells_before,
    genes_before_create = genes_before,
    cells_after_create = ncol(obj),
    median_nCount_RNA = stats::median(meta$nCount_RNA),
    median_nFeature_RNA = stats::median(meta$nFeature_RNA),
    median_percent_mt = stats::median(meta$percent.mt),
    median_percent_ribo = stats::median(meta$percent.ribo),
    response_harmonized = row$response_harmonized,
    stringsAsFactors = FALSE
  )
  objects[[sample_id]] <- obj
}

summary_df <- do.call(rbind, summary_rows)
summary_csv <- file.path(out_dir, "GSE301741_RAW_ROUTE_PILOT_QC.csv")
utils::write.csv(summary_df, summary_csv, row.names = FALSE)

object_path <- file.path(out_dir, "GSE301741_RAW_ROUTE_PILOT_OBJECTS.rds")
saveRDS(objects, object_path)

tcr <- qc[
  qc$modality == "scTCR" &
    qc$extracted_exists == TRUE &
    qc$tcr_status == "ok",
  ,
  drop = FALSE
]
tcr_match <- tcr[
  tcr$patient_id %in% unique(pilot$patient_id) &
    tcr$timepoint_harmonized %in% unique(pilot$timepoint_harmonized),
  ,
  drop = FALSE
]
tcr_csv <- file.path(out_dir, "GSE301741_RAW_ROUTE_PILOT_TCR_CANDIDATES.csv")
utils::write.csv(tcr_match, tcr_csv, row.names = FALSE)

report <- c(
  "# GSE301741 Raw Route Pilot",
  "",
  paste0("Created: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  paste0("- Requested max scRNA files: ", max_scrna),
  paste0("- Loaded scRNA files: ", nrow(summary_df)),
  paste0("- Total cells before CreateSeuratObject: ", sum(summary_df$cells_before_create)),
  paste0("- Total cells after CreateSeuratObject: ", sum(summary_df$cells_after_create)),
  paste0("- Candidate matched TCR files: ", nrow(tcr_match)),
  "",
  "## Files Written",
  "",
  "- `03_rebuild/validation/GSE301741_pilot/GSE301741_RAW_ROUTE_PILOT_QC.csv`",
  "- `03_rebuild/validation/GSE301741_pilot/GSE301741_RAW_ROUTE_PILOT_OBJECTS.rds`",
  "- `03_rebuild/validation/GSE301741_pilot/GSE301741_RAW_ROUTE_PILOT_TCR_CANDIDATES.csv`"
)
report_md <- file.path(out_dir, "GSE301741_RAW_ROUTE_PILOT_REPORT.md")
writeLines(report, report_md)
cat(paste(report, collapse = "\n"))
cat("\n")
