workspace <- Sys.getenv("GSE200996_WORKSPACE")
if (workspace == "") {
  workspace <- normalizePath(getwd(), mustWork = FALSE)
}

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
})

manifest_path <- file.path(
  workspace, "03_rebuild", "manifests", "GSE301741_RAW_ROUTE_SAMPLE_QC.csv"
)
module_path <- file.path(
  workspace, "03_rebuild", "results", "external_validation",
  "GSE123813_gene_set_manifest.csv"
)
out_dir <- Sys.getenv("GSE301741_REBUILD_OUTPUT_DIR")
if (out_dir == "") {
  out_dir <- file.path(
    workspace, "03_rebuild", "validation", "GSE301741_raw_reconstruction"
  )
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_csv <- function(path) {
  if (!file.exists(path)) stop("Missing input: ", path)
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

manifest <- read_csv(manifest_path)
modules <- read_csv(module_path)
scrna <- manifest[
  manifest$modality == "scRNA" &
    manifest$h5_status == "ok" &
    manifest$extracted_exists %in% c(TRUE, "TRUE", "True", "true"),
  ,
  drop = FALSE
]
tcr_manifest <- manifest[
  manifest$modality == "scTCR" &
    manifest$extracted_exists %in% c(TRUE, "TRUE", "True", "true"),
  ,
  drop = FALSE
]
if (nrow(scrna) == 0) stop("No readable GSE301741 scRNA H5 files.")
max_scrna <- suppressWarnings(as.integer(Sys.getenv("GSE301741_MAX_SCRNA")))
if (is.finite(max_scrna) && max_scrna > 0) {
  scrna <- head(scrna, max_scrna)
}

lineage_markers <- list(
  T_cell = c("CD3D", "CD3E", "TRAC", "CD247", "LCK"),
  NK = c("NKG7", "GNLY", "KLRD1", "FCGR3A", "PRF1", "TRAC"),
  B_cell = c("MS4A1", "CD79A", "CD37", "CD22", "CD74", "HLA-DRA"),
  Plasma = c("MZB1", "JCHAIN", "SDC1", "XBP1", "IGKC", "DERL3"),
  Myeloid = c("LST1", "TYROBP", "FCER1G", "AIF1", "CTSS", "LILRB1", "CTSD"),
  DC = c("CD1C", "FCER1A", "CLEC10A", "CST3", "HLA-DPA1", "HLA-DRA"),
  Mast = c("TPSAB1", "TPSB2", "CPA3", "KIT", "MS4A2", "HDC"),
  Epithelial = c("EPCAM", "KRT8", "KRT18", "KRT19", "KRT14", "KRT5", "SFN"),
  Fibroblast = c("COL1A1", "COL1A2", "DCN", "COL3A1", "LUM", "COL6A1"),
  Endothelial = c("PECAM1", "VWF", "EMCN", "KDR", "RAMP2", "PLVAP")
)
lineage_core <- list(
  T_cell = c("CD3D", "CD3E", "TRAC", "CD247", "LCK"),
  NK = c("NKG7", "GNLY", "KLRD1", "FCGR3A"),
  B_cell = c("MS4A1", "CD79A", "CD37", "CD22"),
  Plasma = c("MZB1", "JCHAIN", "SDC1", "XBP1", "DERL3"),
  Myeloid = c("LST1", "TYROBP", "FCER1G", "AIF1", "CTSS", "LILRB1"),
  DC = c("CD1C", "FCER1A", "CLEC10A"),
  Mast = c("TPSAB1", "TPSB2", "CPA3", "KIT"),
  Epithelial = c("EPCAM", "KRT8", "KRT18", "KRT19", "KRT14", "KRT5"),
  Fibroblast = c("COL1A1", "COL1A2", "DCN", "COL3A1", "LUM"),
  Endothelial = c("PECAM1", "VWF", "EMCN", "KDR", "RAMP2", "PLVAP")
)

immune_lineages <- c("T_cell", "NK", "B_cell", "Plasma", "Myeloid", "DC", "Mast")
nonimmune_lineages <- c("Epithelial", "Fibroblast", "Endothelial")
lymphoid_lineages <- c("T_cell", "NK", "B_cell", "Plasma")

allowed_lineages <- function(fraction) {
  if (fraction == "CD3") return(c("T_cell", "NK"))
  if (fraction == "CD45pos") return(immune_lineages)
  if (fraction == "CD45neg") return(nonimmune_lineages)
  c(immune_lineages, nonimmune_lineages)
}

marker_score <- function(norm, genes) {
  present <- intersect(genes, rownames(norm))
  if (!length(present)) return(rep(0, ncol(norm)))
  as.numeric(Matrix::colMeans(norm[present, , drop = FALSE]))
}

marker_detected <- function(counts, genes) {
  present <- intersect(genes, rownames(counts))
  if (!length(present)) return(rep(0L, ncol(counts)))
  as.integer(Matrix::colSums(counts[present, , drop = FALSE] > 0))
}

find_tcr_barcodes <- function(sample_row) {
  if (sample_row$assay_or_library != "5prime") return(character())
  hit <- tcr_manifest[
    tcr_manifest$patient_id == sample_row$patient_id &
      tcr_manifest$timepoint_harmonized == sample_row$timepoint_harmonized &
      tcr_manifest$sort_or_cell_fraction == sample_row$sort_or_cell_fraction,
    ,
    drop = FALSE
  ]
  if (!nrow(hit)) return(character())
  path <- hit$extracted_path[[1]]
  if (!file.exists(path)) return(character())
  tab <- tryCatch(
    read.csv(gzfile(path), stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) NULL
  )
  if (is.null(tab) || !nrow(tab) || !"barcode" %in% names(tab)) return(character())
  keep <- rep(TRUE, nrow(tab))
  if ("is_cell" %in% names(tab)) {
    keep <- keep & tab$is_cell %in% c(TRUE, "True", "TRUE", "true")
  }
  if ("high_confidence" %in% names(tab)) {
    keep <- keep & tab$high_confidence %in% c(TRUE, "True", "TRUE", "true")
  }
  if ("productive" %in% names(tab)) {
    keep <- keep & tab$productive %in% c(TRUE, "True", "TRUE", "true")
  }
  unique(tab$barcode[keep & !is.na(tab$barcode)])
}

module_genes <- unique(unlist(strsplit(
  paste(modules$genes_defined, collapse = ";"), ";", fixed = TRUE
)))
module_genes <- sort(unique(trimws(module_genes[nzchar(module_genes)])))

meta_final <- file.path(out_dir, "GSE301741_RAW_REBUILT_CELL_METADATA.csv.gz")
meta_partial <- paste0(meta_final, ".partial")
if (file.exists(meta_partial)) unlink(meta_partial)
meta_con <- gzfile(meta_partial, open = "wt", compression = 6)
meta_header <- TRUE

qc_rows <- list()
pseudobulk_rows <- list()

for (idx in seq_len(nrow(scrna))) {
  row <- scrna[idx, , drop = FALSE]
  message(
    sprintf(
      "[%02d/%02d] %s %s %s %s",
      idx, nrow(scrna), row$geo_accession, row$patient_id,
      row$timepoint_harmonized, row$sort_or_cell_fraction
    )
  )
  counts <- Read10X_h5(row$extracted_path, use.names = TRUE, unique.features = TRUE)
  if (is.list(counts)) {
    if ("Gene Expression" %in% names(counts)) {
      counts <- counts[["Gene Expression"]]
    } else {
      counts <- counts[[1]]
    }
  }
  counts <- as(counts, "dgCMatrix")
  lib_size <- Matrix::colSums(counts)
  keep_nonzero <- lib_size > 0
  counts <- counts[, keep_nonzero, drop = FALSE]
  lib_size <- lib_size[keep_nonzero]
  n_feature <- Matrix::colSums(counts > 0)
  mt_genes <- grep("^MT-", rownames(counts), ignore.case = TRUE, value = TRUE)
  mt_count <- if (length(mt_genes)) {
    Matrix::colSums(counts[mt_genes, , drop = FALSE])
  } else {
    rep(0, ncol(counts))
  }
  percent_mt <- 100 * mt_count / lib_size

  norm <- Matrix::t(Matrix::t(counts) / lib_size * 10000)
  norm@x <- log1p(norm@x)

  score_mat <- vapply(
    lineage_markers,
    function(genes) marker_score(norm, genes),
    numeric(ncol(norm))
  )
  detected_mat <- vapply(
    lineage_core,
    function(genes) marker_detected(counts, genes),
    integer(ncol(counts))
  )
  colnames(score_mat) <- names(lineage_markers)
  colnames(detected_mat) <- names(lineage_core)

  eligible <- matrix(
    FALSE, nrow = ncol(counts), ncol = length(lineage_markers),
    dimnames = list(colnames(counts), names(lineage_markers))
  )
  eligible[, "T_cell"] <- detected_mat[, "T_cell"] >= 1
  if (nrow(counts[intersect(c("GNLY", "KLRD1", "FCGR3A"), rownames(counts)), , drop = FALSE])) {
    nk_specific <- Matrix::colSums(
      counts[intersect(c("GNLY", "KLRD1", "FCGR3A"), rownames(counts)), , drop = FALSE] > 0
    )
    eligible[, "NK"] <- detected_mat[, "NK"] >= 2 & nk_specific >= 1
  } else {
    eligible[, "NK"] <- FALSE
  }
  eligible[, "B_cell"] <- detected_mat[, "B_cell"] >= 1
  eligible[, "Plasma"] <- detected_mat[, "Plasma"] >= 1
  eligible[, "Myeloid"] <- detected_mat[, "Myeloid"] >= 2
  dc_support <- marker_detected(counts, c("LST1", "CST3", "HLA-DRA"))
  eligible[, "DC"] <- detected_mat[, "DC"] >= 1 & dc_support >= 1
  eligible[, "Mast"] <- detected_mat[, "Mast"] >= 1
  eligible[, "Epithelial"] <- detected_mat[, "Epithelial"] >= 2
  eligible[, "Fibroblast"] <- detected_mat[, "Fibroblast"] >= 2
  eligible[, "Endothelial"] <- detected_mat[, "Endothelial"] >= 2

  allowed <- allowed_lineages(row$sort_or_cell_fraction)
  eligible[, setdiff(colnames(eligible), allowed)] <- FALSE

  tcr_barcodes <- find_tcr_barcodes(row)
  tcr_positive <- colnames(counts) %in% tcr_barcodes
  if ("T_cell" %in% allowed) eligible[tcr_positive, "T_cell"] <- TRUE

  candidate_scores <- score_mat
  candidate_scores[!eligible] <- -Inf
  best_index <- max.col(candidate_scores, ties.method = "first")
  top_score <- candidate_scores[cbind(seq_len(nrow(candidate_scores)), best_index)]
  second_score <- apply(candidate_scores, 1, function(x) {
    finite <- sort(x[is.finite(x)], decreasing = TRUE)
    if (length(finite) >= 2) finite[[2]] else 0
  })
  top_label <- colnames(candidate_scores)[best_index]
  top_label[!is.finite(top_score)] <- "Unresolved"
  top_score[!is.finite(top_score)] <- 0
  score_ratio <- (top_score + 1e-8) / (second_score + 1e-8)
  confident <- top_label != "Unresolved" &
    (score_ratio >= 1.15 | second_score <= 0)

  # Productive TCR is the strongest available discriminator between T and NK.
  top_label[tcr_positive & "T_cell" %in% allowed] <- "T_cell"
  confident[tcr_positive & "T_cell" %in% allowed] <- TRUE
  top_label[!confident] <- "Unresolved"

  lymphoid <- top_label %in% lymphoid_lineages
  feature_gate <- ifelse(lymphoid, 500L, 1000L)
  qc_pass <- n_feature >= feature_gate & percent_mt <= 20 &
    top_label != "Unresolved"
  final_label <- ifelse(qc_pass, top_label, "Filtered_or_unresolved")
  target_lineage <- ifelse(
    final_label == "T_cell", "T_cell",
    ifelse(final_label %in% c("Myeloid", "DC"), "Myeloid", "Other")
  )

  metadata <- data.frame(
    cell_id = paste(row$geo_accession, colnames(counts), sep = "__"),
    barcode = colnames(counts),
    geo_accession = row$geo_accession,
    patient_id = row$patient_id,
    timepoint = row$timepoint_harmonized,
    fraction = row$sort_or_cell_fraction,
    library = row$assay_or_library,
    n_count = as.numeric(lib_size),
    n_feature = as.integer(n_feature),
    percent_mt = as.numeric(percent_mt),
    tcr_productive_high_confidence = tcr_positive,
    preliminary_lineage = top_label,
    final_lineage = final_label,
    target_lineage = target_lineage,
    top_marker_score = as.numeric(top_score),
    second_marker_score = as.numeric(second_score),
    marker_score_ratio = as.numeric(score_ratio),
    t_marker_score = as.numeric(score_mat[, "T_cell"]),
    myeloid_marker_score = as.numeric(score_mat[, "Myeloid"]),
    annotation_method = "RAW_marker_reconstruction",
    stringsAsFactors = FALSE
  )
  write.table(
    metadata, meta_con, sep = ",", row.names = FALSE,
    col.names = meta_header, quote = TRUE, append = !meta_header
  )
  meta_header <- FALSE

  sample_counts <- table(factor(final_label, levels = c(
    names(lineage_markers), "Filtered_or_unresolved"
  )))
  for (label in names(sample_counts)) {
    qc_rows[[length(qc_rows) + 1L]] <- data.frame(
      geo_accession = row$geo_accession,
      patient_id = row$patient_id,
      timepoint = row$timepoint_harmonized,
      fraction = row$sort_or_cell_fraction,
      library = row$assay_or_library,
      lineage = label,
      n_cells = as.integer(sample_counts[[label]]),
      total_cells_loaded = ncol(counts),
      median_n_count = median(lib_size),
      median_n_feature = median(n_feature),
      median_percent_mt = median(percent_mt),
      stringsAsFactors = FALSE
    )
  }

  for (target in c("T_cell", "Myeloid")) {
    cells <- which(target_lineage == target)
    present <- intersect(module_genes, rownames(norm))
    values <- rep(NA_real_, length(module_genes))
    names(values) <- module_genes
    if (length(cells) && length(present)) {
      values[present] <- as.numeric(Matrix::rowMeans(norm[present, cells, drop = FALSE]))
    }
    pseudobulk_rows[[length(pseudobulk_rows) + 1L]] <- data.frame(
      geo_accession = row$geo_accession,
      patient_id = row$patient_id,
      timepoint = row$timepoint_harmonized,
      fraction = row$sort_or_cell_fraction,
      library = row$assay_or_library,
      target_lineage = target,
      n_target_cells = length(cells),
      gene = names(values),
      mean_log1p_cp10k = as.numeric(values),
      stringsAsFactors = FALSE
    )
  }

  rm(counts, norm, score_mat, detected_mat, eligible, candidate_scores, metadata)
  invisible(gc())
}

close(meta_con)
if (file.exists(meta_final)) unlink(meta_final)
if (!file.rename(meta_partial, meta_final)) {
  stop("Could not promote cell metadata output: ", meta_partial)
}

qc <- do.call(rbind, qc_rows)
pseudobulk <- do.call(rbind, pseudobulk_rows)
write.csv(
  qc,
  file.path(out_dir, "GSE301741_RAW_REBUILT_SAMPLE_LINEAGE_QC.csv"),
  row.names = FALSE
)
write.csv(
  pseudobulk,
  gzfile(file.path(out_dir, "GSE301741_RAW_TARGET_LINEAGE_PSEUDOBULK.csv.gz")),
  row.names = FALSE
)

summary_by_lineage <- aggregate(n_cells ~ lineage, qc, sum)
summary_by_lineage <- summary_by_lineage[
  order(summary_by_lineage$n_cells, decreasing = TRUE), ,
  drop = FALSE
]
report <- c(
  "# GSE301741 RAW Cell-Metadata Reconstruction",
  "",
  paste0("Created: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Scope",
  "",
  paste0("- scRNA H5 files processed independently: ", nrow(scrna)),
  paste0("- Raw Cell Ranger-filtered barcodes loaded: ", sum(qc$n_cells)),
  paste0("- Conservatively retained cells: ",
         sum(qc$n_cells[qc$lineage != "Filtered_or_unresolved"])),
  "- Peak memory is bounded by one H5 sample rather than the deposited 12.35 GB RDS.",
  "- Broad labels are reconstructed from prespecified canonical markers and sort-fraction constraints.",
  "- These are not the unavailable original author labels and must remain provenance-labelled as reconstructed.",
  "- QC follows the paper's 500-feature lymphoid, 1,000-feature non-lymphoid, and 20% mitochondrial thresholds.",
  "- Ambiguous marker assignments with top/second score ratio below 1.15 are excluded.",
  "- Productive high-confidence TCR calls are used to distinguish T cells from NK cells for matching 5-prime libraries.",
  "",
  "## Reconstructed Cell Counts",
  "",
  "| Lineage | Cells |",
  "|---|---:|"
)
for (i in seq_len(nrow(summary_by_lineage))) {
  report <- c(
    report,
    sprintf("| %s | %s |", summary_by_lineage$lineage[i],
            format(summary_by_lineage$n_cells[i], big.mark = ",", scientific = FALSE))
  )
}
report <- c(
  report,
  "",
  "## Interpretation Gate",
  "",
  "The reconstructed labels support lineage-aware external validation and sensitivity analysis. They do not recreate the authors' exact UMAP/Louvain assignments, CellBender correction, Scrublet calls, 3-prime/5-prime centering, or CNA-based malignant labels.",
  "",
  "## Outputs",
  "",
  "- `GSE301741_RAW_REBUILT_CELL_METADATA.csv.gz`",
  "- `GSE301741_RAW_REBUILT_SAMPLE_LINEAGE_QC.csv`",
  "- `GSE301741_RAW_TARGET_LINEAGE_PSEUDOBULK.csv.gz`"
)
writeLines(
  report,
  file.path(out_dir, "GSE301741_RAW_RECONSTRUCTION_REPORT.md")
)
writeLines(
  c(
    paste0("completed=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
    paste0("scrna_files=", nrow(scrna)),
    paste0("cells=", sum(qc$n_cells))
  ),
  file.path(out_dir, "GSE301741_RAW_RECONSTRUCTION_SUCCESS.txt")
)

cat(paste(report, collapse = "\n"), "\n")
