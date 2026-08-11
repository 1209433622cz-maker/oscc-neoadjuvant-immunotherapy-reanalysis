args <- commandArgs(trailingOnly = TRUE)
extract_archive <- "--extract" %in% args

workspace <- Sys.getenv("GSE200996_WORKSPACE")
if (workspace == "") {
  workspace <- normalizePath(getwd(), mustWork = FALSE)
}

tar_path <- file.path(
  workspace,
  "00_raw_data",
  "external_validation",
  "GSE301741",
  "GSE301741_RAW.tar"
)
extract_dir <- file.path(
  workspace,
  "00_raw_data",
  "external_validation",
  "GSE301741",
  "RAW_extracted"
)
manifest_path <- file.path(
  workspace,
  "03_rebuild",
  "manifests",
  "EXTERNAL_COHORT_SAMPLE_MANIFEST.csv"
)
out_dir <- file.path(workspace, "03_rebuild", "manifests")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(tar_path)) {
  stop("GSE301741 raw tar not found: ", tar_path)
}
if (!file.exists(manifest_path)) {
  stop("External sample manifest not found: ", manifest_path)
}

tar_members <- utils::untar(tar_path, list = TRUE)
tar_df <- data.frame(
  raw_archive_file_name = basename(tar_members),
  raw_archive_member = tar_members,
  stringsAsFactors = FALSE
)

if (extract_archive) {
  dir.create(extract_dir, recursive = TRUE, showWarnings = FALSE)
  message("Extracting GSE301741_RAW.tar to: ", extract_dir)
  utils::untar(tar_path, exdir = extract_dir)
} else if (!dir.exists(extract_dir)) {
  message("RAW extraction directory does not exist yet: ", extract_dir)
  message("Run with --extract, or use the PowerShell wrapper with -Extract, when the machine is idle.")
}

sample_manifest <- utils::read.csv(manifest_path, stringsAsFactors = FALSE, check.names = FALSE)
sample_manifest <- sample_manifest[sample_manifest$accession == "GSE301741", , drop = FALSE]
sample_manifest <- merge(
  sample_manifest,
  tar_df,
  by = "raw_archive_file_name",
  all.x = TRUE,
  sort = FALSE
)

sample_manifest$extracted_path <- ifelse(
  !is.na(sample_manifest$raw_archive_member),
  file.path(extract_dir, sample_manifest$raw_archive_member),
  ""
)
sample_manifest$extracted_exists <- file.exists(sample_manifest$extracted_path)

summarize_h5 <- function(path) {
  empty <- data.frame(
    n_features = NA_integer_,
    n_barcodes = NA_integer_,
    nnz = NA_integer_,
    h5_status = "not_checked",
    stringsAsFactors = FALSE
  )
  if (!file.exists(path)) {
    empty$h5_status <- "missing"
    return(empty)
  }
  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    empty$h5_status <- "hdf5r_missing"
    return(empty)
  }
  out <- tryCatch({
    h5 <- hdf5r::H5File$new(path, mode = "r")
    on.exit(h5$close_all(), add = TRUE)
    shape <- h5[["matrix/shape"]]$read()
    data_dims <- h5[["matrix/data"]]$dims
    data.frame(
      n_features = as.integer(shape[1]),
      n_barcodes = as.integer(shape[2]),
      nnz = as.integer(data_dims[1]),
      h5_status = "ok",
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    data.frame(
      n_features = NA_integer_,
      n_barcodes = NA_integer_,
      nnz = NA_integer_,
      h5_status = paste0("error: ", conditionMessage(e)),
      stringsAsFactors = FALSE
    )
  })
  out
}

count_gz_lines <- function(path) {
  if (!file.exists(path)) {
    return(data.frame(tcr_rows = NA_integer_, tcr_status = "missing", stringsAsFactors = FALSE))
  }
  out <- tryCatch({
    con <- gzfile(path, open = "rt")
    on.exit(close(con), add = TRUE)
    n <- 0L
    repeat {
      lines <- readLines(con, n = 100000L, warn = FALSE)
      if (!length(lines)) break
      n <- n + length(lines)
    }
    data.frame(
      tcr_rows = max(0L, n - 1L),
      tcr_status = "ok",
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    data.frame(
      tcr_rows = NA_integer_,
      tcr_status = paste0("error: ", conditionMessage(e)),
      stringsAsFactors = FALSE
    )
  })
  out
}

qc_rows <- list()
for (i in seq_len(nrow(sample_manifest))) {
  row <- sample_manifest[i, , drop = FALSE]
  path <- row$extracted_path
  if (!isTRUE(row$extracted_exists)) {
    qc <- data.frame(
      n_features = NA_integer_,
      n_barcodes = NA_integer_,
      nnz = NA_integer_,
      h5_status = ifelse(row$modality == "scRNA", "not_extracted", "not_applicable"),
      tcr_rows = NA_integer_,
      tcr_status = ifelse(row$modality == "scTCR", "not_extracted", "not_applicable"),
      stringsAsFactors = FALSE
    )
  } else if (grepl("\\.h5$", path, ignore.case = TRUE)) {
    h5_qc <- summarize_h5(path)
    qc <- cbind(
      h5_qc,
      data.frame(tcr_rows = NA_integer_, tcr_status = "not_applicable", stringsAsFactors = FALSE)
    )
  } else if (grepl("\\.csv\\.gz$", path, ignore.case = TRUE)) {
    tcr_qc <- count_gz_lines(path)
    qc <- cbind(
      data.frame(n_features = NA_integer_, n_barcodes = NA_integer_, nnz = NA_integer_, h5_status = "not_applicable", stringsAsFactors = FALSE),
      tcr_qc
    )
  } else {
    qc <- data.frame(
      n_features = NA_integer_,
      n_barcodes = NA_integer_,
      nnz = NA_integer_,
      h5_status = "unknown_file_type",
      tcr_rows = NA_integer_,
      tcr_status = "unknown_file_type",
      stringsAsFactors = FALSE
    )
  }
  qc_rows[[i]] <- cbind(row, qc)
}

qc_df <- do.call(rbind, qc_rows)
out_csv <- file.path(out_dir, "GSE301741_RAW_ROUTE_SAMPLE_QC.csv")
utils::write.csv(qc_df, out_csv, row.names = FALSE)

summary_df <- aggregate(
  raw_archive_file_name ~ modality + timepoint_harmonized + sort_or_cell_fraction + extracted_exists,
  data = qc_df,
  FUN = length
)
names(summary_df)[names(summary_df) == "raw_archive_file_name"] <- "n_samples"
out_summary <- file.path(out_dir, "GSE301741_RAW_ROUTE_SUMMARY.csv")
utils::write.csv(summary_df, out_summary, row.names = FALSE)

report <- c(
  "# GSE301741 Raw Route Preparation",
  "",
  paste0("Created: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  paste0("- Raw tar: `", tar_path, "`"),
  paste0("- Extract directory: `", extract_dir, "`"),
  paste0("- Extract requested: `", extract_archive, "`"),
  paste0("- Manifest rows: ", nrow(qc_df)),
  paste0("- Extracted files found: ", sum(qc_df$extracted_exists)),
  paste0("- scRNA H5 rows: ", sum(qc_df$modality == "scRNA")),
  paste0("- scTCR rows: ", sum(qc_df$modality == "scTCR")),
  "",
  "## Files Written",
  "",
  "- `03_rebuild/manifests/GSE301741_RAW_ROUTE_SAMPLE_QC.csv`",
  "- `03_rebuild/manifests/GSE301741_RAW_ROUTE_SUMMARY.csv`"
)
out_md <- file.path(out_dir, "GSE301741_RAW_ROUTE_PREP_REPORT.md")
writeLines(report, out_md)
cat(paste(report, collapse = "\n"))
cat("\n")
