#!/usr/bin/env Rscript

rm(list = ls())
gc()
options(warn = 1)

workspace <- Sys.getenv("GSE200996_WORKSPACE")
if (workspace == "") {
  workspace <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

rds_path <- Sys.getenv("GSE301741_RDS_PATH")
if (rds_path == "") {
  rds_path <- file.path(
    workspace,
    "00_raw_data",
    "external_validation",
    "GSE301741",
    "GSE301741_Seurat_Object_QCpass_137020cells_withMetaData.rds"
  )
}
out_dir <- Sys.getenv("GSE301741_METADATA_OUTPUT_DIR")
if (out_dir == "") out_dir <- file.path(workspace, "03_rebuild", "manifests")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(rds_path)) {
  stop("GSE301741 RDS not found: ", rds_path)
}

atomic_save_rds <- function(value, path) {
  temporary <- paste0(path, ".partial")
  if (file.exists(temporary)) unlink(temporary)
  saveRDS(value, temporary, compress = "gzip")
  if (file.exists(path)) unlink(path)
  if (!file.rename(temporary, path)) stop("Could not finalize: ", path)
}

atomic_write_csv <- function(value, path, compress = FALSE) {
  temporary <- paste0(path, ".partial")
  if (file.exists(temporary)) unlink(temporary)
  if (compress) {
    connection <- gzfile(temporary, open = "wt", compression = 6)
    on.exit(close(connection), add = TRUE)
    utils::write.csv(value, connection, row.names = FALSE, na = "")
    close(connection)
    on.exit(NULL, add = FALSE)
  } else {
    utils::write.csv(value, temporary, row.names = FALSE, na = "")
  }
  if (file.exists(path)) unlink(path)
  if (!file.rename(temporary, path)) stop("Could not finalize: ", path)
}

timestamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
message("[", timestamp(), "] Input: ", rds_path)
message("[", timestamp(), "] File size GB: ", round(file.info(rds_path)$size / 1024^3, 3))

header_connection <- file(rds_path, open = "rb")
header <- readBin(header_connection, what = "raw", n = 16)
close(header_connection)
header_hex <- paste(sprintf("%02X", as.integer(header)), collapse = "-")
message("[", timestamp(), "] Serialization header: ", header_hex)
message("[", timestamp(), "] Loading the RDS in an isolated process.")

suppressPackageStartupMessages(library(SeuratObject))
load_started <- Sys.time()
obj <- readRDS(rds_path)
load_seconds <- as.numeric(difftime(Sys.time(), load_started, units = "secs"))
message("[", timestamp(), "] RDS load complete in ", round(load_seconds, 1), " seconds.")
message("[", timestamp(), "] Object class: ", paste(class(obj), collapse = ";"))
message("[", timestamp(), "] Reported object.size GB: ", round(as.numeric(object.size(obj)) / 1024^3, 3))

if (inherits(obj, "Seurat")) {
  meta <- obj@meta.data
  object_summary <- data.frame(
    metric = c(
      "class",
      "cells",
      "metadata_columns",
      "assays",
      "reductions",
      "object_size_bytes",
      "rds_size_bytes",
      "rds_header_hex",
      "load_seconds"
    ),
    value = c(
      paste(class(obj), collapse = ";"),
      nrow(meta),
      ncol(meta),
      paste(names(obj@assays), collapse = ";"),
      paste(names(obj@reductions), collapse = ";"),
      as.character(as.numeric(object.size(obj))),
      as.character(file.info(rds_path)$size),
      header_hex,
      as.character(load_seconds)
    ),
    stringsAsFactors = FALSE
  )
} else if (is.list(obj) && !is.null(obj$meta.data)) {
  meta <- obj$meta.data
  object_summary <- data.frame(
    metric = c("class", "cells", "metadata_columns", "object_size_bytes", "rds_size_bytes", "rds_header_hex", "load_seconds"),
    value = c(
      paste(class(obj), collapse = ";"),
      nrow(meta),
      ncol(meta),
      as.character(as.numeric(object.size(obj))),
      as.character(file.info(rds_path)$size),
      header_hex,
      as.character(load_seconds)
    ),
    stringsAsFactors = FALSE
  )
} else if (is.data.frame(obj)) {
  meta <- obj
  object_summary <- data.frame(
    metric = c("class", "cells", "metadata_columns", "object_size_bytes", "rds_size_bytes", "rds_header_hex", "load_seconds"),
    value = c(
      paste(class(obj), collapse = ";"),
      nrow(meta),
      ncol(meta),
      as.character(as.numeric(object.size(obj))),
      as.character(file.info(rds_path)$size),
      header_hex,
      as.character(load_seconds)
    ),
    stringsAsFactors = FALSE
  )
} else {
  stop("Could not find a metadata table in the loaded object.")
}

if (is.null(rownames(meta)) || any(!nzchar(rownames(meta)))) {
  stop("The metadata table does not have complete cell row names.")
}

meta_export <- data.frame(
  cell_barcode = rownames(meta),
  meta,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

field_summary <- do.call(
  rbind,
  lapply(names(meta_export), function(column) {
    values <- meta_export[[column]]
    data.frame(
      field = column,
      storage_class = paste(class(values), collapse = ";"),
      n_non_missing = sum(!is.na(values) & nzchar(as.character(values))),
      n_unique_non_missing = length(unique(as.character(values[!is.na(values) & nzchar(as.character(values))]))),
      stringsAsFactors = FALSE
    )
  })
)

candidate_pattern <- paste(
  c(
    "patient", "sample", "response", "outcome", "ptr", "path",
    "treat", "time", "site", "tumou?r", "orig.ident",
    "cell.?type", "annotation", "cluster", "histology", "disease"
  ),
  collapse = "|"
)
candidate_fields <- grep(candidate_pattern, names(meta), ignore.case = TRUE, value = TRUE)
compact_fields <- candidate_fields[
  vapply(
    candidate_fields,
    function(column) {
      n_unique <- field_summary$n_unique_non_missing[field_summary$field == column]
      length(n_unique) == 1 && n_unique <= 500
    },
    logical(1)
  )
]

value_count_rows <- list()
for (column in compact_fields) {
  values <- as.character(meta[[column]])
  values[is.na(values) | !nzchar(values)] <- "<missing>"
  counts <- sort(table(values, useNA = "ifany"), decreasing = TRUE)
  value_count_rows[[length(value_count_rows) + 1]] <- data.frame(
    field = column,
    value = names(counts),
    n_cells = as.integer(counts),
    stringsAsFactors = FALSE
  )
}
value_counts <- if (length(value_count_rows)) {
  do.call(rbind, value_count_rows)
} else {
  data.frame(field = character(), value = character(), n_cells = integer())
}

compact_metadata <- if (length(compact_fields)) {
  unique(meta[, compact_fields, drop = FALSE])
} else {
  data.frame(note = "No compact candidate metadata fields detected.", stringsAsFactors = FALSE)
}

metadata_csv <- file.path(out_dir, "GSE301741_RDS_CELL_METADATA.csv.gz")
metadata_rds <- file.path(out_dir, "GSE301741_RDS_CELL_METADATA.rds")
field_csv <- file.path(out_dir, "GSE301741_RDS_METADATA_FIELD_SUMMARY.csv")
value_csv <- file.path(out_dir, "GSE301741_RDS_METADATA_VALUE_COUNTS.csv")
compact_csv <- file.path(out_dir, "GSE301741_RDS_COMPACT_CLINICAL_CELLTYPE_FIELDS.csv")
object_csv <- file.path(out_dir, "GSE301741_RDS_OBJECT_STRUCTURE.csv")
report_path <- file.path(out_dir, "GSE301741_RDS_METADATA_EXTRACTION_REPORT.md")
success_path <- file.path(out_dir, "GSE301741_RDS_METADATA_EXTRACTION_SUCCESS.txt")

message("[", timestamp(), "] Writing metadata outputs.")
atomic_write_csv(meta_export, metadata_csv, compress = TRUE)
atomic_save_rds(meta_export, metadata_rds)
atomic_write_csv(field_summary, field_csv)
atomic_write_csv(value_counts, value_csv)
atomic_write_csv(compact_metadata, compact_csv)
atomic_write_csv(object_summary, object_csv)

response_fields <- grep("response|outcome|ptr|path", names(meta), ignore.case = TRUE, value = TRUE)
report <- c(
  "# GSE301741 RDS Metadata Extraction",
  "",
  paste0("Completed: ", timestamp()),
  "",
  "## Input",
  "",
  paste0("- RDS: `", rds_path, "`"),
  paste0("- Size: ", round(file.info(rds_path)$size / 1024^3, 3), " GB"),
  paste0("- Serialization header: `", header_hex, "`"),
  paste0("- Load time: ", round(load_seconds, 1), " seconds"),
  "",
  "## Object",
  "",
  paste0("- Class: ", paste(class(obj), collapse = ";")),
  paste0("- Cells in metadata: ", nrow(meta)),
  paste0("- Metadata columns: ", ncol(meta)),
  paste0("- Reported object.size: ", round(as.numeric(object.size(obj)) / 1024^3, 3), " GB"),
  paste0("- Candidate response/outcome fields: ", if (length(response_fields)) paste(response_fields, collapse = "; ") else "none"),
  "",
  "## Outputs",
  "",
  paste0("- Cell metadata CSV: `", metadata_csv, "`"),
  paste0("- Cell metadata RDS: `", metadata_rds, "`"),
  paste0("- Field summary: `", field_csv, "`"),
  paste0("- Value counts: `", value_csv, "`"),
  paste0("- Compact clinical/cell-type fields: `", compact_csv, "`"),
  paste0("- Object structure: `", object_csv, "`")
)
writeLines(report, report_path, useBytes = TRUE)
writeLines(
  c(
    paste0("completed=", timestamp()),
    paste0("cells=", nrow(meta)),
    paste0("metadata_columns=", ncol(meta)),
    paste0("response_fields=", paste(response_fields, collapse = ";"))
  ),
  success_path,
  useBytes = TRUE
)

message("[", timestamp(), "] Metadata extraction complete.")
rm(meta_export, meta, obj)
gc()
