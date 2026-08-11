#!/usr/bin/env Rscript

# Compatibility wrapper. The canonical supplementary-table builder is the
# Python/openpyxl script because it produces an xlsx that validates cleanly
# with both Excel-style readers and the bundled Python runtime.

get_script_path <- function() {
  cmd <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", cmd, value = TRUE)
  if (length(hit)) return(normalizePath(sub("^--file=", "", hit[1]), winslash = "/", mustWork = TRUE))
  normalizePath("03_rebuild/analysis/07_make_supplementary_tables.R", winslash = "/", mustWork = TRUE)
}

script_path <- get_script_path()
py_script <- file.path(dirname(script_path), "07_make_supplementary_tables.py")
python_candidates <- c(
  file.path(Sys.getenv("USERPROFILE"), ".cache", "codex-runtimes", "codex-primary-runtime", "dependencies", "python", "python.exe"),
  Sys.which("python"),
  Sys.which("py")
)
python_candidates <- python_candidates[nzchar(python_candidates)]
python_exe <- python_candidates[file.exists(python_candidates)][1]

if (is.na(python_exe) || !nzchar(python_exe)) {
  stop("Could not find Python to run: ", py_script, call. = FALSE)
}

status <- system2(python_exe, shQuote(py_script))
quit(status = status)

if (FALSE) {

get_script_path <- function() {
  cmd <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", cmd, value = TRUE)
  if (length(hit)) return(normalizePath(sub("^--file=", "", hit[1]), winslash = "/", mustWork = TRUE))
  normalizePath("03_rebuild/analysis/07_make_supplementary_tables.R", winslash = "/", mustWork = TRUE)
}

script_path <- get_script_path()
rebuild_dir <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
workspace_dir <- normalizePath(file.path(rebuild_dir, ".."), winslash = "/", mustWork = TRUE)
results_dir <- file.path(rebuild_dir, "results")
tables_dir <- file.path(rebuild_dir, "tables", "submission")
csv_dir <- file.path(tables_dir, "csv")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)

required_packages <- c("readr", "dplyr", "tibble", "openxlsx")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(openxlsx)
})

path <- function(...) file.path(...)

read_csv_required <- function(rel_path, label) {
  f <- file.path(results_dir, rel_path)
  if (!file.exists(f)) stop("Required file not found for ", label, ": ", f, call. = FALSE)
  readr::read_csv(f, show_col_types = FALSE, na = c("", "NA", "NaN"))
}

read_matrix_as_table <- function(rel_path, label) {
  f <- file.path(results_dir, rel_path)
  if (!file.exists(f)) stop("Required file not found for ", label, ": ", f, call. = FALSE)
  read.csv(f, check.names = FALSE, row.names = 1) %>%
    tibble::rownames_to_column("gene")
}

sheet_name <- function(x) {
  out <- gsub("[\\[\\]\\*\\?/\\\\:]", "_", x)
  substr(out, 1, 31)
}

safe_file_stem <- function(x) {
  gsub("[^A-Za-z0-9_]+", "_", x)
}

table_specs <- tibble::tribble(
  ~sheet, ~rel_path, ~description,
  "S1_raw_h5_manifest", "data_audit/raw_h5_manifest.csv", "Raw h5 manifest retained after workspace cleanup.",
  "S1_patient_presence", "data_audit/cd45_tumor_patient_timepoint_presence.csv", "CD45 tumor patient pre/post availability and response labels.",
  "S1_patient_response", "data_audit/cd45_tumor_patient_response_table.csv", "Patient-level response and cohort labels.",
  "S1_celltype_counts", "data_audit/cd45_tumor_celltype_counts.csv", "Overall CD45 tumor immune-cell type counts.",
  "S1_patient_cell_counts", "data_audit/cd45_tumor_patient_celltype_counts.csv", "Patient/timepoint/cell-type counts.",
  "S2_baseline_comp", "pre_baseline/baseline_pre_composition_trend_respOrd_limma_logit.csv", "Baseline pre-treatment composition ordinal response model.",
  "S2_dynamic_comp_ord", "dynamic_paired/Fig4A_composition_delta_logit_limma_respOrd_trend.csv", "Paired post-pre composition ordinal response model.",
  "S2_dynamic_comp_RvsNR", "dynamic_paired/Fig4A_composition_delta_logit_limma_RvsNR.csv", "Paired post-pre composition binary R/NR sensitivity model.",
  "S3_baseline_T_DE", "pre_baseline/Baseline_T_cell_DE_trend_respOrd.csv", "Baseline T-cell pseudobulk response-trend DE.",
  "S3_baseline_T_GSEA", "pre_baseline/Baseline_T_cell_GSEA_Hallmark.csv", "Baseline T-cell Hallmark GSEA.",
  "S3_baseline_M_DE", "pre_baseline/Baseline_Myeloid_DE_trend_respOrd.csv", "Baseline myeloid pseudobulk response-trend DE.",
  "S3_baseline_M_GSEA", "pre_baseline/Baseline_Myeloid_GSEA_Hallmark.csv", "Baseline myeloid Hallmark GSEA.",
  "S4_dynamic_T_DE", "dynamic_paired/Fig4B_T_cell_interaction_DE_trend.csv", "Dynamic T-cell post-by-response interaction DE.",
  "S4_dynamic_T_GSEA", "dynamic_paired/Fig4B_T_cell_GSEA_Hallmark.csv", "Dynamic T-cell Hallmark GSEA.",
  "S4_T_signature_delta", "dynamic_paired/Fig4B_T_cell_SignatureDelta_source.csv", "T-cell signature deltas in paired patients.",
  "S5_dynamic_M_DE", "dynamic_paired/Fig4B_Myeloid_interaction_DE_trend.csv", "Dynamic myeloid post-by-response interaction DE.",
  "S5_dynamic_M_GSEA", "dynamic_paired/Fig4B_Myeloid_GSEA_Hallmark.csv", "Dynamic myeloid Hallmark GSEA.",
  "S5_M_signature_delta", "dynamic_paired/Fig4B_Myeloid_SignatureDelta_source.csv", "Myeloid signature deltas in paired patients.",
  "S6_LOO_stability", "sensitivity_leave_one_out/LOO_model_stability_summary.csv", "Leave-one-patient model stability summary.",
  "S6_LOO_key_pathways", "sensitivity_leave_one_out/LOO_key_pathway_NES.csv", "Leave-one-patient key pathway NES.",
  "S6_LOO_direction", "sensitivity_leave_one_out/LOO_top25_gene_direction_concordance.csv", "Leave-one-patient top-25 direction concordance."
)

matrix_specs <- tibble::tribble(
  ~sheet, ~rel_path, ~description,
  "S4_T_delta_matrix", "dynamic_paired/Fig4B_Tcell_delta_matrix.csv", "T-cell top-gene post-pre delta matrix.",
  "S5_M_delta_matrix", "dynamic_paired/Fig4B_Myeloid_delta_matrix.csv", "Myeloid top-gene post-pre delta matrix."
)

message("Workspace: ", workspace_dir)
message("Results:   ", results_dir)
message("Tables:    ", tables_dir)

tables <- list()
manifest_rows <- list()

for (i in seq_len(nrow(table_specs))) {
  spec <- table_specs[i, ]
  dat <- read_csv_required(spec$rel_path, spec$sheet)
  tables[[spec$sheet]] <- dat
  manifest_rows[[length(manifest_rows) + 1]] <- tibble(
    sheet = spec$sheet,
    source_file = file.path(results_dir, spec$rel_path),
    n_rows = nrow(dat),
    n_cols = ncol(dat),
    description = spec$description
  )
}

for (i in seq_len(nrow(matrix_specs))) {
  spec <- matrix_specs[i, ]
  dat <- read_matrix_as_table(spec$rel_path, spec$sheet)
  tables[[spec$sheet]] <- dat
  manifest_rows[[length(manifest_rows) + 1]] <- tibble(
    sheet = spec$sheet,
    source_file = file.path(results_dir, spec$rel_path),
    n_rows = nrow(dat),
    n_cols = ncol(dat),
    description = spec$description
  )
}

manifest <- bind_rows(manifest_rows) %>%
  arrange(sheet)

tables <- c(list(README = manifest), tables)

for (nm in names(tables)) {
  readr::write_csv(tables[[nm]], file.path(csv_dir, paste0(safe_file_stem(nm), ".csv")))
}

wb <- openxlsx::createWorkbook()
for (nm in names(tables)) {
  s <- sheet_name(nm)
  openxlsx::addWorksheet(wb, s)
  dat <- tables[[nm]]
  openxlsx::writeDataTable(wb, s, dat, tableStyle = "TableStyleMedium2")
  openxlsx::freezePane(wb, s, firstRow = TRUE)
  openxlsx::setColWidths(wb, s, cols = seq_len(ncol(dat)), widths = "auto")
}

out_xlsx <- file.path(tables_dir, "Supplementary_Tables.xlsx")
openxlsx::saveWorkbook(wb, out_xlsx, overwrite = TRUE)

readr::write_csv(
  tibble(
    generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    workbook = out_xlsx,
    csv_directory = csv_dir,
    n_tables = length(tables) - 1
  ),
  file.path(tables_dir, "SUPPLEMENTARY_TABLES_MANIFEST.csv")
)

message("Done. Supplementary tables written to: ", out_xlsx)
}
