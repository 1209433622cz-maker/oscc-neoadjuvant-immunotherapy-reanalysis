#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
  library(ggplot2)
})

workspace <- Sys.getenv("GSE200996_WORKSPACE")
if (workspace == "") {
  workspace <- "H:/SCI2/OSCC-GSE200996-2025.12"
}
workspace <- normalizePath(workspace, winslash = "/", mustWork = TRUE)
raw_dir <- file.path(workspace, "00_raw_data", "GSE200996_RAW")
meta_dir <- file.path(workspace, "00_raw_data", "GSE200996_metadata")
out_dir <- file.path(workspace, "03_rebuild", "results", "data_audit")
fig_dir <- file.path(workspace, "03_rebuild", "figures", "data_audit")
log_dir <- file.path(workspace, "03_rebuild", "logs")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
sink(file.path(log_dir, paste0("01_data_audit_", stamp, ".log")), split = TRUE)
on.exit(sink(), add = TRUE)

message("Data audit started: ", Sys.time())
message("Workspace: ", workspace)

parse_h5 <- function(path) {
  name <- basename(path)
  type <- fifelse(grepl("tumor", name, ignore.case = TRUE), "tumor",
                  fifelse(grepl("PBMC", name, ignore.case = TRUE), "PBMC",
                          fifelse(grepl("sorted", name, ignore.case = TRUE), "sorted", "other")))
  sample_tag <- str_match(name, "raw_feature_bc_matrix_(.*?)_GEX_sc_")[, 2]
  patient_tokens <- str_extract_all(sample_tag %||% "", "P\\d+")[[1]]
  timepoint <- fifelse(grepl("pre-Tx", name, ignore.case = TRUE), "pre",
                       fifelse(grepl("post-Tx", name, ignore.case = TRUE), "post", NA_character_))
  data.table(
    file = name,
    path = normalizePath(path, winslash = "/", mustWork = TRUE),
    file_type = type,
    sample_tag = sample_tag,
    patient_tokens = paste(patient_tokens, collapse = ";"),
    n_patient_tokens = length(patient_tokens),
    timepoint = timepoint,
    size_mb = round(file.info(path)$size / 1024^2, 3)
  )
}

`%||%` <- function(a, b) if (!is.na(a) && length(a) > 0) a else b

h5_files <- list.files(raw_dir, pattern = "\\.h5$", full.names = TRUE)
h5_manifest <- rbindlist(lapply(h5_files, parse_h5), fill = TRUE)
setorder(h5_manifest, file_type, file)
fwrite(h5_manifest, file.path(out_dir, "raw_h5_manifest.csv"))

tumor_h5 <- h5_manifest[file_type == "tumor"]
tumor_patient_long <- tumor_h5[, .(
  patient_id = unlist(strsplit(patient_tokens, ";", fixed = TRUE))
), by = .(file, timepoint)]
tumor_patient_long <- tumor_patient_long[nchar(patient_id) > 0]
tumor_patient_presence <- dcast(unique(tumor_patient_long), patient_id ~ timepoint, value.var = "file", fun.aggregate = length)
if (!"pre" %in% names(tumor_patient_presence)) tumor_patient_presence[, pre := 0L]
if (!"post" %in% names(tumor_patient_presence)) tumor_patient_presence[, post := 0L]
tumor_patient_presence[, paired := pre > 0 & post > 0]
setorder(tumor_patient_presence, patient_id)
fwrite(tumor_patient_presence, file.path(out_dir, "raw_tumor_patient_timepoint_presence.csv"))

meta_cd45_tumor <- file.path(meta_dir, "GSE200996_CD45.tumor.single.cell.meta.data.txt.gz")
dt <- fread(meta_cd45_tumor)
setnames(
  dt,
  c("V1", "Patient_ID", "Stage", "Cohort", "Path_response", "CellType_ID", "UMAP_1", "UMAP_2"),
  c("cell_id", "patient_id_raw", "sample_stage", "cohort", "path_response", "celltype_id", "umap_1", "umap_2")
)

dt[, sample_tag := sub("^[^_]+_", "", cell_id)]
dt[, patient_id := str_extract(sample_tag, "P\\d+")]
dt[, timepoint := fifelse(grepl("Pre", sample_tag, ignore.case = TRUE), "pre",
                          fifelse(grepl("Post", sample_tag, ignore.case = TRUE), "post", NA_character_))]
dt[, timepoint := factor(timepoint, levels = c("pre", "post"), ordered = TRUE)]
dt[, response_ord := factor(path_response, levels = c("Low", "Medium", "High"), ordered = TRUE)]
dt[, response_bin := fifelse(path_response == "High", "R",
                             fifelse(path_response == "Low", "NR", NA_character_))]

patient_response <- unique(dt[, .(patient_id, cohort, path_response, response_bin)])
patient_response[, response_ord_num := match(path_response, c("Low", "Medium", "High"))]
setorder(patient_response, response_ord_num, patient_id)
fwrite(patient_response, file.path(out_dir, "cd45_tumor_patient_response_table.csv"))

patient_timepoint_counts <- dt[, .(n_cells = .N), by = .(patient_id, cohort, path_response, response_bin, timepoint)]
patient_timepoint_counts[, path_response := factor(path_response, levels = c("Low", "Medium", "High"), ordered = TRUE)]
setorder(patient_timepoint_counts, patient_id, timepoint)
fwrite(patient_timepoint_counts, file.path(out_dir, "cd45_tumor_patient_timepoint_cell_counts.csv"))

celltype_counts <- dt[, .(n_cells = .N), by = .(celltype_id)]
celltype_counts[, fraction := n_cells / sum(n_cells)]
setorder(celltype_counts, -n_cells)
fwrite(celltype_counts, file.path(out_dir, "cd45_tumor_celltype_counts.csv"))

patient_celltype_counts <- dt[, .(n_cells = .N), by = .(patient_id, timepoint, path_response, celltype_id)]
setorder(patient_celltype_counts, patient_id, timepoint, -n_cells)
fwrite(patient_celltype_counts, file.path(out_dir, "cd45_tumor_patient_celltype_counts.csv"))

meta_presence <- dcast(unique(dt[, .(patient_id, timepoint)]), patient_id ~ timepoint, value.var = "timepoint", fun.aggregate = length)
if (!"pre" %in% names(meta_presence)) meta_presence[, pre := 0L]
if (!"post" %in% names(meta_presence)) meta_presence[, post := 0L]
meta_presence[, paired := pre > 0 & post > 0]
meta_presence <- merge(meta_presence, patient_response[, .(patient_id, cohort, path_response, response_bin, response_ord_num)], by = "patient_id", all.x = TRUE)
setorder(meta_presence, response_ord_num, patient_id)
fwrite(meta_presence, file.path(out_dir, "cd45_tumor_patient_timepoint_presence.csv"))

pair_response_summary <- meta_presence[paired == TRUE, .N, by = path_response]
pair_response_summary[, response_ord_num := match(path_response, c("Low", "Medium", "High"))]
setorder(pair_response_summary, response_ord_num)
pair_response_summary[, response_ord_num := NULL]
fwrite(pair_response_summary, file.path(out_dir, "paired_patient_response_summary.csv"))

mismatch_raw_not_meta <- setdiff(tumor_patient_presence$patient_id, meta_presence$patient_id)
mismatch_meta_not_raw <- setdiff(meta_presence$patient_id, tumor_patient_presence$patient_id)
fwrite(data.table(raw_tumor_patient_not_in_meta = mismatch_raw_not_meta), file.path(out_dir, "mismatch_raw_tumor_patients_not_in_cd45_meta.csv"))
fwrite(data.table(cd45_meta_patient_not_in_raw_tumor = mismatch_meta_not_raw), file.path(out_dir, "mismatch_cd45_meta_patients_not_in_raw_tumor.csv"))

p_counts <- ggplot(patient_timepoint_counts, aes(x = patient_id, y = n_cells, fill = timepoint)) +
  geom_col(position = "dodge", width = 0.75) +
  facet_grid(. ~ path_response, scales = "free_x", space = "free_x") +
  labs(x = "Patient", y = "CD45 tumor cells", fill = "Timepoint") +
  theme_classic(base_size = 8) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(file.path(fig_dir, "cd45_tumor_patient_timepoint_cell_counts.pdf"), p_counts, width = 7.2, height = 3.2, useDingbats = FALSE)
ggsave(file.path(fig_dir, "cd45_tumor_patient_timepoint_cell_counts.png"), p_counts, width = 7.2, height = 3.2, dpi = 300)

p_celltype <- ggplot(celltype_counts, aes(x = reorder(celltype_id, n_cells), y = n_cells)) +
  geom_col(width = 0.7, fill = "#4C78A8") +
  coord_flip() +
  labs(x = "Cell type", y = "CD45 tumor cells") +
  theme_classic(base_size = 8)
ggsave(file.path(fig_dir, "cd45_tumor_celltype_counts.pdf"), p_celltype, width = 4.5, height = 2.8, useDingbats = FALSE)
ggsave(file.path(fig_dir, "cd45_tumor_celltype_counts.png"), p_celltype, width = 4.5, height = 2.8, dpi = 300)

summary_lines <- c(
  "# Data Audit Summary",
  "",
  paste0("- Audit time: ", Sys.time()),
  paste0("- Raw h5 files retained: ", nrow(h5_manifest)),
  paste0("- Tumor h5 files: ", nrow(tumor_h5)),
  paste0("- PBMC h5 files: ", nrow(h5_manifest[file_type == "PBMC"])),
  paste0("- CD45 tumor immune cells in metadata: ", nrow(dt)),
  paste0("- Patients in CD45 tumor metadata: ", uniqueN(dt$patient_id)),
  paste0("- Patients with pre tumor metadata: ", nrow(meta_presence[pre > 0])),
  paste0("- Patients with post tumor metadata: ", nrow(meta_presence[post > 0])),
  paste0("- Strict pre/post paired patients in CD45 tumor metadata: ", nrow(meta_presence[paired == TRUE])),
  paste0("- Paired response distribution: ", paste(pair_response_summary[, paste0(path_response, "=", N)], collapse = "; ")),
  paste0("- Raw tumor patients not in CD45 metadata: ", ifelse(length(mismatch_raw_not_meta), paste(mismatch_raw_not_meta, collapse = ", "), "none")),
  paste0("- CD45 metadata patients not in raw tumor h5 manifest: ", ifelse(length(mismatch_meta_not_raw), paste(mismatch_meta_not_raw, collapse = ", "), "none")),
  "",
  "## Interpretation Guardrails",
  "",
  "- Patient is the inferential unit.",
  "- Strict longitudinal dynamic analyses should use paired patients only.",
  "- Ordinal response should be primary; binary response should remain a sensitivity analysis.",
  "- No old manuscript statistics should be treated as final until recomputed."
)
writeLines(summary_lines, file.path(out_dir, "DATA_AUDIT_SUMMARY.md"))

message("Data audit complete.")
print(summary_lines)
