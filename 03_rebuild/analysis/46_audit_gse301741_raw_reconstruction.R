workspace <- Sys.getenv("GSE200996_WORKSPACE")
if (workspace == "") {
  workspace <- normalizePath(getwd(), mustWork = FALSE)
}

manifest_path <- file.path(
  workspace, "03_rebuild", "manifests", "GSE301741_RAW_ROUTE_SAMPLE_QC.csv"
)
rebuild_dir <- file.path(
  workspace, "03_rebuild", "validation", "GSE301741_raw_reconstruction"
)
validation_dir <- file.path(
  workspace, "03_rebuild", "validation", "GSE301741_lineage_aware_validation"
)
out_dir <- file.path(workspace, "03_rebuild", "audit")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_csv <- function(path) {
  if (!file.exists(path)) stop("Missing input: ", path)
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

manifest <- read_csv(manifest_path)
manifest <- manifest[
  manifest$modality == "scRNA" &
    manifest$h5_status == "ok" &
    manifest$extracted_exists %in% c(TRUE, "TRUE", "True", "true"),
  ,
  drop = FALSE
]
metadata <- read_csv(file.path(
  rebuild_dir, "GSE301741_RAW_REBUILT_CELL_METADATA.csv.gz"
))
qc <- read_csv(file.path(
  rebuild_dir, "GSE301741_RAW_REBUILT_SAMPLE_LINEAGE_QC.csv"
))
pseudobulk <- read_csv(file.path(
  rebuild_dir, "GSE301741_RAW_TARGET_LINEAGE_PSEUDOBULK.csv.gz"
))
pair_qc <- read_csv(file.path(
  validation_dir, "GSE301741_LINEAGE_AWARE_PAIR_CELL_QC.csv"
))
tests <- read_csv(file.path(
  validation_dir, "GSE301741_LINEAGE_AWARE_RESPONSE_TESTS.csv"
))

checks <- list()
add_check <- function(id, passed, observed, expected) {
  checks[[length(checks) + 1L]] <<- data.frame(
    check_id = id,
    status = if (isTRUE(passed)) "PASS" else "FAIL",
    observed = as.character(observed),
    expected = as.character(expected),
    stringsAsFactors = FALSE
  )
}

expected_cells <- sum(as.numeric(manifest$n_barcodes))
add_check(
  "metadata_row_count",
  nrow(metadata) == expected_cells,
  nrow(metadata),
  expected_cells
)
add_check(
  "unique_cell_id",
  !anyDuplicated(metadata$cell_id),
  length(unique(metadata$cell_id)),
  nrow(metadata)
)
add_check(
  "sample_accession_count",
  length(unique(metadata$geo_accession)) == nrow(manifest),
  length(unique(metadata$geo_accession)),
  nrow(manifest)
)
qc_total <- aggregate(n_cells ~ geo_accession, qc, sum)
manifest_counts <- manifest[c("geo_accession", "n_barcodes")]
manifest_counts$n_barcodes <- as.numeric(manifest_counts$n_barcodes)
reconcile <- merge(manifest_counts, qc_total, by = "geo_accession", all = TRUE)
add_check(
  "sample_count_reconciliation",
  all(reconcile$n_barcodes == reconcile$n_cells),
  sum(reconcile$n_barcodes == reconcile$n_cells),
  nrow(reconcile)
)

retained <- metadata$final_lineage != "Filtered_or_unresolved"
lymphoid <- metadata$final_lineage %in% c("T_cell", "NK", "B_cell", "Plasma")
add_check(
  "retained_mitochondrial_gate",
  all(metadata$percent_mt[retained] <= 20 + 1e-8),
  max(metadata$percent_mt[retained]),
  "<=20"
)
add_check(
  "retained_lymphoid_feature_gate",
  all(metadata$n_feature[retained & lymphoid] >= 500),
  min(metadata$n_feature[retained & lymphoid]),
  ">=500"
)
add_check(
  "retained_nonlymphoid_feature_gate",
  all(metadata$n_feature[retained & !lymphoid] >= 1000),
  min(metadata$n_feature[retained & !lymphoid]),
  ">=1000"
)
add_check(
  "target_lineage_mapping",
  all(
    metadata$target_lineage[metadata$final_lineage == "T_cell"] == "T_cell"
  ) && all(
    metadata$target_lineage[
      metadata$final_lineage %in% c("Myeloid", "DC")
    ] == "Myeloid"
  ),
  "consistent",
  "T_cell and Myeloid/DC only"
)

cd3_bad <- metadata$fraction == "CD3" &
  !metadata$final_lineage %in% c("T_cell", "NK", "Filtered_or_unresolved")
cd45neg_bad <- metadata$fraction == "CD45neg" &
  !metadata$final_lineage %in% c(
    "Epithelial", "Fibroblast", "Endothelial", "Filtered_or_unresolved"
  )
add_check("cd3_fraction_constraint", !any(cd3_bad), sum(cd3_bad), 0)
add_check("cd45neg_fraction_constraint", !any(cd45neg_bad), sum(cd45neg_bad), 0)

genes_per_sample_lineage <- aggregate(
  gene ~ geo_accession + target_lineage,
  pseudobulk,
  function(x) length(unique(x))
)
add_check(
  "pseudobulk_complete_grid",
  nrow(genes_per_sample_lineage) == nrow(manifest) * 2 &&
    length(unique(genes_per_sample_lineage$gene)) == 1,
  paste(nrow(genes_per_sample_lineage),
        paste(unique(genes_per_sample_lineage$gene), collapse = ",")),
  paste(nrow(manifest) * 2, "constant gene count")
)

t_pairs <- pair_qc[
  pair_qc$target_lineage == "T_cell" &
    pair_qc$eligible_min30 %in% c(TRUE, "TRUE", "True", "true"),
  ,
  drop = FALSE
]
m_pairs <- pair_qc[
  pair_qc$target_lineage == "Myeloid" &
    pair_qc$eligible_min30 %in% c(TRUE, "TRUE", "True", "true"),
  ,
  drop = FALSE
]
add_check(
  "tcell_pair_adequacy",
  nrow(t_pairs) == 7 &&
    sum(t_pairs$response_label == "responder") == 3 &&
    sum(t_pairs$response_label == "non_responder") == 4,
  paste(nrow(t_pairs),
        sum(t_pairs$response_label == "responder"),
        sum(t_pairs$response_label == "non_responder"), sep = "/"),
  "7/3/4"
)
add_check(
  "myeloid_inference_blocked",
  nrow(m_pairs) == 2 &&
    all(m_pairs$response_label == "non_responder") &&
    all(tests$inference_status[tests$target_lineage == "Myeloid"] ==
          "descriptive_only_insufficient_response_groups"),
  paste(nrow(m_pairs), paste(unique(m_pairs$response_label), collapse = ",")),
  "2 non-responders; descriptive only"
)
add_check(
  "tcell_exact_p_range",
  all(
    is.finite(tests$exact_permutation_p[tests$target_lineage == "T_cell"]) &
      tests$exact_permutation_p[tests$target_lineage == "T_cell"] >= 0 &
      tests$exact_permutation_p[tests$target_lineage == "T_cell"] <= 1
  ),
  paste(range(tests$exact_permutation_p[tests$target_lineage == "T_cell"]),
        collapse = "-"),
  "finite within 0-1"
)

checks <- do.call(rbind, checks)
write.csv(
  checks,
  file.path(out_dir, "GSE301741_RAW_RECONSTRUCTION_AUDIT.csv"),
  row.names = FALSE
)

report <- c(
  "# GSE301741 RAW Reconstruction Audit",
  "",
  paste0("Created: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  paste0("- Checks: ", nrow(checks)),
  paste0("- PASS: ", sum(checks$status == "PASS")),
  paste0("- FAIL: ", sum(checks$status == "FAIL")),
  "",
  "| Check | Status | Observed | Expected |",
  "|---|---|---|---|"
)
for (i in seq_len(nrow(checks))) {
  report <- c(
    report,
    sprintf(
      "| %s | %s | %s | %s |",
      checks$check_id[i], checks$status[i],
      checks$observed[i], checks$expected[i]
    )
  )
}
writeLines(
  report,
  file.path(out_dir, "GSE301741_RAW_RECONSTRUCTION_AUDIT.md")
)
if (any(checks$status == "FAIL")) {
  stop("GSE301741 RAW reconstruction audit failed.")
}
cat(paste(report, collapse = "\n"), "\n")
