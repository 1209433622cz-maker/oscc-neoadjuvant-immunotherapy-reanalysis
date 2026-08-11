workspace <- Sys.getenv("GSE200996_WORKSPACE")
if (workspace == "") {
  workspace <- normalizePath(getwd(), mustWork = FALSE)
}

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(ggplot2)
})

raw_manifest_path <- file.path(workspace, "03_rebuild", "manifests", "GSE301741_RAW_ROUTE_SAMPLE_QC.csv")
strata_path <- file.path(
  workspace,
  "03_rebuild",
  "validation",
  "GSE301741_response_recovery",
  "GSE301741_VALIDATION_STRATA_WITH_RESPONSE.csv"
)
module_path <- file.path(workspace, "03_rebuild", "results", "external_validation", "GSE123813_gene_set_manifest.csv")
out_dir <- file.path(workspace, "03_rebuild", "validation", "GSE301741_tierA_module_validation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_csv_safe <- function(path) {
  if (!file.exists(path)) stop("Missing input: ", path)
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

manifest <- read_csv_safe(raw_manifest_path)
strata <- read_csv_safe(strata_path)
modules <- read_csv_safe(module_path)

tier_a <- strata[
  strata$validation_tier == "A_same_fraction_same_library" &
    strata$response_label %in% c("responder", "non_responder"),
]
if (nrow(tier_a) == 0) stop("No Tier A response-labelled patients available.")

scrna <- manifest[
  manifest$modality == "scRNA" &
    manifest$h5_status == "ok" &
    manifest$extracted_exists %in% c(TRUE, "TRUE", "True", "true"),
]

pair_rows <- list()
for (idx in seq_len(nrow(tier_a))) {
  patient <- tier_a$patient_id[idx]
  shared <- strsplit(tier_a$shared_scrna_fraction_library[idx], ";", fixed = TRUE)[[1]]
  shared <- shared[nzchar(shared)]
  for (combo in shared) {
    parts <- strsplit(combo, "/", fixed = TRUE)[[1]]
    if (length(parts) != 2) next
    fraction <- parts[[1]]
    library <- parts[[2]]
    subset <- scrna[
      scrna$patient_id == patient &
        scrna$sort_or_cell_fraction == fraction &
        scrna$assay_or_library == library &
        scrna$timepoint_harmonized %in% c("pre", "post"),
    ]
    if (length(unique(subset$timepoint_harmonized)) != 2) next
    for (tp in c("pre", "post")) {
      row <- subset[subset$timepoint_harmonized == tp, ][1, ]
      pair_rows[[length(pair_rows) + 1]] <- data.frame(
        pair_id = paste(patient, fraction, library, sep = "__"),
        patient_id = patient,
        response_label = tier_a$response_label[idx],
        response_status = tier_a$response_status[idx],
        response_source = tier_a$response_source[idx],
        provenance_status = tier_a$provenance_status[idx],
        timepoint = tp,
        sort_or_cell_fraction = fraction,
        assay_or_library = library,
        geo_accession = row$geo_accession,
        extracted_path = row$extracted_path,
        n_barcodes_manifest = row$n_barcodes,
        stringsAsFactors = FALSE
      )
    }
  }
}
pair_table <- do.call(rbind, pair_rows)
if (is.null(pair_table) || nrow(pair_table) == 0) stop("No matched Tier A h5 pairs found.")

valid_pairs <- names(which(table(pair_table$pair_id) == 2))
pair_table <- pair_table[pair_table$pair_id %in% valid_pairs, ]
pair_table <- pair_table[order(pair_table$patient_id, pair_table$pair_id, pair_table$timepoint), ]
write.csv(
  pair_table,
  file.path(out_dir, "GSE301741_TIERA_SELECTED_SAMPLE_PAIRS.csv"),
  row.names = FALSE
)

message("Selected h5 samples: ", nrow(pair_table), " across pair IDs: ", length(unique(pair_table$pair_id)))

profile_list <- list()
qc_rows <- list()
for (idx in seq_len(nrow(pair_table))) {
  sample_row <- pair_table[idx, ]
  h5_path <- sample_row$extracted_path
  message("Reading ", sample_row$geo_accession, ": ", h5_path)
  counts <- Read10X_h5(h5_path, use.names = TRUE)
  if (is.list(counts)) {
    counts <- counts[[1]]
  }
  counts <- as(counts, "dgCMatrix")
  gene_symbols <- rownames(counts)
  lib_size <- Matrix::colSums(counts)
  keep_cells <- lib_size > 0
  counts <- counts[, keep_cells, drop = FALSE]
  lib_size <- lib_size[keep_cells]
  norm <- Matrix::t(Matrix::t(counts) / lib_size * 10000)
  norm@x <- log1p(norm@x)
  gene_mean <- Matrix::rowMeans(norm)
  profile_list[[sample_row$geo_accession]] <- data.frame(
    gene = gene_symbols,
    value = as.numeric(gene_mean),
    stringsAsFactors = FALSE
  )
  qc_rows[[length(qc_rows) + 1]] <- data.frame(
    geo_accession = sample_row$geo_accession,
    patient_id = sample_row$patient_id,
    timepoint = sample_row$timepoint,
    pair_id = sample_row$pair_id,
    cells_loaded = ncol(counts),
    genes_loaded = nrow(counts),
    median_library_size = median(lib_size),
    stringsAsFactors = FALSE
  )
  rm(counts, norm)
  gc()
}

all_genes <- Reduce(union, lapply(profile_list, function(x) x$gene))
profile_mat <- matrix(NA_real_, nrow = length(all_genes), ncol = length(profile_list))
rownames(profile_mat) <- all_genes
colnames(profile_mat) <- names(profile_list)
for (sample_id in names(profile_list)) {
  prof <- profile_list[[sample_id]]
  profile_mat[prof$gene, sample_id] <- prof$value
}

qc <- do.call(rbind, qc_rows)
write.csv(qc, file.path(out_dir, "GSE301741_TIERA_SAMPLE_QC.csv"), row.names = FALSE)

z_mat <- t(scale(t(profile_mat)))
z_mat[is.nan(z_mat)] <- NA_real_

module_rows <- list()
score_rows <- list()
for (idx in seq_len(nrow(modules))) {
  genes <- trimws(strsplit(modules$genes_defined[idx], ";", fixed = TRUE)[[1]])
  genes <- genes[nzchar(genes)]
  present <- intersect(genes, rownames(z_mat))
  module_rows[[length(module_rows) + 1]] <- data.frame(
    signature = modules$signature[idx],
    target_lineage = modules$target_lineage[idx],
    source = modules$source[idx],
    n_genes_defined = length(genes),
    n_genes_present_in_GSE301741_TierA = length(present),
    genes_present = paste(present, collapse = ";"),
    genes_missing = paste(setdiff(genes, present), collapse = ";"),
    stringsAsFactors = FALSE
  )
  if (length(present) == 0) {
    values <- rep(NA_real_, ncol(z_mat))
  } else {
    values <- colMeans(z_mat[present, , drop = FALSE], na.rm = TRUE)
  }
  for (sample_id in names(values)) {
    meta <- pair_table[pair_table$geo_accession == sample_id, ][1, ]
    score_rows[[length(score_rows) + 1]] <- data.frame(
      geo_accession = sample_id,
      patient_id = meta$patient_id,
      pair_id = meta$pair_id,
      timepoint = meta$timepoint,
      response_label = meta$response_label,
      sort_or_cell_fraction = meta$sort_or_cell_fraction,
      assay_or_library = meta$assay_or_library,
      signature = modules$signature[idx],
      target_lineage = modules$target_lineage[idx],
      score = as.numeric(values[[sample_id]]),
      stringsAsFactors = FALSE
    )
  }
}

module_manifest <- do.call(rbind, module_rows)
score_table <- do.call(rbind, score_rows)
write.csv(module_manifest, file.path(out_dir, "GSE301741_TIERA_MODULE_MANIFEST.csv"), row.names = FALSE)
write.csv(score_table, file.path(out_dir, "GSE301741_TIERA_SAMPLE_MODULE_SCORES.csv"), row.names = FALSE)

delta_rows <- list()
for (sig in unique(score_table$signature)) {
  for (pair_id in unique(score_table$pair_id)) {
    sub <- score_table[score_table$signature == sig & score_table$pair_id == pair_id, ]
    if (!all(c("pre", "post") %in% sub$timepoint)) next
    pre <- sub$score[sub$timepoint == "pre"][1]
    post <- sub$score[sub$timepoint == "post"][1]
    meta <- sub[1, ]
    delta_rows[[length(delta_rows) + 1]] <- data.frame(
      patient_id = meta$patient_id,
      pair_id = pair_id,
      response_label = meta$response_label,
      sort_or_cell_fraction = meta$sort_or_cell_fraction,
      assay_or_library = meta$assay_or_library,
      signature = sig,
      target_lineage = meta$target_lineage,
      pre_score = pre,
      post_score = post,
      delta = post - pre,
      stringsAsFactors = FALSE
    )
  }
}
pair_delta <- do.call(rbind, delta_rows)
write.csv(pair_delta, file.path(out_dir, "GSE301741_TIERA_PAIR_MODULE_DELTAS.csv"), row.names = FALSE)

patient_delta <- aggregate(
  delta ~ patient_id + response_label + signature + target_lineage,
  data = pair_delta,
  FUN = mean
)
write.csv(patient_delta, file.path(out_dir, "GSE301741_TIERA_PATIENT_MODULE_DELTAS.csv"), row.names = FALSE)

test_rows <- list()
for (sig in unique(patient_delta$signature)) {
  sub <- patient_delta[patient_delta$signature == sig, ]
  r <- sub$delta[sub$response_label == "responder"]
  nr <- sub$delta[sub$response_label == "non_responder"]
  if (length(r) >= 2 && length(nr) >= 2) {
    t_p <- tryCatch(t.test(r, nr)$p.value, error = function(e) NA_real_)
    w_p <- tryCatch(wilcox.test(r, nr, exact = FALSE)$p.value, error = function(e) NA_real_)
  } else {
    t_p <- NA_real_
    w_p <- NA_real_
  }
  test_rows[[length(test_rows) + 1]] <- data.frame(
    signature = sig,
    target_lineage = sub$target_lineage[1],
    n_responder = length(r),
    n_non_responder = length(nr),
    mean_delta_responder = mean(r, na.rm = TRUE),
    mean_delta_non_responder = mean(nr, na.rm = TRUE),
    diff_responder_minus_non = mean(r, na.rm = TRUE) - mean(nr, na.rm = TRUE),
    t_p = t_p,
    wilcox_p = w_p,
    stringsAsFactors = FALSE
  )
}
tests <- do.call(rbind, test_rows)
tests$t_fdr <- p.adjust(tests$t_p, method = "BH")
tests$wilcox_fdr <- p.adjust(tests$wilcox_p, method = "BH")
tests <- tests[order(tests$t_p, na.last = TRUE), ]
write.csv(tests, file.path(out_dir, "GSE301741_TIERA_RESPONSE_MODULE_DELTA_TESTS.csv"), row.names = FALSE)

plot_modules <- head(tests$signature[!is.na(tests$t_p)], 12)
plot_df <- patient_delta[patient_delta$signature %in% plot_modules, ]
plot_df$signature <- factor(plot_df$signature, levels = rev(plot_modules))
p <- ggplot(plot_df, aes(x = response_label, y = delta, color = response_label)) +
  geom_hline(yintercept = 0, linewidth = 0.25, color = "grey75") +
  geom_point(position = position_jitter(width = 0.08, height = 0), size = 2.3, alpha = 0.9) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.45, color = "black", linewidth = 0.25) +
  facet_wrap(~signature, scales = "free_y", ncol = 3) +
  scale_color_manual(values = c(non_responder = "#255E0F", responder = "#FD65B3")) +
  labs(x = NULL, y = "Post-pre module-score delta", color = NULL) +
  theme_classic(base_size = 8) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 6.5),
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "bottom"
  )
ggsave(file.path(out_dir, "GSE301741_TIERA_MODULE_DELTA_TOP12.png"), p, width = 7.2, height = 5.2, dpi = 450)
ggsave(file.path(out_dir, "GSE301741_TIERA_MODULE_DELTA_TOP12.pdf"), p, width = 7.2, height = 5.2)

top <- head(tests, 8)
report <- c(
  "# GSE301741 Tier A Provisional Module Validation",
  "",
  paste0("Created: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Design",
  "",
  paste0("- Tier A patients: ", length(unique(patient_delta$patient_id))),
  paste0("- Matched fraction/library pair IDs: ", length(unique(pair_delta$pair_id))),
  "- Response labels: publication-supplement-figure-derived; not yet RDS/table-confirmed.",
  "- Scoring unit: sample-level gene profiles from h5 matrices; no cells are treated as independent response observations.",
  "",
  "## Top Module Delta Tests",
  "",
  "| Signature | Lineage | Mean delta R | Mean delta NR | Difference R-NR | P | FDR |",
  "|---|---|---:|---:|---:|---:|---:|"
)
for (i in seq_len(nrow(top))) {
  report <- c(
    report,
    sprintf(
      "| %s | %s | %.3f | %.3f | %.3f | %.3g | %.3g |",
      top$signature[i],
      top$target_lineage[i],
      top$mean_delta_responder[i],
      top$mean_delta_non_responder[i],
      top$diff_responder_minus_non[i],
      top$t_p[i],
      top$t_fdr[i]
    )
  )
}
report <- c(
  report,
  "",
  "## Interpretation Gate",
  "",
  "This is a provisional validation layer. It can guide manuscript refinement and figure planning, but final response-validation claims require RDS metadata or tabular supplement cross-check of patient-level labels.",
  "",
  "## Files Written",
  "",
  "- `GSE301741_TIERA_SELECTED_SAMPLE_PAIRS.csv`",
  "- `GSE301741_TIERA_SAMPLE_QC.csv`",
  "- `GSE301741_TIERA_MODULE_MANIFEST.csv`",
  "- `GSE301741_TIERA_SAMPLE_MODULE_SCORES.csv`",
  "- `GSE301741_TIERA_PAIR_MODULE_DELTAS.csv`",
  "- `GSE301741_TIERA_PATIENT_MODULE_DELTAS.csv`",
  "- `GSE301741_TIERA_RESPONSE_MODULE_DELTA_TESTS.csv`",
  "- `GSE301741_TIERA_MODULE_DELTA_TOP12.png/pdf`"
)
writeLines(report, file.path(out_dir, "GSE301741_TIERA_MODULE_VALIDATION_REPORT.md"))

message("Wrote report: ", file.path(out_dir, "GSE301741_TIERA_MODULE_VALIDATION_REPORT.md"))
