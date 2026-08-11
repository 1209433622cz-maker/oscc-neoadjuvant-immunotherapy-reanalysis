workspace <- Sys.getenv("GSE200996_WORKSPACE")
if (workspace == "") {
  workspace <- normalizePath(getwd(), mustWork = FALSE)
}

suppressPackageStartupMessages(library(ggplot2))

rebuild_dir <- file.path(
  workspace, "03_rebuild", "validation", "GSE301741_raw_reconstruction"
)
old_validation_dir <- file.path(
  workspace, "03_rebuild", "validation", "GSE301741_tierA_module_validation"
)
module_path <- file.path(
  workspace, "03_rebuild", "results", "external_validation",
  "GSE123813_gene_set_manifest.csv"
)
out_dir <- file.path(
  workspace, "03_rebuild", "validation", "GSE301741_lineage_aware_validation"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_csv <- function(path) {
  if (!file.exists(path)) stop("Missing input: ", path)
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

pseudobulk <- read_csv(file.path(
  rebuild_dir, "GSE301741_RAW_TARGET_LINEAGE_PSEUDOBULK.csv.gz"
))
sample_qc <- read_csv(file.path(
  rebuild_dir, "GSE301741_RAW_REBUILT_SAMPLE_LINEAGE_QC.csv"
))
pairs <- read_csv(file.path(
  old_validation_dir, "GSE301741_TIERA_SELECTED_SAMPLE_PAIRS.csv"
))
modules <- read_csv(module_path)

fraction_allowed <- list(
  T_cell = c("CD3", "CD45pos", "CD45ratio", "Unsorted"),
  Myeloid = c("CD45pos", "CD45ratio", "Unsorted")
)

pb_cells <- unique(pseudobulk[c(
  "geo_accession", "target_lineage", "n_target_cells"
)])
pair_lineages <- list()
for (lineage in names(fraction_allowed)) {
  sub <- pairs[pairs$sort_or_cell_fraction %in% fraction_allowed[[lineage]], , drop = FALSE]
  sub$target_lineage <- lineage
  sub <- merge(
    sub, pb_cells,
    by = c("geo_accession", "target_lineage"),
    all.x = TRUE, sort = FALSE
  )
  pair_lineages[[lineage]] <- sub
}
eligible_samples <- do.call(rbind, pair_lineages)
eligible_samples$n_target_cells[is.na(eligible_samples$n_target_cells)] <- 0

pair_cell_qc <- aggregate(
  n_target_cells ~ pair_id + patient_id + response_label +
    target_lineage + sort_or_cell_fraction + assay_or_library,
  eligible_samples,
  function(x) min(x, na.rm = TRUE)
)
names(pair_cell_qc)[names(pair_cell_qc) == "n_target_cells"] <-
  "minimum_cells_across_pair"
pair_cell_qc$eligible_min30 <- pair_cell_qc$minimum_cells_across_pair >= 30
write.csv(
  pair_cell_qc,
  file.path(out_dir, "GSE301741_LINEAGE_AWARE_PAIR_CELL_QC.csv"),
  row.names = FALSE
)

score_rows <- list()
module_coverage <- list()
for (i in seq_len(nrow(modules))) {
  lineage <- modules$target_lineage[i]
  if (!lineage %in% names(fraction_allowed)) next
  genes <- trimws(strsplit(modules$genes_defined[i], ";", fixed = TRUE)[[1]])
  genes <- genes[nzchar(genes)]
  pb <- pseudobulk[
    pseudobulk$target_lineage == lineage &
      pseudobulk$gene %in% genes,
    ,
    drop = FALSE
  ]
  present <- unique(pb$gene[!is.na(pb$mean_log1p_cp10k)])
  module_coverage[[length(module_coverage) + 1L]] <- data.frame(
    signature = modules$signature[i],
    target_lineage = lineage,
    n_genes_defined = length(genes),
    n_genes_present = length(present),
    genes_present = paste(sort(present), collapse = ";"),
    genes_missing = paste(sort(setdiff(genes, present)), collapse = ";"),
    stringsAsFactors = FALSE
  )
  if (!length(present)) next

  # Gene-wise standardization is fitted across all reconstructed samples
  # within the target lineage, then frozen for the matched-pair comparison.
  pb$gene_z <- ave(
    pb$mean_log1p_cp10k, pb$gene,
    FUN = function(x) {
      s <- sd(x, na.rm = TRUE)
      if (!is.finite(s) || s == 0) return(rep(0, length(x)))
      (x - mean(x, na.rm = TRUE)) / s
    }
  )
  sample_score <- aggregate(
    gene_z ~ geo_accession + patient_id + timepoint + fraction +
      library + target_lineage + n_target_cells,
    pb,
    mean,
    na.rm = TRUE
  )
  sample_score$signature <- modules$signature[i]
  score_rows[[length(score_rows) + 1L]] <- sample_score
}
coverage <- do.call(rbind, module_coverage)
scores <- do.call(rbind, score_rows)
write.csv(
  coverage,
  file.path(out_dir, "GSE301741_LINEAGE_AWARE_MODULE_COVERAGE.csv"),
  row.names = FALSE
)
write.csv(
  scores,
  file.path(out_dir, "GSE301741_LINEAGE_AWARE_SAMPLE_MODULE_SCORES.csv"),
  row.names = FALSE
)

analysis_samples <- merge(
  eligible_samples,
  scores,
  by = c("geo_accession", "patient_id", "target_lineage"),
  all = FALSE,
  sort = FALSE,
  suffixes = c("_pair", "_score")
)
analysis_samples <- analysis_samples[analysis_samples$n_target_cells_score >= 30, , drop = FALSE]

delta_rows <- list()
keys <- unique(analysis_samples[c("pair_id", "signature", "target_lineage")])
for (i in seq_len(nrow(keys))) {
  sub <- analysis_samples[
    analysis_samples$pair_id == keys$pair_id[i] &
      analysis_samples$signature == keys$signature[i] &
      analysis_samples$target_lineage == keys$target_lineage[i],
    ,
    drop = FALSE
  ]
  if (!all(c("pre", "post") %in% sub$timepoint_pair)) next
  pre <- sub$gene_z[sub$timepoint_pair == "pre"][[1]]
  post <- sub$gene_z[sub$timepoint_pair == "post"][[1]]
  meta <- sub[1, , drop = FALSE]
  delta_rows[[length(delta_rows) + 1L]] <- data.frame(
    patient_id = meta$patient_id,
    pair_id = meta$pair_id,
    response_label = meta$response_label,
    fraction = meta$sort_or_cell_fraction,
    library = meta$assay_or_library,
    signature = meta$signature,
    target_lineage = meta$target_lineage,
    pre_score = pre,
    post_score = post,
    delta = post - pre,
    minimum_cells_across_pair = min(sub$n_target_cells_score),
    stringsAsFactors = FALSE
  )
}
if (!length(delta_rows)) stop("No lineage-aware matched pairs passed the 30-cell gate.")
pair_delta <- do.call(rbind, delta_rows)
write.csv(
  pair_delta,
  file.path(out_dir, "GSE301741_LINEAGE_AWARE_PAIR_MODULE_DELTAS.csv"),
  row.names = FALSE
)

patient_delta <- aggregate(
  delta ~ patient_id + response_label + signature + target_lineage,
  pair_delta,
  mean
)
write.csv(
  patient_delta,
  file.path(out_dir, "GSE301741_LINEAGE_AWARE_PATIENT_MODULE_DELTAS.csv"),
  row.names = FALSE
)

exact_permutation_p <- function(values, groups) {
  n_r <- sum(groups == "responder")
  n <- length(values)
  if (n_r < 2 || n - n_r < 2) return(NA_real_)
  observed <- mean(values[groups == "responder"]) -
    mean(values[groups == "non_responder"])
  combinations <- combn(seq_len(n), n_r)
  null <- apply(combinations, 2, function(ix) {
    mean(values[ix]) - mean(values[-ix])
  })
  mean(abs(null) >= abs(observed) - 1e-12)
}

test_rows <- list()
for (sig in unique(patient_delta$signature)) {
  sub <- patient_delta[patient_delta$signature == sig, , drop = FALSE]
  r <- sub$delta[sub$response_label == "responder"]
  nr <- sub$delta[sub$response_label == "non_responder"]
  inferential <- length(r) >= 2 && length(nr) >= 2
  test_rows[[length(test_rows) + 1L]] <- data.frame(
    signature = sig,
    target_lineage = sub$target_lineage[[1]],
    n_responder = length(r),
    n_non_responder = length(nr),
    mean_delta_responder = if (length(r)) mean(r) else NA_real_,
    mean_delta_non_responder = if (length(nr)) mean(nr) else NA_real_,
    diff_responder_minus_non = if (length(r) && length(nr)) mean(r) - mean(nr) else NA_real_,
    t_p = if (inferential) tryCatch(t.test(r, nr)$p.value, error = function(e) NA_real_) else NA_real_,
    wilcox_p = if (inferential) tryCatch(
      wilcox.test(r, nr, exact = FALSE)$p.value,
      error = function(e) NA_real_
    ) else NA_real_,
    exact_permutation_p = if (inferential) {
      exact_permutation_p(sub$delta, sub$response_label)
    } else {
      NA_real_
    },
    inference_status = if (inferential) {
      "eligible_small_n"
    } else {
      "descriptive_only_insufficient_response_groups"
    },
    stringsAsFactors = FALSE
  )
}
tests <- do.call(rbind, test_rows)
tests$permutation_fdr <- ave(
  tests$exact_permutation_p, tests$target_lineage,
  FUN = function(x) p.adjust(x, method = "BH")
)
tests <- tests[order(tests$target_lineage, tests$exact_permutation_p, na.last = TRUE), ]
write.csv(
  tests,
  file.path(out_dir, "GSE301741_LINEAGE_AWARE_RESPONSE_TESTS.csv"),
  row.names = FALSE
)

plot_sigs <- tests$signature[
  tests$target_lineage == "T_cell" &
    tests$inference_status == "eligible_small_n"
]
plot_df <- patient_delta[patient_delta$signature %in% plot_sigs, , drop = FALSE]
short_labels <- c(
  T_LE_INTERFERON_ALPHA_RESPONSE = "IFN-alpha",
  T_LE_INTERFERON_GAMMA_RESPONSE = "IFN-gamma",
  T_LE_union_core = "Union core",
  T_LE_TNFA_SIGNALING_VIA_NFKB = "TNF-NF-kB",
  T_DE_FDR05_positive = "DE FDR < 0.05",
  T_LE_MTORC1_SIGNALING = "mTORC1",
  T_LE_P53_PATHWAY = "p53"
)
plot_df$signature_short <- unname(short_labels[plot_df$signature])
plot_order <- unname(short_labels[plot_sigs])
plot_df$signature_short <- factor(plot_df$signature_short, levels = plot_order)
plot_df$response_label <- factor(
  plot_df$response_label,
  levels = c("non_responder", "responder"),
  labels = c("Non-responder", "Responder")
)
p <- ggplot(plot_df, aes(response_label, delta, color = response_label)) +
  geom_hline(yintercept = 0, linewidth = 0.25, color = "#9A9A9A") +
  geom_point(
    position = position_jitter(width = 0.06, height = 0, seed = 301741),
    size = 2.0, alpha = 0.9
  ) +
  stat_summary(
    fun = mean, geom = "crossbar", width = 0.42,
    color = "black", linewidth = 0.3
  ) +
  facet_wrap(~signature_short, scales = "free_y", ncol = 4) +
  scale_color_manual(values = c(
    "Non-responder" = "#3C5488",
    "Responder" = "#E64B35"
  )) +
  scale_x_discrete(labels = c(
    "Non-responder" = "NR",
    "Responder" = "R"
  )) +
  labs(x = NULL, y = "Delta module score (z)", color = NULL) +
  theme_classic(base_size = 8) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 7),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "none",
    panel.spacing = grid::unit(5, "pt")
  )
ggsave(
  file.path(out_dir, "GSE301741_LINEAGE_AWARE_TCELL_DELTAS.pdf"),
  p, width = 7.2, height = 5.2
)
ggsave(
  file.path(out_dir, "GSE301741_LINEAGE_AWARE_TCELL_DELTAS.png"),
  p, width = 7.2, height = 5.2, dpi = 600
)

top <- head(tests, 12)
report <- c(
  "# GSE301741 Lineage-Aware External Validation",
  "",
  paste0("Created: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Design",
  "",
  "- Unit of inference: patient, using post-treatment minus baseline change.",
  "- Only same-fraction, same-library Tier A pairs are retained.",
  "- T-cell modules exclude CD45-negative pairs.",
  "- Myeloid modules exclude CD3 and CD45-negative pairs.",
  "- Both timepoints must contain at least 30 reconstructed target-lineage cells.",
  "- Exact two-sided label permutation is the primary small-sample test.",
  "- Response labels remain supplement-figure-derived because the deposited RDS stream is unreadable.",
  "",
  "## Cohort Adequacy",
  ""
)
adequacy <- unique(tests[c(
  "target_lineage", "n_responder", "n_non_responder", "inference_status"
)])
for (i in seq_len(nrow(adequacy))) {
  report <- c(
    report,
    sprintf(
      "- %s: %d responders, %d non-responders; %s.",
      adequacy$target_lineage[i], adequacy$n_responder[i],
      adequacy$n_non_responder[i], adequacy$inference_status[i]
    )
  )
}
report <- c(
  report,
  "",
  "## Module Results",
  "",
  "| Signature | Lineage | n R | n NR | Difference R-NR | Exact P | FDR | Status |",
  "|---|---|---:|---:|---:|---:|---:|---|"
)
for (i in seq_len(nrow(top))) {
  report <- c(
    report,
    sprintf(
      "| %s | %s | %d | %d | %s | %s | %s | %s |",
      top$signature[i], top$target_lineage[i],
      top$n_responder[i], top$n_non_responder[i],
      ifelse(is.na(top$diff_responder_minus_non[i]), "NA",
             sprintf("%.3f", top$diff_responder_minus_non[i])),
      ifelse(is.na(top$exact_permutation_p[i]), "NA",
             sprintf("%.3g", top$exact_permutation_p[i])),
      ifelse(is.na(top$permutation_fdr[i]), "NA",
             sprintf("%.3g", top$permutation_fdr[i])),
      top$inference_status[i]
    )
  )
}
report <- c(
  report,
  "",
  "## Interpretation Gate",
  "",
  "This analysis is an independent, lineage-aware sensitivity layer reconstructed from deposited raw H5 matrices. Small patient counts, figure-derived response labels, and non-identical author cell labels preclude confirmatory wording. Myeloid results remain descriptive whenever either response group has fewer than two patients.",
  "",
  "## Superseded Output",
  "",
  "The earlier whole-sample Tier A module analysis is retained for provenance but must not support manuscript claims because it mixed target and non-target lineages."
)
writeLines(
  report,
  file.path(out_dir, "GSE301741_LINEAGE_AWARE_VALIDATION_REPORT.md")
)
writeLines(
  c(
    paste0("completed=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
    paste0("pair_deltas=", nrow(pair_delta)),
    paste0("patients=", length(unique(patient_delta$patient_id)))
  ),
  file.path(out_dir, "GSE301741_LINEAGE_AWARE_VALIDATION_SUCCESS.txt")
)
cat(paste(report, collapse = "\n"), "\n")
