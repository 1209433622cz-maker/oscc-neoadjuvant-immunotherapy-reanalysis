#!/usr/bin/env Rscript

# Build Extended Data figures that support the main abundance/state-remodeling story.

get_script_path <- function() {
  cmd <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", cmd, value = TRUE)
  if (length(hit)) return(normalizePath(sub("^--file=", "", hit[1]), winslash = "/", mustWork = TRUE))
  normalizePath("03_rebuild/analysis/06_make_extended_data_figures.R", winslash = "/", mustWork = TRUE)
}

script_path <- get_script_path()
rebuild_dir <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
workspace_dir <- normalizePath(file.path(rebuild_dir, ".."), winslash = "/", mustWork = TRUE)
results_dir <- file.path(rebuild_dir, "results")
figure_dir <- file.path(rebuild_dir, "figures", "submission")
source_dir <- file.path(figure_dir, "source_data")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

required_packages <- c(
  "ggplot2", "dplyr", "tidyr", "readr", "stringr", "forcats",
  "patchwork", "scales", "grid", "tibble"
)
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(forcats)
  library(patchwork)
  library(scales)
})

theme_submission <- function(base_size = 7, base_family = "") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.3, colour = "black"),
      axis.text = element_text(colour = "black"),
      plot.title = element_text(face = "bold", size = rel(1.05), margin = margin(b = 4)),
      plot.subtitle = element_text(size = rel(0.9), colour = "grey25", margin = margin(b = 4)),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", colour = "black"),
      legend.title = element_text(face = "bold"),
      legend.key.size = grid::unit(0.35, "cm"),
      plot.tag = element_text(face = "bold", size = rel(1.25)),
      plot.margin = margin(5, 5, 5, 5)
    )
}
theme_set(theme_submission())

pal_response <- c(Low = "#0072B2", Medium = "#E69F00", High = "#CC79A7")
pal_timepoint <- c(pre = "#6BAED6", post = "#253494")
pal_celltype <- c(
  `T cell` = "#4C78A8",
  Myeloid = "#F58518",
  `B cell` = "#54A24B",
  `NK cell` = "#B279A2",
  Cycling = "#E45756",
  Mast = "#72B7B2"
)
pal_sig <- c(`Not significant` = "#8F969E", `FDR < 0.05` = "#B22222")

path <- function(...) file.path(...)

read_csv_required <- function(...) {
  f <- path(...)
  if (!file.exists(f)) stop("Required file not found: ", f, call. = FALSE)
  readr::read_csv(f, show_col_types = FALSE, na = c("", "NA", "NaN"))
}

fmt_p <- function(x) {
  ifelse(is.na(x), "NA", ifelse(x < 0.001, formatC(x, format = "e", digits = 2), sprintf("%.3f", x)))
}

clean_pathway <- function(x) {
  recode(
    x,
    HALLMARK_TNFA_SIGNALING_VIA_NFKB = "TNFA/NF-kB",
    HALLMARK_MTORC1_SIGNALING = "mTORC1",
    HALLMARK_P53_PATHWAY = "p53",
    HALLMARK_INTERFERON_ALPHA_RESPONSE = "IFN-alpha",
    HALLMARK_INTERFERON_GAMMA_RESPONSE = "IFN-gamma",
    HALLMARK_INFLAMMATORY_RESPONSE = "Inflammatory",
    HALLMARK_KRAS_SIGNALING_UP = "KRAS up",
    HALLMARK_MYC_TARGETS = "MYC targets",
    HALLMARK_COMPLEMENT = "Complement",
    HALLMARK_HYPOXIA = "Hypoxia",
    HALLMARK_UV_RESPONSE_UP = "UV response up",
    HALLMARK_ESTROGEN_RESPONSE_LATE = "Estrogen late",
    HALLMARK_ALLOGRAFT_REJECTION = "Allograft",
    HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION = "EMT",
    HALLMARK_OXIDATIVE_PHOSPHORYLATION = "OxPhos",
    HALLMARK_FATTY_ACID_METABOLISM = "Fatty acid",
    HALLMARK_IL6_JAK_STAT3_SIGNALING = "IL6/JAK/STAT3",
    HALLMARK_IL2_STAT5_SIGNALING = "IL2/STAT5",
    HALLMARK_DNA_REPAIR = "DNA repair",
    HALLMARK_APOPTOSIS = "Apoptosis",
    .default = str_to_sentence(str_replace_all(str_remove(x, "^HALLMARK_"), "_", " "))
  )
}

add_ci <- function(df) {
  df %>%
    mutate(
      se = ifelse(is.finite(t) & abs(t) > 1e-8, abs(logFC / t), NA_real_),
      ci_low = logFC - 1.96 * se,
      ci_high = logFC + 1.96 * se,
      sig = ifelse(adj.P.Val < 0.05, "FDR < 0.05", "Not significant")
    )
}

save_plot_set <- function(plot, stem, width, height, dpi = 360) {
  png_file <- file.path(figure_dir, paste0(stem, ".png"))
  pdf_file <- file.path(figure_dir, paste0(stem, ".pdf"))
  svg_file <- file.path(figure_dir, paste0(stem, ".svg"))
  if (requireNamespace("ragg", quietly = TRUE)) {
    ggsave(png_file, plot, width = width, height = height, dpi = dpi, bg = "white", device = ragg::agg_png)
  } else {
    ggsave(png_file, plot, width = width, height = height, dpi = dpi, bg = "white")
  }
  ggsave(pdf_file, plot, width = width, height = height, bg = "white", device = grDevices::pdf)
  if (requireNamespace("svglite", quietly = TRUE)) {
    ggsave(svg_file, plot, width = width, height = height, bg = "white", device = svglite::svglite)
  }
  invisible(c(png = png_file, pdf = pdf_file, svg = if (file.exists(svg_file)) svg_file else NA_character_))
}

write_source_data <- function(stem, sheets) {
  for (nm in names(sheets)) {
    safe_nm <- gsub("[^A-Za-z0-9_]+", "_", nm)
    readr::write_csv(sheets[[nm]], file.path(source_dir, paste0(stem, "_", safe_nm, ".csv")))
  }
  fixer <- file.path(rebuild_dir, "analysis", "11_rebuild_figure_source_workbooks.py")
  python_candidates <- unique(c(
    Sys.getenv("GSE200996_PYTHON", unset = ""),
    Sys.getenv("CODEX_PYTHON", unset = ""),
    "C:/Users/Administrator/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe",
    unname(Sys.which("python"))
  ))
  python_candidates <- python_candidates[nzchar(python_candidates)]
  python_candidates <- python_candidates[file.exists(python_candidates)]
  if (file.exists(fixer) && length(python_candidates)) {
    status <- system2(
      python_candidates[1],
      args = c(normalizePath(fixer, winslash = "/", mustWork = TRUE), "--stem", stem),
      stdout = TRUE,
      stderr = TRUE
    )
    exit_code <- attr(status, "status")
    if (!is.null(exit_code) && exit_code != 0) {
      warning("Failed to rebuild source-data XLSX for ", stem, ": ", paste(status, collapse = "\n"))
    }
  } else {
    warning("Skipped source-data XLSX rebuild for ", stem, "; per-panel CSV files were written.")
  }
}

top_gsea <- function(df, n = 8) {
  df %>%
    filter(!is.na(NES), !is.na(padj)) %>%
    arrange(padj, desc(abs(NES))) %>%
    slice_head(n = n) %>%
    mutate(
      pathway_label = clean_pathway(pathway),
      pathway_label = factor(pathway_label, levels = rev(pathway_label)),
      neg_log10_fdr = pmin(-log10(padj), 12)
    )
}

fit_signature_slopes <- function(df) {
  df %>%
    group_by(signature) %>%
    summarise(
      slope_per_response_step = coef(lm(delta ~ resp_num))[["resp_num"]],
      p_value = summary(lm(delta ~ resp_num))$coefficients["resp_num", "Pr(>|t|)"],
      mean_delta = mean(delta, na.rm = TRUE),
      .groups = "drop"
    )
}

message("Workspace: ", workspace_dir)
message("Results:   ", results_dir)
message("Figures:   ", figure_dir)

presence <- read_csv_required(results_dir, "data_audit", "cd45_tumor_patient_timepoint_presence.csv")
sample_counts <- read_csv_required(results_dir, "data_audit", "cd45_tumor_patient_timepoint_cell_counts.csv")
celltype_counts <- read_csv_required(results_dir, "data_audit", "cd45_tumor_celltype_counts.csv")
patient_cell_counts <- read_csv_required(results_dir, "data_audit", "cd45_tumor_patient_celltype_counts.csv")
paired_summary <- read_csv_required(results_dir, "data_audit", "paired_patient_response_summary.csv")

baseline_comp <- read_csv_required(results_dir, "pre_baseline", "baseline_pre_composition_trend_respOrd_limma_logit.csv") %>%
  add_ci()
t_base_de <- read_csv_required(results_dir, "pre_baseline", "Baseline_T_cell_DE_trend_respOrd.csv")
m_base_de <- read_csv_required(results_dir, "pre_baseline", "Baseline_Myeloid_DE_trend_respOrd.csv")
t_base_gsea <- read_csv_required(results_dir, "pre_baseline", "Baseline_T_cell_GSEA_Hallmark.csv")
m_base_gsea <- read_csv_required(results_dir, "pre_baseline", "Baseline_Myeloid_GSEA_Hallmark.csv")

t_sig <- read_csv_required(results_dir, "dynamic_paired", "Fig4B_T_cell_SignatureDelta_source.csv") %>%
  mutate(celltype = "T cell")
m_sig <- read_csv_required(results_dir, "dynamic_paired", "Fig4B_Myeloid_SignatureDelta_source.csv") %>%
  mutate(celltype = "Myeloid")
sig_df <- bind_rows(t_sig, m_sig) %>%
  mutate(response_ord = factor(response_ord, levels = c("Low", "Medium", "High")))
sig_slopes <- sig_df %>%
  group_by(celltype) %>%
  group_modify(~ fit_signature_slopes(.x)) %>%
  ungroup()

loo_stability <- read_csv_required(results_dir, "sensitivity_leave_one_out", "LOO_model_stability_summary.csv")
loo_pathways <- read_csv_required(results_dir, "sensitivity_leave_one_out", "LOO_key_pathway_NES.csv")
loo_direction <- read_csv_required(results_dir, "sensitivity_leave_one_out", "LOO_top25_gene_direction_concordance.csv")

## Extended Data Fig. 1: data audit

patient_order <- presence %>%
  arrange(response_ord_num, cohort, patient_id) %>%
  pull(patient_id)

sample_counts_plot <- sample_counts %>%
  mutate(
    patient_id = factor(patient_id, levels = patient_order),
    path_response = factor(path_response, levels = c("Low", "Medium", "High")),
    timepoint = factor(timepoint, levels = c("pre", "post"))
  )

ed1a <- ggplot(sample_counts_plot, aes(patient_id, n_cells, fill = timepoint)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.62) +
  facet_grid(path_response ~ ., scales = "free_x", space = "free_x") +
  scale_fill_manual(values = pal_timepoint, name = "Timepoint", labels = c(pre = "Pre", post = "Post")) +
  scale_y_continuous(labels = comma) +
  labs(title = "CD45+ tumor immune-cell counts", x = NULL, y = "Cells") +
  theme_submission() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")

ed1b <- celltype_counts %>%
  mutate(celltype_id = factor(celltype_id, levels = rev(celltype_id))) %>%
  ggplot(aes(fraction, celltype_id, fill = celltype_id)) +
  geom_col(width = 0.7) +
  scale_x_continuous(labels = percent_format(accuracy = 1)) +
  scale_fill_manual(values = pal_celltype, guide = "none") +
  labs(title = "Overall CD45+ immune composition", x = "Fraction of cells", y = NULL) +
  theme_submission()

patient_comp <- patient_cell_counts %>%
  group_by(patient_id, timepoint, path_response) %>%
  mutate(fraction = n_cells / sum(n_cells)) %>%
  ungroup() %>%
  mutate(
    patient_id = factor(patient_id, levels = patient_order),
    timepoint = factor(timepoint, levels = c("pre", "post"), labels = c("Pre", "Post")),
    path_response = factor(path_response, levels = c("Low", "Medium", "High")),
    celltype_id = factor(celltype_id, levels = names(pal_celltype))
  )

ed1c <- ggplot(patient_comp, aes(patient_id, fraction, fill = celltype_id)) +
  geom_col(width = 0.78) +
  facet_grid(timepoint ~ path_response, scales = "free_x", space = "free_x") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_fill_manual(values = pal_celltype, name = "Cell type") +
  labs(title = "Patient-level immune composition", x = NULL, y = "Fraction within sample") +
  theme_submission() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")

ed1d <- paired_summary %>%
  mutate(path_response = factor(path_response, levels = c("Low", "Medium", "High"))) %>%
  ggplot(aes(path_response, N, fill = path_response)) +
  geom_col(width = 0.62) +
  geom_text(aes(label = N), vjust = -0.4, size = 2.2) +
  scale_fill_manual(values = pal_response, guide = "none") +
  scale_y_continuous(breaks = pretty_breaks(n = 4), expand = expansion(mult = c(0, 0.12))) +
  labs(title = "Strict paired-patient response distribution", x = "Response depth", y = "Patients") +
  theme_submission()

ed1 <- (ed1a | ed1b) / (ed1c | ed1d) +
  plot_layout(heights = c(1, 1.15), widths = c(1.35, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 8))

save_plot_set(ed1, "ExtendedData1_data_audit", width = 7.2, height = 5.56, dpi = 300)
write_source_data(
  "ExtendedData1",
  list(
    sample_counts = sample_counts,
    celltype_counts = celltype_counts,
    patient_celltype_counts = patient_cell_counts,
    paired_response_summary = paired_summary
  )
)

## Extended Data Fig. 2: baseline context

cell_order <- c("T cell", "Myeloid", "Mast", "B cell", "NK cell", "Cycling")

ed2a <- baseline_comp %>%
  mutate(
    celltype = factor(celltype, levels = rev(cell_order)),
    sig = factor(sig, levels = c("Not significant", "FDR < 0.05"))
  ) %>%
  ggplot(aes(logFC, celltype)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.35, colour = "grey45") +
  geom_segment(aes(x = ci_low, xend = ci_high, y = celltype, yend = celltype), linewidth = 0.42, colour = "grey45", na.rm = TRUE) +
  geom_point(aes(fill = sig), shape = 21, size = 2.3, stroke = 0.35, colour = "black") +
  scale_fill_manual(values = pal_sig, name = NULL) +
  labs(title = "Baseline abundance trends", x = "logit abundance trend per response step", y = NULL) +
  theme_submission() +
  theme(legend.position = "none")

make_gsea_panel <- function(df, title) {
  plot_df <- top_gsea(df, n = 8)
  ggplot(plot_df, aes(NES, pathway_label)) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.35, colour = "grey50") +
    geom_point(aes(size = neg_log10_fdr, colour = NES > 0), alpha = 0.95) +
    scale_colour_manual(values = c(`FALSE` = "#2C7BB6", `TRUE` = "#B22222"), labels = c(`FALSE` = "Negative", `TRUE` = "Positive"), name = "NES sign") +
    scale_size_continuous(
      range = c(1.4, 4.6),
      limits = c(0, 12),
      breaks = c(3, 6, 9),
      name = "-log10 FDR"
    ) +
    labs(title = title, x = "NES", y = NULL) +
    theme_submission() +
    theme(legend.position = "right")
}

ed2b <- make_gsea_panel(t_base_gsea, "Baseline T-cell Hallmark context") +
  theme(legend.position = "none")
ed2c <- make_gsea_panel(m_base_gsea, "Baseline myeloid Hallmark context")

top_base_de <- bind_rows(t_base_de, m_base_de) %>%
  group_by(celltype) %>%
  arrange(padj, desc(abs(stat)), .by_group = TRUE) %>%
  slice_head(n = 8) %>%
  ungroup() %>%
  mutate(
    gene_label = paste(gene, celltype, sep = " | "),
    gene_label = fct_reorder(gene_label, stat),
    sig = ifelse(padj < 0.05, "FDR < 0.05", "Not significant")
  )

ed2d <- ggplot(top_base_de, aes(stat, gene_label)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.35, colour = "grey50") +
  geom_segment(aes(x = 0, xend = stat, y = gene_label, yend = gene_label), colour = "grey55", linewidth = 0.45) +
  geom_point(aes(fill = sig), shape = 21, size = 2.1, stroke = 0.35, colour = "black") +
  facet_wrap(~ celltype, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = pal_sig, name = NULL) +
  scale_y_discrete(labels = function(x) sub(" \\| .*", "", x)) +
  labs(title = "Top baseline response-trend genes", x = "DESeq2 Wald statistic", y = NULL) +
  theme_submission() +
  theme(legend.position = "bottom")

ed2 <- (ed2a | ed2d) / (ed2b | ed2c) +
  plot_layout(heights = c(1.05, 1), widths = c(1, 1.1)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 8))

save_plot_set(ed2, "ExtendedData2_baseline_context", width = 7.2, height = 5.67, dpi = 300)
write_source_data(
  "ExtendedData2",
  list(
    baseline_composition = baseline_comp,
    baseline_Tcell_DE_top8 = top_base_de %>% filter(celltype == "T cell"),
    baseline_Myeloid_DE_top8 = top_base_de %>% filter(celltype == "Myeloid"),
    baseline_Tcell_GSEA_top8 = top_gsea(t_base_gsea, n = 8),
    baseline_Myeloid_GSEA_top8 = top_gsea(m_base_gsea, n = 8)
  )
)

## Extended Data Fig. 3: signature deltas

make_sig_panel <- function(df, title) {
  ggplot(df, aes(resp_num, delta, colour = response_ord)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.35, colour = "grey55") +
    geom_smooth(aes(group = signature), method = "lm", formula = y ~ x, se = FALSE, colour = "grey35", linewidth = 0.45) +
    geom_point(size = 1.8, alpha = 0.95) +
    facet_wrap(~ signature, scales = "free_y", nrow = 2) +
    scale_colour_manual(values = pal_response, name = "Response") +
    scale_x_continuous(breaks = c(1, 2, 3), labels = c("Low", "Medium", "High")) +
    labs(title = title, x = "Response depth", y = "Signature delta (post-pre)") +
    theme_submission() +
    theme(legend.position = "bottom")
}

ed3a <- make_sig_panel(sig_df %>% filter(celltype == "T cell"), "T-cell signature deltas")
ed3b <- make_sig_panel(sig_df %>% filter(celltype == "Myeloid"), "Myeloid signature deltas")
ed3b <- ed3b + theme(legend.position = "none")

ed3c <- sig_slopes %>%
  mutate(
    signature = fct_reorder(signature, slope_per_response_step),
    sig = ifelse(p_value < 0.05, "P < 0.05", "Not significant")
  ) %>%
  ggplot(aes(slope_per_response_step, signature)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.35, colour = "grey55") +
  geom_point(aes(fill = sig), shape = 21, size = 2.2, stroke = 0.35, colour = "black") +
  facet_wrap(~ celltype, scales = "free") +
  scale_fill_manual(values = c(`Not significant` = "#8F969E", `P < 0.05` = "#B22222"), guide = "none") +
  labs(title = "Signature slope summaries", x = "Slope per response step", y = NULL) +
  theme_submission() +
  theme(legend.position = "bottom")

ed3 <- (ed3a / ed3b / ed3c) +
  plot_layout(heights = c(1, 1, 0.72)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 8))

save_plot_set(ed3, "ExtendedData3_signature_deltas", width = 7.2, height = 6.5, dpi = 300)
write_source_data(
  "ExtendedData3",
  list(
    signature_deltas = sig_df,
    signature_slope_summaries = sig_slopes
  )
)

## Extended Data Fig. 4: leave-one-patient sensitivity

loo_plot_df <- loo_stability %>%
  filter(leave_out != "NONE") %>%
  mutate(
    leave_out = factor(leave_out, levels = c("P18", "P23", "P24", "P27", "P29", "P32")),
    high_influence = ifelse(leave_out == "P32", "P32 removed", "Other patient removed")
  )

ed4a <- ggplot(loo_plot_df, aes(leave_out, stat_spearman, fill = high_influence)) +
  geom_col(width = 0.7) +
  facet_wrap(~ celltype, nrow = 1) +
  scale_fill_manual(values = c(`Other patient removed` = "#8F969E", `P32 removed` = "#B22222"), name = NULL) +
  scale_y_continuous(limits = c(0, 1), labels = number_format(accuracy = 0.1)) +
  labs(title = "Gene-rank similarity to full model", x = "Left-out patient", y = "Spearman correlation") +
  theme_submission() +
  theme(legend.position = "bottom")

ed4b <- ggplot(loo_plot_df, aes(leave_out, top100_overlap, fill = high_influence)) +
  geom_col(width = 0.7) +
  facet_wrap(~ celltype, nrow = 1) +
  scale_fill_manual(values = c(`Other patient removed` = "#8F969E", `P32 removed` = "#B22222"), guide = "none") +
  scale_y_continuous(limits = c(0, 1), labels = percent_format(accuracy = 1)) +
  labs(title = "Top-100 gene overlap", x = "Left-out patient", y = "Overlap with full model") +
  theme_submission()

key_pathways <- c(
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_MTORC1_SIGNALING",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_COMPLEMENT",
  "HALLMARK_P53_PATHWAY"
)

loo_pathway_plot <- loo_pathways %>%
  filter(pathway %in% key_pathways) %>%
  mutate(
    leave_out = factor(leave_out, levels = c("NONE", "P18", "P23", "P24", "P27", "P29", "P32")),
    pathway_label = factor(clean_pathway(pathway), levels = rev(unique(clean_pathway(key_pathways)))),
    fdr_sig = ifelse(padj < 0.05, "FDR < 0.05", "Not significant")
  )

ed4c <- ggplot(loo_pathway_plot, aes(leave_out, pathway_label, fill = NES)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  geom_point(aes(shape = fdr_sig), size = 1.2, colour = "black") +
  facet_wrap(~ celltype, nrow = 1, scales = "free_y") +
  scale_fill_gradient2(low = "#2C7BB6", mid = "#F7F7F7", high = "#B22222", midpoint = 0, name = "NES") +
  scale_shape_manual(values = c(`FDR < 0.05` = 16, `Not significant` = 1), name = NULL) +
  labs(title = "Key pathway NES after leave-one-patient refits", x = "Left-out patient", y = NULL) +
  theme_submission() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")

direction_summary <- loo_direction %>%
  filter(leave_out != "NONE") %>%
  group_by(celltype, leave_out) %>%
  summarise(direction_concordance = mean(same_direction, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    leave_out = factor(leave_out, levels = c("P18", "P23", "P24", "P27", "P29", "P32")),
    high_influence = ifelse(leave_out == "P32", "P32 removed", "Other patient removed")
  )

ed4d <- ggplot(direction_summary, aes(leave_out, direction_concordance, fill = high_influence)) +
  geom_col(width = 0.7) +
  facet_wrap(~ celltype, nrow = 1) +
  scale_fill_manual(values = c(`Other patient removed` = "#8F969E", `P32 removed` = "#B22222"), guide = "none") +
  scale_y_continuous(limits = c(0, 1), labels = percent_format(accuracy = 1)) +
  labs(title = "Full top-25 direction concordance", x = "Left-out patient", y = "Same direction") +
  theme_submission()

ed4 <- (ed4a | ed4b) / ed4c / ed4d +
  plot_layout(heights = c(0.9, 1.25, 0.85), widths = c(1, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 8))

save_plot_set(ed4, "ExtendedData4_leave_one_out_sensitivity", width = 13.2, height = 12.0)
write_source_data(
  "ExtendedData4",
  list(
    loo_model_stability = loo_stability,
    loo_key_pathway_NES = loo_pathways,
    loo_top25_direction_long = loo_direction,
    loo_top25_direction_summary = direction_summary
  )
)

manifest <- tibble::tibble(
  generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  figure = paste0("ExtendedData", 1:4),
  stem = c(
    "ExtendedData1_data_audit",
    "ExtendedData2_baseline_context",
    "ExtendedData3_signature_deltas",
    "ExtendedData4_leave_one_out_sensitivity"
  ),
  primary_message = c(
    "Audit of patient/timepoint coverage and CD45+ immune composition.",
    "Baseline-only abundance, gene and pathway context for the main negative abundance framing.",
    "Signature deltas are supportive and not statistically strong enough to carry the main claim.",
    "Leave-one-patient diagnostics support pathway-level claims but caution against fixed-gene signatures."
  )
)
readr::write_csv(manifest, file.path(figure_dir, "EXTENDED_DATA_FIGURE_MANIFEST.csv"))

message("Done. Extended Data figures and source data written to: ", figure_dir)
