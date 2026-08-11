#!/usr/bin/env Rscript

# Build polished submission-oriented figures from the recomputed result tables.

get_script_path <- function() {
  cmd <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", cmd, value = TRUE)
  if (length(hit)) return(normalizePath(sub("^--file=", "", hit[1]), winslash = "/", mustWork = TRUE))
  normalizePath("03_rebuild/analysis/05_make_core_figures.R", winslash = "/", mustWork = TRUE)
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
  "patchwork", "scales", "grid"
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
pal_binary <- c(`Not significant` = "#8F969E", `FDR < 0.05` = "#B22222")
pal_exact <- c(`Not significant` = "#8F969E", `Exact FDR < 0.05` = "#B22222")
pal_present <- c(`No sample` = "#F2F2F2", `Sample available` = "#2F4858")
pal_cohort <- c(Mono = "#6C6CB3", Combo = "#C75D4A")
heat_cols <- c("#2C7BB6", "#F7F7F7", "#D7191C")

path <- function(...) file.path(...)

read_csv_required <- function(...) {
  f <- path(...)
  if (!file.exists(f)) stop("Required file not found: ", f, call. = FALSE)
  readr::read_csv(f, show_col_types = FALSE, na = c("", "NA", "NaN"))
}

read_matrix_required <- function(...) {
  f <- path(...)
  if (!file.exists(f)) stop("Required file not found: ", f, call. = FALSE)
  x <- read.csv(f, check.names = FALSE, row.names = 1)
  as.matrix(x)
}

read_annotation_required <- function(...) {
  f <- path(...)
  if (!file.exists(f)) stop("Required file not found: ", f, call. = FALSE)
  x <- read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
  names(x)[1] <- "patient"
  tibble::as_tibble(x)
}

fmt_p <- function(x) {
  ifelse(is.na(x), "NA", ifelse(x < 0.001, formatC(x, format = "e", digits = 2), sprintf("%.3f", x)))
}

clean_pathway <- function(x) {
  recode(
    x,
    HALLMARK_TNFA_SIGNALING_VIA_NFKB = "TNFA/NF-kB signaling",
    HALLMARK_MTORC1_SIGNALING = "mTORC1 signaling",
    HALLMARK_P53_PATHWAY = "p53 pathway",
    HALLMARK_INTERFERON_ALPHA_RESPONSE = "IFN-alpha response",
    HALLMARK_INTERFERON_GAMMA_RESPONSE = "IFN-gamma response",
    HALLMARK_INFLAMMATORY_RESPONSE = "Inflammatory response",
    HALLMARK_KRAS_SIGNALING_UP = "KRAS signaling up",
    HALLMARK_MYC_TARGETS = "MYC targets",
    HALLMARK_COMPLEMENT = "Complement",
    HALLMARK_HYPOXIA = "Hypoxia",
    HALLMARK_UV_RESPONSE_UP = "UV response up",
    HALLMARK_ESTROGEN_RESPONSE_LATE = "Estrogen response late",
    HALLMARK_ALLOGRAFT_REJECTION = "Allograft rejection",
    HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION = "EMT",
    HALLMARK_OXIDATIVE_PHOSPHORYLATION = "Oxidative phosphorylation",
    HALLMARK_FATTY_ACID_METABOLISM = "Fatty acid metabolism",
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

row_z <- function(mat) {
  z <- t(scale(t(mat)))
  z[!is.finite(z)] <- 0
  pmax(pmin(z, 2.5), -2.5)
}

matrix_to_long <- function(mat, ann, top_n = 18) {
  keep <- seq_len(min(nrow(mat), top_n))
  mat2 <- mat[keep, , drop = FALSE]
  z <- row_z(mat2)
  df <- as.data.frame(z, check.names = FALSE) %>%
    tibble::rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to = "patient", values_to = "z_delta") %>%
    left_join(ann %>% select(patient, response_ord, response_bin, cohort), by = "patient")

  patient_order <- ann %>%
    mutate(
      response_ord = factor(response_ord, levels = c("Low", "Medium", "High")),
      patient = factor(patient, levels = patient)
    ) %>%
    arrange(response_ord, patient) %>%
    pull(patient) %>%
    as.character()

  df %>%
    mutate(
      patient = factor(patient, levels = patient_order),
      gene = factor(gene, levels = rev(rownames(mat2))),
      response_ord = factor(response_ord, levels = c("Low", "Medium", "High"))
    )
}

save_plot_set <- function(plot, stem, width, height, dpi = 360) {
  png_file <- file.path(figure_dir, paste0(stem, ".png"))
  pdf_file <- file.path(figure_dir, paste0(stem, ".pdf"))
  svg_file <- file.path(figure_dir, paste0(stem, ".svg"))

  if (requireNamespace("ragg", quietly = TRUE)) {
    ggplot2::ggsave(png_file, plot, width = width, height = height, dpi = dpi, bg = "white", device = ragg::agg_png)
  } else {
    ggplot2::ggsave(png_file, plot, width = width, height = height, dpi = dpi, bg = "white")
  }
  ggplot2::ggsave(pdf_file, plot, width = width, height = height, bg = "white", device = grDevices::pdf)
  if (requireNamespace("svglite", quietly = TRUE)) {
    ggplot2::ggsave(svg_file, plot, width = width, height = height, bg = "white", device = svglite::svglite)
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

message("Workspace: ", workspace_dir)
message("Results:   ", results_dir)
message("Figures:   ", figure_dir)

presence <- read_csv_required(results_dir, "data_audit", "cd45_tumor_patient_timepoint_presence.csv")
response_table <- read_csv_required(results_dir, "data_audit", "cd45_tumor_patient_response_table.csv")
paired_summary <- read_csv_required(results_dir, "data_audit", "paired_patient_response_summary.csv")
baseline_comp <- read_csv_required(results_dir, "pre_baseline", "baseline_pre_composition_trend_respOrd_limma_logit.csv") %>%
  add_ci()
dyn_comp_ord <- read_csv_required(results_dir, "dynamic_paired", "Fig4A_composition_delta_logit_limma_respOrd_trend.csv") %>%
  add_ci() %>%
  mutate(model = "Ordinal response trend")
dyn_comp_bin <- read_csv_required(results_dir, "dynamic_paired", "Fig4A_composition_delta_logit_limma_RvsNR.csv") %>%
  add_ci() %>%
  mutate(model = "High versus Low sensitivity")
exact_abundance <- read_csv_required(
  results_dir,
  "sensitivity_exact_permutation",
  "ABUNDANCE_EXACT_PERMUTATION_RESULTS.csv"
)
exact_ordinal <- exact_abundance %>%
  filter(analysis == "paired_delta_ordinal_unadjusted") %>%
  select(celltype, exact_p_two_sided, exact_bh_fdr)
exact_binary <- exact_abundance %>%
  filter(analysis == "paired_delta_High_vs_Low") %>%
  select(celltype, exact_p_two_sided, exact_bh_fdr)
dyn_comp_ord <- dyn_comp_ord %>% left_join(exact_ordinal, by = "celltype")
dyn_comp_bin <- dyn_comp_bin %>% left_join(exact_binary, by = "celltype")
paired_prop <- read_csv_required(results_dir, "dynamic_paired", "Fig4A_patient_timepoint_celltype_prop_long_allpaired.csv")

t_gsea <- read_csv_required(results_dir, "dynamic_paired", "Fig4B_T_cell_GSEA_Hallmark.csv")
m_gsea <- read_csv_required(results_dir, "dynamic_paired", "Fig4B_Myeloid_GSEA_Hallmark.csv")
t_mat <- read_matrix_required(results_dir, "dynamic_paired", "Fig4B_Tcell_delta_matrix.csv")
m_mat <- read_matrix_required(results_dir, "dynamic_paired", "Fig4B_Myeloid_delta_matrix.csv")
t_ann <- read_annotation_required(results_dir, "dynamic_paired", "Fig4B_Tcell_sample_annotation.csv")
m_ann <- read_annotation_required(results_dir, "dynamic_paired", "Fig4B_Myeloid_sample_annotation.csv")

study_counts <- tibble::tibble(
  item = c(
    "Tumor scRNA-seq libraries",
    "CD45+ immune cells",
    "Patients with tumor immune metadata",
    "Strict pre/post paired patients",
    "Primary model"
  ),
  value = c(
    "25 tumor h5",
    "74,557 cells",
    "19 patients",
    paste0(sum(presence$paired, na.rm = TRUE), " patients: ",
           paste0(paired_summary$path_response, "=", paired_summary$N, collapse = ", ")),
    "patient-level abundance and pseudobulk interaction"
  )
)

panel1a_nodes <- tibble::tibble(
  x = 1:5,
  y = 0,
  step = c("Dataset", "Tumor libraries", "Immune cells", "Patient boundary", "Inference"),
  value = c("GSE200996\nOSCC tumors", "25 tumor h5", "74,557 CD45+\nimmune cells", "19 patients\n6 strict pairs", "Abundance +\nstate remodeling"),
  group = c("input", "input", "boundary", "boundary", "model")
)

p1a <- ggplot(panel1a_nodes, aes(x, y)) +
  geom_segment(
    data = tibble::tibble(x = 1:4, xend = 2:5, y = 1, yend = 1),
    aes(x = x + 0.12, xend = xend - 0.12, y = 0, yend = 0),
    inherit.aes = FALSE,
    linewidth = 0.42,
    colour = "grey35",
    arrow = grid::arrow(length = grid::unit(0.13, "cm"), type = "closed")
  ) +
  geom_point(aes(fill = group), shape = 21, size = 4.8, stroke = 0.4, colour = "black") +
  geom_text(aes(label = step), y = 0.22, fontface = "bold", size = 2.1, vjust = 0) +
  geom_text(aes(label = value), y = -0.22, size = 1.9, lineheight = 0.9, vjust = 1) +
  scale_fill_manual(values = c(input = "#D7E7F3", boundary = "#DDEED7", model = "#F4DEC8"), guide = "none") +
  coord_cartesian(xlim = c(0.55, 5.58), ylim = c(-0.50, 0.48), clip = "off") +
  labs(title = "Study design") +
  theme_void(base_size = 7) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0, size = 7),
    plot.margin = margin(4, 5, 1, 5)
  )

presence_long <- presence %>%
  mutate(
    path_response = factor(path_response, levels = c("Low", "Medium", "High")),
    patient_id = factor(patient_id, levels = presence %>% arrange(response_ord_num, cohort, patient_id) %>% pull(patient_id)),
    cohort = factor(cohort, levels = c("Mono", "Combo"))
  ) %>%
  pivot_longer(c(pre, post), names_to = "timepoint", values_to = "available") %>%
  mutate(
    timepoint = factor(timepoint, levels = c("pre", "post"), labels = c("Pre", "Post")),
    available_label = ifelse(available == 1, "Sample available", "No sample")
  )

p1b <- ggplot(presence_long, aes(timepoint, patient_id)) +
  geom_tile(aes(fill = available_label), colour = "white", linewidth = 0.4, width = 0.92, height = 0.82) +
  geom_point(
    data = presence_long %>% distinct(patient_id, path_response, cohort, paired) %>% mutate(x = 2.72),
    aes(x = x, y = patient_id, colour = cohort, shape = paired),
    inherit.aes = FALSE,
    size = 1.6,
    stroke = 0.4
  ) +
  facet_grid(path_response ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_manual(values = pal_present, name = NULL) +
  scale_colour_manual(values = pal_cohort, name = "Cohort") +
  scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16), labels = c(`FALSE` = "No", `TRUE` = "Yes"), name = "Paired") +
  coord_cartesian(xlim = c(0.55, 3.0), clip = "off") +
  labs(title = "Patient-level sample availability", x = NULL, y = NULL) +
  theme_submission() +
  theme(
    legend.position = "bottom",
    panel.spacing.y = grid::unit(0.1, "lines"),
    axis.ticks.y = element_blank(),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0)
  )

cell_order <- c("T cell", "Myeloid", "Mast", "B cell", "NK cell", "Cycling")

baseline_plot_df <- baseline_comp %>%
  mutate(
    celltype = factor(celltype, levels = rev(cell_order)),
    fdr_label = paste0("FDR ", fmt_p(adj.P.Val))
  )

p1c <- ggplot(baseline_plot_df, aes(logFC, celltype)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.35, colour = "grey45") +
  geom_segment(aes(x = ci_low, xend = ci_high, y = celltype, yend = celltype), linewidth = 0.45, colour = "grey45", na.rm = TRUE) +
  geom_point(aes(fill = sig), shape = 21, size = 2.2, stroke = 0.35, colour = "black") +
  scale_fill_manual(values = pal_binary, guide = "none") +
  labs(
    title = "Baseline abundance trends",
    x = "logit abundance trend per response step",
    y = NULL
  ) +
  theme_submission()

paired_focus <- paired_prop %>%
  filter(celltype %in% c("T cell", "Myeloid", "Mast")) %>%
  mutate(
    response_ord = factor(response_ord, levels = c("Low", "Medium", "High")),
    celltype = factor(celltype, levels = c("T cell", "Myeloid", "Mast")),
    timepoint = factor(timepoint, levels = c("pre", "post"), labels = c("Pre", "Post"))
  )

p1d <- ggplot(paired_focus, aes(timepoint, prop, group = patient, colour = response_ord, shape = response_ord, linetype = response_ord)) +
  geom_line(linewidth = 0.5, alpha = 0.85) +
  geom_point(size = 1.8, alpha = 0.95) +
  facet_wrap(~ celltype, scales = "free_y", nrow = 1) +
  scale_colour_manual(values = pal_response, name = "Response") +
  scale_shape_manual(values = c(Low = 16, Medium = 17, High = 15), name = "Response") +
  scale_linetype_manual(values = c(Low = "solid", Medium = "dashed", High = "dotted"), name = "Response") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    title = "Paired immune-abundance trajectories",
    x = NULL,
    y = "Immune-cell fraction"
  ) +
  theme_submission() +
  theme(legend.position = "bottom")

delta_model_df <- bind_rows(dyn_comp_ord, dyn_comp_bin) %>%
  mutate(
    model = factor(model, levels = c("Ordinal response trend", "High versus Low sensitivity")),
    celltype = factor(celltype, levels = rev(cell_order)),
    exact_sig = factor(
      ifelse(exact_bh_fdr < 0.05, "Exact FDR < 0.05", "Not significant"),
      levels = c("Not significant", "Exact FDR < 0.05")
    )
  )

p1e <- ggplot(delta_model_df, aes(logFC, celltype)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.35, colour = "grey45") +
  geom_segment(aes(x = ci_low, xend = ci_high, y = celltype, yend = celltype), linewidth = 0.42, colour = "grey45", na.rm = TRUE) +
  geom_point(aes(fill = exact_sig), shape = 21, size = 2.25, stroke = 0.35, colour = "black") +
  facet_wrap(~ model, nrow = 1, scales = "free_x") +
  scale_fill_manual(values = pal_exact, guide = "none", drop = TRUE) +
  labs(
    title = "Exact-permutation abundance effects",
    subtitle = NULL,
    x = "post-pre logit abundance effect",
    y = NULL
  ) +
  theme_submission() +
  theme(legend.position = "bottom")

fig1 <- p1a / ((p1b | p1c) + plot_layout(widths = c(1.1, 1))) / ((p1d | p1e) + plot_layout(widths = c(1.2, 1))) +
  plot_layout(heights = c(0.42, 1.38, 1.25)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 8))

save_plot_set(fig1, "Figure1_submission_abundance_and_design", width = 7.2, height = 5.56, dpi = 300)

write_source_data(
  "Figure1",
  list(
    study_design_counts = study_counts,
    patient_availability = presence_long,
    baseline_composition = baseline_comp,
    paired_celltype_proportions = paired_prop,
    paired_delta_ordinal = dyn_comp_ord,
    paired_delta_High_vs_Low = dyn_comp_bin,
    paired_delta_exact_permutation = exact_abundance
  )
)

make_heatmap_plot <- function(mat, ann, title, top_n = 18, show_legend = FALSE) {
  ann2 <- ann %>%
    mutate(
      response_ord = factor(response_ord, levels = c("Low", "Medium", "High")),
      patient = factor(patient, levels = patient)
    ) %>%
    arrange(response_ord, patient)
  df <- matrix_to_long(mat, ann2, top_n = top_n)

  ggplot(df, aes(patient, gene, fill = z_delta)) +
    geom_tile(colour = "white", linewidth = 0.25) +
    scale_fill_gradient2(
      low = heat_cols[1], mid = heat_cols[2], high = heat_cols[3],
      midpoint = 0, limits = c(-2.5, 2.5), oob = squish, name = "Row z",
      guide = if (show_legend) guide_colorbar(barheight = grid::unit(2.2, "cm")) else "none"
    ) +
    labs(title = title, subtitle = "Columns ordered by response: Low, Medium, High", x = NULL, y = NULL) +
    theme_submission() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      axis.ticks = element_blank(),
      legend.position = "right"
    )
}

prepare_gsea <- function(df, top_n = 10) {
  df %>%
    filter(!is.na(NES), !is.na(padj)) %>%
    arrange(padj, desc(abs(NES))) %>%
    slice_head(n = top_n) %>%
    mutate(
      pathway_label = clean_pathway(pathway),
      pathway_label = factor(pathway_label, levels = rev(pathway_label)),
      neg_log10_fdr = pmin(-log10(padj), 25)
    )
}

make_gsea_plot <- function(df, title) {
  plot_df <- prepare_gsea(df, top_n = 10)
  x_low <- max(0, floor((min(plot_df$NES, na.rm = TRUE) - 0.08) * 10) / 10)
  x_high <- ceiling((max(plot_df$NES, na.rm = TRUE) + 0.08) * 10) / 10
  ggplot(plot_df, aes(NES, pathway_label)) +
    geom_point(aes(size = neg_log10_fdr), colour = "#B22222", alpha = 0.95) +
    scale_size_continuous(range = c(1.6, 4.7), name = "-log10 FDR") +
    scale_x_continuous(limits = c(x_low, x_high), breaks = pretty_breaks(n = 4)) +
    labs(title = title, x = "Normalized enrichment score", y = NULL) +
    theme_submission() +
    theme(
      legend.position = "right",
      legend.key.size = grid::unit(0.28, "cm"),
      plot.title = element_text(size = 8.8, face = "bold", margin = margin(b = 2)),
      axis.text.y = element_text(size = 7.4),
      plot.margin = margin(2, 4, 2, 4)
    )
}

p2a <- make_heatmap_plot(t_mat, t_ann, "T-cell dynamic interaction genes", top_n = 18, show_legend = FALSE)
p2b <- make_heatmap_plot(m_mat, m_ann, "Myeloid dynamic interaction genes", top_n = 18, show_legend = TRUE)
p2c <- make_gsea_plot(t_gsea, "T-cell Hallmark enrichment")
p2d <- make_gsea_plot(m_gsea, "Myeloid Hallmark enrichment")

fig2 <- (p2a | p2b) / (p2c | p2d) +
  plot_layout(heights = c(1.16, 0.88), guides = "keep") +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 8))

save_plot_set(fig2, "Figure2_submission_state_remodeling", width = 13.2, height = 10.4)

write_source_data(
  "Figure2",
  list(
    Tcell_delta_matrix_top18 = tibble::rownames_to_column(as.data.frame(t_mat[seq_len(min(18, nrow(t_mat))), , drop = FALSE]), "gene"),
    Tcell_sample_annotation = t_ann,
    Tcell_GSEA_top10 = prepare_gsea(t_gsea, top_n = 10),
    Myeloid_delta_matrix_top18 = tibble::rownames_to_column(as.data.frame(m_mat[seq_len(min(18, nrow(m_mat))), , drop = FALSE]), "gene"),
    Myeloid_sample_annotation = m_ann,
    Myeloid_GSEA_top10 = prepare_gsea(m_gsea, top_n = 10)
  )
)

manifest <- tibble::tibble(
  generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  figure = c("Figure1", "Figure2"),
  stem = c("Figure1_submission_abundance_and_design", "Figure2_submission_state_remodeling"),
  primary_message = c(
    "Baseline abundance is weak, while paired post-pre deltas nominate treatment-linked compartment shifts.",
    "T-cell and myeloid pseudobulk interactions converge on inflammatory, interferon and mTORC1 pathway remodeling."
  )
)
readr::write_csv(manifest, file.path(figure_dir, "SUBMISSION_FIGURE_MANIFEST.csv"))

session_file <- file.path(figure_dir, "SESSION_INFO.txt")
sink(session_file)
cat("Generated at: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n", sep = "")
cat("Workspace: ", workspace_dir, "\n", sep = "")
cat("Script: ", script_path, "\n\n", sep = "")
sessionInfo()
sink()

message("Done. Submission figures and source data written to: ", figure_dir)
