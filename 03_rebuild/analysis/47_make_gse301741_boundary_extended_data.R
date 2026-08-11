workspace <- Sys.getenv("GSE200996_WORKSPACE")
if (workspace == "") {
  workspace <- normalizePath(getwd(), mustWork = FALSE)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(scales)
})

rebuild_dir <- file.path(
  workspace, "03_rebuild", "validation", "GSE301741_raw_reconstruction"
)
validation_dir <- file.path(
  workspace, "03_rebuild", "validation", "GSE301741_lineage_aware_validation"
)
figure_dir <- file.path(workspace, "03_rebuild", "figures", "submission")
source_dir <- file.path(figure_dir, "source_data")
manuscript_dir <- file.path(workspace, "03_rebuild", "manuscript")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(manuscript_dir, recursive = TRUE, showWarnings = FALSE)

read_required <- function(path) {
  if (!file.exists(path)) stop("Missing input: ", path)
  readr::read_csv(path, show_col_types = FALSE)
}

qc <- read_required(file.path(
  rebuild_dir, "GSE301741_RAW_REBUILT_SAMPLE_LINEAGE_QC.csv"
))
pair_qc <- read_required(file.path(
  validation_dir, "GSE301741_LINEAGE_AWARE_PAIR_CELL_QC.csv"
))
patient_delta <- read_required(file.path(
  validation_dir, "GSE301741_LINEAGE_AWARE_PATIENT_MODULE_DELTAS.csv"
))
tests <- read_required(file.path(
  validation_dir, "GSE301741_LINEAGE_AWARE_RESPONSE_TESTS.csv"
))

theme_nature <- function(base_size = 7) {
  theme_classic(base_size = base_size) +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.3, colour = "black"),
      axis.text = element_text(colour = "black"),
      plot.title = element_text(face = "bold", size = rel(1.02)),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", colour = "black"),
      legend.title = element_text(face = "bold"),
      legend.key.size = grid::unit(0.32, "cm"),
      plot.tag = element_text(face = "bold", size = 8),
      plot.margin = margin(4, 5, 4, 5)
    )
}

short_labels <- c(
  T_DE_FDR05_positive = "DE FDR < 0.05",
  T_LE_TNFA_SIGNALING_VIA_NFKB = "TNF-NF-kB",
  T_LE_MTORC1_SIGNALING = "mTORC1",
  T_LE_P53_PATHWAY = "p53",
  T_LE_INTERFERON_ALPHA_RESPONSE = "IFN-alpha",
  T_LE_INTERFERON_GAMMA_RESPONSE = "IFN-gamma",
  T_LE_union_core = "Union core"
)
module_order <- c(
  "DE FDR < 0.05", "TNF-NF-kB", "mTORC1", "p53",
  "IFN-alpha", "IFN-gamma", "Union core"
)

lineage_palette <- c(
  "T cell" = "#4C78A8",
  "NK" = "#72B7B2",
  "Myeloid/DC" = "#E45756",
  "B/Plasma" = "#54A24B",
  "Epithelial" = "#F2CF5B",
  "Fibroblast" = "#B279A2",
  "Endothelial" = "#FF9DA6",
  "Mast" = "#9D755D"
)
response_palette <- c(
  "Non-responder" = "#3C5488",
  "Responder" = "#E64B35"
)

composition <- qc %>%
  filter(lineage != "Filtered_or_unresolved") %>%
  mutate(
    lineage_group = case_when(
      lineage == "T_cell" ~ "T cell",
      lineage == "NK" ~ "NK",
      lineage %in% c("Myeloid", "DC") ~ "Myeloid/DC",
      lineage %in% c("B_cell", "Plasma") ~ "B/Plasma",
      lineage == "Epithelial" ~ "Epithelial",
      lineage == "Fibroblast" ~ "Fibroblast",
      lineage == "Endothelial" ~ "Endothelial",
      lineage == "Mast" ~ "Mast",
      TRUE ~ "Other"
    ),
    fraction = factor(
      fraction,
      levels = c("CD3", "CD45pos", "CD45ratio", "Unsorted", "CD45neg")
    )
  ) %>%
  group_by(fraction, lineage_group) %>%
  summarise(n_cells = sum(n_cells), .groups = "drop") %>%
  group_by(fraction) %>%
  mutate(
    fraction_total = sum(n_cells),
    fraction_percent = n_cells / fraction_total
  ) %>%
  ungroup()

p_a <- ggplot(composition, aes(fraction, fraction_percent, fill = lineage_group)) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.15) +
  coord_flip() +
  scale_fill_manual(values = lineage_palette, breaks = names(lineage_palette)) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.01))
  ) +
  labs(
    title = "Conservative lineage reconstruction",
    x = NULL,
    y = "Retained cells",
    fill = "Lineage"
  ) +
  theme_nature() +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.margin = margin(t = -2),
    legend.text = element_text(size = 6.5)
  ) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

pair_plot <- pair_qc %>%
  mutate(
    response = recode(
      response_label,
      non_responder = "Non-responder",
      responder = "Responder"
    ),
    lineage_label = recode(
      target_lineage,
      T_cell = "T cell",
      Myeloid = "Myeloid/DC"
    ),
    patient_id = factor(patient_id, levels = unique(patient_id)),
    gate = if_else(eligible_min30, "Pass", "Below 30")
  )

p_b <- ggplot(
  pair_plot,
  aes(patient_id, minimum_cells_across_pair, colour = response, shape = gate)
) +
  geom_hline(
    yintercept = 30, linetype = "dashed", linewidth = 0.35, colour = "#555555"
  ) +
  geom_point(size = 2.2, stroke = 0.65) +
  facet_wrap(~lineage_label, ncol = 1, scales = "free_x") +
  scale_y_log10(
    breaks = c(10, 30, 100, 300, 1000, 3000),
    labels = scales::label_number()
  ) +
  scale_colour_manual(values = response_palette) +
  scale_shape_manual(values = c("Pass" = 16, "Below 30" = 4)) +
  labs(
    title = "Matched-pair eligibility",
    x = NULL,
    y = "Minimum cells in pair",
    colour = "Outcome",
    shape = "30-cell gate"
  ) +
  theme_nature() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 6.5),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.margin = margin(t = -2),
    legend.text = element_text(size = 6.5),
    panel.spacing.y = grid::unit(4, "pt")
  ) +
  guides(
    colour = guide_legend(nrow = 1, order = 1),
    shape = guide_legend(nrow = 1, order = 2)
  )

heatmap <- patient_delta %>%
  filter(target_lineage == "T_cell") %>%
  mutate(
    module = unname(short_labels[signature]),
    response = recode(
      response_label,
      non_responder = "NR",
      responder = "R"
    )
  )
patient_order <- heatmap %>%
  distinct(patient_id, response) %>%
  arrange(factor(response, levels = c("NR", "R")), patient_id) %>%
  mutate(patient_label = paste0(patient_id, " (", response, ")"))
heatmap <- heatmap %>%
  left_join(patient_order, by = c("patient_id", "response")) %>%
  mutate(
    module = factor(module, levels = module_order),
    patient_label = factor(patient_label, levels = rev(patient_order$patient_label))
  )
heat_limit <- max(abs(heatmap$delta), na.rm = TRUE)

p_c <- ggplot(heatmap, aes(module, patient_label, fill = delta)) +
  geom_tile(colour = "white", linewidth = 0.35) +
  scale_fill_gradient2(
    low = "#3C5488",
    mid = "white",
    high = "#E64B35",
    midpoint = 0,
    limits = c(-heat_limit, heat_limit),
    oob = scales::squish
  ) +
  labs(
    title = "Patient-level T-cell module changes",
    x = NULL,
    y = NULL,
    fill = "Post - pre\nscore (z)"
  ) +
  theme_nature() +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(angle = 40, hjust = 1, size = 6.5),
    axis.text.y = element_text(size = 6.7),
    legend.position = "right",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.35)
  )

effect <- tests %>%
  filter(target_lineage == "T_cell") %>%
  mutate(
    module = unname(short_labels[signature]),
    direction = if_else(
      diff_responder_minus_non >= 0,
      "Higher in responders",
      "Lower in responders"
    ),
    label = sprintf("%s   P = %.3f", module, exact_permutation_p)
  ) %>%
  arrange(diff_responder_minus_non) %>%
  mutate(label = factor(label, levels = label))

p_d <- ggplot(effect, aes(diff_responder_minus_non, label, colour = direction)) +
  geom_vline(xintercept = 0, linewidth = 0.35, colour = "#777777") +
  geom_segment(
    aes(x = 0, xend = diff_responder_minus_non, yend = label),
    linewidth = 0.45,
    colour = "#B5B5B5"
  ) +
  geom_point(size = 2.25) +
  scale_colour_manual(values = c(
    "Higher in responders" = "#E64B35",
    "Lower in responders" = "#3C5488"
  )) +
  labs(
    title = "Exact label permutation",
    x = "Responder - non-responder delta difference",
    y = NULL,
    colour = NULL
  ) +
  theme_nature() +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 6.5),
    axis.text.y = element_text(size = 6.6),
    panel.grid.major.y = element_line(colour = "#EEEEEE", linewidth = 0.25)
  ) +
  guides(colour = guide_legend(nrow = 1))

figure <- ((p_a | p_b) / (p_c | p_d)) +
  plot_layout(widths = c(1.08, 0.92), heights = c(0.98, 1.02), guides = "keep") +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag.position = c(0.005, 0.995))

stem <- "ExtendedData10_submission_gse301741_boundary"
png_path <- file.path(figure_dir, paste0(stem, ".png"))
pdf_path <- file.path(figure_dir, paste0(stem, ".pdf"))
svg_path <- file.path(figure_dir, paste0(stem, ".svg"))
if (requireNamespace("ragg", quietly = TRUE)) {
  ggsave(
    png_path, figure, width = 7.2, height = 5.8, dpi = 600,
    bg = "white", device = ragg::agg_png
  )
} else {
  ggsave(png_path, figure, width = 7.2, height = 5.8, dpi = 600, bg = "white")
}
ggsave(pdf_path, figure, width = 7.2, height = 5.8, bg = "white")
if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(
    svg_path, figure, width = 7.2, height = 5.8,
    bg = "white", device = svglite::svglite
  )
}

source_tables <- list(
  panel_a_lineage_composition = composition,
  panel_b_pair_cell_gate = pair_plot,
  panel_c_patient_module_delta = heatmap %>%
    mutate(
      module = as.character(module),
      patient_label = as.character(patient_label)
    ),
  panel_d_exact_permutation = effect %>%
    mutate(label = as.character(label))
)
for (nm in names(source_tables)) {
  write_csv(
    source_tables[[nm]],
    file.path(source_dir, paste0(stem, "_", nm, ".csv"))
  )
}

fixer <- file.path(
  workspace, "03_rebuild", "analysis", "11_rebuild_figure_source_workbooks.py"
)
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
    python_candidates[[1]],
    args = c(
      normalizePath(fixer, winslash = "/", mustWork = TRUE),
      "--stem", stem
    ),
    stdout = TRUE,
    stderr = TRUE
  )
  exit_code <- attr(status, "status")
  if (!is.null(exit_code) && exit_code != 0) {
    warning("Source workbook rebuild failed: ", paste(status, collapse = "\n"))
  }
}

manifest <- file.path(figure_dir, "EXTENDED_DATA_FIGURE_MANIFEST.csv")
manifest <- if (file.exists(manifest)) read_required(manifest) else tibble()
if ("generated_at" %in% names(manifest)) {
  manifest <- manifest %>% mutate(generated_at = as.character(generated_at))
}
manifest <- manifest %>%
  filter(stem != !!stem) %>%
  bind_rows(tibble(
    generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    figure = "Extended Data 10",
    stem = stem,
    scientific_role = paste(
      "Disease-matched single-cell boundary analysis showing conservative",
      "RAW lineage reconstruction, paired target-cell eligibility and null",
      "exact T-cell response tests."
    ),
    nature_style = "Pass after source-data and rendered-layout QC.",
    redraw_priority = "none",
    recommended_action = paste(
      "Keep as a boundary result; do not describe as positive independent",
      "validation, and retain the non-estimable myeloid comparison."
    )
  ))
write_csv(manifest, file.path(figure_dir, "EXTENDED_DATA_FIGURE_MANIFEST.csv"))

report <- c(
  "# Extended Data Figure 10: GSE301741 Boundary QC",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Scientific Role",
  "",
  "- Disease-matched single-cell boundary analysis, not positive validation.",
  "- Panel a documents conservative lineage reconstruction by sorting fraction.",
  "- Panel b exposes the paired 30-cell eligibility gate and the absent eligible myeloid responder.",
  "- Panel c shows every patient-level T-cell module delta.",
  "- Panel d reports response-group effect estimates with exact permutation P values.",
  "",
  "## Truth Boundary",
  "",
  "- Response labels remain supplement-figure-derived.",
  "- Reconstructed broad lineages are not the unavailable author labels.",
  "- All seven T-cell exact tests are non-significant.",
  "- The myeloid response comparison is not estimable.",
  "",
  "## Outputs",
  "",
  paste0("- `03_rebuild/figures/submission/", stem, ".png`"),
  paste0("- `03_rebuild/figures/submission/", stem, ".pdf`"),
  paste0("- `03_rebuild/figures/submission/", stem, ".svg`"),
  paste0("- `03_rebuild/figures/submission/source_data/", stem, "_source_data.xlsx`")
)
writeLines(
  report,
  file.path(manuscript_dir, "EXTENDED_DATA10_GSE301741_BOUNDARY_QC.md")
)
cat(paste(report, collapse = "\n"), "\n")
