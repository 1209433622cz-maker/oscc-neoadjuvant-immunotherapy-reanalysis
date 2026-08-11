#!/usr/bin/env Rscript

rm(list = ls())
gc()
set.seed(1234)

workspace <- Sys.getenv("GSE200996_WORKSPACE")
if (workspace == "") workspace <- normalizePath(".", winslash = "/", mustWork = TRUE)
rebuild <- file.path(workspace, "03_rebuild")
result_dir <- file.path(rebuild, "results", "discovery_lineage_composition_sensitivity")
figure_dir <- file.path(rebuild, "figures", "submission")
source_dir <- file.path(figure_dir, "source_data")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(stringr)
  library(scales)
})
if (requireNamespace("BiocParallel", quietly = TRUE)) {
  BiocParallel::register(BiocParallel::SerialParam())
}

strict_patients <- c("P24", "P27", "P29", "P18", "P23", "P32")
response_map <- c(P24 = "Low", P27 = "Low", P29 = "Low",
                  P18 = "Medium", P23 = "Medium", P32 = "High")
cohort_map <- c(P24 = "Combo", P27 = "Mono", P29 = "Mono",
                P18 = "Combo", P23 = "Mono", P32 = "Combo")
heat_colours <- c("#2C7BB6", "#F7F7F7", "#D7191C")

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

theme_nature <- function(base_size = 7) {
  theme_classic(base_size = base_size) +
    theme(
      axis.line = element_line(linewidth = 0.3, colour = "black"),
      axis.ticks = element_line(linewidth = 0.25, colour = "black"),
      axis.text = element_text(colour = "black"),
      plot.title = element_text(face = "bold", size = 7, margin = margin(b = 3)),
      plot.subtitle = element_text(size = 7.0, colour = "grey30", margin = margin(b = 3)),
      plot.tag = element_text(face = "bold", size = 8),
      legend.title = element_text(face = "bold", size = 6.8),
      legend.text = element_text(size = 6.6),
      legend.key.height = grid::unit(3.5, "mm"),
      plot.margin = margin(3, 3, 3, 3)
    )
}

row_z <- function(matrix) {
  z <- t(scale(t(matrix)))
  z[!is.finite(z)] <- 0
  pmax(pmin(z, 2.5), -2.5)
}

fit_and_delta <- function(lineage) {
  stem <- gsub(" ", "_", lineage)
  pb <- readRDS(file.path(result_dir, paste0(stem, "_REFINED_PSEUDOBULK_INPUT.rds")))
  md <- pb$meta
  rownames(md) <- md$pb_id
  md$patient_id <- factor(md$patient_id)
  md$cohort <- factor(md$cohort)
  design <- ~ patient_id + post + post:resp_num
  dds <- DESeqDataSetFromMatrix(
    round(pb$counts[rowSums(pb$counts) >= 10, , drop = FALSE]),
    md,
    design
  )
  dds <- DESeq(dds, quiet = TRUE)
  transformed <- assay(varianceStabilizingTransformation(dds, blind = FALSE))
  delta <- sapply(strict_patients, function(patient) {
    transformed[, paste0(patient, "|post")] - transformed[, paste0(patient, "|pre")]
  })
  result <- read.csv(
    file.path(result_dir, paste0(stem, "_primary_GENE_RESULTS.csv")),
    stringsAsFactors = FALSE
  )
  result <- result[is.finite(result$pvalue), , drop = FALSE]
  result <- result[order(result$padj, result$pvalue), , drop = FALSE]
  ordered_genes <- result$gene[result$gene %in% rownames(delta)]
  delta <- delta[ordered_genes, strict_patients, drop = FALSE]
  annotation <- data.frame(
    patient = strict_patients,
    response_ord = unname(response_map[strict_patients]),
    cohort = unname(cohort_map[strict_patients]),
    stringsAsFactors = FALSE
  )
  list(delta = delta, annotation = annotation, result = result)
}

make_heatmap <- function(delta, annotation, title, show_legend) {
  top <- delta[seq_len(min(18, nrow(delta))), , drop = FALSE]
  z <- row_z(top)
  long <- as.data.frame(z, check.names = FALSE) |>
    tibble::rownames_to_column("gene") |>
    pivot_longer(-gene, names_to = "patient", values_to = "z_delta") |>
    left_join(annotation, by = "patient") |>
    mutate(
      patient = factor(patient, levels = strict_patients),
      gene = factor(gene, levels = rev(rownames(top)))
    )
  ggplot(long, aes(patient, gene, fill = z_delta)) +
    geom_tile(colour = "white", linewidth = 0.20) +
    scale_fill_gradient2(
      low = heat_colours[1], mid = heat_colours[2], high = heat_colours[3],
      midpoint = 0, limits = c(-2.5, 2.5), oob = squish,
      name = "Row z",
      guide = if (show_legend) guide_colorbar(barheight = grid::unit(17, "mm")) else "none"
    ) +
    labs(
      title = title,
      subtitle = "Low                 Medium          High",
      x = NULL,
      y = NULL
    ) +
    theme_nature() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 6.7),
      axis.text.y = element_text(size = 6.4),
      axis.ticks = element_blank(),
      legend.position = "right"
    )
}

prepare_gsea <- function(lineage, top_n = 10) {
  stem <- gsub(" ", "_", lineage)
  primary <- read.csv(
    file.path(result_dir, paste0(stem, "_primary_HALLMARK_GSEA.csv")),
    stringsAsFactors = FALSE
  )
  adjusted <- read.csv(
    file.path(result_dir, paste0(stem, "_composition_adjusted_HALLMARK_GSEA.csv")),
    stringsAsFactors = FALSE
  ) |>
    select(pathway, adjusted_NES = NES, adjusted_FDR = padj)
  primary |>
    filter(is.finite(NES), is.finite(padj)) |>
    arrange(padj, desc(abs(NES))) |>
    slice_head(n = top_n) |>
    left_join(adjusted, by = "pathway") |>
    mutate(
      robustness = ifelse(
        is.finite(adjusted_FDR) &
          adjusted_FDR < 0.05 &
          sign(NES) == sign(adjusted_NES),
        "Composition-robust",
        "Composition-sensitive"
      ),
      pathway_label = clean_pathway(pathway),
      pathway_label = factor(pathway_label, levels = rev(pathway_label)),
      neg_log10_fdr = pmin(-log10(padj), 25),
      lineage = lineage
    )
}

make_gsea <- function(data, title) {
  x_low <- floor((min(data$NES) - 0.08) * 10) / 10
  x_high <- ceiling((max(data$NES) + 0.08) * 10) / 10
  ggplot(data, aes(NES, pathway_label)) +
    geom_point(
      aes(size = neg_log10_fdr, shape = robustness),
      colour = "#B22222",
      stroke = 0.65
    ) +
    scale_shape_manual(
      values = c("Composition-robust" = 16, "Composition-sensitive" = 1)
    ) +
    scale_size_continuous(range = c(1.5, 4.2), name = "-log10 FDR") +
    scale_x_continuous(limits = c(x_low, x_high), breaks = pretty_breaks(n = 4)) +
    labs(title = title, x = "Primary normalized enrichment score", y = NULL) +
    guides(
      shape = guide_legend(title = NULL, order = 1),
      size = guide_legend(title = "-log10 FDR", order = 2)
    ) +
    theme_nature() +
    theme(
      axis.text.y = element_text(size = 6.5),
      legend.position = "right",
      legend.key.height = grid::unit(3.0, "mm")
    )
}

t_cell <- fit_and_delta("T cell")
myeloid <- fit_and_delta("Myeloid")
t_gsea <- prepare_gsea("T cell")
m_gsea <- prepare_gsea("Myeloid")

p_a <- make_heatmap(t_cell$delta, t_cell$annotation, "Refined T-cell interaction genes", FALSE)
p_b <- make_heatmap(myeloid$delta, myeloid$annotation, "Refined myeloid interaction genes", TRUE)
p_c <- make_gsea(t_gsea, "T-cell Hallmark enrichment")
p_d <- make_gsea(m_gsea, "Myeloid Hallmark enrichment")

figure <- (p_a | p_b) / (p_c | p_d) +
  plot_layout(heights = c(1.1, 0.9), guides = "keep") +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 8))

stem <- file.path(figure_dir, "Figure2_submission_state_remodeling")
ggsave(paste0(stem, ".png"), figure, width = 7.2, height = 6.1, dpi = 600, bg = "white")
ggsave(paste0(stem, ".pdf"), figure, width = 7.2, height = 6.1, bg = "white")
ggsave(paste0(stem, ".svg"), figure, width = 7.2, height = 6.1, bg = "white")

write.csv(
  tibble::rownames_to_column(
    as.data.frame(t_cell$delta[seq_len(min(18, nrow(t_cell$delta))), , drop = FALSE]),
    "gene"
  ),
  file.path(source_dir, "Figure2_Tcell_delta_matrix_top18.csv"),
  row.names = FALSE
)
write.csv(
  tibble::rownames_to_column(
    as.data.frame(myeloid$delta[seq_len(min(18, nrow(myeloid$delta))), , drop = FALSE]),
    "gene"
  ),
  file.path(source_dir, "Figure2_Myeloid_delta_matrix_top18.csv"),
  row.names = FALSE
)
write.csv(t_cell$annotation, file.path(source_dir, "Figure2_Tcell_annotation.csv"), row.names = FALSE)
write.csv(myeloid$annotation, file.path(source_dir, "Figure2_Myeloid_annotation.csv"), row.names = FALSE)
write.csv(t_gsea, file.path(source_dir, "Figure2_Tcell_GSEA_top10.csv"), row.names = FALSE)
write.csv(m_gsea, file.path(source_dir, "Figure2_Myeloid_GSEA_top10.csv"), row.names = FALSE)

report <- c(
  "# Refined Discovery Figure 2",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "- Heatmaps use author-resolved CD4/CD8 and myeloid lineage definitions.",
  "- Display values are VST post-minus-pre deltas, row-standardized only for visualization.",
  "- Genes are ordered by the refined primary DESeq2 interaction result.",
  "- Filled GSEA points retain direction and FDR < 0.05 after response-blind composition-PC1 adjustment.",
  "- Open GSEA points are composition-sensitive and cannot carry the robust mechanistic claim.",
  "",
  paste0("T-cell primary FDR < 0.05 genes: ", sum(t_cell$result$padj < 0.05, na.rm = TRUE)),
  paste0("Myeloid primary FDR < 0.05 genes: ", sum(myeloid$result$padj < 0.05, na.rm = TRUE)),
  paste0("T-cell robust displayed pathways: ", sum(t_gsea$robustness == "Composition-robust"), "/", nrow(t_gsea)),
  paste0("Myeloid robust displayed pathways: ", sum(m_gsea$robustness == "Composition-robust"), "/", nrow(m_gsea))
)
writeLines(report, file.path(result_dir, "REFINED_FIGURE2_REPORT.md"))
message("Completed refined Figure 2: ", stem)
