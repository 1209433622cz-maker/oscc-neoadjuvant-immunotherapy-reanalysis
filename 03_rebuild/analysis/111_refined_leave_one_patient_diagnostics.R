rm(list = ls())
gc()
set.seed(20260727)

base_dir <- Sys.getenv("GSE200996_BASEDIR")
if (base_dir == "") base_dir <- "H:/SCI2/OSCC-GSE200996-2025.12/03_rebuild"

input_dir <- file.path(base_dir, "results", "discovery_lineage_composition_sensitivity")
out_dir <- file.path(base_dir, "results", "refined_leave_one_out")
figure_dir <- file.path(base_dir, "figures", "submission")
source_dir <- file.path(figure_dir, "source_data")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(DESeq2)
  library(fgsea)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

if (requireNamespace("BiocParallel", quietly = TRUE)) {
  BiocParallel::register(BiocParallel::SerialParam())
}

load_hallmark_sets <- function() {
  cache_file <- file.path(Sys.getenv("LOCALAPPDATA"), "R", "cache", "R", "msigdbr", "msigdb.2025.1.Hs.rds")
  if (!file.exists(cache_file)) stop("Local MSigDB cache not found: ", cache_file)
  msig <- readRDS(cache_file)
  msig <- msig[
    msig$gs_collection == "H" &
      msig$db_target_species == "HS" &
      !is.na(msig$db_gene_symbol),
  ]
  split(msig$db_gene_symbol, msig$gs_name)
}

hallmark <- load_hallmark_sets()

key_pathways <- list(
  "T cell" = c(
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    "HALLMARK_P53_PATHWAY",
    "HALLMARK_E2F_TARGETS",
    "HALLMARK_G2M_CHECKPOINT",
    "HALLMARK_MITOTIC_SPINDLE",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE"
  ),
  "Myeloid" = c(
    "HALLMARK_MTORC1_SIGNALING",
    "HALLMARK_INTERFERON_ALPHA_RESPONSE",
    "HALLMARK_INFLAMMATORY_RESPONSE",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    "HALLMARK_COMPLEMENT"
  )
)

fit_interaction <- function(pb, leave_out = NA_character_, seed = 20260727L) {
  md <- pb$meta
  counts <- pb$counts
  if (!is.na(leave_out)) {
    md <- md[as.character(md$patient_id) != leave_out, , drop = FALSE]
    counts <- counts[, md$pb_id, drop = FALSE]
  }

  md$patient_id <- droplevels(factor(as.character(md$patient_id)))
  md$timepoint <- factor(as.character(md$timepoint), levels = c("pre", "post"))
  md$post <- as.integer(md$timepoint == "post")
  md$resp_num <- as.numeric(md$resp_num)
  rownames(md) <- md$pb_id
  counts <- round(counts[rowSums(counts) >= 10, rownames(md), drop = FALSE])

  if (length(unique(md$patient_id)) < 4) stop("Too few paired patients after leave-out")
  if (length(unique(md$resp_num)) < 2) stop("Too few response levels after leave-out")

  dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = md,
    design = ~ patient_id + post + post:resp_num
  )
  dds <- DESeq(dds, quiet = TRUE)
  int_name <- grep("post.*resp_num|resp_num.*post", resultsNames(dds), value = TRUE)[1]
  if (is.na(int_name)) stop("Interaction coefficient not found")

  result <- as.data.frame(results(dds, name = int_name))
  result$gene <- rownames(result)
  result <- result[is.finite(result$stat), , drop = FALSE]
  result <- result[
    order(ifelse(is.na(result$padj), Inf, result$padj), ifelse(is.na(result$pvalue), Inf, result$pvalue)),
    ,
    drop = FALSE
  ]
  ranks <- setNames(result$stat, result$gene)
  ranks <- sort(ranks[is.finite(ranks)], decreasing = TRUE)
  set.seed(seed)
  gsea <- as.data.frame(
    fgseaMultilevel(
      pathways = hallmark,
      stats = ranks,
      minSize = 10,
      maxSize = 500,
      BPPARAM = BiocParallel::SerialParam()
    )
  )
  gsea <- gsea[order(gsea$padj), , drop = FALSE]
  list(
    result = result,
    gsea = gsea,
    n_patients = length(unique(md$patient_id)),
    n_samples = nrow(md),
    residual_df = nrow(md) - ncol(model.matrix(~ patient_id + post + post:resp_num, data = md))
  )
}

top_genes <- function(result, n = 100) {
  head(
    result$gene[
      order(ifelse(is.na(result$padj), Inf, result$padj), ifelse(is.na(result$pvalue), Inf, result$pvalue))
    ],
    n
  )
}

summary_list <- list()
pathway_list <- list()
direction_list <- list()
full_results <- list()

inputs <- c(
  "T cell" = "T_cell_REFINED_PSEUDOBULK_INPUT.rds",
  "Myeloid" = "Myeloid_REFINED_PSEUDOBULK_INPUT.rds"
)

for (lineage in names(inputs)) {
  message("Refined leave-one-out: ", lineage)
  pb <- readRDS(file.path(input_dir, inputs[[lineage]]))
  patients <- sort(unique(as.character(pb$meta$patient_id)))
  full <- fit_interaction(pb, seed = 20260727L + match(lineage, names(inputs)))
  full_results[[lineage]] <- full$result
  full_stat <- setNames(full$result$stat, full$result$gene)
  full_top25 <- top_genes(full$result, 25)
  full_top50 <- top_genes(full$result, 50)
  full_top100 <- top_genes(full$result, 100)

  for (index in seq_along(c("NONE", patients))) {
    leave_out <- c("NONE", patients)[index]
    message("  leave_out = ", leave_out)
    fit <- if (leave_out == "NONE") {
      full
    } else {
      tryCatch(
        fit_interaction(
          pb,
          leave_out = leave_out,
          seed = 20260727L + 100L * match(lineage, names(inputs)) + index
        ),
        error = function(error) error
      )
    }

    if (inherits(fit, "error")) {
      summary_list[[length(summary_list) + 1]] <- data.frame(
        lineage = lineage,
        leave_out = leave_out,
        status = "FAIL",
        n_patients = NA_integer_,
        n_samples = NA_integer_,
        residual_df = NA_integer_,
        n_common_genes = NA_integer_,
        stat_spearman = NA_real_,
        top50_overlap = NA_real_,
        top100_overlap = NA_real_,
        full_top25_direction_concordance = NA_real_,
        note = conditionMessage(fit)
      )
      next
    }

    loo_stat <- setNames(fit$result$stat, fit$result$gene)
    common <- intersect(names(full_stat), names(loo_stat))
    top25_common <- intersect(full_top25, names(loo_stat))
    summary_list[[length(summary_list) + 1]] <- data.frame(
      lineage = lineage,
      leave_out = leave_out,
      status = "OK",
      n_patients = fit$n_patients,
      n_samples = fit$n_samples,
      residual_df = fit$residual_df,
      n_common_genes = length(common),
      stat_spearman = suppressWarnings(cor(full_stat[common], loo_stat[common], method = "spearman")),
      top50_overlap = length(intersect(full_top50, top_genes(fit$result, 50))) / 50,
      top100_overlap = length(intersect(full_top100, top_genes(fit$result, 100))) / 100,
      full_top25_direction_concordance = mean(sign(full_stat[top25_common]) == sign(loo_stat[top25_common])),
      note = ""
    )

    direction_list[[length(direction_list) + 1]] <- data.frame(
      lineage = lineage,
      leave_out = leave_out,
      gene = full_top25,
      full_stat = as.numeric(full_stat[full_top25]),
      loo_stat = as.numeric(loo_stat[full_top25]),
      same_direction = sign(full_stat[full_top25]) == sign(loo_stat[full_top25])
    )

    pathways <- fit$gsea[
      fit$gsea$pathway %in% key_pathways[[lineage]],
      c("pathway", "NES", "pval", "padj", "size"),
      drop = FALSE
    ]
    pathways$lineage <- lineage
    pathways$leave_out <- leave_out
    pathway_list[[length(pathway_list) + 1]] <- pathways
  }
}

summary_df <- bind_rows(summary_list)
pathway_df <- bind_rows(pathway_list)
direction_df <- bind_rows(direction_list)

write.csv(summary_df, file.path(out_dir, "REFINED_LOO_MODEL_STABILITY.csv"), row.names = FALSE)
write.csv(pathway_df, file.path(out_dir, "REFINED_LOO_KEY_PATHWAY_GSEA.csv"), row.names = FALSE)
write.csv(direction_df, file.path(out_dir, "REFINED_LOO_TOP25_DIRECTION.csv"), row.names = FALSE)

clean_pathway <- function(pathway) {
  pathway |>
    sub("^HALLMARK_", "", x = _) |>
    gsub("_", " ", x = _) |>
    tools::toTitleCase()
}

compact_theme <- theme_classic(base_size = 7) +
  theme(
    plot.title = element_text(face = "bold", size = 7),
    strip.text = element_text(face = "bold", size = 7),
    plot.tag = element_text(face = "bold", size = 8),
    legend.key.height = grid::unit(3.5, "mm")
  )

loo_plot <- summary_df |>
  filter(status == "OK", leave_out != "NONE") |>
  mutate(
    leave_out = factor(leave_out, levels = c("P18", "P23", "P24", "P27", "P29", "P32")),
    influence = ifelse(leave_out == "P32", "P32 removed", "Other patient removed")
  )

p_a <- ggplot(loo_plot, aes(leave_out, stat_spearman, fill = influence)) +
  geom_col(width = 0.68) +
  facet_wrap(~ lineage, nrow = 1) +
  scale_fill_manual(values = c("Other patient removed" = "#858D96", "P32 removed" = "#B4473F"), name = NULL) +
  scale_y_continuous(limits = c(0, 1), labels = number_format(accuracy = 0.1)) +
  compact_theme +
  theme(legend.position = "bottom") +
  labs(title = "Genome-wide rank stability", x = "Left-out patient", y = "Spearman correlation")

p_b <- ggplot(loo_plot, aes(leave_out, top100_overlap, fill = influence)) +
  geom_col(width = 0.68) +
  facet_wrap(~ lineage, nrow = 1) +
  scale_fill_manual(values = c("Other patient removed" = "#858D96", "P32 removed" = "#B4473F"), guide = "none") +
  scale_y_continuous(limits = c(0, 1), labels = percent_format(accuracy = 1)) +
  compact_theme +
  labs(title = "Top-100 gene overlap", x = "Left-out patient", y = "Overlap")

path_plot <- pathway_df |>
  mutate(
    leave_out = factor(leave_out, levels = c("NONE", "P18", "P23", "P24", "P27", "P29", "P32")),
    pathway_label = clean_pathway(pathway),
    fdr_status = ifelse(padj < 0.05, "FDR < 0.05", "Not significant")
  )

p_c <- ggplot(path_plot, aes(leave_out, pathway_label, fill = NES)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  geom_point(aes(shape = fdr_status), size = 1.15, colour = "black") +
  facet_wrap(~ lineage, nrow = 1, scales = "free_y") +
  scale_fill_gradient2(low = "#3E73B9", mid = "#F7F7F7", high = "#B4473F", midpoint = 0, name = "NES") +
  scale_shape_manual(values = c("FDR < 0.05" = 16, "Not significant" = 1), name = NULL) +
  compact_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right") +
  labs(title = "Fixed pathway enrichment after refined leave-one-out refits", x = "Left-out patient", y = NULL)

direction_plot <- direction_df |>
  filter(leave_out != "NONE") |>
  group_by(lineage, leave_out) |>
  summarise(direction_concordance = mean(same_direction, na.rm = TRUE), .groups = "drop") |>
  mutate(
    leave_out = factor(leave_out, levels = c("P18", "P23", "P24", "P27", "P29", "P32")),
    influence = ifelse(leave_out == "P32", "P32 removed", "Other patient removed")
  )

p_d <- ggplot(direction_plot, aes(leave_out, direction_concordance, fill = influence)) +
  geom_col(width = 0.68) +
  facet_wrap(~ lineage, nrow = 1) +
  scale_fill_manual(values = c("Other patient removed" = "#858D96", "P32 removed" = "#B4473F"), guide = "none") +
  scale_y_continuous(limits = c(0, 1), labels = percent_format(accuracy = 1)) +
  compact_theme +
  labs(title = "Direction concordance of refined full-model top 25", x = "Left-out patient", y = "Same direction")

figure <- (p_a | p_b) / p_c / p_d +
  plot_layout(heights = c(0.9, 1.25, 0.85)) +
  plot_annotation(tag_levels = "a")

stem <- file.path(figure_dir, "ExtendedData4_refined_leave_one_out_sensitivity")
ggsave(paste0(stem, ".png"), figure, width = 7.2, height = 6.5, dpi = 300)
ggsave(paste0(stem, ".pdf"), figure, width = 7.2, height = 6.5)
ggsave(paste0(stem, ".svg"), figure, width = 7.2, height = 6.5)

write.csv(summary_df, file.path(source_dir, "ExtendedData4_refined_loo_stability.csv"), row.names = FALSE)
write.csv(pathway_df, file.path(source_dir, "ExtendedData4_refined_loo_pathways.csv"), row.names = FALSE)
write.csv(direction_df, file.path(source_dir, "ExtendedData4_refined_loo_top25_direction.csv"), row.names = FALSE)
write.csv(direction_plot, file.path(source_dir, "ExtendedData4_refined_loo_direction_summary.csv"), row.names = FALSE)

p32 <- summary_df |> filter(leave_out == "P32", status == "OK")
report <- c(
  "# Refined-lineage leave-one-patient diagnostics",
  "",
  paste0("- Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "- Input: author-resolved refined T-cell and myeloid pseudobulk objects.",
  "- Model: `~ patient_id + post + post:resp_num`.",
  "- GSEA: MSigDB 2025.1.Hs Hallmark, deterministic seeds, serial execution.",
  "",
  "## P32 removal",
  "",
  paste(capture.output(print(p32)), collapse = "\n"),
  "",
  "## Interpretation",
  "",
  "- Leave-one-patient diagnostics are descriptive influence analyses in a six-patient discovery cohort.",
  "- P32 removal eliminates the only High response level; the remaining Low-versus-Medium gradient is not equivalent to the full ordinal contrast.",
  "- Genome-wide rank, top-gene overlap, gene-direction and pathway-level results must be reported separately."
)
writeLines(report, file.path(out_dir, "REFINED_LEAVE_ONE_OUT_REPORT.md"))
message("Refined leave-one-patient diagnostics complete: ", out_dir)
