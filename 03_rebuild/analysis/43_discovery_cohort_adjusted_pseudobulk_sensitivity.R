#!/usr/bin/env Rscript

rm(list = ls())
gc()
set.seed(1234)

workspace <- Sys.getenv("GSE200996_WORKSPACE")
if (workspace == "") workspace <- normalizePath(".", winslash = "/", mustWork = TRUE)
rebuild <- file.path(workspace, "03_rebuild")
object_path <- file.path(rebuild, "obj_full_QC_logUMAP.rds")
out_dir <- file.path(rebuild, "results", "sensitivity_cohort_adjusted_pseudobulk")
figure_dir <- file.path(rebuild, "figures", "submission")
source_dir <- file.path(figure_dir, "source_data")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(DESeq2)
  library(fgsea)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(stringr)
})
if (requireNamespace("BiocParallel", quietly = TRUE)) {
  BiocParallel::register(BiocParallel::SerialParam())
}

load_hallmark <- function() {
  cache <- file.path(
    Sys.getenv("LOCALAPPDATA"), "R", "cache", "R", "msigdbr",
    "msigdb.2025.1.Hs.rds"
  )
  if (!file.exists(cache)) stop("Hallmark cache not found: ", cache)
  x <- readRDS(cache)
  x <- x[
    x$gs_collection == "H" &
      x$db_target_species == "HS" &
      !is.na(x$db_gene_symbol),
  ]
  split(x$db_gene_symbol, x$gs_name)
}

safe_fgsea <- function(result, hallmark, seed) {
  set.seed(seed)
  ranks <- result$stat
  names(ranks) <- result$gene
  ranks <- sort(ranks[is.finite(ranks)], decreasing = TRUE)
  fg <- as.data.frame(
    fgseaMultilevel(
      pathways = hallmark,
      stats = ranks,
      minSize = 10,
      maxSize = 500,
      nproc = 1,
      BPPARAM = BiocParallel::SerialParam(RNGseed = seed)
    )
  )
  if ("leadingEdge" %in% names(fg)) {
    fg$leadingEdge <- vapply(fg$leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))
  }
  fg[order(fg$padj, -abs(fg$NES)), ]
}

make_pb <- function(obj, celltype, min_cells = 30) {
  meta <- obj@meta.data
  cells <- rownames(meta)[
    meta$cd45_celltype_id == celltype &
      meta$timepoint %in% c("pre", "post") &
      !is.na(meta$response_ord)
  ]
  x <- subset(obj, cells = cells)
  x$pb_id <- paste(x$patient_id, x$timepoint, sep = "|")
  n_by_pb <- table(x$pb_id)
  keep <- names(n_by_pb)[n_by_pb >= min_cells]
  x <- subset(x, cells = Cells(x)[x$pb_id %in% keep])
  counts <- round(
    AggregateExpression(
      x, assays = "RNA", slot = "counts", group.by = "pb_id",
      return.seurat = FALSE
    )$RNA
  )
  md <- data.frame(pb_id = colnames(counts), stringsAsFactors = FALSE)
  split_id <- str_split_fixed(md$pb_id, "\\|", 2)
  md$patient_id <- split_id[, 1]
  md$timepoint <- factor(split_id[, 2], levels = c("pre", "post"))
  patient_meta <- obj@meta.data |>
    tibble::rownames_to_column("cell") |>
    filter(patient_id %in% md$patient_id) |>
    distinct(patient_id, response_ord, cohort)
  md <- left_join(md, patient_meta, by = "patient_id")
  md$resp_num <- as.numeric(factor(md$response_ord, levels = c("Low", "Medium", "High"), ordered = TRUE))
  md$post <- as.integer(md$timepoint == "post")
  md$patient_id <- factor(md$patient_id)
  md$cohort <- factor(md$cohort, levels = c("Mono", "Combo"))
  rownames(md) <- md$pb_id
  counts <- counts[, rownames(md), drop = FALSE]
  rm(x)
  gc()
  list(counts = counts, meta = md)
}

restrict_strict_pairs <- function(pb) {
  md <- pb$meta
  pair_status <- md |>
    mutate(patient_id_character = as.character(patient_id)) |>
    group_by(patient_id_character) |>
    summarise(
      has_pre = any(post == 0L),
      has_post = any(post == 1L),
      .groups = "drop"
    )
  paired_ids <- pair_status$patient_id_character[pair_status$has_pre & pair_status$has_post]
  keep <- as.character(md$patient_id) %in% paired_ids
  md <- droplevels(md[keep, , drop = FALSE])
  counts <- pb$counts[, rownames(md), drop = FALSE]
  list(counts = counts, meta = md)
}

fit_model <- function(pb, adjusted = FALSE) {
  counts <- pb$counts[rowSums(pb$counts) >= 10, , drop = FALSE]
  md <- pb$meta
  formula <- if (adjusted) {
    ~ patient_id + post + post:cohort + post:resp_num
  } else {
    ~ patient_id + post + post:resp_num
  }
  design_matrix <- model.matrix(formula, md)
  if (qr(design_matrix)$rank != ncol(design_matrix)) {
    stop("Rank-deficient design: ", deparse(formula))
  }
  dds <- DESeqDataSetFromMatrix(counts, md, design = formula)
  dds <- DESeq(dds, quiet = TRUE)
  response_name <- grep("post.*resp_num|resp_num.*post", resultsNames(dds), value = TRUE)[1]
  if (is.na(response_name)) stop("Response interaction coefficient not found")
  result <- as.data.frame(results(dds, name = response_name))
  result$gene <- rownames(result)
  result <- result[order(result$padj, result$pvalue), ]
  list(
    dds = dds,
    result = result,
    coefficient = response_name,
    formula = deparse(formula),
    residual_df = nrow(md) - ncol(design_matrix)
  )
}

unique_permutations <- function(values) {
  values <- as.numeric(values)
  levels <- sort(unique(values))
  counts <- table(factor(values, levels = levels))
  output <- vector("list", factorial(length(values)) / prod(factorial(as.numeric(counts))))
  index <- 1L
  build <- function(prefix, remaining) {
    if (length(prefix) == length(values)) {
      output[[index]] <<- prefix
      index <<- index + 1L
      return(invisible(NULL))
    }
    for (i in seq_along(levels)) {
      if (remaining[i] > 0L) {
        next_remaining <- remaining
        next_remaining[i] <- next_remaining[i] - 1L
        build(c(prefix, levels[i]), next_remaining)
      }
    }
    invisible(NULL)
  }
  build(numeric(0), as.integer(counts))
  do.call(rbind, output)
}

stratified_assignments <- function(patient_info) {
  mono <- patient_info$patient_id[patient_info$cohort == "Mono"]
  combo <- patient_info$patient_id[patient_info$cohort == "Combo"]
  mono_values <- patient_info$resp_num[match(mono, patient_info$patient_id)]
  combo_values <- patient_info$resp_num[match(combo, patient_info$patient_id)]
  mono_perm <- unique_permutations(mono_values)
  combo_perm <- unique_permutations(combo_values)
  output <- list()
  index <- 1
  for (i in seq_len(nrow(mono_perm))) {
    for (j in seq_len(nrow(combo_perm))) {
      values <- setNames(rep(NA_real_, nrow(patient_info)), patient_info$patient_id)
      values[mono] <- mono_perm[i, ]
      values[combo] <- combo_perm[j, ]
      output[[index]] <- values
      index <- index + 1
    }
  }
  output
}

pathway_patient_deltas <- function(dds, pathways, hallmark) {
  transformed <- assay(varianceStabilizingTransformation(dds, blind = FALSE))
  md <- as.data.frame(colData(dds))
  md$pb_id <- rownames(md)
  md$patient_id <- as.character(md$patient_id)
  md$timepoint <- ifelse(md$post == 1, "post", "pre")
  bind_rows(lapply(pathways, function(pathway) {
    genes <- intersect(hallmark[[pathway]], rownames(transformed))
    if (!length(genes)) stop("No genes available for pathway: ", pathway)
    scores <- colMeans(transformed[genes, , drop = FALSE])
    data.frame(
      pb_id = names(scores),
      pathway = pathway,
      score = as.numeric(scores),
      n_genes = length(genes),
      stringsAsFactors = FALSE
    ) |>
      left_join(md[, c("pb_id", "patient_id", "timepoint", "cohort", "resp_num")], by = "pb_id") |>
      select(patient_id, cohort, resp_num, pathway, n_genes, timepoint, score) |>
      tidyr::pivot_wider(names_from = timepoint, values_from = score) |>
      filter(is.finite(pre), is.finite(post)) |>
      mutate(delta = post - pre)
  }))
}

response_coef <- function(delta, cohort, response) {
  fit <- lm(delta ~ cohort + response)
  unname(coef(fit)[["response"]])
}

key_pathways <- list(
  "T cell" = c(
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    "HALLMARK_MTORC1_SIGNALING",
    "HALLMARK_P53_PATHWAY",
    "HALLMARK_INTERFERON_ALPHA_RESPONSE",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "HALLMARK_INFLAMMATORY_RESPONSE"
  ),
  "Myeloid" = c(
    "HALLMARK_MTORC1_SIGNALING",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    "HALLMARK_INTERFERON_ALPHA_RESPONSE",
    "HALLMARK_INFLAMMATORY_RESPONSE",
    "HALLMARK_COMPLEMENT"
  )
)

hallmark <- load_hallmark()
pb_paths <- setNames(
  file.path(out_dir, paste0(gsub(" ", "_", c("T cell", "Myeloid")), "_PSEUDOBULK_INPUT.rds")),
  c("T cell", "Myeloid")
)
need_object <- any(!file.exists(pb_paths))
if (need_object) {
  message("Loading cached discovery object: ", object_path)
  obj <- readRDS(object_path)
  missing_celltypes <- names(pb_paths)[!file.exists(pb_paths)]
  for (celltype in missing_celltypes) {
    message("Pre-aggregating ", celltype, " pseudobulk before model fitting")
    pb_missing <- make_pb(obj, celltype)
    saveRDS(pb_missing, pb_paths[[celltype]])
    rm(pb_missing)
    gc()
  }
  rm(obj)
  gc()
  message("Released the Seurat object after pseudobulk aggregation.")
} else {
  message("Both pseudobulk cache files exist; the Seurat object will not be loaded.")
}

gene_outputs <- list()
gsea_outputs <- list()
key_outputs <- list()
permutation_outputs <- list()
null_outputs <- list()

for (celltype in c("T cell", "Myeloid")) {
  if (!file.exists(pb_paths[[celltype]])) stop("Pseudobulk cache missing: ", pb_paths[[celltype]])
  message("Reusing cached ", celltype, " pseudobulk")
  pb_cached <- readRDS(pb_paths[[celltype]])
  pb <- restrict_strict_pairs(pb_cached)
  if (ncol(pb$counts) != ncol(pb_cached$counts)) {
    message(
      "Restricting ", celltype, " pseudobulk from ", ncol(pb_cached$counts),
      " samples to ", ncol(pb$counts), " strict-pair samples"
    )
    saveRDS(pb, pb_paths[[celltype]])
  }
  rm(pb_cached)
  write.csv(pb$meta, file.path(out_dir, paste0(gsub(" ", "_", celltype), "_PSEUDOBULK_METADATA.csv")), row.names = FALSE)

  unadjusted <- fit_model(pb, adjusted = FALSE)
  adjusted <- fit_model(pb, adjusted = TRUE)
  for (model_name in c("unadjusted", "cohort_adjusted")) {
    fit <- if (model_name == "unadjusted") unadjusted else adjusted
    gene <- fit$result
    gene$celltype <- celltype
    gene$model <- model_name
    gene$formula <- fit$formula
    gene$residual_df <- fit$residual_df
    gene_outputs[[length(gene_outputs) + 1]] <- gene
    write.csv(
      gene,
      file.path(out_dir, paste0(gsub(" ", "_", celltype), "_", model_name, "_GENE_RESULTS.csv")),
      row.names = FALSE
    )
    seed_offset <- if (celltype == "T cell") 100L else 200L
    seed_offset <- seed_offset + if (model_name == "unadjusted") 1L else 2L
    fg <- safe_fgsea(gene, hallmark, seed = 1234L + seed_offset)
    fg$celltype <- celltype
    fg$model <- model_name
    gsea_outputs[[length(gsea_outputs) + 1]] <- fg
    write.csv(
      fg,
      file.path(out_dir, paste0(gsub(" ", "_", celltype), "_", model_name, "_HALLMARK_GSEA.csv")),
      row.names = FALSE
    )
  }

  fg_u <- gsea_outputs[[length(gsea_outputs) - 1]]
  fg_a <- gsea_outputs[[length(gsea_outputs)]]
  key <- full_join(
    fg_u |> select(pathway, unadjusted_NES = NES, unadjusted_p = pval, unadjusted_FDR = padj),
    fg_a |> select(pathway, adjusted_NES = NES, adjusted_p = pval, adjusted_FDR = padj),
    by = "pathway"
  ) |>
    filter(pathway %in% key_pathways[[celltype]]) |>
    mutate(celltype = celltype)
  key_outputs[[length(key_outputs) + 1]] <- key

  pathway_deltas <- pathway_patient_deltas(
    adjusted$dds, key_pathways[[celltype]], hallmark
  )
  patient_info <- pathway_deltas |>
    distinct(patient_id, cohort, resp_num) |>
    mutate(patient_id = as.character(patient_id), cohort = as.character(cohort))
  assignments <- stratified_assignments(patient_info)
  observed_scores <- setNames(
    vapply(key_pathways[[celltype]], function(pathway) {
      x <- pathway_deltas[pathway_deltas$pathway == pathway, ]
      response_coef(x$delta, x$cohort, x$resp_num)
    }, numeric(1)),
    key_pathways[[celltype]]
  )
  null_matrix <- matrix(
    NA_real_, nrow = length(assignments), ncol = length(observed_scores),
    dimnames = list(NULL, names(observed_scores))
  )
  for (assignment_index in seq_along(assignments)) {
    map <- assignments[[assignment_index]]
    null_matrix[assignment_index, ] <- vapply(
      names(observed_scores),
      function(pathway) {
        x <- pathway_deltas[pathway_deltas$pathway == pathway, ]
        response_coef(
          x$delta,
          x$cohort,
          as.numeric(map[x$patient_id])
        )
      },
      numeric(1)
    )
  }
  permutation <- data.frame(
    celltype = celltype,
    pathway = names(observed_scores),
    observed_pathway_delta_response_coef = as.numeric(observed_scores),
    n_paired_patients = nrow(patient_info),
    n_unique_stratified_assignments = nrow(null_matrix),
    exact_p_two_sided = vapply(
      seq_along(observed_scores),
      function(column) mean(abs(null_matrix[, column]) >= abs(observed_scores[column]) - 1e-12),
      numeric(1)
    )
  )
  permutation$exact_bh_fdr <- p.adjust(permutation$exact_p_two_sided, method = "BH")
  permutation_outputs[[length(permutation_outputs) + 1]] <- permutation
  null_outputs[[length(null_outputs) + 1]] <- data.frame(
    celltype = celltype,
    assignment_id = rep(seq_len(nrow(null_matrix)), each = ncol(null_matrix)),
    pathway = rep(colnames(null_matrix), times = nrow(null_matrix)),
    permuted_pathway_delta_response_coef = as.vector(t(null_matrix))
  )
  rm(pb, unadjusted, adjusted)
  gc()
}

genes <- bind_rows(gene_outputs)
gsea <- bind_rows(gsea_outputs)
key <- bind_rows(key_outputs)
permutation <- bind_rows(permutation_outputs)
null <- bind_rows(null_outputs)
write.csv(key, file.path(out_dir, "KEY_PATHWAY_MODEL_COMPARISON.csv"), row.names = FALSE)
write.csv(permutation, file.path(out_dir, "KEY_PATHWAY_STRATIFIED_EXACT_PERMUTATION.csv"), row.names = FALSE)
write.csv(null, file.path(out_dir, "KEY_PATHWAY_STRATIFIED_NULL.csv"), row.names = FALSE)
write.csv(key, file.path(source_dir, "ExtendedData9_key_pathway_model_comparison.csv"), row.names = FALSE)
write.csv(permutation, file.path(source_dir, "ExtendedData9_key_pathway_exact_permutation.csv"), row.names = FALSE)

clean_pathway <- function(x) {
  str_to_sentence(str_replace_all(str_remove(x, "^HALLMARK_"), "_", " ")) |>
    str_replace_all(c("Tnfa" = "TNFA", "nfkb" = "NF-kB", "Mtorc1" = "mTORC1"))
}
gene_compare <- genes |>
  select(celltype, model, gene, stat) |>
  tidyr::pivot_wider(names_from = model, values_from = stat) |>
  filter(is.finite(unadjusted), is.finite(cohort_adjusted))
rank_summary <- gene_compare |>
  group_by(celltype) |>
  summarise(spearman = cor(unadjusted, cohort_adjusted, method = "spearman"), .groups = "drop")
rho_t <- rank_summary$spearman[rank_summary$celltype == "T cell"]
rho_m <- rank_summary$spearman[rank_summary$celltype == "Myeloid"]
compact_theme <- theme_classic(base_size = 7) +
  theme(
    plot.title = element_text(face = "bold", size = 9),
    plot.tag = element_text(face = "bold", size = 8),
    legend.title = element_blank(),
    legend.key.height = grid::unit(3.5, "mm")
  )
p_a <- ggplot(filter(gene_compare, celltype == "T cell"), aes(unadjusted, cohort_adjusted)) +
  geom_point(size = 0.35, alpha = 0.20, colour = "#3E73B9") +
  geom_abline(slope = 1, intercept = 0, colour = "grey45", linewidth = 0.35) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.15, vjust = 1.25, label = sprintf("Spearman rho = %.3f", rho_t), size = 2.4) +
  compact_theme + labs(title = "T-cell gene statistics", x = "Unadjusted Wald statistic", y = "Cohort-adjusted Wald statistic")
p_b <- ggplot(filter(gene_compare, celltype == "Myeloid"), aes(unadjusted, cohort_adjusted)) +
  geom_point(size = 0.35, alpha = 0.20, colour = "#B85A4A") +
  geom_abline(slope = 1, intercept = 0, colour = "grey45", linewidth = 0.35) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.15, vjust = 1.25, label = sprintf("Spearman rho = %.3f", rho_m), size = 2.4) +
  compact_theme + labs(title = "Myeloid gene statistics", x = "Unadjusted Wald statistic", y = "Cohort-adjusted Wald statistic")
p_c <- key |>
  mutate(pathway_label = clean_pathway(pathway)) |>
  ggplot(aes(unadjusted_NES, adjusted_NES, colour = celltype)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey45", linewidth = 0.35) +
  geom_point(size = 1.7) +
  scale_colour_manual(values = c("T cell" = "#3E73B9", "Myeloid" = "#B85A4A")) +
  compact_theme + labs(title = "Key-pathway enrichment", x = "Unadjusted NES", y = "Cohort-adjusted NES")
p_d <- permutation |>
  mutate(pathway_label = clean_pathway(pathway)) |>
  ggplot(aes(exact_p_two_sided, reorder(pathway_label, exact_p_two_sided), colour = celltype)) +
  geom_vline(xintercept = 0.05, linetype = "dashed", colour = "grey45", linewidth = 0.35) +
  geom_point(size = 1.7, position = position_dodge(width = 0.45)) +
  scale_colour_manual(values = c("T cell" = "#3E73B9", "Myeloid" = "#B85A4A"), guide = "none") +
  scale_x_continuous(limits = c(0, 0.82), breaks = c(0, 0.2, 0.4, 0.6, 0.8)) +
  compact_theme + theme(axis.text.y = element_text(size = 6.5)) +
  labs(title = "Within-cohort exact permutation", x = "Two-sided exact P", y = NULL)
figure <- (p_a | p_b) / (p_c | p_d) +
  plot_annotation(tag_levels = "a") +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
stem <- file.path(figure_dir, "ExtendedData9_submission_discovery_cohort_sensitivity")
ggsave(paste0(stem, ".png"), figure, width = 7.2, height = 5.6, dpi = 600)
ggsave(paste0(stem, ".pdf"), figure, width = 7.2, height = 5.6)
ggsave(paste0(stem, ".svg"), figure, width = 7.2, height = 5.6)
report <- c(
  "# Discovery Cohort-Adjusted Pseudobulk Sensitivity",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Design",
  "",
  "- Original model: `~ patient_id + post + post:resp_num`.",
  "- Sensitivity model: `~ patient_id + post + post:cohort + post:resp_num`.",
  "- Exact pathway sensitivity: all 18 unique response-label assignments within Mono/Combo strata.",
  "- The exact test uses a fixed pathway panel selected from the primary discovery GSEA and is a sensitivity analysis, not independent validation.",
  "",
  "## Gene-statistic concordance",
  "",
  capture.output(print(as.data.frame(rank_summary), row.names = FALSE)),
  "",
  "## Key pathways",
  "",
  capture.output(print(key, row.names = FALSE)),
  "",
  "## Stratified exact permutation",
  "",
  capture.output(print(permutation, row.names = FALSE)),
  "",
  "Interpret the primary discovery mechanism as cohort-robust only if key pathway direction and rank structure remain stable after the `post:cohort` term. Exact P values are resolution-limited by 18 assignments."
)
writeLines(report, file.path(out_dir, "DISCOVERY_COHORT_ADJUSTED_PSEUDOBULK_REPORT.md"))
message("Completed cohort-adjusted discovery sensitivity: ", out_dir)
