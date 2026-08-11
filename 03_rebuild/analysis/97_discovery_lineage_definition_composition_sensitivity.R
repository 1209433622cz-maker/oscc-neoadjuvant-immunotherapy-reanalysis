#!/usr/bin/env Rscript

rm(list = ls())
gc()
set.seed(1234)

workspace <- Sys.getenv("GSE200996_WORKSPACE")
if (workspace == "") workspace <- normalizePath(".", winslash = "/", mustWork = TRUE)
rebuild <- file.path(workspace, "03_rebuild")
metadata_dir <- file.path(workspace, "00_raw_data", "GSE200996_metadata")
object_path <- file.path(rebuild, "obj_full_QC_logUMAP.rds")
out_dir <- file.path(rebuild, "results", "discovery_lineage_composition_sensitivity")
figure_dir <- file.path(rebuild, "figures", "submission")
source_dir <- file.path(figure_dir, "source_data")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(DESeq2)
  library(fgsea)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(stringr)
})
if (requireNamespace("BiocParallel", quietly = TRUE)) {
  BiocParallel::register(BiocParallel::SerialParam())
}

strict_patients <- c("P18", "P23", "P24", "P27", "P29", "P32")
response_map <- c(P18 = 2, P23 = 2, P24 = 1, P27 = 1, P29 = 1, P32 = 3)
response_label_map <- c(P18 = "Medium", P23 = "Medium", P24 = "Low",
                        P27 = "Low", P29 = "Low", P32 = "High")
cohort_map <- c(P18 = "Combo", P23 = "Mono", P24 = "Combo",
                P27 = "Mono", P29 = "Mono", P32 = "Combo")

read_meta <- function(filename) {
  x <- read.delim(
    gzfile(file.path(metadata_dir, filename)),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  names(x)[1] <- "cell_id"
  x$timepoint <- ifelse(x$Stage == "Pre-Tx", "pre", "post")
  x
}

cd45 <- read_meta("GSE200996_CD45.tumor.single.cell.meta.data.txt.gz")
cd4 <- read_meta("GSE200996_CD4.tumor.single.cell.meta.data.txt.gz")
cd8 <- read_meta("GSE200996_CD8.tumor.single.cell.meta.data.txt.gz")
myeloid <- read_meta("GSE200996_Myeloid.tumor.single.cell.meta.data.txt.gz")

cd45_type <- setNames(cd45$CellType_ID, cd45$cell_id)
lineage_membership <- bind_rows(
  data.frame(
    refined_lineage = "T cell",
    refined_source = "CD4",
    broad_label = unname(cd45_type[cd4$cell_id]),
    stringsAsFactors = FALSE
  ),
  data.frame(
    refined_lineage = "T cell",
    refined_source = "CD8",
    broad_label = unname(cd45_type[cd8$cell_id]),
    stringsAsFactors = FALSE
  ),
  data.frame(
    refined_lineage = "Myeloid",
    refined_source = "Myeloid",
    broad_label = unname(cd45_type[myeloid$cell_id]),
    stringsAsFactors = FALSE
  )
) |>
  mutate(broad_label = ifelse(is.na(broad_label), "Not in CD45 metadata", broad_label)) |>
  count(refined_lineage, refined_source, broad_label, name = "n_cells")

broad_unresolved <- bind_rows(
  data.frame(
    refined_lineage = "T cell",
    refined_source = "Broad-only",
    broad_label = "T cell not in CD4/CD8 refined metadata",
    n_cells = sum(
      cd45$CellType_ID == "T cell" &
        !cd45$cell_id %in% union(cd4$cell_id, cd8$cell_id)
    )
  ),
  data.frame(
    refined_lineage = "Myeloid",
    refined_source = "Broad-only",
    broad_label = "Myeloid not in refined metadata",
    n_cells = sum(
      cd45$CellType_ID == "Myeloid" &
        !cd45$cell_id %in% myeloid$cell_id
    )
  )
)
lineage_membership <- bind_rows(lineage_membership, broad_unresolved)
write.csv(
  lineage_membership,
  file.path(out_dir, "LINEAGE_DEFINITION_OVERLAP_AUDIT.csv"),
  row.names = FALSE
)

t_refined_meta <- bind_rows(
  cd4 |>
    transmute(
      cell_id, patient_id = Patient_ID, timepoint,
      cluster = paste0("CD4_C", CellType_ID)
    ),
  cd8 |>
    transmute(
      cell_id, patient_id = Patient_ID, timepoint,
      cluster = paste0("CD8_C", CellType_ID)
    )
)
m_refined_meta <- myeloid |>
  transmute(
    cell_id, patient_id = Patient_ID, timepoint,
    cluster = paste0("M_C", CellType_ID)
  )

cluster_labels <- c(
  CD4_C1 = "CD4 C1", CD4_C2 = "CD4 C2", CD4_C3 = "CD4 C3",
  CD4_C4 = "CD4 C4 (CXCL13+ effector)",
  CD4_C5 = "CD4 C5 (cycling)",
  CD4_C6 = "CD4 C6", CD4_C7 = "CD4 C7",
  CD8_C1 = "CD8 C1 (GZMK+)", CD8_C2 = "CD8 C2 (GZMK+)",
  CD8_C3 = "CD8 C3 (ITGAE+ Trm)",
  CD8_C4 = "CD8 C4 (cycling ITGAE+ Trm)",
  CD8_C5 = "CD8 C5 (IL7R+ root)",
  CD8_C6 = "CD8 C6", CD8_C7 = "CD8 C7",
  M_C1 = "Myeloid C1 (CD14+ monocyte)",
  M_C2 = "Myeloid C2", M_C3 = "Myeloid C3", M_C4 = "Myeloid C4",
  M_C5 = "Myeloid C5", M_C6 = "Myeloid C6", M_C7 = "Myeloid C7",
  M_C8 = "Myeloid C8 (LAMP3+ mregDC)", M_C9 = "Myeloid C9"
)

make_sample_metadata <- function() {
  expand.grid(
    patient_id = strict_patients,
    timepoint = c("pre", "post"),
    stringsAsFactors = FALSE
  ) |>
    mutate(
      pb_id = paste(patient_id, timepoint, sep = "|"),
      response_ord = unname(response_label_map[patient_id]),
      cohort = unname(cohort_map[patient_id]),
      resp_num = unname(response_map[patient_id]),
      post = as.integer(timepoint == "post")
    )
}

make_refined_pb <- function(obj, cell_meta, label) {
  cell_meta <- cell_meta |>
    filter(patient_id %in% strict_patients, timepoint %in% c("pre", "post")) |>
    distinct(cell_id, .keep_all = TRUE)
  metadata_cells <- nrow(cell_meta)
  pb_levels <- make_sample_metadata()$pb_id
  count_layers <- Layers(obj[["RNA"]], search = "^counts")
  if (!length(count_layers)) stop("No RNA count layers found")
  aggregated <- NULL
  matched_cells <- character()
  sample_cells <- setNames(integer(length(pb_levels)), pb_levels)
  for (layer_name in count_layers) {
    layer_counts <- LayerData(obj[["RNA"]], layer = layer_name)
    layer_cells <- intersect(cell_meta$cell_id, colnames(layer_counts))
    if (!length(layer_cells)) {
      rm(layer_counts)
      next
    }
    layer_meta <- cell_meta[match(layer_cells, cell_meta$cell_id), , drop = FALSE]
    pb_id <- paste(layer_meta$patient_id, layer_meta$timepoint, sep = "|")
    group <- factor(pb_id, levels = pb_levels)
    group_matrix <- sparse.model.matrix(~ 0 + group)
    colnames(group_matrix) <- sub("^group", "", colnames(group_matrix))
    partial <- layer_counts[, layer_cells, drop = FALSE] %*% group_matrix
    if (is.null(aggregated)) {
      aggregated <- partial
    } else {
      if (!identical(rownames(aggregated), rownames(partial))) {
        partial <- partial[match(rownames(aggregated), rownames(partial)), , drop = FALSE]
      }
      aggregated <- aggregated + partial
    }
    matched_cells <- c(matched_cells, layer_cells)
    sample_cells <- sample_cells + table(factor(pb_id, levels = pb_levels))
    rm(layer_counts, partial, group_matrix)
    gc()
  }
  if (anyDuplicated(matched_cells)) stop(label, " cells appeared in more than one count layer")
  if (length(matched_cells) < 0.98 * metadata_cells) {
    stop(
      label, " refined-cell match below 98%: ",
      length(matched_cells), "/", metadata_cells
    )
  }
  if (any(sample_cells < 30)) {
    stop(label, " has a strict-pair sample below the 30-cell gate")
  }
  aggregated <- aggregated[, pb_levels, drop = FALSE]
  md <- make_sample_metadata()
  rownames(md) <- md$pb_id
  md$patient_id <- factor(md$patient_id, levels = strict_patients)
  md$timepoint <- factor(md$timepoint, levels = c("pre", "post"))
  md$cohort <- factor(md$cohort, levels = c("Mono", "Combo"))
  list(
    counts = aggregated,
    meta = md,
    cell_summary = data.frame(
      lineage = label,
      metadata_cells = metadata_cells,
      object_matched_cells = length(matched_cells),
      minimum_sample_cells = min(sample_cells),
      maximum_sample_cells = max(sample_cells)
    )
  )
}

pb_paths <- c(
  "T cell" = file.path(out_dir, "T_cell_REFINED_PSEUDOBULK_INPUT.rds"),
  "Myeloid" = file.path(out_dir, "Myeloid_REFINED_PSEUDOBULK_INPUT.rds")
)
if (any(!file.exists(pb_paths))) {
  message("Loading discovery Seurat object once for refined-lineage aggregation")
  obj <- readRDS(object_path)
  if (!file.exists(pb_paths[["T cell"]])) {
    saveRDS(make_refined_pb(obj, t_refined_meta, "T cell"), pb_paths[["T cell"]])
  }
  if (!file.exists(pb_paths[["Myeloid"]])) {
    saveRDS(make_refined_pb(obj, m_refined_meta, "Myeloid"), pb_paths[["Myeloid"]])
  }
  rm(obj)
  gc()
  message("Released discovery Seurat object")
}

pb_cell_summary <- bind_rows(lapply(pb_paths, function(path) readRDS(path)$cell_summary))
write.csv(
  pb_cell_summary,
  file.path(out_dir, "REFINED_PSEUDOBULK_CELL_AUDIT.csv"),
  row.names = FALSE
)

make_composition <- function(cell_meta, lineage) {
  expected_samples <- make_sample_metadata()
  expected_clusters <- sort(unique(cell_meta$cluster))
  counts <- cell_meta |>
    filter(patient_id %in% strict_patients, timepoint %in% c("pre", "post")) |>
    count(patient_id, timepoint, cluster, name = "n_cells") |>
    complete(
      patient_id = strict_patients,
      timepoint = c("pre", "post"),
      cluster = expected_clusters,
      fill = list(n_cells = 0L)
    ) |>
    group_by(patient_id, timepoint) |>
    mutate(
      lineage = lineage,
      total_cells = sum(n_cells),
      fraction = n_cells / total_cells,
      adjusted_fraction = (n_cells + 0.5) / (total_cells + 0.5 * n_distinct(cluster)),
      logit_fraction = qlogis(adjusted_fraction),
      log_count = log(n_cells + 0.5),
      clr = log_count - mean(log_count)
    ) |>
    ungroup() |>
    mutate(
      pb_id = paste(patient_id, timepoint, sep = "|"),
      cluster_label = unname(cluster_labels[cluster])
    )
  if (any(counts$total_cells < 30)) stop(lineage, " composition sample below 30 cells")

  delta_wide <- counts |>
    select(patient_id, timepoint, cluster, clr) |>
    pivot_wider(names_from = timepoint, values_from = clr) |>
    mutate(delta_clr = post - pre) |>
    select(patient_id, cluster, delta_clr) |>
    pivot_wider(names_from = cluster, values_from = delta_clr) |>
    arrange(match(patient_id, strict_patients))
  delta_matrix <- as.matrix(delta_wide[, -1, drop = FALSE])
  rownames(delta_matrix) <- delta_wide$patient_id
  pca <- prcomp(delta_matrix, center = TRUE, scale. = FALSE)
  for (component in seq_len(ncol(pca$rotation))) {
    anchor <- which.max(abs(pca$rotation[, component]))
    if (pca$rotation[anchor, component] < 0) {
      pca$rotation[, component] <- -pca$rotation[, component]
      pca$x[, component] <- -pca$x[, component]
    }
  }
  variance <- pca$sdev^2 / sum(pca$sdev^2)
  patient_scores <- data.frame(
    lineage = lineage,
    patient_id = rownames(pca$x),
    comp_pc1 = pca$x[, 1],
    comp_pc2 = if (ncol(pca$x) >= 2) pca$x[, 2] else NA_real_,
    pc1_variance = variance[1],
    pc2_variance = if (length(variance) >= 2) variance[2] else NA_real_,
    response_ord = unname(response_label_map[rownames(pca$x)]),
    resp_num = unname(response_map[rownames(pca$x)]),
    cohort = unname(cohort_map[rownames(pca$x)]),
    stringsAsFactors = FALSE
  )
  loadings <- data.frame(
    lineage = lineage,
    cluster = rownames(pca$rotation),
    cluster_label = unname(cluster_labels[rownames(pca$rotation)]),
    pc1_loading = pca$rotation[, 1],
    pc2_loading = if (ncol(pca$rotation) >= 2) pca$rotation[, 2] else NA_real_,
    stringsAsFactors = FALSE
  )
  list(counts = counts, scores = patient_scores, loadings = loadings)
}

composition <- list(
  "T cell" = make_composition(t_refined_meta, "T cell"),
  "Myeloid" = make_composition(m_refined_meta, "Myeloid")
)
composition_counts <- bind_rows(lapply(composition, `[[`, "counts"))
composition_scores <- bind_rows(lapply(composition, `[[`, "scores"))
composition_loadings <- bind_rows(lapply(composition, `[[`, "loadings"))
write.csv(composition_counts, file.path(out_dir, "REFINED_CLUSTER_COMPOSITION.csv"), row.names = FALSE)
write.csv(composition_scores, file.path(out_dir, "COMPOSITION_PC_SCORES.csv"), row.names = FALSE)
write.csv(composition_loadings, file.path(out_dir, "COMPOSITION_PC_LOADINGS.csv"), row.names = FALSE)

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

all_assignments <- function(patient_info) {
  permutations <- unique_permutations(patient_info$resp_num)
  lapply(seq_len(nrow(permutations)), function(index) {
    setNames(permutations[index, ], patient_info$patient_id)
  })
}

stratified_assignments <- function(patient_info) {
  mono <- patient_info$patient_id[patient_info$cohort == "Mono"]
  combo <- patient_info$patient_id[patient_info$cohort == "Combo"]
  mono_values <- patient_info$resp_num[match(mono, patient_info$patient_id)]
  combo_values <- patient_info$resp_num[match(combo, patient_info$patient_id)]
  mono_perm <- unique_permutations(mono_values)
  combo_perm <- unique_permutations(combo_values)
  output <- list()
  index <- 1L
  for (i in seq_len(nrow(mono_perm))) {
    for (j in seq_len(nrow(combo_perm))) {
      values <- setNames(rep(NA_real_, nrow(patient_info)), patient_info$patient_id)
      values[mono] <- mono_perm[i, ]
      values[combo] <- combo_perm[j, ]
      output[[index]] <- values
      index <- index + 1L
    }
  }
  output
}

response_coefficient <- function(values, cohort, response, composition = NULL) {
  if (is.null(composition)) {
    fit <- lm(values ~ cohort + response)
  } else {
    fit <- lm(values ~ cohort + composition + response)
  }
  unname(coef(fit)[["response"]])
}

exact_test <- function(values, patient_info, assignments, composition = NULL) {
  observed <- response_coefficient(
    values,
    patient_info$cohort,
    patient_info$resp_num,
    composition
  )
  null <- vapply(assignments, function(assignment) {
    response_coefficient(
      values,
      patient_info$cohort,
      as.numeric(assignment[patient_info$patient_id]),
      composition
    )
  }, numeric(1))
  c(
    response_coefficient = observed,
    exact_p_two_sided = mean(abs(null) >= abs(observed) - 1e-12),
    n_assignments = length(null)
  )
}

cluster_delta <- composition_counts |>
  select(lineage, patient_id, timepoint, cluster, cluster_label, logit_fraction) |>
  pivot_wider(names_from = timepoint, values_from = logit_fraction) |>
  mutate(
    delta_logit_fraction = post - pre,
    cohort = unname(cohort_map[patient_id]),
    resp_num = unname(response_map[patient_id])
  )

cluster_tests <- bind_rows(lapply(
  split(cluster_delta, interaction(cluster_delta$lineage, cluster_delta$cluster, drop = TRUE)),
  function(x) {
    x <- x[match(strict_patients, x$patient_id), , drop = FALSE]
    info <- x[, c("patient_id", "cohort", "resp_num")]
    assignments <- stratified_assignments(info)
    result <- exact_test(x$delta_logit_fraction, info, assignments)
    data.frame(
      lineage = x$lineage[1],
      cluster = x$cluster[1],
      cluster_label = x$cluster_label[1],
      response_coefficient = result[["response_coefficient"]],
      exact_p_two_sided = result[["exact_p_two_sided"]],
      n_assignments = result[["n_assignments"]]
    )
  }
)) |>
  group_by(lineage) |>
  mutate(exact_bh_fdr = p.adjust(exact_p_two_sided, method = "BH")) |>
  ungroup()
write.csv(cluster_tests, file.path(out_dir, "CLUSTER_COMPOSITION_EXACT_TESTS.csv"), row.names = FALSE)

pc_tests <- bind_rows(lapply(split(composition_scores, composition_scores$lineage), function(x) {
  x <- x[match(strict_patients, x$patient_id), , drop = FALSE]
  info <- x[, c("patient_id", "cohort", "resp_num")]
  result <- exact_test(x$comp_pc1, info, stratified_assignments(info))
  data.frame(
    lineage = x$lineage[1],
    component = "PC1",
    variance_explained = x$pc1_variance[1],
    response_coefficient = result[["response_coefficient"]],
    exact_p_two_sided = result[["exact_p_two_sided"]],
    n_assignments = result[["n_assignments"]]
  )
}))
write.csv(pc_tests, file.path(out_dir, "COMPOSITION_PC_EXACT_TESTS.csv"), row.names = FALSE)

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

model_formulas <- list(
  primary = ~ patient_id + post + post:resp_num,
  cohort_adjusted = ~ patient_id + post + post:cohort + post:resp_num,
  composition_adjusted = ~ patient_id + post + post:comp_pc1 + post:resp_num,
  cohort_composition_adjusted = ~ patient_id + post + post:cohort + post:comp_pc1 + post:resp_num
)

fit_model <- function(pb, formula) {
  counts <- round(pb$counts[rowSums(pb$counts) >= 10, , drop = FALSE])
  md <- pb$meta
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
    formula = paste(deparse(formula), collapse = ""),
    residual_df = nrow(md) - ncol(design_matrix)
  )
}

hallmark <- load_hallmark()
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

gene_outputs <- list()
gsea_outputs <- list()
fit_objects <- list()
model_index <- 0L
for (lineage in c("T cell", "Myeloid")) {
  pb <- readRDS(pb_paths[[lineage]])
  score <- composition_scores |>
    filter(lineage == !!lineage) |>
    select(patient_id, comp_pc1)
  pb$meta <- pb$meta |>
    left_join(score, by = "patient_id") |>
    mutate(comp_pc1 = as.numeric(scale(comp_pc1)))
  rownames(pb$meta) <- pb$meta$pb_id
  pb$counts <- pb$counts[, rownames(pb$meta), drop = FALSE]
  fit_objects[[lineage]] <- list()
  for (model_name in names(model_formulas)) {
    message("Fitting ", lineage, ": ", model_name)
    fit <- fit_model(pb, model_formulas[[model_name]])
    fit_objects[[lineage]][[model_name]] <- fit
    gene <- fit$result
    gene$lineage <- lineage
    gene$model <- model_name
    gene$formula <- fit$formula
    gene$residual_df <- fit$residual_df
    gene_outputs[[length(gene_outputs) + 1L]] <- gene
    write.csv(
      gene,
      file.path(
        out_dir,
        paste0(gsub(" ", "_", lineage), "_", model_name, "_GENE_RESULTS.csv")
      ),
      row.names = FALSE
    )
    model_index <- model_index + 1L
    fg <- safe_fgsea(gene, hallmark, seed = 1234L + model_index)
    fg$lineage <- lineage
    fg$model <- model_name
    gsea_outputs[[length(gsea_outputs) + 1L]] <- fg
    write.csv(
      fg,
      file.path(
        out_dir,
        paste0(gsub(" ", "_", lineage), "_", model_name, "_HALLMARK_GSEA.csv")
      ),
      row.names = FALSE
    )
  }
  rm(pb)
  gc()
}

genes <- bind_rows(gene_outputs)
gsea <- bind_rows(gsea_outputs)
key_gsea <- gsea |>
  filter(pathway %in% unlist(key_pathways)) |>
  select(lineage, model, pathway, NES, pval, padj, size, leadingEdge)
write.csv(key_gsea, file.path(out_dir, "KEY_PATHWAY_MODEL_COMPARISON.csv"), row.names = FALSE)

comparison_metrics <- bind_rows(lapply(c("T cell", "Myeloid"), function(lineage) {
  broad_path <- file.path(
    rebuild, "results", "sensitivity_cohort_adjusted_pseudobulk",
    paste0(gsub(" ", "_", lineage), "_unadjusted_GENE_RESULTS.csv")
  )
  broad <- read.csv(broad_path, stringsAsFactors = FALSE) |>
    select(gene, broad_stat = stat, broad_padj = padj)
  x <- genes |>
    filter(lineage == !!lineage) |>
    select(gene, model, stat, padj) |>
    pivot_wider(names_from = model, values_from = c(stat, padj))
  metric <- function(reference, target, label) {
    finite <- is.finite(reference$stat) & is.finite(target$stat)
    top_reference <- head(reference$gene[order(-abs(reference$stat))], 100)
    top_target <- head(target$gene[order(-abs(target$stat))], 100)
    selected <- is.finite(reference$padj) & reference$padj < 0.05 & finite
    data.frame(
      lineage = lineage,
      comparison = label,
      n_common_genes = sum(finite),
      spearman_stat = cor(reference$stat[finite], target$stat[finite], method = "spearman"),
      top100_overlap = length(intersect(top_reference, top_target)),
      primary_fdr05_genes = sum(is.finite(reference$padj) & reference$padj < 0.05),
      target_fdr05_genes = sum(is.finite(target$padj) & target$padj < 0.05),
      primary_fdr05_sign_concordance = if (sum(selected)) {
        mean(sign(reference$stat[selected]) == sign(target$stat[selected]))
      } else {
        NA_real_
      }
    )
  }
  refined_primary <- data.frame(
    gene = x$gene, stat = x$stat_primary, padj = x$padj_primary
  )
  comp <- data.frame(
    gene = x$gene, stat = x$stat_composition_adjusted,
    padj = x$padj_composition_adjusted
  )
  combined <- data.frame(
    gene = x$gene, stat = x$stat_cohort_composition_adjusted,
    padj = x$padj_cohort_composition_adjusted
  )
  broad_join <- inner_join(
    refined_primary,
    broad,
    by = "gene"
  )
  broad_reference <- data.frame(
    gene = broad_join$gene, stat = broad_join$broad_stat, padj = broad_join$broad_padj
  )
  refined_common <- data.frame(
    gene = broad_join$gene, stat = broad_join$stat, padj = broad_join$padj
  )
  bind_rows(
    metric(broad_reference, refined_common, "broad_definition_vs_refined_definition"),
    metric(refined_primary, comp, "refined_primary_vs_composition_adjusted"),
    metric(refined_primary, combined, "refined_primary_vs_cohort_composition_adjusted")
  )
}))
write.csv(
  comparison_metrics,
  file.path(out_dir, "GENE_STATISTIC_CONCORDANCE.csv"),
  row.names = FALSE
)

pathway_patient_deltas <- function(dds, pathways, hallmark, lineage) {
  transformed <- assay(varianceStabilizingTransformation(dds, blind = FALSE))
  md <- as.data.frame(colData(dds))
  md$pb_id <- rownames(md)
  md$patient_id <- as.character(md$patient_id)
  md$timepoint <- ifelse(md$post == 1, "post", "pre")
  bind_rows(lapply(pathways, function(pathway) {
    pathway_genes <- intersect(hallmark[[pathway]], rownames(transformed))
    scores <- colMeans(transformed[pathway_genes, , drop = FALSE])
    data.frame(
      pb_id = names(scores),
      pathway = pathway,
      score = as.numeric(scores),
      n_genes = length(pathway_genes),
      stringsAsFactors = FALSE
    ) |>
      left_join(
        md[, c("pb_id", "patient_id", "timepoint", "cohort", "resp_num", "comp_pc1")],
        by = "pb_id"
      ) |>
      select(patient_id, cohort, resp_num, comp_pc1, pathway, n_genes, timepoint, score) |>
      pivot_wider(names_from = timepoint, values_from = score) |>
      filter(is.finite(pre), is.finite(post)) |>
      mutate(lineage = lineage, delta = post - pre)
  }))
}

pathway_exact <- bind_rows(lapply(c("T cell", "Myeloid"), function(lineage) {
  delta <- pathway_patient_deltas(
    fit_objects[[lineage]][["primary"]]$dds,
    key_pathways[[lineage]],
    hallmark,
    lineage
  )
  bind_rows(lapply(split(delta, delta$pathway), function(x) {
    x <- x[match(strict_patients, x$patient_id), , drop = FALSE]
    info <- x[, c("patient_id", "cohort", "resp_num")]
    assignments <- stratified_assignments(info)
    base <- exact_test(x$delta, info, assignments)
    adjusted <- exact_test(x$delta, info, assignments, composition = x$comp_pc1)
    bind_rows(
      data.frame(
        lineage = lineage, pathway = x$pathway[1],
        model = "cohort_adjusted",
        response_coefficient = base[["response_coefficient"]],
        exact_p_two_sided = base[["exact_p_two_sided"]],
        n_assignments = base[["n_assignments"]]
      ),
      data.frame(
        lineage = lineage, pathway = x$pathway[1],
        model = "cohort_composition_adjusted",
        response_coefficient = adjusted[["response_coefficient"]],
        exact_p_two_sided = adjusted[["exact_p_two_sided"]],
        n_assignments = adjusted[["n_assignments"]]
      )
    )
  }))
})) |>
  group_by(lineage, model) |>
  mutate(exact_bh_fdr = p.adjust(exact_p_two_sided, method = "BH")) |>
  ungroup()
write.csv(
  pathway_exact,
  file.path(out_dir, "KEY_PATHWAY_COMPOSITION_ADJUSTED_EXACT_TESTS.csv"),
  row.names = FALSE
)

write.csv(
  composition_counts,
  file.path(source_dir, "ExtendedData14_refined_cluster_composition.csv"),
  row.names = FALSE
)
write.csv(
  comparison_metrics,
  file.path(source_dir, "ExtendedData14_gene_statistic_concordance.csv"),
  row.names = FALSE
)
write.csv(
  key_gsea,
  file.path(source_dir, "ExtendedData14_key_pathway_models.csv"),
  row.names = FALSE
)
write.csv(
  pc_tests,
  file.path(source_dir, "ExtendedData14_composition_pc_exact_tests.csv"),
  row.names = FALSE
)

clean_pathway <- function(x) {
  str_to_sentence(str_replace_all(str_remove(x, "^HALLMARK_"), "_", " ")) |>
    str_replace_all(c("Tnfa" = "TNFA", "nfkb" = "NF-kB", "Mtorc1" = "mTORC1"))
}
compact_theme <- theme_classic(base_size = 7) +
  theme(
    plot.title = element_text(face = "bold", size = 9),
    plot.tag = element_text(face = "bold", size = 8),
    legend.title = element_blank(),
    legend.key.height = grid::unit(3.5, "mm")
  )
palette <- c("Low" = "#3E73B9", "Medium" = "#E5A33D", "High" = "#B44E52")

p_a <- composition_scores |>
  mutate(response_ord = factor(response_ord, levels = c("Low", "Medium", "High"))) |>
  ggplot(aes(response_ord, comp_pc1, colour = response_ord, shape = cohort)) +
  geom_hline(yintercept = 0, colour = "grey75", linewidth = 0.3) +
  geom_point(size = 2.2, position = position_jitter(width = 0.08, height = 0)) +
  facet_wrap(~ lineage, scales = "free_y") +
  scale_colour_manual(values = palette, guide = "none") +
  compact_theme +
  labs(
    title = "Response-blind subcluster composition PC1",
    x = "Pathological response",
    y = "Paired-delta PC1 score"
  )

gene_plot <- genes |>
  filter(model %in% c("primary", "composition_adjusted")) |>
  select(lineage, model, gene, stat) |>
  pivot_wider(names_from = model, values_from = stat) |>
  filter(is.finite(primary), is.finite(composition_adjusted)) |>
  left_join(
    comparison_metrics |>
      filter(comparison == "refined_primary_vs_composition_adjusted") |>
      select(lineage, spearman_stat),
    by = "lineage"
  )
p_b <- ggplot(gene_plot, aes(primary, composition_adjusted)) +
  geom_point(size = 0.28, alpha = 0.16, colour = "#4B6F8A") +
  geom_abline(slope = 1, intercept = 0, colour = "grey45", linewidth = 0.35) +
  geom_text(
    data = gene_plot |> distinct(lineage, spearman_stat),
    aes(
      x = -Inf, y = Inf,
      label = sprintf("Spearman rho = %.3f", spearman_stat)
    ),
    inherit.aes = FALSE,
    hjust = -0.08, vjust = 1.25, size = 2.3
  ) +
  facet_wrap(~ lineage, scales = "free") +
  compact_theme +
  labs(
    title = "Gene-statistic stability after composition adjustment",
    x = "Refined-lineage primary Wald statistic",
    y = "Composition-adjusted Wald statistic"
  )

key_plot <- key_gsea |>
  filter(model %in% c("primary", "composition_adjusted")) |>
  select(lineage, model, pathway, NES) |>
  pivot_wider(names_from = model, values_from = NES) |>
  mutate(pathway_label = clean_pathway(pathway))
p_c <- ggplot(key_plot, aes(primary, composition_adjusted, colour = lineage)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey45", linewidth = 0.35) +
  geom_point(size = 1.8) +
  scale_colour_manual(values = c("T cell" = "#3E73B9", "Myeloid" = "#B85A4A")) +
  compact_theme +
  labs(
    title = "Key-pathway enrichment stability",
    x = "Refined-lineage primary NES",
    y = "Composition-adjusted NES"
  )

loading_plot <- composition_loadings |>
  group_by(lineage) |>
  slice_max(abs(pc1_loading), n = 6, with_ties = FALSE) |>
  ungroup() |>
  mutate(cluster_label = str_wrap(cluster_label, width = 25)) |>
  ggplot(aes(pc1_loading, reorder(cluster_label, pc1_loading), fill = lineage)) +
  geom_col(width = 0.72) +
  facet_wrap(~ lineage, scales = "free_y") +
  scale_fill_manual(values = c("T cell" = "#3E73B9", "Myeloid" = "#B85A4A"), guide = "none") +
  compact_theme +
  theme(axis.text.y = element_text(size = 6.5)) +
  labs(title = "Largest composition-PC1 loadings", x = "PC1 loading", y = NULL)

figure <- (p_a | p_b) / (p_c | loading_plot) +
  plot_annotation(tag_levels = "a") +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
stem <- file.path(figure_dir, "ExtendedData14_submission_discovery_lineage_composition_sensitivity")
ggsave(paste0(stem, ".png"), figure, width = 7.2, height = 6.0, dpi = 600)
ggsave(paste0(stem, ".pdf"), figure, width = 7.2, height = 6.0)
ggsave(paste0(stem, ".svg"), figure, width = 7.2, height = 6.0)

summary_lines <- c(
  "# Discovery Refined-Lineage and Composition Sensitivity",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Frozen design",
  "",
  "- T-cell refined definition: the union of author-resolved CD4 and CD8 tumor-cell metadata.",
  "- Myeloid refined definition: author-resolved tumor-myeloid metadata.",
  "- Composition: author cluster fractions, 0.5-cell pseudocount, CLR transform and paired post-minus-pre delta.",
  "- Composition PC1 was derived without response labels and oriented by its largest absolute loading.",
  "- Gene models: primary, treatment-cohort adjusted, composition-PC1 adjusted and jointly adjusted.",
  "- Exact composition and pathway tests enumerate all 18 response assignments within treatment cohorts.",
  "- Joint models have two residual degrees of freedom and are interpreted as stress tests.",
  "",
  "## Refined pseudobulk cell audit",
  "",
  capture.output(print(pb_cell_summary, row.names = FALSE)),
  "",
  "## Composition PC tests",
  "",
  capture.output(print(pc_tests, row.names = FALSE)),
  "",
  "## Gene-statistic concordance",
  "",
  capture.output(print(comparison_metrics, row.names = FALSE)),
  "",
  "## Cluster composition tests",
  "",
  capture.output(print(cluster_tests[order(cluster_tests$exact_p_two_sided), ], row.names = FALSE)),
  "",
  "## Key pathway exact tests",
  "",
  capture.output(print(pathway_exact, row.names = FALSE)),
  "",
  "Primary claims should use the refined lineage only if its model is full-rank and its gene/pathway direction is not materially dependent on composition PC1. Broad-lineage results are retained solely as annotation-definition sensitivity."
)
writeLines(
  summary_lines,
  file.path(out_dir, "DISCOVERY_LINEAGE_COMPOSITION_SENSITIVITY_REPORT.md")
)
message("Completed refined-lineage composition sensitivity: ", out_dir)
