base_dir <- Sys.getenv("GSE200996_BASEDIR")
if (base_dir == "") {
  base_dir <- "H:/SCI2/OSCC-GSE200996-2025.12/03_rebuild"
}

result_dir <- file.path(base_dir, "results")
out_file <- file.path(result_dir, "REANALYSIS_RESULT_READOUT.md")

read_csv <- function(...) {
  read.csv(file.path(...), check.names = FALSE, stringsAsFactors = FALSE)
}

fmt_num <- function(x, digits = 3) {
  ifelse(is.na(x), "NA", formatC(x, digits = digits, format = "fg", flag = "#"))
}

fmt_p <- function(x) {
  ifelse(is.na(x), "NA", formatC(x, digits = 3, format = "e"))
}

sort_by_padj <- function(x) {
  if (!("padj" %in% names(x))) return(x)
  x[order(ifelse(is.na(x$padj), Inf, x$padj),
          ifelse(is.na(x$pvalue), Inf, x$pvalue)), , drop = FALSE]
}

sort_gsea <- function(x) {
  x[order(ifelse(is.na(x$padj), Inf, x$padj),
          ifelse(is.na(x$pval), Inf, x$pval)), , drop = FALSE]
}

md_table <- function(df) {
  if (nrow(df) == 0) return("_No rows._")
  df[] <- lapply(df, as.character)
  header <- paste0("| ", paste(names(df), collapse = " | "), " |")
  sep <- paste0("| ", paste(rep("---", ncol(df)), collapse = " | "), " |")
  rows <- apply(df, 1, function(z) paste0("| ", paste(z, collapse = " | "), " |"))
  paste(c(header, sep, rows), collapse = "\n")
}

top_de_table <- function(path, n = 10) {
  x <- sort_by_padj(read.csv(path, check.names = FALSE, stringsAsFactors = FALSE))
  cols <- intersect(c("gene", "log2FoldChange", "stat", "pvalue", "padj"), names(x))
  y <- head(x[, cols, drop = FALSE], n)
  if ("log2FoldChange" %in% names(y)) y$log2FoldChange <- fmt_num(y$log2FoldChange)
  if ("stat" %in% names(y)) y$stat <- fmt_num(y$stat)
  if ("pvalue" %in% names(y)) y$pvalue <- fmt_p(y$pvalue)
  if ("padj" %in% names(y)) y$padj <- fmt_p(y$padj)
  list(
    table = y,
    n_fdr_05 = sum(!is.na(x$padj) & x$padj < 0.05),
    n_fdr_10 = sum(!is.na(x$padj) & x$padj < 0.10)
  )
}

top_gsea_table <- function(path, n = 10) {
  x <- sort_gsea(read.csv(path, check.names = FALSE, stringsAsFactors = FALSE))
  y <- head(x[, c("pathway", "NES", "pval", "padj", "size"), drop = FALSE], n)
  y$NES <- fmt_num(y$NES)
  y$pval <- fmt_p(y$pval)
  y$padj <- fmt_p(y$padj)
  y
}

signature_slopes <- function(path) {
  x <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  out <- do.call(rbind, lapply(split(x, x$signature), function(d) {
    fit <- summary(lm(delta ~ resp_num, data = d))
    data.frame(
      signature = d$signature[1],
      slope_per_response_step = coef(fit)[2, 1],
      p_value = coef(fit)[2, 4],
      mean_delta = mean(d$delta, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  out <- out[order(out$p_value), , drop = FALSE]
  out$slope_per_response_step <- fmt_num(out$slope_per_response_step)
  out$p_value <- fmt_p(out$p_value)
  out$mean_delta <- fmt_num(out$mean_delta)
  out
}

main_story <- read_csv(result_dir, "main_story", "MainStory_Alignment_Table.csv")
baseline_comp <- read_csv(result_dir, "pre_baseline", "baseline_pre_composition_trend_respOrd_limma_logit.csv")
dyn_comp_ord <- read_csv(result_dir, "dynamic_paired", "Fig4A_composition_delta_logit_limma_respOrd_trend.csv")
dyn_comp <- read_csv(result_dir, "dynamic_paired", "Fig4A_composition_delta_logit_limma_RvsNR.csv")

baseline_comp_out <- baseline_comp[, c("celltype", "logFC", "P.Value", "adj.P.Val")]
baseline_comp_out$logFC <- fmt_num(baseline_comp_out$logFC)
baseline_comp_out$P.Value <- fmt_p(baseline_comp_out$P.Value)
baseline_comp_out$adj.P.Val <- fmt_p(baseline_comp_out$adj.P.Val)

dyn_comp_out <- dyn_comp[, c("celltype", "logFC", "P.Value", "adj.P.Val")]
dyn_comp_out$logFC <- fmt_num(dyn_comp_out$logFC)
dyn_comp_out$P.Value <- fmt_p(dyn_comp_out$P.Value)
dyn_comp_out$adj.P.Val <- fmt_p(dyn_comp_out$adj.P.Val)

dyn_comp_ord_out <- dyn_comp_ord[, c("celltype", "logFC", "P.Value", "adj.P.Val")]
dyn_comp_ord_out$logFC <- fmt_num(dyn_comp_ord_out$logFC)
dyn_comp_ord_out$P.Value <- fmt_p(dyn_comp_ord_out$P.Value)
dyn_comp_ord_out$adj.P.Val <- fmt_p(dyn_comp_ord_out$adj.P.Val)

t_de <- top_de_table(file.path(result_dir, "dynamic_paired", "Fig4B_T_cell_interaction_DE_trend.csv"))
m_de <- top_de_table(file.path(result_dir, "dynamic_paired", "Fig4B_Myeloid_interaction_DE_trend.csv"))

t_gsea <- top_gsea_table(file.path(result_dir, "dynamic_paired", "Fig4B_T_cell_GSEA_Hallmark.csv"))
m_gsea <- top_gsea_table(file.path(result_dir, "dynamic_paired", "Fig4B_Myeloid_GSEA_Hallmark.csv"))

t_sig <- signature_slopes(file.path(result_dir, "dynamic_paired", "Fig4B_T_cell_SignatureDelta_source.csv"))
m_sig <- signature_slopes(file.path(result_dir, "dynamic_paired", "Fig4B_Myeloid_SignatureDelta_source.csv"))

lines <- c(
  "# Reanalysis Result Readout",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Core Interpretation",
  "",
  "- Baseline immune composition trends are weak after patient-level modelling.",
  "- In the all-paired ordinal composition model, mast-cell post-pre change increases with response depth; myeloid change trends downward but is not FDR-significant.",
  "- In the binary paired sensitivity model, responder status is associated with lower myeloid post-pre delta and higher mast-cell delta.",
  "- The main manuscript should therefore emphasize treatment-induced immune-state remodeling, with myeloid and T-cell pathway activation as the mechanistic layer.",
  "- Claims should remain hypothesis-generating because the strict paired cohort has six patients and the binary R/NR paired subset is smaller.",
  "",
  "## Main Story Alignment Table",
  "",
  md_table(within(main_story, {
                baseline_trend_logFC <- fmt_num(baseline_trend_logFC)
                baseline_trend_P <- fmt_p(baseline_trend_P)
                baseline_trend_FDR <- fmt_p(baseline_trend_FDR)
                dyn_delta_ord_logFC <- fmt_num(dyn_delta_ord_logFC)
                dyn_delta_ord_P <- fmt_p(dyn_delta_ord_P)
                dyn_delta_ord_FDR <- fmt_p(dyn_delta_ord_FDR)
                dyn_delta_RvsNR_logFC <- fmt_num(dyn_delta_RvsNR_logFC)
                dyn_delta_RvsNR_P <- fmt_p(dyn_delta_RvsNR_P)
                dyn_delta_RvsNR_FDR <- fmt_p(dyn_delta_RvsNR_FDR)
  })),
  "",
  "## Composition Results",
  "",
  "Baseline ordinal-response trend:",
  "",
  md_table(baseline_comp_out),
  "",
  "Paired post-pre ordinal-response delta:",
  "",
  md_table(dyn_comp_ord_out),
  "",
  "Paired post-pre binary R/NR delta (sensitivity):",
  "",
  md_table(dyn_comp_out),
  "",
  "## Dynamic Pseudobulk Interaction: Top Genes",
  "",
  paste0("T cell: FDR<0.05 genes = ", t_de$n_fdr_05, "; FDR<0.10 genes = ", t_de$n_fdr_10),
  "",
  md_table(t_de$table),
  "",
  paste0("Myeloid: FDR<0.05 genes = ", m_de$n_fdr_05, "; FDR<0.10 genes = ", m_de$n_fdr_10),
  "",
  md_table(m_de$table),
  "",
  "## Dynamic Hallmark GSEA",
  "",
  "T cell:",
  "",
  md_table(t_gsea),
  "",
  "Myeloid:",
  "",
  md_table(m_gsea),
  "",
  "## Signature Delta Trend Summaries",
  "",
  "T cell signature deltas:",
  "",
  md_table(t_sig),
  "",
  "Myeloid signature deltas:",
  "",
  md_table(m_sig),
  "",
  "## Output Package Notes",
  "",
  "- Key dynamic panels now have matching source CSV/RDS files.",
  "- Key ggplot and pheatmap panels now have PDF outputs in addition to PNG.",
  "- Final figure styling still needs a Nature-style polish pass before manuscript submission."
)

writeLines(lines, out_file)
writeLines(lines)
