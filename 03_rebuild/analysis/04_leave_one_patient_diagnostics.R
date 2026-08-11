rm(list = ls())
gc()
set.seed(1234)

base_dir <- Sys.getenv("GSE200996_BASEDIR")
if (base_dir == "") base_dir <- "H:/SCI2/OSCC-GSE200996-2025.12/03_rebuild"

meta_cd45_tumor_path <- Sys.getenv("GSE200996_META")
if (meta_cd45_tumor_path == "") {
  meta_cd45_tumor_path <- "H:/SCI2/OSCC-GSE200996-2025.12/00_raw_data/GSE200996_metadata/GSE200996_CD45.tumor.single.cell.meta.data.txt.gz"
}

out_dir <- file.path(base_dir, "results", "sensitivity_leave_one_out")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(DESeq2)
  library(fgsea)
})

if (requireNamespace("BiocParallel", quietly = TRUE)) {
  BiocParallel::register(BiocParallel::SerialParam())
}

message("Reading metadata...")
dt0 <- fread(meta_cd45_tumor_path)
setnames(
  dt0,
  c("V1","Patient_ID","Stage","Cohort","Path_response","CellType_ID","UMAP_1","UMAP_2"),
  c("cell_id","patient_id_raw","sample_stage","cohort","path_response","celltype_id","umap_1","umap_2")
)
dt0[, sample_tag := sub("^[^_]+_", "", cell_id)]
dt0[, patient_id := str_extract(sample_tag, "P\\d+")]
dt0[, timepoint := fifelse(grepl("Pre", sample_tag, ignore.case=TRUE), "pre",
                           fifelse(grepl("Post", sample_tag, ignore.case=TRUE), "post", NA_character_))]
dt0[, response_bin := fifelse(path_response == "High", "R",
                              fifelse(path_response == "Low", "NR", NA_character_))]
dt0[, response_ord := factor(path_response, levels=c("Low","Medium","High"), ordered=TRUE)]

patient_map <- dt0[!is.na(patient_id),
                   .(cohort=unique(cohort)[1],
                     path_response=unique(path_response)[1],
                     response_bin=unique(response_bin)[1],
                     response_ord=unique(response_ord)[1]),
                   by=patient_id]

load_hallmark_sets <- function() {
  cache_file <- file.path(Sys.getenv("LOCALAPPDATA"), "R", "cache", "R", "msigdbr", "msigdb.2025.1.Hs.rds")
  if (file.exists(cache_file)) {
    msig_local <- readRDS(cache_file)
    msig_local <- msig_local[
      msig_local$gs_collection == "H" &
        msig_local$db_target_species == "HS" &
        !is.na(msig_local$db_gene_symbol),
    ]
    return(split(msig_local$db_gene_symbol, msig_local$gs_name))
  }
  stop("Local MSigDB cache not found: ", cache_file)
}

hallmark <- load_hallmark_sets()

message("Loading cached Seurat object...")
obj_rds <- file.path(base_dir, "obj_full_QC_logUMAP.rds")
obj_full <- readRDS(obj_rds)
meta0 <- obj_full@meta.data

paired_patients <- intersect(unique(meta0$patient_id[meta0$timepoint == "pre"]),
                             unique(meta0$patient_id[meta0$timepoint == "post"]))
cells_dyn_all <- rownames(meta0)[
  meta0$patient_id %in% paired_patients &
    meta0$timepoint %in% c("pre","post") &
    !is.na(meta0$cd45_celltype_id) &
    !is.na(meta0$response_ord)
]
obj_dyn <- subset(obj_full, cells = cells_dyn_all)

make_pb_one_ct <- function(obj, ct, min_cells = 30) {
  sub_cells <- Cells(obj)[obj$cd45_celltype_id == ct]
  if (length(sub_cells) == 0) stop("No cells for celltype: ", ct)
  x <- subset(obj, cells = sub_cells)
  x$pb_id <- paste(x$patient_id, x$timepoint, sep="|")
  pb_n <- table(x$pb_id)
  keep_pb <- names(pb_n)[pb_n >= min_cells]
  x <- subset(x, cells = Cells(x)[x$pb_id %in% keep_pb])

  cnt <- AggregateExpression(x, assays="RNA", slot="counts",
                             group.by="pb_id", return.seurat=FALSE)$RNA
  cnt <- round(cnt)

  md <- data.frame(pb_id=colnames(cnt), stringsAsFactors=FALSE)
  tmp <- stringr::str_split_fixed(md$pb_id, "\\|", 2)
  md$patient_id <- tmp[,1]
  md$timepoint <- factor(tmp[,2], levels=c("pre","post"))

  pm <- as.data.frame(patient_map)[, c("patient_id","cohort","response_bin","response_ord")]
  pm$resp_num <- as.numeric(pm$response_ord)
  md <- dplyr::left_join(md, pm, by="patient_id")
  md$post <- as.integer(md$timepoint == "post")

  tp_tab <- table(md$patient_id, md$timepoint)
  keep_pat <- rownames(tp_tab)[apply(tp_tab > 0, 1, all)]
  md <- md[md$patient_id %in% keep_pat, , drop=FALSE]
  cnt <- cnt[, md$pb_id, drop=FALSE]

  rownames(md) <- md$pb_id
  md$patient_id <- factor(md$patient_id)
  list(counts=cnt, meta=md)
}

fit_interaction <- function(pb, leave_out = NA_character_) {
  cnt <- pb$counts
  md <- pb$meta

  if (!is.na(leave_out)) {
    md <- md[as.character(md$patient_id) != leave_out, , drop=FALSE]
    cnt <- cnt[, rownames(md), drop=FALSE]
  }

  md$patient_id <- droplevels(factor(as.character(md$patient_id)))
  md$timepoint <- factor(as.character(md$timepoint), levels=c("pre","post"))
  md$post <- as.integer(md$timepoint == "post")
  md$resp_num <- as.numeric(md$resp_num)

  if (length(unique(md$patient_id)) < 4) stop("Too few paired patients after leave-out.")
  if (length(unique(md$resp_num)) < 2) stop("Too few response levels after leave-out.")

  cnt <- cnt[rowSums(cnt) >= 10, , drop=FALSE]
  dds <- DESeqDataSetFromMatrix(countData=cnt, colData=md,
                                design=~ patient_id + post + post:resp_num)
  dds <- DESeq(dds, quiet=TRUE)
  rn <- resultsNames(dds)
  int_name <- grep("post.*resp_num|resp_num.*post", rn, value=TRUE)[1]
  if (is.na(int_name)) stop("Cannot find interaction term: ", paste(rn, collapse=" | "))

  res <- as.data.frame(results(dds, name=int_name))
  res$gene <- rownames(res)
  res <- res[is.finite(res$stat), , drop=FALSE]
  res <- res[order(ifelse(is.na(res$padj), Inf, res$padj),
                   ifelse(is.na(res$pvalue), Inf, res$pvalue)), , drop=FALSE]

  ranks <- res$stat
  names(ranks) <- res$gene
  ranks <- sort(ranks[is.finite(ranks)], decreasing=TRUE)
  fg <- fgseaMultilevel(pathways=hallmark, stats=ranks, minSize=10, maxSize=500)
  fg <- fg[order(fg$padj), ]

  list(dds=dds, res=res, fg=fg, n_patients=length(unique(md$patient_id)), n_samples=nrow(md))
}

top_genes <- function(res, n = 100) {
  head(res$gene[order(ifelse(is.na(res$padj), Inf, res$padj),
                      ifelse(is.na(res$pvalue), Inf, res$pvalue))], n)
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

all_summary <- list()
all_paths <- list()
all_gene_conc <- list()

for (ct in c("T cell", "Myeloid")) {
  message("Processing ", ct)
  pb <- make_pb_one_ct(obj_dyn, ct, min_cells=30)
  pts <- sort(unique(as.character(pb$meta$patient_id)))

  full <- fit_interaction(pb, leave_out=NA_character_)
  full_res <- full$res
  full_stat <- setNames(full_res$stat, full_res$gene)
  full_top50 <- top_genes(full_res, 50)
  full_top100 <- top_genes(full_res, 100)
  full_top25 <- top_genes(full_res, 25)

  for (leave in c("NONE", pts)) {
    message("  leave_out = ", leave)
    fit <- if (leave == "NONE") full else tryCatch(
      fit_interaction(pb, leave_out=leave),
      error=function(e) e
    )

    if (inherits(fit, "error")) {
      all_summary[[length(all_summary)+1]] <- data.frame(
        celltype=ct, leave_out=leave, status="FAIL",
        n_patients=NA_integer_, n_samples=NA_integer_,
        stat_spearman=NA_real_, top50_overlap=NA_real_, top100_overlap=NA_real_,
        full_top25_direction_concordance=NA_real_,
        note=conditionMessage(fit)
      )
      next
    }

    res <- fit$res
    stat <- setNames(res$stat, res$gene)
    common <- intersect(names(full_stat), names(stat))
    full_common <- full_stat[common]
    stat_common <- stat[common]

    stat_cor <- suppressWarnings(cor(full_common, stat_common, method="spearman", use="complete.obs"))
    this_top50 <- top_genes(res, 50)
    this_top100 <- top_genes(res, 100)
    top50_overlap <- length(intersect(full_top50, this_top50)) / length(full_top50)
    top100_overlap <- length(intersect(full_top100, this_top100)) / length(full_top100)

    top25_common <- intersect(full_top25, names(stat))
    direction_conc <- mean(sign(full_stat[top25_common]) == sign(stat[top25_common]), na.rm=TRUE)

    all_summary[[length(all_summary)+1]] <- data.frame(
      celltype=ct, leave_out=leave, status="OK",
      n_patients=fit$n_patients, n_samples=fit$n_samples,
      stat_spearman=stat_cor,
      top50_overlap=top50_overlap,
      top100_overlap=top100_overlap,
      full_top25_direction_concordance=direction_conc,
      note=""
    )

    gene_tab <- data.frame(
      celltype=ct,
      leave_out=leave,
      gene=full_top25,
      full_stat=as.numeric(full_stat[full_top25]),
      loo_stat=as.numeric(stat[full_top25]),
      same_direction=sign(full_stat[full_top25]) == sign(stat[full_top25]),
      stringsAsFactors=FALSE
    )
    all_gene_conc[[length(all_gene_conc)+1]] <- gene_tab

    fg <- as.data.frame(fit$fg)
    fg$key_rank <- match(fg$pathway, key_pathways[[ct]])
    path_tab <- fg[fg$pathway %in% key_pathways[[ct]],
                   c("pathway","NES","pval","padj","size"), drop=FALSE]
    path_tab$celltype <- ct
    path_tab$leave_out <- leave
    all_paths[[length(all_paths)+1]] <- path_tab
  }
}

summary_df <- bind_rows(all_summary)
paths_df <- bind_rows(all_paths)
gene_conc_df <- bind_rows(all_gene_conc)

write.csv(summary_df, file.path(out_dir, "LOO_model_stability_summary.csv"), row.names=FALSE)
write.csv(paths_df, file.path(out_dir, "LOO_key_pathway_NES.csv"), row.names=FALSE)
write.csv(gene_conc_df, file.path(out_dir, "LOO_top25_gene_direction_concordance.csv"), row.names=FALSE)

summary_ok <- summary_df %>% filter(status == "OK", leave_out != "NONE")

p1 <- ggplot(summary_ok, aes(x=leave_out, y=stat_spearman, fill=celltype)) +
  geom_col(position=position_dodge(width=0.75), width=0.7) +
  coord_cartesian(ylim=c(0, 1)) +
  theme_bw() +
  labs(x="Left-out patient", y="Spearman correlation with full-model gene statistics",
       title="Leave-one-patient stability of interaction statistics")
ggsave(file.path(out_dir, "LOO_stat_correlation.png"), p1, width=8, height=4.5, dpi=300)
ggsave(file.path(out_dir, "LOO_stat_correlation.pdf"), p1, width=8, height=4.5)

p2 <- ggplot(summary_ok, aes(x=leave_out, y=top100_overlap, fill=celltype)) +
  geom_col(position=position_dodge(width=0.75), width=0.7) +
  coord_cartesian(ylim=c(0, 1)) +
  theme_bw() +
  labs(x="Left-out patient", y="Top-100 gene overlap with full model",
       title="Leave-one-patient top-gene stability")
ggsave(file.path(out_dir, "LOO_top100_overlap.png"), p2, width=8, height=4.5, dpi=300)
ggsave(file.path(out_dir, "LOO_top100_overlap.pdf"), p2, width=8, height=4.5)

paths_plot <- paths_df %>%
  filter(leave_out != "NONE") %>%
  mutate(pathway_short = sub("^HALLMARK_", "", pathway))

p3 <- ggplot(paths_plot, aes(x=leave_out, y=NES, group=pathway_short, colour=pathway_short)) +
  geom_hline(yintercept=0, linewidth=0.2) +
  geom_point(size=1.8) +
  geom_line(linewidth=0.4) +
  facet_wrap(~ celltype, scales="free_y") +
  theme_bw() +
  theme(legend.position="bottom") +
  labs(x="Left-out patient", y="NES", colour="Pathway",
       title="Leave-one-patient stability of key Hallmark programs")
ggsave(file.path(out_dir, "LOO_key_pathway_NES.png"), p3, width=10, height=6, dpi=300)
ggsave(file.path(out_dir, "LOO_key_pathway_NES.pdf"), p3, width=10, height=6)

fmt <- function(x) ifelse(is.na(x), "NA", formatC(x, digits=3, format="fg", flag="#"))

readout <- c(
  "# Leave-One-Patient Diagnostic Readout",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Purpose",
  "",
  "Assess whether the central T-cell and myeloid dynamic interaction signals are dominated by a single paired patient.",
  "",
  "## Model",
  "",
  "`~ patient_id + post + post:resp_num` was refit for T cell and Myeloid after leaving out each paired patient.",
  "",
  "## Stability Summary",
  "",
  paste(capture.output(print(summary_df)), collapse="\n"),
  "",
  "## Interpretation Notes",
  "",
  "- High Spearman correlation indicates broad gene-rank stability against the full model.",
  "- Top-gene overlap is expected to be more sensitive in a six-patient design.",
  "- Pathway-level stability is more important than isolated-gene stability for the manuscript claim.",
  "- If leaving out P32 collapses the signal, the final manuscript must state that the High-responder contrast is influential."
)

writeLines(readout, file.path(out_dir, "LEAVE_ONE_OUT_READOUT.md"))
message("Leave-one-patient diagnostics complete: ", out_dir)
