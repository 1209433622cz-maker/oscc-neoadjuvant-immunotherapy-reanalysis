############################################################
# GSE200996 | Full pipeline (Windows local, Seurat v5 friendly)
# Author: Chen Zhi
# Last update: 2025-12-28
#
# 主线目标：
# 1) 读取 tumor CD45 meta -> patient_map（cohort/response）
# 2) 读取 tumor scRNA raw_feature .h5 -> 合并 obj_full
# 3) QC + log流程 -> PCA/UMAP/cluster + 免疫注释挂回
# 4) Baseline(pre-only) 做强：
#    4.1 组成差异：R/NR（二分类） + response_ord（趋势）
#    4.2 pseudobulk + DESeq2：R/NR（二分类） + response_ord（趋势）
#    4.3 baseline “主文风格”热图：T/Myeloid Top30 genes（按 response_ord 排序）
# 5) Dynamic (paired) Fig4 做机制增强：
#    5.1 Fig4A：组成变化 Δ(post-pre) 的 R/NR 差异（样本少，支持性）
#    5.2 Fig4B：T/Myeloid 做趋势交互模型 + signatureΔ + Top30 Δ热图（按 response_ord 排序）
#    5.3 Fig4B：GSEA Hallmark（修复 write.csv list 列导致的报错）
############################################################

## =========================
## 0. Setup
## =========================
rm(list=ls()); gc()
set.seed(1234)

# -------- Paths (cross-platform) --------
repo_dir <- normalizePath(getwd())
base_dir <- Sys.getenv("GSE200996_BASEDIR")
if (base_dir == "") base_dir <- repo_dir
raw_dir  <- Sys.getenv("GSE200996_RAWDIR")
if (raw_dir == "") raw_dir <- file.path(base_dir, "data", "raw")
meta_cd45_tumor_path <- Sys.getenv("GSE200996_META")
if (meta_cd45_tumor_path == "") meta_cd45_tumor_path <- file.path(raw_dir, "GSE200996_CD45.tumor.single.cell.meta.data.txt.gz")

# Outputs
out_dir_pre <- file.path(base_dir, "results", "pre_baseline")
out_dir_dyn <- file.path(base_dir, "results", "dynamic_paired")
out_dir_main <- file.path(base_dir, "results", "main_story")
dir.create(out_dir_pre, recursive=TRUE, showWarnings=FALSE)
dir.create(out_dir_dyn, recursive=TRUE, showWarnings=FALSE)
dir.create(out_dir_main, recursive=TRUE, showWarnings=FALSE)

# -------- Rebuild parameters --------
min_counts_prefilter   <- 200
min_features_prefilter <- 100

qc_min_features <- 200
qc_max_features <- 6000
qc_max_mt       <- 20

npcs_use   <- 50
dims_use   <- 1:30
resolution <- 0.6

min_cells_per_patient_celltype <- 30
min_cells_pb_dyn <- 30

# -------- Packages --------
suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
  library(Matrix)
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(limma)
  library(DESeq2)
  # 全矢量热图（可选，不装也不影响导出 delta/ann）
    if (requireNamespace("ComplexHeatmap", quietly=TRUE)) library(ComplexHeatmap)
    if (requireNamespace("circlize", quietly=TRUE)) library(circlize)
  # 机制/通路/热图
  # install.packages(c("msigdbr","fgsea","pheatmap"))
  library(msigdbr)
  library(fgsea)
  library(pheatmap)
})

## =========================
## 1. 读取 tumor CD45 meta -> patient_map
## =========================
message("Step1: Reading tumor CD45 meta...")

dt0 <- fread(meta_cd45_tumor_path)

setnames(dt0,
         c("V1","Patient_ID","Stage","Cohort","Path_response","CellType_ID","UMAP_1","UMAP_2"),
         c("cell_id","patient_id_raw","sample_stage","cohort","path_response","celltype_id","umap_1","umap_2"))

# 从 cell_id 解析 patient 与 pre/post
dt0[, sample_tag := sub("^[^_]+_", "", cell_id)]  # P13_Post-Tx
dt0[, patient_id := str_extract(sample_tag, "P\\d+")]
dt0[, timepoint  := fifelse(grepl("Pre",  sample_tag, ignore.case=TRUE), "pre",
                            fifelse(grepl("Post", sample_tag, ignore.case=TRUE), "post", NA_character_))]

# response：主文二分类 High vs Low（Medium -> NA）
dt0[, response_bin := fifelse(path_response=="High","R",
                              fifelse(path_response=="Low","NR", NA_character_))]

# response：机制/趋势用有序 Low<Medium<High
dt0[, response_ord := factor(path_response, levels=c("Low","Medium","High"), ordered=TRUE)]

# 患者级映射（一人一行）
patient_map <- dt0[!is.na(patient_id),
                   .(cohort        = unique(cohort)[1],
                     path_response  = unique(path_response)[1],
                     response_bin   = unique(response_bin)[1],
                     response_ord   = unique(response_ord)[1]),
                   by = patient_id]

# 检查唯一性
dup_pat <- patient_map[, .N, by=patient_id][N>1]
print(dup_pat)
stopifnot(nrow(dup_pat) == 0)

message("Step1 done: patient_map rows = ", nrow(patient_map))

obj_rds <- file.path(base_dir, "obj_full_QC_logUMAP.rds")
use_cached_obj <- file.exists(obj_rds) &&
  tolower(Sys.getenv("GSE200996_USE_CACHED_OBJ", "TRUE")) != "false"

if (use_cached_obj) {
  message("Loading cached Seurat object: ", obj_rds)
  obj_full <- readRDS(obj_rds)
  message("Loaded cached obj_full: cells = ", ncol(obj_full), " | genes = ", nrow(obj_full))
} else {

## =========================
## 2. 列出 tumor h5
## =========================
message("Step2: Listing tumor h5 files...")

tumor_h5 <- list.files(raw_dir, pattern="GEX_sc_tumor\\.h5$", full.names=TRUE)
message("tumor_h5 n = ", length(tumor_h5))
print(head(basename(tumor_h5), 5))
stopifnot(length(tumor_h5) > 0)

## =========================
## 3. 从文件名解析 sample_tag
## =========================
get_sample_tag <- function(fn){
  x <- basename(fn)
  tag <- str_match(x, "_(P\\d+_pre-Tx|P\\d+_post-Tx)_GEX_sc_tumor\\.h5$")[,2]
  if (is.na(tag)) stop("Cannot parse sample_tag from filename: ", x)
  tag <- str_replace(tag, "_pre-Tx$",  "_Pre-Tx")
  tag <- str_replace(tag, "_post-Tx$", "_Post-Tx")
  tag
}

## =========================
## 4. 读单个 tumor h5 + 空液滴粗过滤 + CreateSeuratObject
## =========================
read_one_tumor_prefilter <- function(h5,
                                     min_counts=200,
                                     min_features=100){
  tag <- get_sample_tag(h5)  # e.g. P13_Post-Tx
  
  mat <- Read10X_h5(h5)
  if (is.list(mat)) mat <- mat[["Gene Expression"]] %||% mat[[1]]
  
  # 统一 cell_id 格式：BARCODE_Pxx_Pre/Post-Tx
  bc <- gsub("-1$", "", colnames(mat))
  colnames(mat) <- paste0(bc, "_", tag)
  
  # 空液滴粗过滤（否则会出现上千万 cell）
  nCount   <- Matrix::colSums(mat)
  nFeature <- Matrix::colSums(mat != 0)
  keep <- (nCount >= min_counts) & (nFeature >= min_features)
  mat <- mat[, keep, drop=FALSE]
  
  obj <- CreateSeuratObject(mat, project="GSE200996_fulltumor")
  obj$sample_tag <- tag
  obj$patient_id <- str_extract(tag, "P\\d+")
  obj$timepoint  <- ifelse(grepl("Pre-Tx", tag), "pre",
                           ifelse(grepl("Post-Tx", tag), "post", NA))
  obj
}

## =========================
## 5. 合并为 obj_full
## =========================
message("Step3-5: Reading all tumor h5 -> merging...")

objs <- lapply(tumor_h5, read_one_tumor_prefilter,
               min_counts=min_counts_prefilter,
               min_features=min_features_prefilter)

obj_full <- Reduce(function(x,y) merge(x,y), objs)
rm(objs); gc()
message("Merged obj_full: cells = ", ncol(obj_full), " | genes = ", nrow(obj_full))

## =========================
## 6. 回填 cohort/response 到 obj_full
## =========================
message("Step6: Filling cohort/response from patient_map...")

pm <- as.data.frame(patient_map)
rownames(pm) <- pm$patient_id

obj_full$cohort        <- pm[obj_full$patient_id, "cohort"]
obj_full$path_response <- pm[obj_full$patient_id, "path_response"]
obj_full$response_bin  <- pm[obj_full$patient_id, "response_bin"]
obj_full$response_ord  <- pm[obj_full$patient_id, "response_ord"]

print(table(obj_full$cohort, useNA="ifany"))
print(table(obj_full$response_bin, useNA="ifany"))
print(table(obj_full$response_ord, useNA="ifany"))

## =========================
## 7. 挂 CD45 免疫细胞类型（cd45_celltype_id）
## =========================
message("Step7: Attaching CD45 celltype_id ...")

cd45_map <- dt0[, .(cell_id, cd45_celltype_id = celltype_id)]
cd45_map <- as.data.frame(cd45_map)
rownames(cd45_map) <- cd45_map$cell_id
obj_full$cd45_celltype_id <- cd45_map[Cells(obj_full), "cd45_celltype_id"]

message("CD45 annotated cells = ", sum(!is.na(obj_full$cd45_celltype_id)),
        " / total cells = ", ncol(obj_full))
print(table(is.na(obj_full$cd45_celltype_id)))  # TRUE=非免疫

## =========================
## 8. QC 指标 + QC 过滤
## =========================
message("Step8: QC metrics & filtering...")

DefaultAssay(obj_full) <- "RNA"
obj_full[["percent.mt"]]   <- PercentageFeatureSet(obj_full, pattern="^MT-")
obj_full[["percent.ribo"]] <- PercentageFeatureSet(obj_full, pattern="^RPS|^RPL")

obj_full <- subset(obj_full,
                   subset = nFeature_RNA > qc_min_features &
                     nFeature_RNA < qc_max_features &
                     percent.mt < qc_max_mt)
message("After QC: cells = ", ncol(obj_full))

## =========================
## 9. 轻量 log 流程：PCA/UMAP/聚类
## =========================
message("Step9: LogNormalize -> PCA/UMAP ...")

obj_full <- NormalizeData(obj_full, verbose=FALSE)
obj_full <- FindVariableFeatures(obj_full, selection.method="vst", nfeatures=3000, verbose=FALSE)
obj_full <- ScaleData(obj_full, features=VariableFeatures(obj_full), verbose=FALSE)
obj_full <- RunPCA(obj_full, features=VariableFeatures(obj_full), npcs=npcs_use, verbose=FALSE)

obj_full <- FindNeighbors(obj_full, dims=dims_use, verbose=FALSE)
obj_full <- FindClusters(obj_full, resolution=resolution, verbose=FALSE)
obj_full <- RunUMAP(obj_full, dims=dims_use, verbose=FALSE)

## =========================
## 10. marker 快速检查（免疫/上皮/基质）
## =========================
message("Step10: FeaturePlot marker check...")

FeaturePlot(obj_full,
            features = c("PTPRC","EPCAM","KRT8","KRT18","COL1A1","DCN","PECAM1","VWF"),
            ncol = 4, order = TRUE, reduction = "umap")

## =========================
## 11. Compartment 打分（可选）
## =========================
message("Step11: Compartment scoring...")

obj_full <- AddModuleScore(obj_full, features=list(c("PTPRC","LST1","CD3D","NKG7")), name="ImmuneScore_")
obj_full <- AddModuleScore(obj_full, features=list(c("EPCAM","KRT8","KRT18","KRT19")), name="EpiScore_")
obj_full <- AddModuleScore(obj_full, features=list(c("COL1A1","DCN","LUM","COL1A2")), name="FibScore_")
obj_full <- AddModuleScore(obj_full, features=list(c("PECAM1","VWF","KDR")), name="EndoScore_")

scores <- obj_full@meta.data[, c("ImmuneScore_1","EpiScore_1","FibScore_1","EndoScore_1")]
lab <- colnames(scores)[max.col(scores, ties.method="first")]

obj_full$compartment2 <- dplyr::recode(lab,
                                       "ImmuneScore_1"="Immune",
                                       "EpiScore_1"="Epithelial/Tumor",
                                       "FibScore_1"="Fibroblast",
                                       "EndoScore_1"="Endothelial")

DimPlot(obj_full, group.by="compartment2", label=TRUE, repel=TRUE)

## =========================
## 12. 保存对象（强烈建议）
## =========================
saveRDS(obj_full, obj_rds)
message("Saved obj_full_QC_logUMAP.rds")

}

############################################################
## 16. BASELINE (pre-only) 做强
## 16A: 组成差异（propeller-like：logit + limma）
##      A1) R/NR 二分类
##      A2) response_ord 趋势（Low/Medium/High -> 1/2/3）
## 16B: pseudobulk + DESeq2
##      B1) R/NR 二分类
##      B2) response_ord 趋势
## 16C: baseline 主文风格热图（T/Myeloid Top30，按 response_ord 排序）
############################################################
message("Step16: BASELINE (pre-only) 강화 ...")

meta0 <- obj_full@meta.data

# baseline：pre + immune only
meta_pre_all <- meta0 %>%
  filter(timepoint=="pre", !is.na(cd45_celltype_id))

# baseline: 二分类用
meta_pre_bin <- meta_pre_all %>% filter(response_bin %in% c("NR","R"))
meta_pre_bin$response_bin <- factor(meta_pre_bin$response_bin, levels=c("NR","R"))

# baseline: 趋势用（把 Medium 也吃进来增强 power）
meta_pre_trend <- meta_pre_all %>% filter(!is.na(response_ord))
meta_pre_trend$resp_num <- as.numeric(meta_pre_trend$response_ord)  # Low=1 Med=2 High=3

message("Baseline pre immune cells (bin)   = ", nrow(meta_pre_bin),
        " | patients = ", length(unique(meta_pre_bin$patient_id)))
message("Baseline pre immune cells (trend) = ", nrow(meta_pre_trend),
        " | patients = ", length(unique(meta_pre_trend$patient_id)))

## ---------- 16A1) 组成差异：R/NR ----------
baseline_comp_limma <- function(meta_df, mode=c("bin","trend")){
  mode <- match.arg(mode)
  
  ct_counts <- table(meta_df$cd45_celltype_id, meta_df$patient_id)
  ct_counts <- as.matrix(ct_counts)
  ct_props  <- sweep(ct_counts, 2, colSums(ct_counts), "/")
  
  eps <- 1e-4
  ct_logit <- log((ct_props + eps)/(1-ct_props + eps))
  
  # patient-level info（每人一行）
  if (mode=="bin"){
    si <- unique(meta_df[, c("patient_id","response_bin","cohort")])
    si <- si[match(colnames(ct_logit), si$patient_id), ]
    stopifnot(all(si$patient_id == colnames(ct_logit)))
    si$cohort <- factor(si$cohort)
    # cohort 可能只有一类，避免非满秩
    design <- if (length(unique(si$cohort))>1) model.matrix(~ cohort + response_bin, si)
    else model.matrix(~ response_bin, si)
    fit <- eBayes(lmFit(ct_logit, design))
    out <- topTable(fit, coef=grep("response_binR", colnames(design), value=TRUE)[1],
                    number=Inf, sort.by="P")
    out$celltype <- rownames(out)
    return(out)
  } else {
    si <- unique(meta_df[, c("patient_id","resp_num","response_ord","cohort")])
    si <- si[match(colnames(ct_logit), si$patient_id), ]
    stopifnot(all(si$patient_id == colnames(ct_logit)))
    si$cohort <- factor(si$cohort)
    # cohort 可能只有一类，避免非满秩
    design <- if (length(unique(si$cohort))>1) model.matrix(~ cohort + resp_num, si)
    else model.matrix(~ resp_num, si)
    fit <- eBayes(lmFit(ct_logit, design))
    out <- topTable(fit, coef=grep("resp_num", colnames(design), value=TRUE)[1],
                    number=Inf, sort.by="P")
    out$celltype <- rownames(out)
    return(out)
  }
}

# A1: R/NR
comp_pre_bin <- baseline_comp_limma(meta_pre_bin, mode="bin")
write.csv(comp_pre_bin, file.path(out_dir_pre, "baseline_pre_composition_RvsNR_limma_logit.csv"),
          row.names=FALSE)

# A2: trend
comp_pre_trend <- baseline_comp_limma(meta_pre_trend, mode="trend")
write.csv(comp_pre_trend, file.path(out_dir_pre, "baseline_pre_composition_trend_respOrd_limma_logit.csv"),
          row.names=FALSE)

## （可选）propeller 真·组成检验：有 speckle 就跑，没有就跳过
tryCatch({
  if (requireNamespace("speckle", quietly=TRUE)) {
    library(speckle)
    # bin
    prop_bin <- speckle::propeller(clusters = meta_pre_bin$cd45_celltype_id,
                                   sample   = meta_pre_bin$patient_id,
                                   group    = meta_pre_bin$response_bin,
                                   transform="logit")
    write.csv(prop_bin, file.path(out_dir_pre, "baseline_pre_composition_RvsNR_propeller.csv"),
              row.names=FALSE)
  } else {
    message("speckle 未安装：跳过 propeller（limma-logit 结果已足够替代）")
  }
}, error=function(e){
  message("[propeller ERROR] ", conditionMessage(e))
})

## ---------- 16B) baseline pseudobulk（共用函数） ----------
make_pb_baseline_pre <- function(obj, min_cells=30, use_bin=TRUE, use_trend=TRUE){
  # 只用 pre + immune
  cells_pre_imm <- rownames(obj@meta.data)[
    obj$timepoint=="pre" & !is.na(obj$cd45_celltype_id)
  ]
  x <- subset(obj, cells=cells_pre_imm)
  
  # pb_id：patient|celltype
  x$pb_id <- paste(x$patient_id, x$cd45_celltype_id, sep="|")
  
  pb_n <- table(x$pb_id)
  keep_pb <- names(pb_n)[pb_n >= min_cells]
  x <- subset(x, cells = Cells(x)[x$pb_id %in% keep_pb])
  
  cnt <- AggregateExpression(x, assays="RNA", slot="counts", group.by="pb_id",
                             return.seurat=FALSE)$RNA
  cnt <- round(cnt)
  
  md <- data.frame(pb_id=colnames(cnt), stringsAsFactors=FALSE)
  tmp <- str_split_fixed(md$pb_id, "\\|", 2)
  md$patient_id <- tmp[,1]
  md$celltype   <- tmp[,2]
  
  pm <- as.data.frame(patient_map)[, c("patient_id","cohort","response_bin","response_ord")]
  pm$resp_num <- as.numeric(pm$response_ord)
  md <- left_join(md, pm, by="patient_id")
  md$cohort <- factor(md$cohort)
  md$response_bin <- factor(md$response_bin, levels=c("NR","R"))
  
  rownames(md) <- md$pb_id
  list(counts=cnt, meta=md)
}

pb_pre <- make_pb_baseline_pre(obj_full, min_cells=min_cells_per_patient_celltype)

run_deseq2_baseline_one_ct <- function(pb, ct, mode=c("bin","trend")){
  mode <- match.arg(mode)
  
  md <- pb$meta
  cnt <- pb$counts
  
  idx <- which(md$celltype == ct)
  if (length(idx) < 4) return(NULL)
  
  md2  <- md[idx, , drop=FALSE]
  cnt2 <- cnt[, idx, drop=FALSE]
  
  # gene filter
  cnt2 <- cnt2[rowSums(cnt2) >= 10, , drop=FALSE]
  if (nrow(cnt2) < 200) return(NULL)
  
  if (mode=="bin"){
    md2 <- md2[!is.na(md2$response_bin), , drop=FALSE]
    cnt2 <- cnt2[, rownames(md2), drop=FALSE]
    tab <- table(md2$response_bin)
    if (length(tab)<2 || any(tab < 2)) return(NULL)
    
    design <- if (length(unique(md2$cohort))>1) ~ cohort + response_bin else ~ response_bin
    dds <- DESeqDataSetFromMatrix(cnt2, md2, design=design)
    dds <- DESeq(dds, quiet=TRUE)
    res <- as.data.frame(results(dds, contrast=c("response_bin","R","NR")))
    res$gene <- rownames(res)
    res$celltype <- ct
    res <- res[order(res$padj), ]
    return(list(dds=dds, res=res, coef_name=NULL))
  } else {
    md2 <- md2[!is.na(md2$resp_num), , drop=FALSE]
    cnt2 <- cnt2[, rownames(md2), drop=FALSE]
    if (length(unique(md2$resp_num)) < 2) return(NULL)
    
    design <- if (length(unique(md2$cohort))>1) ~ cohort + resp_num else ~ resp_num
    dds <- DESeqDataSetFromMatrix(cnt2, md2, design=design)
    dds <- DESeq(dds, quiet=TRUE)
    rn <- resultsNames(dds)
    coef_name <- grep("resp_num", rn, value=TRUE)[1]
    res <- as.data.frame(results(dds, name=coef_name))
    res$gene <- rownames(res)
    res$celltype <- ct
    res <- res[order(res$padj), ]
    return(list(dds=dds, res=res, coef_name=coef_name))
  }
}

## ---------- 16C) baseline 主文热图（T/Myeloid Top30, 按 response_ord 排序） ----------
# Hallmark gene sets（✅修复 '.' 问题）
load_hallmark_sets <- function() {
  cache_file <- file.path(
    Sys.getenv("LOCALAPPDATA"),
    "R", "cache", "R", "msigdbr", "msigdb.2025.1.Hs.rds"
  )
  if (file.exists(cache_file)) {
    message("Loading local MSigDB cache: ", cache_file)
    msig_local <- readRDS(cache_file)
    msig_local <- msig_local[
      msig_local$gs_collection == "H" &
        msig_local$db_target_species == "HS" &
        !is.na(msig_local$db_gene_symbol),
    ]
    return(split(msig_local$db_gene_symbol, msig_local$gs_name))
  }
  message("Local MSigDB cache not found; falling back to msigdbr download.")
  msig <- msigdbr(species="Homo sapiens", collection="H")
  split(msig$gene_symbol, msig$gs_name)
}
hallmark <- load_hallmark_sets()

# 写 csv 前把 list 列转成字符（✅修复 EncodeElement list 报错）
safe_write_fgsea <- function(fg, file){
  if ("leadingEdge" %in% colnames(fg)) {
    fg$leadingEdge <- vapply(fg$leadingEdge, function(x) paste(x, collapse=";"), FUN.VALUE=character(1))
  }
  write.csv(fg, file, row.names=FALSE)
}

plot_baseline_heatmap_top30 <- function(dds, coef_name, ct, out_png){
  vst_mat <- assay(vst(dds, blind=TRUE))
  md <- as.data.frame(colData(dds))
  md$pb_id <- rownames(md)
  
  res <- as.data.frame(results(dds, name=coef_name))
  res$gene <- rownames(res)
  # Top30
  if (all(is.na(res$padj))) res <- res[order(abs(res$stat), decreasing=TRUE), ]
  else res <- res[order(res$padj), ]
  top_genes <- head(res$gene[is.finite(res$stat)], 30)
  top_genes <- intersect(top_genes, rownames(vst_mat))
  if (length(top_genes) < 10) return(NULL)
  
  # 样本=patient（pre-only），按 response_ord 排序
  md2 <- md %>% arrange(response_ord)
  mat <- vst_mat[top_genes, rownames(md2), drop=FALSE]
  mat <- t(scale(t(mat))); mat[is.na(mat)] <- 0
  
  ann_col <- md2 %>% select(cohort, response_bin, response_ord)
  rownames(ann_col) <- rownames(md2)
  
  pheatmap(mat,
           cluster_rows=TRUE, cluster_cols=FALSE,
           annotation_col=ann_col,
           show_colnames=TRUE, fontsize_col=9,
           main=paste0("Baseline(pre) | ", ct, " | Top30 genes (ordered by response_ord)"),
           filename=out_png, width=9, height=10)
  if (grepl("\\.png$", out_png, ignore.case=TRUE)) {
    pheatmap(mat,
             cluster_rows=TRUE, cluster_cols=FALSE,
             annotation_col=ann_col,
             show_colnames=TRUE, fontsize_col=9,
             main=paste0("Baseline(pre) | ", ct, " | Top30 genes (ordered by response_ord)"),
             filename=sub("\\.png$", ".pdf", out_png, ignore.case=TRUE),
             width=9, height=10)
  }
  invisible(top_genes)
}

run_baseline_focus_ct <- function(ct){
  # trend 模式更强（把 Medium 也用上）
  rr <- run_deseq2_baseline_one_ct(pb_pre, ct, mode="trend")
  if (is.null(rr)) {
    message("[Baseline skip] ", ct, " trend: insufficient samples")
    return(NULL)
  }
  
  out_prefix <- file.path(out_dir_pre, paste0("Baseline_", gsub("[ /]","_",ct)))
  
  # 保存 DE
  write.csv(rr$res, paste0(out_prefix, "_DE_trend_respOrd.csv"), row.names=FALSE)
  
  # GSEA（用 stat 排序）
  ranks <- rr$res$stat; names(ranks) <- rr$res$gene
  ranks <- ranks[is.finite(ranks)]
  ranks <- sort(ranks, decreasing=TRUE)
  
  fg <- fgseaMultilevel(pathways=hallmark, stats=ranks, minSize=10, maxSize=500)
  fg <- fg[order(fg$padj), ]
  safe_write_fgsea(fg, paste0(out_prefix, "_GSEA_Hallmark.csv"))
  
  # bubble（不需要 leadingEdge）
  fg_plot <- fg %>% as.data.frame() %>%
    filter(!is.na(padj)) %>% arrange(padj) %>% slice_head(n=30) %>% arrange(NES)
  
  p_bub <- ggplot(fg_plot, aes(x=NES, y=reorder(pathway, NES), size=-log10(padj))) +
    geom_point() + theme_bw() +
    labs(title=paste0("Baseline(pre) | ", ct, " | GSEA Hallmark (response_ord trend)"),
         x="NES", y="Pathway")
  ggsave(paste0(out_prefix, "_GSEA_bubble.png"), p_bub, width=9, height=6, dpi=300)
  ggsave(paste0(out_prefix, "_GSEA_bubble.pdf"), p_bub, width=9, height=6)
  
  # baseline Top30 heatmap
  plot_baseline_heatmap_top30(rr$dds, rr$coef_name, ct, paste0(out_prefix, "_Top30_Heatmap.png"))
  
  rr
}

base_T <- run_baseline_focus_ct("T cell")
base_M <- run_baseline_focus_ct("Myeloid")

############################################################
## 17. Dynamic paired analysis (Fig4)
## Fig4A: 组成变化（post-pre）R/NR 差异（支持性）
## Fig4B: 核心（T/Myeloid）趋势交互模型 + signatureΔ + Top30 Δ热图
############################################################
message("Step17: Dynamic paired analysis (Fig4) ...")

# paired patients
paired_patients <- intersect(unique(meta0$patient_id[meta0$timepoint=="pre"]),
                             unique(meta0$patient_id[meta0$timepoint=="post"]))
message("Paired patients (pre & post) = ", length(paired_patients))
print(paired_patients)

## ---------- 17A Fig4A：组成变化 Δ(post-pre) 的 R/NR 差异 ----------
meta_dyn <- meta0 %>%
  filter(patient_id %in% paired_patients,
         timepoint %in% c("pre","post"),
         response_bin %in% c("R","NR"),
         !is.na(cd45_celltype_id))

meta_dyn$response_bin <- factor(meta_dyn$response_bin, levels=c("NR","R"))
meta_dyn$timepoint <- factor(meta_dyn$timepoint, levels=c("pre","post"))

meta_dyn_all <- meta0 %>%
  filter(patient_id %in% paired_patients,
         timepoint %in% c("pre","post"),
         !is.na(response_ord),
         !is.na(cd45_celltype_id))
meta_dyn_all$response_ord <- factor(meta_dyn_all$response_ord, levels=c("Low","Medium","High"), ordered=TRUE)
meta_dyn_all$response_bin <- factor(meta_dyn_all$response_bin, levels=c("NR","R"))
meta_dyn_all$timepoint <- factor(meta_dyn_all$timepoint, levels=c("pre","post"))
meta_dyn_all$resp_num <- as.numeric(meta_dyn_all$response_ord)

if (nrow(meta_dyn_all) > 0 && length(unique(meta_dyn_all$patient_id)) >= 3) {
  dt_counts_all <- as.data.table(meta_dyn_all)
  ct_pt_all <- dt_counts_all[, .N, by=.(cd45_celltype_id, patient_id, timepoint)]
  setnames(ct_pt_all, "N", "n_celltype")
  tot_pt_all <- dt_counts_all[, .N, by=.(patient_id, timepoint)]
  setnames(tot_pt_all, "N", "n_total")
  ct_pt_all <- merge(ct_pt_all, tot_pt_all, by=c("patient_id","timepoint"), all.x=TRUE)
  ct_pt_all[, prop := n_celltype / n_total]

  eps <- 1e-4
  ct_pt_all[, logit_prop := log((prop+eps)/(1-prop+eps))]
  wide_lp_all <- dcast(ct_pt_all, cd45_celltype_id + patient_id ~ timepoint,
                       value.var="logit_prop", fill=log((0+eps)/(1-0+eps)))
  wide_lp_all[, delta := post - pre]

  pt_info_ord <- unique(dt_counts_all[, .(patient_id, response_ord, response_bin, cohort)])
  pt_info_ord[, resp_num := as.numeric(factor(response_ord, levels=c("Low","Medium","High"), ordered=TRUE))]
  patients_order_all <- unique(wide_lp_all$patient_id)
  pt_info_ord <- pt_info_ord[match(patients_order_all, pt_info_ord$patient_id), ]

  delta_mat_all <- dcast(wide_lp_all, cd45_celltype_id ~ patient_id, value.var="delta")
  ct_names_all <- delta_mat_all$cd45_celltype_id
  delta_mat_all <- as.matrix(delta_mat_all[, -1, with=FALSE])
  rownames(delta_mat_all) <- ct_names_all

  fit_dyn_ord <- eBayes(lmFit(delta_mat_all, model.matrix(~ resp_num, data=pt_info_ord)))
  comp_dyn_ord <- topTable(fit_dyn_ord, coef="resp_num", number=Inf, sort.by="P")
  comp_dyn_ord$celltype <- rownames(comp_dyn_ord)
  write.csv(comp_dyn_ord, file.path(out_dir_dyn, "Fig4A_composition_delta_logit_limma_respOrd_trend.csv"),
            row.names=FALSE)

  if (all(is.na(comp_dyn_ord$adj.P.Val))) top_ct_ord <- head(comp_dyn_ord$celltype[order(comp_dyn_ord$P.Value)], 6)
  else top_ct_ord <- head(comp_dyn_ord$celltype[order(comp_dyn_ord$adj.P.Val)], 6)

  plot_df_ord <- ct_pt_all[cd45_celltype_id %in% top_ct_ord, .(patient_id, timepoint, cd45_celltype_id, prop)]
  plot_df_ord <- merge(plot_df_ord, pt_info_ord[, .(patient_id, response_ord, response_bin)], by="patient_id", all.x=TRUE)
  plot_df_ord$timepoint <- factor(plot_df_ord$timepoint, levels=c("pre","post"))

  p_fig4A_ord <- ggplot(plot_df_ord, aes(x=timepoint, y=prop, group=patient_id, colour=response_ord)) +
    geom_line(alpha=0.65) +
    geom_point(size=1.7) +
    facet_wrap(~ cd45_celltype_id, scales="free_y", ncol=3) +
    theme_bw() +
    labs(title="Fig4A | Paired composition change, ordinal response trend",
         x="Timepoint", y="Proportion within immune cells", colour="Response")
  ggsave(file.path(out_dir_dyn, "Fig4A_paired_composition_top_celltypes_ord.png"),
         p_fig4A_ord, width=14, height=7.5, dpi=300)
  ggsave(file.path(out_dir_dyn, "Fig4A_paired_composition_top_celltypes_ord.pdf"),
         p_fig4A_ord, width=14, height=7.5)
  saveRDS(plot_df_ord, file.path(out_dir_dyn, "Fig4A_plot_df_topCT_ord.rds"))
} else {
  message("Fig4A ordinal trend skipped: paired ordinal immune samples too few.")
}

if (nrow(meta_dyn) > 0 && length(unique(meta_dyn$patient_id)) >= 3) {
  
  dt_counts <- as.data.table(meta_dyn)
  ct_pt <- dt_counts[, .N, by=.(cd45_celltype_id, patient_id, timepoint)]
  setnames(ct_pt, "N", "n_celltype")
  tot_pt <- dt_counts[, .N, by=.(patient_id, timepoint)]
  setnames(tot_pt, "N", "n_total")
  
  ct_pt <- merge(ct_pt, tot_pt, by=c("patient_id","timepoint"), all.x=TRUE)
  ct_pt[, prop := n_celltype/n_total]
## ===== [NEW] Fig4A 长表导出：patient×timepoint×celltype×prop + response/cohort =====
  pt_info_all <- unique(dt_counts[, .(patient_id, response_bin, response_ord, cohort)])
  fig4A_long <- ct_pt[, .(
    patient   = patient_id,
    timepoint = timepoint,
    celltype  = cd45_celltype_id,
    prop      = prop
  )]
  fig4A_long <- merge(fig4A_long, pt_info_all, by.x="patient", by.y="patient_id", all.x=TRUE)
  
  write.csv(fig4A_long,
            file.path(out_dir_dyn, "Fig4A_patient_timepoint_celltype_prop_long.csv"),
            row.names = FALSE)
# 可选：把你用于作图的 plot_df 也存一份，方便完全复刻 Fig4A
  
  eps <- 1e-4
  ct_pt[, logit_prop := log((prop+eps)/(1-prop+eps))]
  
  wide_lp <- dcast(ct_pt, cd45_celltype_id + patient_id ~ timepoint,
                   value.var="logit_prop", fill=log((0+eps)/(1-0+eps)))
  wide_lp[, delta := post - pre]
  
  pt_info <- unique(dt_counts[, .(patient_id, response_bin)])
  patients_order <- unique(wide_lp$patient_id)
  pt_info <- pt_info[match(patients_order, pt_info$patient_id), ]
  
  delta_mat <- dcast(wide_lp, cd45_celltype_id ~ patient_id, value.var="delta")
  ct_names <- delta_mat$cd45_celltype_id
  delta_mat <- as.matrix(delta_mat[, -1, with=FALSE])
  rownames(delta_mat) <- ct_names
  
  fit_dyn <- eBayes(lmFit(delta_mat, model.matrix(~ response_bin, data=pt_info)))
  comp_dyn <- topTable(fit_dyn, coef="response_binR", number=Inf, sort.by="P")
  comp_dyn$celltype <- rownames(comp_dyn)
  write.csv(comp_dyn, file.path(out_dir_dyn, "Fig4A_composition_delta_logit_limma_RvsNR.csv"),
            row.names=FALSE)
  
  # top 6 celltypes画 paired lines
  if (all(is.na(comp_dyn$adj.P.Val))) top_ct <- head(comp_dyn$celltype[order(comp_dyn$P.Value)], 6)
  else top_ct <- head(comp_dyn$celltype[order(comp_dyn$adj.P.Val)], 6)
  
  plot_df <- ct_pt[cd45_celltype_id %in% top_ct, .(patient_id, timepoint, cd45_celltype_id, prop)]
  plot_df <- merge(plot_df, pt_info, by="patient_id", all.x=TRUE)
  plot_df$timepoint <- factor(plot_df$timepoint, levels=c("pre","post"))
  
  p_fig4A <- ggplot(plot_df, aes(x=timepoint, y=prop, group=patient_id, linetype=response_bin)) +
    geom_line(alpha=0.6) +
    geom_point(size=1.6) +
    facet_wrap(~ cd45_celltype_id, scales="free_y", ncol=3) +
    theme_bw() +
    labs(title="Fig4A | Paired composition change (pre→post), stratified by response_bin",
         x="Timepoint", y="Proportion within immune cells")
  ggsave(file.path(out_dir_dyn, "Fig4A_paired_composition_top_celltypes.png"),
         p_fig4A, width=14, height=7.5, dpi=300)
  ggsave(file.path(out_dir_dyn, "Fig4A_paired_composition_top_celltypes.pdf"),
         p_fig4A, width=14, height=7.5)
  saveRDS(plot_df, file.path(out_dir_dyn, "Fig4A_plot_df_topCT.rds"))
} else {
  message("Fig4A skipped: paired R/NR immune samples too few.")
}
## ---------- ✅ Fig4A 长表导出：patient×timepoint×celltype prop ----------
# ct_pt 里已经有 patient_id/timepoint/cd45_celltype_id/prop（上面算出来的）
# 补齐 response_bin/response_ord/cohort（每个 patient 一行）
pt_anno <- unique(dt_counts_all[, .(patient_id, response_bin, response_ord, cohort)])

fig4A_long <- merge(
  ct_pt_all[, .(patient_id, timepoint, celltype = cd45_celltype_id, prop)],
  pt_anno,
  by = "patient_id",
  all.x = TRUE
)

# 列名对齐你需要的格式：patient/timepoint/celltype/prop/response...
setnames(fig4A_long, "patient_id", "patient")
fig4A_long <- fig4A_long[order(patient, timepoint, celltype)]

write.csv(
  fig4A_long,
  file.path(out_dir_dyn, "Fig4A_patient_timepoint_celltype_prop_long.csv"),
  row.names = FALSE
)
write.csv(
  fig4A_long,
  file.path(out_dir_dyn, "Fig4A_patient_timepoint_celltype_prop_long_allpaired.csv"),
  row.names = FALSE
)
message("[OK] Fig4A long table saved: ", file.path(out_dir_dyn, "Fig4A_patient_timepoint_celltype_prop_long.csv"))


## ---------- 17B Fig4B：T & Myeloid（趋势交互 + signatureΔ + Top30 Δ热图） ----------
sig_T <- list(
  Cytotoxic  = c("NKG7","PRF1","GZMB","GNLY","GZMH"),
  Exhaustion = c("PDCD1","LAG3","HAVCR2","TIGIT","CTLA4"),
  IFNG       = c("IFNG","CXCL9","CXCL10","STAT1","IRF1"),
  Prolif     = c("MKI67","TOP2A","STMN1","TYMS")
)
sig_M <- list(
  Inflammatory = c("S100A8","S100A9","LYZ","FCGR3A","LGALS3"),
  Antigen_Pres = c("HLA-DRA","HLA-DRB1","CD74","HLA-DPA1","HLA-DPB1"),
  ISG          = c("ISG15","IFIT1","IFIT3","MX1","IFI6"),
  M2_like      = c("MRC1","CD163","MSR1","C1QA","C1QB")
)

# 动态用免疫细胞：paired + response_ord 不 NA
cells_dyn_all <- rownames(meta0)[
  meta0$patient_id %in% paired_patients &
    meta0$timepoint %in% c("pre","post") &
    !is.na(meta0$cd45_celltype_id) &
    !is.na(meta0$response_ord)
]
obj_dyn <- subset(obj_full, cells = cells_dyn_all)

# helper：某个 celltype 的 pseudobulk（patient|timepoint），严格配对
make_pb_one_ct <- function(obj, ct, min_cells = 30){
  sub_cells <- Cells(obj)[obj$cd45_celltype_id == ct]
  if (length(sub_cells) == 0) stop("No cells for celltype: ", ct)
  
  x <- subset(obj, cells = sub_cells)
  x$pb_id <- paste(x$patient_id, x$timepoint, sep="|")  # patient|timepoint
  
  pb_n <- table(x$pb_id)
  keep_pb <- names(pb_n)[pb_n >= min_cells]
  x <- subset(x, cells = Cells(x)[x$pb_id %in% keep_pb])
  
  cnt <- AggregateExpression(x, assays="RNA", slot="counts", group.by="pb_id",
                             return.seurat=FALSE)$RNA
  cnt <- round(cnt)
  
  md <- data.frame(pb_id = colnames(cnt), stringsAsFactors = FALSE)
  tmp <- stringr::str_split_fixed(md$pb_id, "\\|", 2)
  md$patient_id <- tmp[,1]
  md$timepoint  <- factor(tmp[,2], levels=c("pre","post"))
  
  pm <- as.data.frame(patient_map)[, c("patient_id","cohort","response_bin","response_ord")]
  pm$resp_num <- as.numeric(pm$response_ord)  # Low=1 Med=2 High=3
  md <- dplyr::left_join(md, pm, by="patient_id")
  
  # ✅用 post(0/1) 保证满秩：design ~ patient_id + post + post:resp_num
  md$post <- as.integer(md$timepoint == "post")
  
  # 只保留 pre/post 都存在的 patient（严格配对）
  tp_tab <- table(md$patient_id, md$timepoint)
  keep_pat <- rownames(tp_tab)[apply(tp_tab > 0, 1, all)]
  md <- md[md$patient_id %in% keep_pat, , drop=FALSE]
  cnt <- cnt[, md$pb_id, drop=FALSE]
  
  rownames(md) <- md$pb_id
  md$patient_id <- factor(md$patient_id)
  
  if (length(unique(md$resp_num)) < 2) stop(ct, " : resp_num 无变化，无法做趋势模型")
  list(counts=cnt, meta=md)
}

plot_delta_heatmap_top30 <- function(dds, int_name, ct, patient_map,
                                     out_png = NULL,
                                     export_prefix = NULL,
                                     top_n = 30){
  
  vst_mat <- assay(vst(dds, blind=TRUE))
  md <- as.data.frame(colData(dds))
  md$pb_id <- rownames(md)
  md$timepoint <- ifelse(md$post==1, "post", "pre")
  md$patient_id <- as.character(md$patient_id)
  
  res_int <- as.data.frame(results(dds, name=int_name))
  res_int$gene <- rownames(res_int)
  res_int <- res_int[is.finite(res_int$stat), , drop=FALSE]
  
  if (all(is.na(res_int$padj))) res_int <- res_int[order(abs(res_int$stat), decreasing=TRUE), ]
  else res_int <- res_int[order(res_int$padj), ]
  
  top_genes <- head(res_int$gene, top_n)
  top_genes <- intersect(top_genes, rownames(vst_mat))
  if (length(top_genes) < 10) {
    message("[Heatmap skip] ", ct, " : top genes < 10")
    return(NULL)
  }
  
  # -------- patient-level Δ(post-pre) --------
  pats <- unique(md$patient_id)
  delta <- matrix(NA_real_, nrow=length(top_genes), ncol=length(pats),
                  dimnames=list(top_genes, pats))
  
  for (p in pats) {
    pre_id  <- md$pb_id[md$patient_id==p & md$timepoint=="pre"]
    post_id <- md$pb_id[md$patient_id==p & md$timepoint=="post"]
    if (length(pre_id)==1 && length(post_id)==1) {
      delta[, p] <- vst_mat[top_genes, post_id] - vst_mat[top_genes, pre_id]
    }
  }
  
  # -------- annotation（列注释）并按 response_ord 排序 --------
  ann <- as.data.frame(patient_map)[, c("patient_id","response_ord","response_bin","cohort")]
  ann$patient_id <- as.character(ann$patient_id)
  ann <- ann[ann$patient_id %in% colnames(delta), , drop=FALSE]
  rownames(ann) <- ann$patient_id
  ann$patient_id <- NULL
  
  # 统一排序：Low < Medium < High（若你这里没有 Medium 也没关系）
  if (!is.factor(ann$response_ord)) {
    ann$response_ord <- factor(ann$response_ord, levels=c("Low","Medium","High"), ordered=TRUE)
  }
  ord <- order(ann$response_ord)
  ann <- ann[ord, , drop=FALSE]
  delta <- delta[, rownames(ann), drop=FALSE]
  
  # -------- 画 PNG（仍然用 pheatmap，保持你现有风格）--------
  if (!is.null(out_png)) {
    mat_plot <- t(scale(t(delta))); mat_plot[is.na(mat_plot)] <- 0
    
    pheatmap::pheatmap(
      mat_plot,
      cluster_rows = TRUE,
      cluster_cols = FALSE,
      annotation_col = ann,
      show_colnames = TRUE,
      fontsize_col = 9,
      main = paste0("Fig4B | ", ct, " | Top", top_n, " genes Δ(post-pre) (ordered by response_ord)"),
      filename = out_png,
      width = 9, height = 10
    )
    if (grepl("\\.png$", out_png, ignore.case=TRUE)) {
      pheatmap::pheatmap(
        mat_plot,
        cluster_rows = TRUE,
        cluster_cols = FALSE,
        annotation_col = ann,
        show_colnames = TRUE,
        fontsize_col = 9,
        main = paste0("Fig4B | ", ct, " | Top", top_n, " genes Δ(post-pre) (ordered by response_ord)"),
        filename = sub("\\.png$", ".pdf", out_png, ignore.case=TRUE),
        width = 9, height = 10
      )
    }
  }
  
  # -------- ✅ 导出 heatmap 输入：delta + ann --------
  if (!is.null(export_prefix)) {
    saveRDS(list(delta=delta, ann=ann), paste0(export_prefix, "_delta_and_ann.rds"))
    write.csv(delta, paste0(export_prefix, "_delta_matrix.csv"))
    write.csv(ann,   paste0(export_prefix, "_sample_annotation.csv"))
    message("[OK] Heatmap inputs saved: ", paste0(export_prefix, "_delta_and_ann.rds"))
  }
  
  invisible(list(top_genes=top_genes, delta=delta, ann=ann))
}


run_fig4B_ct <- function(ct, sig_list){
  message("\n===== Fig4B: ", ct, " =====")
  
  pb <- make_pb_one_ct(obj_dyn, ct, min_cells=min_cells_pb_dyn)
  cnt <- pb$counts
  md  <- pb$meta
  
  message("Pseudobulk samples = ", ncol(cnt), " | paired patients = ", length(unique(md$patient_id)))
  print(table(ifelse(md$post==1,"post","pre")))
  print(table(md$response_ord, useNA="ifany"))
  
  cnt <- cnt[rowSums(cnt) >= 10, , drop=FALSE]
  
  out_prefix <- file.path(out_dir_dyn, paste0("Fig4B_", gsub("[ /]","_",ct)))
  
  # 1) DESeq2（满秩）
  dds <- DESeqDataSetFromMatrix(countData=cnt, colData=md,
                                design = ~ patient_id + post + post:resp_num)
  dds <- DESeq(dds, quiet=TRUE)
  
  rn <- resultsNames(dds)
  int_name <- grep("post.*resp_num|resp_num.*post", rn, value=TRUE)[1]
  if (is.na(int_name)) stop("Cannot find interaction term. resultsNames=", paste(rn, collapse=" | "))
  
  res_int <- as.data.frame(results(dds, name=int_name))
  res_int$gene <- rownames(res_int)
  res_int <- res_int[order(res_int$padj), ]
  write.csv(res_int, paste0(out_prefix, "_interaction_DE_trend.csv"), row.names=FALSE)
  
  # 2) GSEA（✅修复保存 list 列）
  ranks <- res_int$stat; names(ranks) <- res_int$gene
  ranks <- ranks[is.finite(ranks)]
  ranks <- sort(ranks, decreasing=TRUE)
  
  fg <- fgseaMultilevel(pathways=hallmark, stats=ranks, minSize=10, maxSize=500)
  fg <- fg[order(fg$padj), ]
  safe_write_fgsea(fg, paste0(out_prefix, "_GSEA_Hallmark.csv"))
  
  fg_plot <- fg %>% as.data.frame() %>%
    filter(!is.na(padj)) %>% arrange(padj) %>% slice_head(n=30) %>% arrange(NES)
  
  p_bub <- ggplot(fg_plot, aes(x=NES, y=reorder(pathway, NES), size=-log10(padj))) +
    geom_point() + theme_bw() +
    labs(title=paste0("Fig4B | ", ct, " | GSEA Hallmark (post×response_ord trend)"),
         x="NES", y="Pathway")
  ggsave(paste0(out_prefix, "_GSEA_bubble.png"), p_bub, width=9, height=6, dpi=300)
  ggsave(paste0(out_prefix, "_GSEA_bubble.pdf"), p_bub, width=9, height=6)
  
  # 3) signature Δ(post-pre) vs response_ord
  vst_mat <- assay(vst(dds, blind=TRUE))
  cd <- as.data.frame(colData(dds))
  cd$pb_id <- rownames(cd)
  cd$timepoint <- ifelse(cd$post==1,"post","pre")
  
  sig_score <- function(genes){
    genes <- intersect(genes, rownames(vst_mat))
    if (length(genes)==0) return(rep(NA_real_, ncol(vst_mat)))
    colMeans(vst_mat[genes, , drop=FALSE])
  }
  
  sig_df <- lapply(names(sig_list), function(nm){
    sc <- sig_score(sig_list[[nm]])
    data.frame(pb_id=colnames(vst_mat), signature=nm, score=sc, stringsAsFactors=FALSE)
  }) %>% bind_rows() %>% left_join(cd, by="pb_id")
  
  sig_w <- sig_df %>%
    select(patient_id, timepoint, response_ord, resp_num, signature, score) %>%
    pivot_wider(names_from=timepoint, values_from=score) %>%
    mutate(delta = post - pre)
  write.csv(sig_w, paste0(out_prefix, "_SignatureDelta_source.csv"), row.names=FALSE)
  
  p_sig <- ggplot(sig_w, aes(x=resp_num, y=delta, label=as.character(patient_id))) +
    geom_point() +
    geom_text(vjust=-0.4, size=3) +
    geom_smooth(method="lm", se=FALSE) +
    facet_wrap(~ signature, scales="free_y") +
    theme_bw() +
    labs(title=paste0("Fig4B | ", ct, " | Signature Δ(post-pre) vs response_ord trend"),
         x="Response ordinal (Low=1, Medium=2, High=3)", y="Δ signature (VST)")
  ggsave(paste0(out_prefix, "_SignatureDelta.png"), p_sig, width=11, height=6.5, dpi=300)
  ggsave(paste0(out_prefix, "_SignatureDelta.pdf"), p_sig, width=11, height=6.5)
  
  # 4) ✅ Top30 Δ热图 + 同步导出 delta/ann（用于 ComplexHeatmap 全矢量）
  export_prefix <- if (ct == "T cell") {
    file.path(out_dir_dyn, "Fig4B_Tcell")
  } else if (ct == "Myeloid") {
    file.path(out_dir_dyn, "Fig4B_Myeloid")
  } else {
    file.path(out_dir_dyn, paste0("Fig4B_", gsub("[ /]","_",ct)))
  }
  
  plot_delta_heatmap_top30(
    dds = dds,
    int_name = int_name,
    ct = ct,
    patient_map = patient_map,
    out_png = paste0(out_prefix, "_Top30_DeltaHeatmap.png"),
    export_prefix = export_prefix,
    top_n = 30
  )
  
  invisible(list(dds=dds, res_int=res_int, fgsea=fg, sig_delta=sig_w))
}

# 只跑 T cell & Myeloid
res_T <- run_fig4B_ct("T cell",  sig_T)
res_M <- run_fig4B_ct("Myeloid", sig_M)

############################################################
## 18. 主线对齐表：baseline(趋势) + dynamic(Fig4A bin) + Fig4B GSEA top
############################################################
message("Step18: Main story alignment table...")

read_if_exists <- function(f) if (file.exists(f)) read.csv(f) else data.frame()

base_comp_trend <- read_if_exists(file.path(out_dir_pre, "baseline_pre_composition_trend_respOrd_limma_logit.csv"))
dyn_comp_ord    <- read_if_exists(file.path(out_dir_dyn, "Fig4A_composition_delta_logit_limma_respOrd_trend.csv"))
dyn_comp_bin    <- read_if_exists(file.path(out_dir_dyn, "Fig4A_composition_delta_logit_limma_RvsNR.csv"))

focus_ct <- c("T cell","Myeloid")

base_comp2 <- base_comp_trend %>%
  filter(celltype %in% focus_ct) %>%
  transmute(celltype,
            baseline_trend_logFC = logFC,
            baseline_trend_P     = P.Value,
            baseline_trend_FDR   = adj.P.Val)

dyn_comp_ord2 <- dyn_comp_ord %>%
  filter(celltype %in% focus_ct) %>%
  transmute(celltype,
            dyn_delta_ord_logFC = logFC,
            dyn_delta_ord_P     = P.Value,
            dyn_delta_ord_FDR   = adj.P.Val)

dyn_comp2 <- dyn_comp_bin %>%
  filter(celltype %in% focus_ct) %>%
  transmute(celltype,
            dyn_delta_RvsNR_logFC = logFC,
            dyn_delta_RvsNR_P     = P.Value,
            dyn_delta_RvsNR_FDR   = adj.P.Val)

read_top_path <- function(ct){
  f <- file.path(out_dir_dyn, paste0("Fig4B_", gsub("[ /]","_",ct), "_GSEA_Hallmark.csv"))
  if (!file.exists(f)) return(data.frame(celltype=ct, top_pathways=NA))
  x <- read.csv(f)
  x <- x[!is.na(x$padj), ]
  x <- x[order(x$padj), ]
  top <- head(x$pathway, 8)
  data.frame(celltype=ct, top_pathways=paste(top, collapse="; "))
}
path_df <- bind_rows(read_top_path("T cell"), read_top_path("Myeloid"))

story_tab <- base_comp2 %>%
  full_join(dyn_comp_ord2, by="celltype") %>%
  full_join(dyn_comp2, by="celltype") %>%
  left_join(path_df, by="celltype")

write.csv(story_tab, file.path(out_dir_main, "MainStory_Alignment_Table.csv"), row.names=FALSE)
print(story_tab)

message("ALL DONE. Key outputs:")
message(" - Baseline: ", out_dir_pre, " (composition + DE trend + baseline heatmaps + GSEA)")
message(" - Dynamic Fig4: ", out_dir_dyn, " (Fig4A + Fig4B bubble/signature/Δheatmap)")
message(" - Main story table: ", file.path(out_dir_main, "MainStory_Alignment_Table.csv"))

