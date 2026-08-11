repo_dir <- Sys.getenv("GSE200996_WORKSPACE")
if (repo_dir == "") repo_dir <- getwd()
repo_dir <- normalizePath(repo_dir, mustWork = FALSE)

env_dir <- file.path(repo_dir, "03_rebuild", "env")
log_dir <- file.path(repo_dir, "03_rebuild", "logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

pkg_file <- file.path(env_dir, "required_R_packages.csv")
manifest <- read.csv(pkg_file, stringsAsFactors = FALSE)

status <- manifest
status$available <- vapply(status$package, requireNamespace, logical(1), quietly = TRUE)
status$version <- vapply(status$package, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) return(NA_character_)
  as.character(utils::packageVersion(pkg))
}, character(1))

write.csv(status, file.path(log_dir, "R_package_status_latest.csv"), row.names = FALSE)

path_status <- data.frame(
  item = c(
    "workspace",
    "raw_h5_dir",
    "cd45_tumor_metadata",
    "cached_seurat_object",
    "msigdb_cache",
    "main_pipeline",
    "leave_one_out_script",
    "result_readout"
  ),
  path = c(
    repo_dir,
    file.path(repo_dir, "00_raw_data", "GSE200996_RAW"),
    file.path(repo_dir, "00_raw_data", "GSE200996_metadata", "GSE200996_CD45.tumor.single.cell.meta.data.txt.gz"),
    file.path(repo_dir, "03_rebuild", "obj_full_QC_logUMAP.rds"),
    file.path(Sys.getenv("LOCALAPPDATA"), "R", "cache", "R", "msigdbr", "msigdb.2025.1.Hs.rds"),
    file.path(repo_dir, "03_rebuild", "analysis", "run_pipeline_rebuild.R"),
    file.path(repo_dir, "03_rebuild", "analysis", "04_leave_one_patient_diagnostics.R"),
    file.path(repo_dir, "03_rebuild", "results", "REANALYSIS_RESULT_READOUT.md")
  ),
  stringsAsFactors = FALSE
)
path_status$exists <- file.exists(path_status$path) | dir.exists(path_status$path)
path_status$size_mb <- vapply(path_status$path, function(p) {
  if (!file.exists(p)) return(NA_real_)
  round(file.info(p)$size / 1024^2, 2)
}, numeric(1))

write.csv(path_status, file.path(log_dir, "rebuild_path_status_latest.csv"), row.names = FALSE)

session_txt <- capture.output(sessionInfo())
writeLines(session_txt, file.path(log_dir, "R_sessionInfo_latest.txt"))

stamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
md <- c(
  "# Rebuild Environment Status",
  "",
  paste0("Checked: ", stamp),
  "",
  "## R Runtime",
  "",
  paste0("- R home: `", R.home(), "`"),
  paste0("- R version: `", R.version.string, "`"),
  paste0("- Platform: `", R.version$platform, "`"),
  paste0("- Working directory: `", getwd(), "`"),
  "",
  "## Required Paths",
  "",
  "| Item | Exists | Size MB | Path |",
  "|---|---:|---:|---|",
  apply(path_status, 1, function(x) {
    paste0("| ", x[["item"]], " | ", x[["exists"]], " | ", x[["size_mb"]], " | `", x[["path"]], "` |")
  }),
  "",
  "## Package Status",
  "",
  "| Package | Source | Tier | Available | Version |",
  "|---|---|---|---:|---:|",
  apply(status, 1, function(x) {
    paste0("| ", x[["package"]], " | ", x[["source"]], " | ", x[["tier"]], " | ", x[["available"]], " | ", x[["version"]], " |")
  }),
  "",
  "## Files Written",
  "",
  "- `03_rebuild/logs/R_package_status_latest.csv`",
  "- `03_rebuild/logs/rebuild_path_status_latest.csv`",
  "- `03_rebuild/logs/R_sessionInfo_latest.txt`"
)

out_md <- file.path(log_dir, "ENVIRONMENT_STATUS_LATEST.md")
writeLines(md, out_md)
cat(paste(md, collapse = "\n"))
cat("\n")
