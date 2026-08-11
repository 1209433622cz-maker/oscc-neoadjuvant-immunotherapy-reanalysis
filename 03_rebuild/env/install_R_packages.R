args <- commandArgs(trailingOnly = TRUE)
install_optional <- any(args %in% c("--optional", "--all")) ||
  tolower(Sys.getenv("GSE200996_INSTALL_OPTIONAL", "FALSE")) %in% c("1", "true", "yes")

repo_dir <- Sys.getenv("GSE200996_WORKSPACE")
if (repo_dir == "") repo_dir <- getwd()
repo_dir <- normalizePath(repo_dir, mustWork = FALSE)
env_dir <- file.path(repo_dir, "03_rebuild", "env")
log_dir <- file.path(repo_dir, "03_rebuild", "logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

pkg_file <- file.path(env_dir, "required_R_packages.csv")
if (!file.exists(pkg_file)) {
  stop("Package manifest not found: ", pkg_file)
}

options(repos = c(CRAN = "https://cloud.r-project.org"))

manifest <- read.csv(pkg_file, stringsAsFactors = FALSE)
tiers <- c("core", "figure", "manuscript", "reproducibility")
if (install_optional) tiers <- c(tiers, "optional")
manifest <- manifest[manifest$tier %in% tiers, , drop = FALSE]

message("Installing/checking R packages")
message("Install optional packages: ", install_optional)
message("Package count: ", nrow(manifest))

install_cran <- function(pkgs) {
  if (!length(pkgs)) return(invisible(NULL))
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    message("Installing CRAN packages: ", paste(missing, collapse = ", "))
    install.packages(missing, dependencies = TRUE)
  } else {
    message("All CRAN packages already available.")
  }
}

install_bioc <- function(pkgs) {
  if (!length(pkgs)) return(invisible(NULL))
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    message("Installing Bioconductor packages: ", paste(missing, collapse = ", "))
    BiocManager::install(missing, ask = FALSE, update = FALSE)
  } else {
    message("All Bioconductor packages already available.")
  }
}

cran_pkgs <- unique(manifest$package[manifest$source == "CRAN"])
bioc_pkgs <- unique(manifest$package[manifest$source == "Bioc"])

install_cran(cran_pkgs)
install_bioc(bioc_pkgs)

status <- manifest
status$available <- vapply(status$package, requireNamespace, logical(1), quietly = TRUE)
status$version <- vapply(status$package, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) return(NA_character_)
  as.character(utils::packageVersion(pkg))
}, character(1))

out_csv <- file.path(log_dir, "R_package_status_latest.csv")
write.csv(status, out_csv, row.names = FALSE)
print(status[, c("package", "source", "tier", "available", "version")], row.names = FALSE)

if (any(!status$available)) {
  failed <- status$package[!status$available]
  stop("Missing packages after installation: ", paste(failed, collapse = ", "))
}

message("R package installation/check complete.")
message("Status CSV: ", out_csv)
