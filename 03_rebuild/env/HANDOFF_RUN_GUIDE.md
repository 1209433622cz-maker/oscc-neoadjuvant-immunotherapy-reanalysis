# Handoff Run Guide

Date: 2026-06-22

This guide is for running the OSCC GSE200996 rebuild on local compute.

## 1. Bootstrap Environment

Run this first on the target Windows machine:

```powershell
powershell -ExecutionPolicy Bypass -File H:\SCI2\OSCC-GSE200996-2025.12\03_rebuild\env\bootstrap_rebuild_env.ps1
```

To also install optional heavy packages:

```powershell
powershell -ExecutionPolicy Bypass -File H:\SCI2\OSCC-GSE200996-2025.12\03_rebuild\env\bootstrap_rebuild_env.ps1 -Optional
```

Outputs:

- `03_rebuild/logs/ENVIRONMENT_STATUS_LATEST.md`
- `03_rebuild/logs/R_package_status_latest.csv`
- `03_rebuild/logs/PYTHON_ENVIRONMENT_STATUS_LATEST.md`

## 2. Run Lightweight Checks

Data audit:

```powershell
powershell -ExecutionPolicy Bypass -File H:\SCI2\OSCC-GSE200996-2025.12\03_rebuild\env\run_rebuild_job.ps1 -Job audit
```

Result summary refresh:

```powershell
powershell -ExecutionPolicy Bypass -File H:\SCI2\OSCC-GSE200996-2025.12\03_rebuild\env\run_rebuild_job.ps1 -Job summary
```

## 3. Run Full Reanalysis

Use cached Seurat object if `03_rebuild/obj_full_QC_logUMAP.rds` already exists:

```powershell
powershell -ExecutionPolicy Bypass -File H:\SCI2\OSCC-GSE200996-2025.12\03_rebuild\env\run_rebuild_job.ps1 -Job full -UseCachedObj true
```

Rebuild the Seurat object from raw h5 files:

```powershell
powershell -ExecutionPolicy Bypass -File H:\SCI2\OSCC-GSE200996-2025.12\03_rebuild\env\run_rebuild_job.ps1 -Job full -UseCachedObj false
```

Expected heavy inputs:

- `00_raw_data/GSE200996_RAW/*.h5`
- `00_raw_data/GSE200996_metadata/GSE200996_CD45.tumor.single.cell.meta.data.txt.gz`

Expected heavy output:

- `03_rebuild/obj_full_QC_logUMAP.rds`
- `03_rebuild/results/pre_baseline/`
- `03_rebuild/results/dynamic_paired/`
- `03_rebuild/results/main_story/`

## 4. Run Leave-One-Patient Diagnostics

```powershell
powershell -ExecutionPolicy Bypass -File H:\SCI2\OSCC-GSE200996-2025.12\03_rebuild\env\run_rebuild_job.ps1 -Job loo
```

Outputs:

- `03_rebuild/results/sensitivity_leave_one_out/LOO_model_stability_summary.csv`
- `03_rebuild/results/sensitivity_leave_one_out/LOO_key_pathway_NES.csv`
- `03_rebuild/results/sensitivity_leave_one_out/LOO_stat_correlation.png`
- `03_rebuild/results/sensitivity_leave_one_out/LOO_key_pathway_NES.png`

## 5. Run Everything

This is the heavy route:

```powershell
powershell -ExecutionPolicy Bypass -File H:\SCI2\OSCC-GSE200996-2025.12\03_rebuild\env\run_rebuild_job.ps1 -Job all-heavy -UseCachedObj true
```

Use `-UseCachedObj false` only when the raw h5 object must be rebuilt.

## 6. If A Job Fails

Do not delete outputs. Send or inspect:

- the newest `*.out.log` in `03_rebuild/logs/`
- the matching newest `*.err.log` in `03_rebuild/logs/`
- `03_rebuild/logs/ENVIRONMENT_STATUS_LATEST.md`

The script names in the logs identify which job failed.

## 7. Division Of Labor

Local user can handle:

- downloading large raw data files;
- running `all-heavy`;
- rebuilding the Seurat object from h5 files;
- running leave-one-patient diagnostics;
- installing optional heavy packages.

Codex can handle:

- editing and checking scripts;
- interpreting logs;
- summarizing results;
- converting outputs into figures, methods, manuscript text and response-ready evidence.
