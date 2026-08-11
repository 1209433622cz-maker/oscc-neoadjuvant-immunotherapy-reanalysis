# Local Compute Task List

Date: 2026-06-22

This file lists tasks that can be handed to local compute when they are too slow or require large downloads.

## Current Status

Completed on this machine:

- Core full reanalysis with cached Seurat object.
- Ordinal and binary paired composition analyses.
- T-cell and myeloid pseudobulk interaction models.
- Hallmark GSEA.
- Leave-one-patient diagnostics.
- Core, figure and manuscript R packages.
- Python helper packages.

## Task A: Rebuild Seurat Object From Raw h5 Files

Use when:

- `03_rebuild/obj_full_QC_logUMAP.rds` is missing;
- raw h5 files changed;
- QC parameters are changed.

Command:

```powershell
powershell -ExecutionPolicy Bypass -File H:\SCI2\OSCC-GSE200996-2025.12\03_rebuild\env\run_rebuild_job.ps1 -Job full -UseCachedObj false
```

Expected output:

- `03_rebuild/obj_full_QC_logUMAP.rds`
- refreshed `03_rebuild/results/`

Send back if failed:

- newest `02_full_reanalysis_*.out.log`
- newest `02_full_reanalysis_*.err.log`

## Task B: Full Reanalysis With Cached Object

Use when:

- scripts changed but Seurat object does not need rebuilding.

Command:

```powershell
powershell -ExecutionPolicy Bypass -File H:\SCI2\OSCC-GSE200996-2025.12\03_rebuild\env\run_rebuild_job.ps1 -Job full -UseCachedObj true
```

Expected output:

- `03_rebuild/results/pre_baseline/`
- `03_rebuild/results/dynamic_paired/`
- `03_rebuild/results/main_story/`

## Task C: Leave-One-Patient Sensitivity

Use when:

- main model or response coding changes;
- reviewer asks about robustness;
- final Extended Data figures need refreshing.

Command:

```powershell
powershell -ExecutionPolicy Bypass -File H:\SCI2\OSCC-GSE200996-2025.12\03_rebuild\env\run_rebuild_job.ps1 -Job loo
```

Expected output:

- `03_rebuild/results/sensitivity_leave_one_out/LOO_model_stability_summary.csv`
- `03_rebuild/results/sensitivity_leave_one_out/LOO_key_pathway_NES.csv`
- `03_rebuild/results/sensitivity_leave_one_out/LOO_stat_correlation.png`
- `03_rebuild/results/sensitivity_leave_one_out/LOO_key_pathway_NES.png`

## Task D: Optional Heavy R Packages

Use when:

- optional sensitivity analyses need `speckle`;
- enrichment visualizations require `clusterProfiler` or `enrichplot`;
- alternative pseudobulk checks require `edgeR`.

Command:

```powershell
powershell -ExecutionPolicy Bypass -File H:\SCI2\OSCC-GSE200996-2025.12\03_rebuild\env\bootstrap_rebuild_env.ps1 -Optional
```

Send back if failed:

- newest `env_install_R_*.out.log`
- newest `env_install_R_*.err.log`

## Task E: Large External Validation Downloads

Use when:

- adding GSE123813 or another validation cohort;
- downloading supplementary matrices or raw files from GEO/SRA.

Download location convention:

```text
00_raw_data/GSE123813_validation/
```

After download, do not rename files unless a manifest is also updated. Tell Codex the exact filenames and source URLs.

## Task F: Final Figure Rendering

Use when:

- high-resolution final panels need to be regenerated after figure script edits.

Expected future command:

```powershell
powershell -ExecutionPolicy Bypass -File H:\SCI2\OSCC-GSE200996-2025.12\03_rebuild\env\run_rebuild_job.ps1 -Job summary
```

This will be updated once final figure scripts are added.
