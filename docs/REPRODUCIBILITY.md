# Reproducibility guide

## Supported environment

The frozen analysis was executed on Windows with:

- R 4.6.0;
- Seurat 5.5.0;
- DESeq2 1.52.0;
- limma 3.68.4;
- edgeR 4.10.1;
- fgsea 1.38.0;
- Python 3.13.7 and 3.12.13 with recorded numerical dependencies.

Complete package lists are in `03_rebuild/env/required_R_packages.csv` and `03_rebuild/env/requirements-python.txt`.

## Workflow levels

1. `run_rebuild_job.ps1 -Job audit` checks the GSE200996 discovery inputs.
2. `run_rebuild_job.ps1 -Job full -UseCachedObj false` rebuilds the discovery object and primary analyses from raw matrices.
3. `run_rebuild_job.ps1 -Job loo` runs leave-one-patient diagnostics.
4. `run_clean_room_reproducibility.ps1 -Stage all` executes the registered discovery, external-cohort, figure and comparison stages in an isolated workspace.
5. The cohort-specific wrappers run GSE195832, GSE232240 and GSE301741 analyses independently when a full clean-room run is not required.

## Memory-sensitive route

Do not lower the GSE301741 RDS memory guard on a 32-GB workstation. Use the raw-H5 reconstruction route:

```powershell
powershell -ExecutionPolicy Bypass -File .\03_rebuild\env\prepare_gse301741_raw_route_when_idle.ps1 -Workspace $PWD -Extract
powershell -ExecutionPolicy Bypass -File .\03_rebuild\env\run_gse301741_raw_reconstruction.ps1 -Workspace $PWD
```

The workflow processes one H5 matrix at a time and applies prespecified same-fraction, same-library and target-lineage cell-count gates.

## Frozen GSE232240 route

The independent IMCISION test is intentionally separated from the discovery workflow. Download and run it with:

```powershell
powershell -ExecutionPolicy Bypass -File .\03_rebuild\env\download_gse232240_validation.ps1 -Workspace $PWD
powershell -ExecutionPolicy Bypass -File .\03_rebuild\env\run_gse232240_frozen_validation.ps1 -Workspace $PWD
```

The analysis streams the text count matrix gene by gene and immediately aggregates cells into patient/timepoint/broad-lineage pseudobulks. The primary eligibility, score, exact test, global-PC adjustment, overlap-preserving null and leave-one-patient-out gates are defined in `03_rebuild/config/gse232240_validation.json`.

## Statistical boundaries

- Clinical-response inference uses patients, not cells, as replicates.
- Small discovery-cohort abundance and pathway summaries use exhaustive label assignments where feasible.
- The primary discovery pseudobulk design is `~ patient_id + post + post:resp_num`.
- The external 16-module family and positive orientation were frozen before external-cohort scoring.
- External null, opposite, non-estimable and global-shift-sensitive results are retained rather than filtered.

## Expected local outputs

Generated outputs are written under `03_rebuild/results`, `figures`, `tables`, `validation`, `logs` and `targets`. These locations are ignored by Git. Publication source-data tables should be generated from the registered scripts, not manually edited.
