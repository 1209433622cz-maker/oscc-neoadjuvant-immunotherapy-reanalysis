# Discovery Cohort-Sensitivity Run Guide

Created: 2026-07-26

## Purpose

The six paired GSE200996 patients contain three Mono and three Combo patients, and response depth is partially associated with treatment cohort. The primary pseudobulk model contains a patient effect and a response-dependent treatment interaction but does not contain a treatment-cohort-dependent change term. This job recomputes T-cell and myeloid pseudobulk models with `post:cohort` included.

## Memory gate

The cached discovery object is approximately 1.17 GB on disk and expands in memory. Run only when at least 18 GB physical memory is free. The wrapper stops before loading the object when this condition is not met.

## PowerShell command

```powershell
powershell -ExecutionPolicy Bypass -File H:\SCI2\OSCC-GSE200996-2025.12\03_rebuild\env\run_discovery_cohort_sensitivity_when_idle.ps1 -MinimumFreeMemoryGB 18
```

Do not use `-MinimumFreeMemoryGB 18` for the 12.35 GB GSE301741 RDS. That separate object is expected to require approximately 40-60 GB working memory and is not suitable for this 31.8 GB machine.

## Models

- Original: `~ patient_id + post + post:resp_num`
- Cohort sensitivity: `~ patient_id + post + post:cohort + post:resp_num`
- Exact pathway sensitivity: all 18 unique response-label assignments obtained by permuting labels within the Mono and Combo strata.

## Expected outputs

- `03_rebuild/results/sensitivity_cohort_adjusted_pseudobulk/`
- `03_rebuild/figures/submission/ExtendedData9_submission_discovery_cohort_sensitivity.png`
- `03_rebuild/figures/submission/ExtendedData9_submission_discovery_cohort_sensitivity.pdf`
- `03_rebuild/figures/submission/ExtendedData9_submission_discovery_cohort_sensitivity.svg`
- `03_rebuild/logs/43_discovery_cohort_sensitivity_*.out.log`
- `03_rebuild/logs/43_discovery_cohort_sensitivity_*.err.log`

## Interpretation gate

The mechanism may remain the primary manuscript result only if:

1. T-cell and myeloid genome-wide Wald-statistic ranks remain directionally concordant after adding `post:cohort`.
2. The principal TNFA/NF-kB, interferon and mTORC1 pathway directions remain positive.
3. No single pathway is promoted solely because of the resolution-limited 18-assignment exact test.
4. Any loss of adjusted significance is reported as a limitation rather than hidden by retaining only the original model.

After the job completes, rerun the registered integration and audit checks before describing the analysis as content locked.
