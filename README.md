# Immune-state remodelling in OSCC immunotherapy

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21862445.svg)](https://doi.org/10.5281/zenodo.21862445)

Reproducible patient-level analysis of longitudinal immune-state remodelling and cross-cohort portability in neoadjuvant oral squamous cell carcinoma (OSCC) immunotherapy.

This repository provides the analysis code, prespecified configuration files, environment setup and reproducibility documentation supporting the study:

> Immune-state remodelling and portability in neoadjuvant oral squamous cell carcinoma immunotherapy

## Study overview

The study asks whether longitudinal immune-state changes observed during neoadjuvant immunotherapy are associated with pathological response and whether those signals remain portable across independent public cohorts.

The discovery analysis uses GSE200996 and treats the patient as the unit of clinical-response inference. It combines broad immune-abundance summaries, lineage-resolved pseudobulk models, exact permutation tests, composition sensitivity analyses and influential-patient diagnostics. External analyses use GSE281729, GSE195832, GSE232240, GSE179730, GSE301741 and GSE123813 to evaluate response association, matched-null behaviour and cross-cohort portability.

The study is an exploratory secondary analysis of public, de-identified datasets. It does not provide a diagnostic classifier, establish causality, validate a clinical biomarker or support treatment-selection claims.

## Evidence structure

- Discovery inference is based on six strictly paired patients, including one High responder.
- Patient-level models and exact permutation procedures are used to avoid treating cells as independent clinical replicates.
- GSE195832 and GSE232240 provide external boundary tests with opposite, nonsignificant directions.
- Null, opposite-direction and non-estimable results are retained as part of the evidence record.
- Independent recomputation and clean-room analyses test numerical reproducibility without being represented as independent clinical adjudication.

## Repository layout

- `03_rebuild/analysis/`: R and Python analysis, figure-generation and audit scripts.
- `03_rebuild/config/`: analysis specifications, response harmonisation rules and data manifests.
- `03_rebuild/env/`: environment bootstrap, download and workflow orchestration scripts.
- `docs/DATA_AVAILABILITY.md`: public datasets and acquisition routes.
- `docs/STUDY_DESIGN.md`: cohort roles, estimands and inferential boundaries.
- `docs/REPRODUCIBILITY.md`: environment and execution guidance.
- `docs/EVIDENCE_STRENGTH.md`: evidence grading and permitted claims.
- `docs/GSE232240_VALIDATION.md`: independent immune-resolved boundary analysis.

Raw data, downloaded GEO/SRA archives, Seurat objects, large matrices, intermediate results, local logs, manuscripts and journal-upload files are intentionally excluded.

## Reproduction

The workflows were developed and audited on Windows PowerShell with R and Python. Commands below assume the repository is the workspace root.

```powershell
# Install and audit the required environment.
powershell -ExecutionPolicy Bypass -File .\03_rebuild\env\bootstrap_rebuild_env.ps1 -Workspace $PWD -Optional

# Download registered public inputs.
powershell -ExecutionPolicy Bypass -File .\03_rebuild\env\download_external_datasets.ps1 -Workspace $PWD

# Reconstruct the discovery analysis from raw inputs.
powershell -ExecutionPolicy Bypass -File .\03_rebuild\env\run_rebuild_job.ps1 -Workspace $PWD -Job full -UseCachedObj false

# Run discovery-cohort sensitivity analysis.
powershell -ExecutionPolicy Bypass -File .\03_rebuild\env\run_discovery_cohort_sensitivity_when_idle.ps1 -Workspace $PWD -MinimumFreeMemoryGB 18

# Recompute the GSE195832 validation layers.
powershell -ExecutionPolicy Bypass -File .\03_rebuild\env\run_gse195832_independent_recompute.ps1 -Workspace $PWD
powershell -ExecutionPolicy Bypass -File .\03_rebuild\env\run_gse195832_locked_family_validation.ps1 -Workspace $PWD

# Run the GSE232240 analysis and its isolated implementation.
powershell -ExecutionPolicy Bypass -File .\03_rebuild\env\download_gse232240_validation.ps1 -Workspace $PWD
powershell -ExecutionPolicy Bypass -File .\03_rebuild\env\run_gse232240_frozen_validation.ps1 -Workspace $PWD
powershell -ExecutionPolicy Bypass -File .\03_rebuild\env\run_gse232240_independent_cleanroom.ps1 -Workspace $PWD
```

Some GSE301741 reconstruction steps are memory guarded and process raw H5 matrices one sample at a time. Consult `docs/REPRODUCIBILITY.md` before running the complete workflow.

## Reproducibility status

- Cached-object and fresh discovery reconstructions completed successfully.
- The isolated clean-room workflow reproduced 185 registered CSV outputs within the prespecified tolerances.
- The independent GSE195832 recomputation reproduced all nine locked primary checks.
- The GSE232240 isolated implementation reproduced all registered model and matched-null layers.
- Figure-source, denominator, statistical-boundary and manuscript-number audits completed without unresolved discrepancies.

These checks establish computational reproducibility within the registered workflow. They do not replace external clinical review or independent statistical sign-off.

## Data and code availability

All analysed datasets are public. Accession numbers and acquisition routes are listed in `docs/DATA_AVAILABILITY.md`. Users must obtain source data from the original repositories and comply with their terms.

The archived code snapshot is available at [10.5281/zenodo.21862445](https://doi.org/10.5281/zenodo.21862445). Citation metadata are provided in `CITATION.cff`.

## Author

- Author and corresponding author: Zhi Chen
- Affiliation: School of Medicine, The Chinese University of Hong Kong, Shenzhen, China
- ORCID: [0009-0001-0072-5576](https://orcid.org/0009-0001-0072-5576)
- Email: zhichen1@link.cuhk.edu.cn

## License

Original analysis code is released under the MIT License. External datasets, third-party dependencies and source-publication materials remain subject to their own licences and terms.
