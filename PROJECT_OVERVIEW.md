# Project overview

## Research question

This project evaluates longitudinal immune-state remodelling during neoadjuvant immunotherapy for oral squamous cell carcinoma and examines whether discovery-cohort signals remain portable across independent public cohorts.

## Analytical design

- Discovery cohort: GSE200996.
- Primary clinical unit: patient.
- Longitudinal unit: strictly paired pretreatment and on-treatment samples where available.
- Discovery layers: immune abundance, lineage-resolved pseudobulk interactions, pathway-level summaries, exact permutation inference, composition sensitivity and influential-patient analysis.
- External boundary cohorts: GSE281729, GSE195832, GSE232240, GSE179730, GSE301741 and GSE123813.
- Reproducibility layers: fresh reconstruction, registered source-table checks, independent recomputation and isolated clean-room implementation.

## Interpretation boundary

The analysis is exploratory and observational. Associations are not interpreted as causal effects. The small strictly paired discovery cohort, response-label heterogeneity and nonsignificant opposite-direction external results limit generalisation. The repository does not support diagnostic, clinically validated biomarker or treatment-selection claims.

## Public release boundary

This repository contains original analysis code, configuration files and reproducibility documentation. It excludes raw data, downloaded archives, large intermediate objects, manuscripts, journal-upload files, local logs and private contact metadata.

The archived code snapshot is available at https://doi.org/10.5281/zenodo.21862445. The DOI identifies a preserved research-code archive; it is not a claim of clinical validation.
