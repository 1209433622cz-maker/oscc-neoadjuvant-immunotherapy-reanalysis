# GSE232240 immune-resolved boundary analysis

## Why this cohort is first

GSE232240 is an independent neoadjuvant HNSCC cohort with matched pretreatment and post-treatment CD45+ single-cell data and deposited patient-level response labels. The downloaded metadata contain 32,399 cells from 18 patients. The frozen primary analysis excludes the source-defined nivolumab-monotherapy patient (Pat04) and requires at least 30 T cells and 30 myeloid cells at both timepoints.

The metadata audit identified 14 patients eligible in both lineages (eight RE and six NR), giving 3,003 exhaustive binary label assignments.

## Download

The complete GEO archive is 61,460,480 bytes (about 58.6 MiB). The command uses IPv4, retries and a `.partial` file, so rerunning resumes an interrupted transfer.

```powershell
powershell -ExecutionPolicy Bypass -File .\03_rebuild\env\download_gse232240_validation.ps1 -Workspace $PWD
```

Do not delete `GSE232240_RAW.tar.partial` after a network interruption. Run the same command again.

## Metadata-only gate

This step is safe while other analyses are running and does not load the count matrix.

```powershell
powershell -ExecutionPolicy Bypass -File .\03_rebuild\env\run_gse232240_frozen_validation.ps1 -Workspace $PWD -MetadataOnly
```

## Full frozen test

After the archive and extracted count file pass integrity checks:

```powershell
powershell -ExecutionPolicy Bypass -File .\03_rebuild\env\run_gse232240_frozen_validation.ps1 -Workspace $PWD
```

The count matrix is streamed one gene at a time. It is not converted to a dense cell-by-gene object, so this route is designed for the current 31.8 GB machine.

## Interpretation lock

Do not alter genes, modules, lineages, cell thresholds or exclusions after viewing the response association. A null, opposite or non-estimable result remains reportable and must not be hidden. Main-text claims may be changed only after the complete gate report is generated.

## Frozen result

The responder-minus-non-responder family effect was -0.230349 (exact two-sided P = 0.110889). The global-PC1-adjusted response effect was -0.201150 (HC3 P = 0.496604), and rank scoring was also negative and null (effect = -0.010128; exact P = 0.456543). T-cell and myeloid lineage effects were both negative. No individual module passed multiplicity correction (minimum FDR = 0.430976), while frozen-gene coverage ranged from 92.96% to 100%.

The observed absolute effect exceeded all 2,000 overlap-preserving matched random families (plus-one empirical P = 0.000500), and all 14 leave-one-patient estimates remained negative. This matched-null extremeness does not establish biological or clinical specificity and cannot override the null primary test or the failure of the prespecified positive orientation, global-PC and rank-direction gates. Two of six gates passed, so the result is classified as `SENSITIVITY_OR_BOUNDARY`, not positive validation.

## Reproducibility gates

```powershell
python .\03_rebuild\analysis\181_audit_validation_upgrade_assets.py --workspace $PWD
python .\03_rebuild\analysis\183_audit_gse232240_content_integration.py --workspace $PWD
powershell -ExecutionPolicy Bypass -File .\03_rebuild\env\run_gse232240_independent_cleanroom.ps1 -Workspace $PWD
```

The frozen asset audit passed 49/49 checks. The isolated implementation was allowed to read only the raw count matrix, raw metadata, frozen configuration and frozen module manifest before unblinding. It independently reimplemented metadata parsing, segmented pseudobulk aggregation, PC1 and HC3 estimation, reproduced all 3,003 exact-null effects and the frozen-seed matched null, and passed 118/118 comparisons. A second matched-null run with an independent seed retained the same extremeness conclusion. Cross-document integration checks then confirmed consistency among the manuscript, Figure 6 and Supplementary Table S19.
