# Evidence strength and reporting boundary

Date: 2026-08-04

## Evidence hierarchy

The frozen analyses support four distinct levels of evidence that must not be conflated:

1. **Computational reproducibility: high.** The discovery clean-room workflow reproduced 185 of 185 registered CSV outputs within the frozen tolerances, the independent GSE195832 route reproduced all nine locked primary checks, and a separately implemented GSE232240 route passed all 118 unblinded comparisons from isolated raw/configuration inputs.
2. **Internal biological association: moderate to low.** In six strictly matched discovery patients, response-associated rankings concentrated in longitudinal T-cell and myeloid state programs rather than broad immune abundance. Interpretation is limited by one High responder, exact patient-level nulls, composition sensitivity and influential-patient diagnostics.
3. **External portability: low and currently unsupported.** GSE281729 was sensitive to treatment timing and global expression structure; GSE195832 and the immune-resolved GSE232240 cohort were directionally opposite and nonsignificant after locked controls; GSE179730 was exactly null; and GSE301741 was negative or non-estimable at the available response-label boundary.
4. **Clinical prediction or mechanism: unsupported.** No locked classifier, decision threshold, calibration analysis, prospective validation or causal experiment was performed.

Computational reproducibility shows that the same inputs and code recover the same results. It does not establish biological transportability, causal mechanism or clinical utility.

## Permitted interpretation

The manuscript may state that longitudinal response-associated signal in the discovery cohort was more evident in lineage-resolved transcriptional state than in broad immune abundance. It may describe the implicated programs as candidates for future longitudinal and mechanistic testing. It must also state that patient-level exact tests, the sole High responder and external cohorts limit pathway specificity and portability.

The following claims are not supported:

- that the programs predict OSCC immunotherapy response;
- that an external cohort validated a clinical biomarker;
- that interferon or myeloid-secretory mechanisms were established;
- that the results can select patients, treatments or surgery timing; or
- that a large single-cell count compensates for the six-patient inferential sample.

## Result-presentation integrity

Results must not be selected according to whether they favor the central claim. The following material findings remain part of the main-text evidentiary boundary:

- six matched discovery patients, including one High responder;
- broad-abundance and patient-level exact-test null results;
- the 12 frozen pathway-family exact-test results;
- P32 leave-one-patient influence diagnostics;
- GSE281729 timing and global-PC sensitivity;
- the opposite or null GSE195832 and GSE179730 results; and
- the negative or non-estimable GSE301741 boundary.

Detailed coefficient tables, scoring variants, random-family distributions, full MCP-counter outputs, per-patient leave-one-out diagnostics, reconstruction QC and complete permutation distributions may be summarized in the main text and retained in supplementary materials. This changes presentation density, not evidentiary meaning.

Contextual analyses may be omitted only as complete prespecified families when they are remote from the central OSCC/HNSCC response question. Favorable fragments from an omitted family must not be retained. Analyses that influenced gene, pathway, score, direction or conclusion selection cannot be removed after the fact.

## Legitimate routes to stronger evidence

Evidence strength can increase only through additional independent information or stronger design, for example:

- a larger fixed-treatment, fixed-sampling-time longitudinal cohort;
- prospectively frozen signatures and analysis plans;
- blinded pathological-response adjudication;
- spatial, protein or TCR measurements linked to matched patients;
- an independent analyst reproducing the registered workflow; and
- functional experiments testing the nominated T-cell and myeloid programs.

Suppressing unfavorable results would increase selective-reporting risk without increasing evidence strength.

## Six-figure presentation

The registered six-main-figure rebuild may prioritize the strongest discovery programs visually, but it does not permit a positive-only narrative. Exact-test limits, influential-patient behavior and independent opposite/null/non-estimable results remain visible in the principal evidentiary sequence. Repetitive diagnostics may move to supplementary materials. The full mapping is recorded in `SIX_FIGURE_REBUILD.md`.
