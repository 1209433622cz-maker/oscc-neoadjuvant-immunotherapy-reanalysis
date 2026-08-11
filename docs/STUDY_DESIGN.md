# Study design

## Objective

To characterise longitudinal immune-state changes during neoadjuvant immunotherapy for oral squamous cell carcinoma and to test whether response-associated discovery signals remain portable across independent public cohorts.

## Cohort roles

| Cohort | Role in the study | Inferential boundary |
| --- | --- | --- |
| GSE200996 | Discovery single-cell cohort | Patient-level inference is limited to six strictly paired patients, including one High responder. |
| GSE281729 | External bulk response analysis | Independent response association; not a replication of the single-cell estimand. |
| GSE195832 | External response-boundary analysis | Opposite, nonsignificant direction retained as boundary evidence. |
| GSE232240 | Independent immune-resolved longitudinal analysis | Fourteen eligible paired patients; opposite, nonsignificant family-level direction. |
| GSE179730 | External portability context | Treatment and sampling differences limit direct effect comparison. |
| GSE301741 | Single-cell portability context | Response-label and treatment heterogeneity are preserved in interpretation. |
| GSE123813 | External biological context | Used as contextual evidence rather than clinical validation. |

## Statistical unit and estimands

The patient is the clinical unit of inference. Cell-level observations are aggregated or modelled within lineage before patient-level contrasts are evaluated. Exact permutation inference is used where the small number of strict pairs makes asymptotic assumptions unreliable. Composition sensitivity and leave-one-patient-out analyses assess whether findings are dependent on lineage abundance or influential individuals.

## Response labels

Response groups follow the source-cohort metadata and documented harmonisation rules. The study does not claim independent pathological reassessment, RECIST adjudication or blinded clinical endpoint review.

## Claim boundary

The results support an exploratory account of heterogeneous immune remodelling and limited cross-cohort portability. They do not establish causality, diagnostic accuracy, prospective clinical utility, a validated biomarker or a treatment-selection rule.
