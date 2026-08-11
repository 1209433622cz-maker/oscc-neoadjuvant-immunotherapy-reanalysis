#!/usr/bin/env python
"""Recover GSE179730 response labels from Table S2 and run exact response tests."""

from __future__ import annotations

import itertools
import json
import math
import os
import re
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests


WORKSPACE = Path(os.environ.get("GSE200996_WORKSPACE", Path.cwd())).resolve()
SOURCE_JSON = (
    WORKSPACE
    / "02_references"
    / "external_supplements"
    / "GSE179730"
    / "PMC8561238_mmc1_bioc.json"
)
VALIDATION_DIR = (
    WORKSPACE / "03_rebuild" / "validation" / "GSE179730_bulk_treatment_direction"
)
PAIRED_PATH = VALIDATION_DIR / "GSE179730_LOCKED_MODULE_PAIRED_DELTAS.csv"
LABEL_OUT = VALIDATION_DIR / "GSE179730_RESPONSE_LABELS_TABLE_S2.csv"
STATS_OUT = VALIDATION_DIR / "GSE179730_LOCKED_MODULE_RESPONSE_EXACT.csv"
REPORT_OUT = VALIDATION_DIR / "GSE179730_RESPONSE_VALIDATION_REPORT.md"
SOURCE_OUT = (
    WORKSPACE
    / "03_rebuild"
    / "figures"
    / "external_validation"
    / "source_data"
    / "GSE179730_response_exact_source.csv"
)


def bh(values: pd.Series) -> np.ndarray:
    result = np.full(len(values), np.nan)
    valid = values.notna().to_numpy()
    if valid.any():
        result[valid] = multipletests(values.loc[valid].astype(float), method="fdr_bh")[1]
    return result


def table_s2_text() -> str:
    payload = json.loads(SOURCE_JSON.read_text(encoding="utf-8"))
    passages = payload[0]["documents"][0]["passages"]
    matches = [
        passage["text"]
        for passage in passages
        if passage.get("text", "").startswith("Table S2. Tumor response")
    ]
    if len(matches) != 1:
        raise RuntimeError(f"Expected one Table S2 passage, found {len(matches)}")
    return matches[0]


def recover_labels() -> pd.DataFrame:
    text = table_s2_text()
    matches = re.findall(r"\bPt(\d{2})\s+(Responder|Progressor|Stable)\b", text)
    if len(matches) != 12:
        raise RuntimeError(f"Expected 12 patient outcomes in Table S2, found {len(matches)}")
    rows = []
    for patient_number, source_outcome in matches:
        response_binary = (
            "non_responder" if source_outcome == "Progressor" else "responder"
        )
        rows.append(
            {
                "article_patient_id": f"Pt{patient_number}",
                "geo_patient_id": f"HN{patient_number.zfill(3)}",
                "source_outcome": source_outcome,
                "response_binary": response_binary,
                "response_definition": (
                    "clinical benefit (Responder or Stable)"
                    if response_binary == "responder"
                    else "progression/no clinical benefit"
                ),
                "source": "PMC8561238 Document S1 Table S2 via NCBI BioC",
                "mapping_rule": "zero-padded numeric patient identifier",
            }
        )
    labels = pd.DataFrame(rows).sort_values("article_patient_id")
    if labels["article_patient_id"].duplicated().any():
        raise RuntimeError("Duplicate Table S2 patient IDs")
    if (labels["response_binary"] == "responder").sum() != 7:
        raise RuntimeError("Table S2 responder count does not match the source article")
    if (labels["response_binary"] == "non_responder").sum() != 5:
        raise RuntimeError("Table S2 non-responder count does not match the source article")
    return labels


def exact_response_models(
    paired: pd.DataFrame, labels: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame]:
    local_patients = sorted(paired["patient_id"].dropna().unique())
    label_map = labels.set_index("geo_patient_id")["response_binary"].to_dict()
    missing = sorted(set(local_patients) - set(label_map))
    if missing:
        raise RuntimeError(f"Local paired patients lack Table S2 labels: {missing}")
    local_labels = labels[labels["geo_patient_id"].isin(local_patients)].copy()
    if len(local_labels) != 11:
        raise RuntimeError(f"Expected 11 RNA-seq paired patients, found {len(local_labels)}")
    if (local_labels["response_binary"] == "responder").sum() != 6:
        raise RuntimeError("Expected six source-defined responders among RNA-seq pairs")
    if (local_labels["response_binary"] == "non_responder").sum() != 5:
        raise RuntimeError("Expected five source-defined non-responders among RNA-seq pairs")

    n_responder = int(
        (local_labels["response_binary"] == "responder").sum()
    )
    assignments = list(itertools.combinations(range(len(local_patients)), n_responder))
    if len(assignments) != math.comb(11, 6):
        raise RuntimeError("Unexpected exact-permutation assignment count")

    rows = []
    for signature, subset in paired.groupby("signature"):
        current = (
            subset.set_index("patient_id")
            .loc[local_patients]
            .reset_index()
            .copy()
        )
        current["response_binary"] = current["patient_id"].map(label_map)
        values = current["delta_post_minus_pre"].astype(float).to_numpy()
        observed_mask = (
            current["response_binary"].to_numpy() == "responder"
        )
        observed = float(
            values[observed_mask].mean() - values[~observed_mask].mean()
        )
        null = []
        for responder_indices in assignments:
            mask = np.zeros(len(local_patients), dtype=bool)
            mask[list(responder_indices)] = True
            null.append(float(values[mask].mean() - values[~mask].mean()))
        null_values = np.asarray(null)
        exact_p = float(
            np.mean(np.abs(null_values) >= abs(observed) - 1e-12)
        )
        responders = values[observed_mask]
        non_responders = values[~observed_mask]
        welch = stats.ttest_ind(
            responders, non_responders, equal_var=False, nan_policy="omit"
        )
        se = math.sqrt(
            responders.var(ddof=1) / len(responders)
            + non_responders.var(ddof=1) / len(non_responders)
        )
        numerator = (
            responders.var(ddof=1) / len(responders)
            + non_responders.var(ddof=1) / len(non_responders)
        ) ** 2
        denominator = (
            (responders.var(ddof=1) / len(responders)) ** 2
            / (len(responders) - 1)
            + (non_responders.var(ddof=1) / len(non_responders)) ** 2
            / (len(non_responders) - 1)
        )
        df = numerator / denominator if denominator > 0 else math.nan
        critical = stats.t.ppf(0.975, df) if pd.notna(df) else math.nan
        rows.append(
            {
                "signature": signature,
                "target_lineage": str(current["target_lineage"].iloc[0]),
                "n_pairs": len(current),
                "n_responder": int(observed_mask.sum()),
                "n_non_responder": int((~observed_mask).sum()),
                "responder_mean_delta": float(responders.mean()),
                "non_responder_mean_delta": float(non_responders.mean()),
                "responder_minus_non_responder": observed,
                "ci95_low": observed - critical * se,
                "ci95_high": observed + critical * se,
                "welch_p": float(welch.pvalue),
                "exact_assignments": len(assignments),
                "exact_p": exact_p,
            }
        )
    output = pd.DataFrame(rows)
    output["welch_fdr"] = bh(output["welch_p"])
    output["exact_fdr"] = bh(output["exact_p"])
    output = output.sort_values(["exact_p", "signature"])
    paired_with_labels = paired.merge(
        local_labels[
            [
                "geo_patient_id",
                "article_patient_id",
                "source_outcome",
                "response_binary",
            ]
        ],
        left_on="patient_id",
        right_on="geo_patient_id",
        how="left",
        validate="many_to_one",
    )
    return output, paired_with_labels


def write_report(labels: pd.DataFrame, stats_df: pd.DataFrame) -> None:
    local = labels[labels["geo_patient_id"].isin(
        pd.read_csv(PAIRED_PATH)["patient_id"].unique()
    )]
    lines = [
        "# GSE179730 Table-S2 Response Validation Report",
        "",
        f"Generated: {datetime.now().astimezone().isoformat(timespec='seconds')}",
        "",
        "## Source Recovery",
        "",
        f"- Source JSON: `{SOURCE_JSON}`",
        "- NCBI BioC passage: PMC8561238 Document S1, Table S2.",
        "- Table S2 contains 12 patient outcomes: seven source-defined responders/clinical-benefit cases and five progressors/non-responders.",
        "- Local RNA-seq identifiers `HN001`, `HN002`, etc. were mapped to article identifiers `Pt01`, `Pt02`, etc. by the shared zero-padded numeric identifier.",
        "- All 11 paired RNA-seq patients mapped to Table S2; Pt06 is the sole Table-S2 responder without a local paired RNA-seq column.",
        f"- RNA-seq response groups: {(local['response_binary'] == 'responder').sum()} responders and {(local['response_binary'] == 'non_responder').sum()} non-responders.",
        "",
        "## Statistical Design",
        "",
        "- Locked module scores and paired post-minus-pre deltas were reused without gene reselection.",
        "- The primary effect is responder minus non-responder mean paired delta.",
        "- All 462 assignments of six responder labels among 11 patients were enumerated for two-sided exact inference.",
        "- Benjamini-Hochberg FDR was calculated across all 16 locked modules.",
        "",
        "## Findings",
        "",
        f"- Positive responder-minus-non-responder deltas: {(stats_df['responder_minus_non_responder'] > 0).sum()}/{len(stats_df)} locked modules.",
        f"- Minimum exact P: {stats_df['exact_p'].min():.4g}.",
        f"- Minimum exact FDR: {stats_df['exact_fdr'].min():.4g}.",
        f"- Exact FDR < 0.05: {(stats_df['exact_fdr'] < 0.05).sum()}/{len(stats_df)} modules.",
        "",
        "## Interpretation",
        "",
        "GSE179730 can now be used as an independent OCSCC bulk response-association check under the source article's broad clinical-benefit definition. The locked module family is directionally higher in responders, especially for myeloid/complement/mTORC1 modules, but no module survives exhaustive patient-label permutation with FDR correction.",
        "This cohort therefore supplies response-annotated directional support and a null exact-inference boundary. It must not be presented as positive statistical validation or as pathological-response-depth replication.",
        "",
        "## Outputs",
        "",
        f"- Frozen response labels: `{LABEL_OUT}`",
        f"- Exact response models: `{STATS_OUT}`",
        f"- Figure source data: `{SOURCE_OUT}`",
    ]
    REPORT_OUT.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    labels = recover_labels()
    paired = pd.read_csv(PAIRED_PATH)
    stats_df, paired_with_labels = exact_response_models(paired, labels)
    LABEL_OUT.parent.mkdir(parents=True, exist_ok=True)
    SOURCE_OUT.parent.mkdir(parents=True, exist_ok=True)
    labels.to_csv(LABEL_OUT, index=False)
    stats_df.to_csv(STATS_OUT, index=False)
    stats_df.to_csv(SOURCE_OUT, index=False)
    paired_with_labels.to_csv(
        VALIDATION_DIR / "GSE179730_LOCKED_MODULE_PAIRED_DELTAS_WITH_RESPONSE.csv",
        index=False,
    )
    write_report(labels, stats_df)
    print(REPORT_OUT)
    print(STATS_OUT)


if __name__ == "__main__":
    main()
