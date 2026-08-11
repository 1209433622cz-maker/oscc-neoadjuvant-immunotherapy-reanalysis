#!/usr/bin/env python
"""Audit GSE281729 response associations against response-adaptive exposure timing."""

from __future__ import annotations

import math
import os
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests


WORKSPACE = Path(os.environ.get("GSE200996_WORKSPACE", Path.cwd())).resolve()
VALIDATION_DIR = WORKSPACE / "03_rebuild" / "validation" / "GSE281729_bulk_module_validation"
PAIRED_PATH = VALIDATION_DIR / "GSE281729_PAIRED_MODULE_DELTAS_ROBUSTNESS.csv"
ANNOTATION_PATH = VALIDATION_DIR / "GSE281729_EMBEDDED_SAMPLE_ANNOTATION.csv"
MODEL_OUT = VALIDATION_DIR / "GSE281729_RESPONSE_ADAPTIVE_TIMING_MODELS.csv"
STRATIFIED_OUT = VALIDATION_DIR / "GSE281729_RESPONSE_ADAPTIVE_TIMING_STRATIFIED.csv"
PATIENT_OUT = VALIDATION_DIR / "GSE281729_RESPONSE_ADAPTIVE_TIMING_PATIENTS.csv"
REPORT_OUT = VALIDATION_DIR / "GSE281729_RESPONSE_ADAPTIVE_TIMING_REPORT.md"
SOURCE_OUT = (
    WORKSPACE
    / "03_rebuild"
    / "figures"
    / "submission"
    / "source_data"
    / "ExtendedData7_response_adaptive_timing_source.csv"
)

SELECTED = [
    "M_LE_INTERFERON_ALPHA_RESPONSE",
    "T_LE_INTERFERON_ALPHA_RESPONSE",
    "M_LE_INTERFERON_GAMMA_RESPONSE",
    "M_LE_union_core",
    "T_LE_INTERFERON_GAMMA_RESPONSE",
    "M_LE_MTORC1_SIGNALING",
    "T_LE_MTORC1_SIGNALING",
    "T_LE_union_core",
]


def bh(values: pd.Series) -> np.ndarray:
    result = np.full(len(values), np.nan)
    valid = values.notna().to_numpy()
    if valid.any():
        result[valid] = multipletests(values.loc[valid].astype(float), method="fdr_bh")[1]
    return result


def load_data() -> tuple[pd.DataFrame, pd.DataFrame]:
    paired = pd.read_csv(PAIRED_PATH)
    annotation = pd.read_csv(ANNOTATION_PATH)
    patients = annotation[
        [
            "patient_id",
            "doses",
            "primary_response_percent",
            "response_harmonized_ordinal",
            "response_ord_num",
            "hpv",
            "second_drug",
        ]
    ].drop_duplicates("patient_id")
    if patients["patient_id"].duplicated().any():
        raise RuntimeError("Patient metadata are not unique")
    paired = paired.drop(columns=["doses", "primary_response_percent"], errors="ignore")
    paired = paired.merge(
        patients[["patient_id", "doses", "primary_response_percent"]],
        on="patient_id",
        how="left",
        validate="many_to_one",
    )
    paired["doses"] = pd.to_numeric(paired["doses"], errors="coerce")
    paired["primary_response_percent"] = pd.to_numeric(
        paired["primary_response_percent"], errors="coerce"
    )
    if paired[["doses", "primary_response_percent"]].isna().any().any():
        raise RuntimeError("Exposure or primary pathologic-response metadata are missing")
    return paired, patients


def categorical_parts(data: pd.DataFrame, columns: list[str]) -> list[pd.DataFrame]:
    parts = []
    for column in columns:
        if data[column].nunique() > 1:
            parts.append(
                pd.get_dummies(
                    data[column].astype(str),
                    prefix=column,
                    drop_first=True,
                    dtype=float,
                )
            )
    return parts


def fit(
    data: pd.DataFrame,
    response_variable: str,
    categorical_covariates: list[str],
    numeric_covariates: list[str],
    interaction: bool = False,
) -> dict[str, float | int]:
    keep = [
        "post_minus_pre",
        response_variable,
        *categorical_covariates,
        *numeric_covariates,
    ]
    model_data = data[keep].replace([np.inf, -np.inf], np.nan).dropna().copy()
    if response_variable == "primary_response_percent":
        model_data["pTR_per_10_percent"] = (
            model_data["primary_response_percent"].astype(float) / 10.0
        )
        response_term = "pTR_per_10_percent"
    else:
        response_term = response_variable

    parts = [model_data[[response_term]].astype(float)]
    parts.extend(categorical_parts(model_data, categorical_covariates))
    for covariate in numeric_covariates:
        parts.append(model_data[[covariate]].astype(float))
    if interaction:
        model_data["dose2"] = (model_data["doses"].astype(float) == 2).astype(float)
        model_data["response_x_dose2"] = (
            model_data[response_term].astype(float) * model_data["dose2"]
        )
        parts = [
            model_data[[response_term, "dose2", "response_x_dose2"]].astype(float)
        ]

    design = sm.add_constant(pd.concat(parts, axis=1), has_constant="add")
    model = sm.OLS(
        model_data["post_minus_pre"].astype(float), design.astype(float)
    ).fit(cov_type="HC3")
    term = "response_x_dose2" if interaction else response_term
    ci = model.conf_int(alpha=0.05).loc[term]
    return {
        "n": int(len(model_data)),
        "coef": float(model.params[term]),
        "std_error": float(model.bse[term]),
        "ci95_low": float(ci.iloc[0]),
        "ci95_high": float(ci.iloc[1]),
        "p_value": float(model.pvalues[term]),
        "condition_number": float(np.linalg.cond(design.astype(float))),
    }


def run_models(paired: pd.DataFrame) -> pd.DataFrame:
    specs = [
        (
            "ordinal_hpv_second_drug_HC3",
            "response_ord_num",
            ["hpv", "second_drug"],
            [],
            False,
        ),
        (
            "ordinal_doses_HC3",
            "response_ord_num",
            [],
            ["doses"],
            False,
        ),
        (
            "ordinal_hpv_second_drug_doses_HC3",
            "response_ord_num",
            ["hpv", "second_drug"],
            ["doses"],
            False,
        ),
        (
            "continuous_pTR_hpv_second_drug_doses_HC3",
            "primary_response_percent",
            ["hpv", "second_drug"],
            ["doses"],
            False,
        ),
        (
            "ordinal_response_by_dose_interaction_HC3",
            "response_ord_num",
            [],
            ["doses"],
            True,
        ),
    ]
    rows = []
    for signature, subset in paired.groupby("signature"):
        lineage = str(subset["target_lineage"].iloc[0])
        for model_name, response_variable, categorical, numeric, interaction in specs:
            result = fit(
                subset,
                response_variable=response_variable,
                categorical_covariates=categorical,
                numeric_covariates=numeric,
                interaction=interaction,
            )
            rows.append(
                {
                    "signature": signature,
                    "target_lineage": lineage,
                    "model": model_name,
                    "coefficient_term": (
                        "response_x_dose2"
                        if interaction
                        else (
                            "pTR_per_10_percent"
                            if response_variable == "primary_response_percent"
                            else response_variable
                        )
                    ),
                    **result,
                }
            )
    output = pd.DataFrame(rows)
    output["fdr"] = np.nan
    for model_name in output["model"].unique():
        mask = output["model"] == model_name
        output.loc[mask, "fdr"] = bh(output.loc[mask, "p_value"])
    return output


def run_stratified(paired: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for signature, subset in paired.groupby("signature"):
        lineage = str(subset["target_lineage"].iloc[0])
        for dose_group in [1, 2]:
            current = subset[subset["doses"] == dose_group].copy()
            result = fit(
                current,
                response_variable="response_ord_num",
                categorical_covariates=[],
                numeric_covariates=[],
            )
            rows.append(
                {
                    "signature": signature,
                    "target_lineage": lineage,
                    "doses": dose_group,
                    **result,
                }
            )
    output = pd.DataFrame(rows)
    output["fdr"] = np.nan
    for dose_group in output["doses"].unique():
        mask = output["doses"] == dose_group
        output.loc[mask, "fdr"] = bh(output.loc[mask, "p_value"])
    return output


def write_report(
    paired: pd.DataFrame,
    patients: pd.DataFrame,
    models: pd.DataFrame,
    stratified: pd.DataFrame,
) -> None:
    paired_patients = (
        paired[
            [
                "patient_id",
                "doses",
                "response_harmonized_ordinal",
                "response_ord_num",
                "hpv",
                "second_drug",
                "primary_response_percent",
            ]
        ]
        .drop_duplicates("patient_id")
        .sort_values("patient_id")
    )
    selected = models[models["signature"].isin(SELECTED)].copy()
    full = selected[
        selected["model"] == "ordinal_hpv_second_drug_doses_HC3"
    ]
    continuous = selected[
        selected["model"] == "continuous_pTR_hpv_second_drug_doses_HC3"
    ]
    interaction = selected[
        selected["model"] == "ordinal_response_by_dose_interaction_HC3"
    ]
    dose_table = (
        paired_patients.groupby(["doses", "response_harmonized_ordinal"])
        .size()
        .unstack(fill_value=0)
    )
    response_levels = ["Low", "Medium", "High"]
    dose_lines = [
        "| Author-provided doses | Low | Medium | High | Total |",
        "|---:|---:|---:|---:|---:|",
    ]
    for dose_group in sorted(dose_table.index):
        counts = [
            int(dose_table.loc[dose_group, level])
            if level in dose_table.columns
            else 0
            for level in response_levels
        ]
        dose_lines.append(
            f"| {int(dose_group)} | {counts[0]} | {counts[1]} | {counts[2]} | {sum(counts)} |"
        )

    lines = [
        "# GSE281729 Response-Adaptive Timing Sensitivity Report",
        "",
        f"Generated: {datetime.now().astimezone().isoformat(timespec='seconds')}",
        "",
        "## Design Issue",
        "",
        "The source trial used response-adaptive surgical timing. The processed expression annotation contains an author-provided `doses` field with values 1 or 2. This field jointly reflects additional treatment exposure and a later surgical sampling time, and it is associated with response category. It is therefore not treated as an ordinary baseline confounder; it is used as a design-sensitivity variable.",
        "",
        "## Analysis Boundary",
        "",
        f"- Response-annotated paired patients: {paired_patients['patient_id'].nunique()}.",
        f"- One-dose/timing group: {int((paired_patients['doses'] == 1).sum())}.",
        f"- Two-dose/timing group: {int((paired_patients['doses'] == 2).sum())}.",
        "- HC3 covariance was used for all sensitivity models.",
        "- FDR was calculated across all 16 locked modules within each model family.",
        "",
        "## Dose/Timing By Response",
        "",
        *dose_lines,
        "",
        "## Leading-Module Findings",
        "",
        f"- Full HPV, second-drug and dose/timing-adjusted ordinal slopes were negative for {int((full['coef'] < 0).sum())}/{len(full)} displayed modules.",
        f"- Full adjusted ordinal FDR < 0.05 was retained for {int((full['fdr'] < 0.05).sum())}/{len(full)} displayed modules; the minimum displayed-module FDR was {full['fdr'].min():.4f}.",
        f"- Continuous primary-pTR slopes were negative for {int((continuous['coef'] < 0).sum())}/{len(continuous)} displayed modules, but FDR < 0.05 was retained for {int((continuous['fdr'] < 0.05).sum())}/{len(continuous)}.",
        f"- Nominal response-by-dose/timing interaction P < 0.05 occurred for {int((interaction['p_value'] < 0.05).sum())}/{len(interaction)} displayed modules; none should be treated as confirmed without model-family FDR.",
        "",
        "## Interpretation",
        "",
        "The negative ordinal response slopes are not removed by adding the author-provided exposure/timing field, so the primary association is not explained solely by the one-versus-two-dose design indicator.",
        "However, within-group estimates are stronger in the two-dose/timing stratum, and continuous pTR sensitivity is weaker after full adjustment. Because the extra exposure and surgical timing were assigned adaptively after early radiographic assessment, these models cannot identify a treatment-independent causal response effect.",
        "GSE281729 should therefore be described as an HNSCC bulk response-association cohort with response-adaptive timing sensitivity, not as an independent validation of a universal longitudinal direction.",
        "",
        "## Outputs",
        "",
        f"- Patient design table: `{PATIENT_OUT}`",
        f"- HC3 timing-sensitivity models: `{MODEL_OUT}`",
        f"- Dose/timing-stratified models: `{STRATIFIED_OUT}`",
        f"- Extended Data source table: `{SOURCE_OUT}`",
    ]
    REPORT_OUT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    paired_patients.to_csv(PATIENT_OUT, index=False)


def main() -> None:
    paired, patients = load_data()
    models = run_models(paired)
    stratified = run_stratified(paired)
    MODEL_OUT.parent.mkdir(parents=True, exist_ok=True)
    SOURCE_OUT.parent.mkdir(parents=True, exist_ok=True)
    models.to_csv(MODEL_OUT, index=False)
    stratified.to_csv(STRATIFIED_OUT, index=False)
    models[models["signature"].isin(SELECTED)].to_csv(SOURCE_OUT, index=False)
    write_report(paired, patients, models, stratified)
    print(REPORT_OUT)
    print(MODEL_OUT)


if __name__ == "__main__":
    main()
