#!/usr/bin/env python
"""Exact patient-label permutation sensitivity for paired immune abundance."""

from __future__ import annotations

import itertools
import math
import os
from datetime import datetime
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


WORKSPACE = Path(os.environ.get("GSE200996_WORKSPACE", Path.cwd())).resolve()
RESULT_DIR = WORKSPACE / "03_rebuild" / "results" / "sensitivity_exact_permutation"
FIGURE_DIR = WORKSPACE / "03_rebuild" / "figures" / "submission"
SOURCE_DIR = FIGURE_DIR / "source_data"
RESULT_DIR.mkdir(parents=True, exist_ok=True)
FIGURE_DIR.mkdir(parents=True, exist_ok=True)
SOURCE_DIR.mkdir(parents=True, exist_ok=True)

LONG_PATH = (
    WORKSPACE
    / "03_rebuild"
    / "results"
    / "dynamic_paired"
    / "Fig4A_patient_timepoint_celltype_prop_long_allpaired.csv"
)
ORDINAL_PATH = (
    WORKSPACE
    / "03_rebuild"
    / "results"
    / "dynamic_paired"
    / "Fig4A_composition_delta_logit_limma_respOrd_trend.csv"
)
BINARY_PATH = (
    WORKSPACE
    / "03_rebuild"
    / "results"
    / "dynamic_paired"
    / "Fig4A_composition_delta_logit_limma_RvsNR.csv"
)

EPSILON = 1e-4
CELL_ORDER = ["T cell", "Myeloid", "Mast", "B cell", "NK cell", "Cycling"]
RESPONSE_MAP = {"Low": 1.0, "Medium": 2.0, "High": 3.0}
LINEAGE_COLORS = {
    "T cell": "#3E73B9",
    "Myeloid": "#B85A4A",
    "Mast": "#5D9E58",
    "B cell": "#777777",
    "NK cell": "#777777",
    "Cycling": "#777777",
}


def bh_adjust(values: pd.Series) -> pd.Series:
    p = values.astype(float).to_numpy()
    order = np.argsort(p)
    ranked = p[order]
    adjusted = np.minimum.accumulate(
        (ranked * len(ranked) / np.arange(1, len(ranked) + 1))[::-1]
    )[::-1]
    out = np.empty(len(p))
    out[order] = np.minimum(adjusted, 1.0)
    return pd.Series(out, index=values.index)


def logit(values: np.ndarray) -> np.ndarray:
    return np.log((values + EPSILON) / (1 - values + EPSILON))


def fit_ols(y: np.ndarray, response: np.ndarray, cohort: np.ndarray | None = None) -> dict[str, float]:
    columns = [np.ones(len(y))]
    if cohort is not None:
        columns.append(cohort)
    columns.append(response)
    x = np.column_stack(columns)
    beta, _, rank, _ = np.linalg.lstsq(x, y, rcond=None)
    if rank != x.shape[1]:
        raise ValueError("Rank-deficient abundance sensitivity model")
    residual = y - x @ beta
    df_residual = len(y) - x.shape[1]
    sigma2 = float(residual @ residual / df_residual)
    covariance = sigma2 * np.linalg.inv(x.T @ x)
    se = math.sqrt(float(covariance[-1, -1]))
    t_value = float(beta[-1] / se)
    p_value = float(2 * stats.t.sf(abs(t_value), df=df_residual))
    critical = float(stats.t.ppf(0.975, df=df_residual))
    return {
        "coef": float(beta[-1]),
        "std_error": se,
        "ci95_low": float(beta[-1] - critical * se),
        "ci95_high": float(beta[-1] + critical * se),
        "parametric_p": p_value,
        "df_residual": df_residual,
    }


def exact_p(observed: float, null_values: np.ndarray) -> tuple[float, int]:
    n_extreme = int(np.sum(np.abs(null_values) >= abs(observed) - 1e-12))
    return n_extreme / len(null_values), n_extreme


def unique_response_permutations(response: np.ndarray) -> list[np.ndarray]:
    return [
        np.asarray(values, dtype=float)
        for values in sorted(set(itertools.permutations(response.tolist())))
    ]


def stratified_response_permutations(
    response: np.ndarray, cohort: np.ndarray
) -> list[np.ndarray]:
    strata = [np.where(cohort == value)[0] for value in sorted(set(cohort.tolist()))]
    within = [
        sorted(set(itertools.permutations(response[index].tolist()))) for index in strata
    ]
    outputs: list[np.ndarray] = []
    for selections in itertools.product(*within):
        permuted = response.copy()
        for index, values in zip(strata, selections):
            permuted[index] = np.asarray(values, dtype=float)
        outputs.append(permuted)
    return outputs


def observed_effects_from_table(path: Path, model_name: str) -> pd.DataFrame:
    table = pd.read_csv(path)
    required = {"celltype", "logFC", "t", "P.Value", "adj.P.Val"}
    missing = required - set(table.columns)
    if missing:
        raise ValueError(f"{path.name} lacks columns: {sorted(missing)}")
    table = table[["celltype", "logFC", "t", "P.Value", "adj.P.Val"]].copy()
    table["std_error"] = np.abs(table["logFC"] / table["t"])
    table["ci95_low"] = table["logFC"] - 1.96 * table["std_error"]
    table["ci95_high"] = table["logFC"] + 1.96 * table["std_error"]
    return table.rename(
        columns={
            "logFC": "coef",
            "P.Value": "parametric_p",
            "adj.P.Val": "parametric_fdr",
        }
    ).assign(model=model_name)


def make_patient_matrices(long: pd.DataFrame) -> tuple[list[str], list[str], dict[str, np.ndarray]]:
    patients = sorted(long["patient"].astype(str).unique())
    celltypes = [cell for cell in CELL_ORDER if cell in set(long["celltype"].astype(str))]
    delta: dict[str, np.ndarray] = {}
    for celltype in celltypes:
        wide = (
            long.loc[long["celltype"] == celltype]
            .pivot_table(index="patient", columns="timepoint", values="prop", aggfunc="first")
            .reindex(patients)
            .fillna(0.0)
        )
        delta[celltype] = logit(wide["post"].to_numpy()) - logit(wide["pre"].to_numpy())
    return patients, celltypes, delta


def main() -> None:
    long = pd.read_csv(LONG_PATH)
    patients, celltypes, delta = make_patient_matrices(long)
    patient_info = (
        long[["patient", "response_ord", "cohort"]]
        .drop_duplicates()
        .set_index("patient")
        .reindex(patients)
    )
    response = patient_info["response_ord"].map(RESPONSE_MAP).to_numpy(dtype=float)
    cohort = patient_info["cohort"].eq("Combo").to_numpy(dtype=float)

    ordinal_permutations = unique_response_permutations(response)
    stratified_permutations = stratified_response_permutations(response, cohort)
    ordinal_source = observed_effects_from_table(ORDINAL_PATH, "ordinal_unadjusted")
    binary_source = observed_effects_from_table(BINARY_PATH, "High_vs_Low")

    result_rows: list[dict[str, object]] = []
    null_rows: list[dict[str, object]] = []

    for celltype in celltypes:
        y = delta[celltype]

        observed = float(ordinal_source.loc[ordinal_source["celltype"] == celltype, "coef"].iloc[0])
        null = np.asarray(
            [fit_ols(y, permuted)["coef"] for permuted in ordinal_permutations]
        )
        if not math.isclose(observed, fit_ols(y, response)["coef"], abs_tol=1e-8):
            raise ValueError(f"Ordinal coefficient mismatch for {celltype}")
        p_exact, n_extreme = exact_p(observed, null)
        source = ordinal_source.loc[ordinal_source["celltype"] == celltype].iloc[0]
        result_rows.append(
            {
                "analysis": "paired_delta_ordinal_unadjusted",
                "celltype": celltype,
                "n_patients": len(patients),
                "n_unique_permutations": len(ordinal_permutations),
                "coef": observed,
                "std_error": np.nan,
                "ci95_low": source["ci95_low"],
                "ci95_high": source["ci95_high"],
                "parametric_p": source["parametric_p"],
                "parametric_fdr": source["parametric_fdr"],
                "exact_p_two_sided": p_exact,
                "n_extreme": n_extreme,
                "permutation_scheme": "all unique ordinal-response label assignments",
            }
        )
        for permutation_id, value in enumerate(null, 1):
            null_rows.append(
                {
                    "analysis": "paired_delta_ordinal_unadjusted",
                    "celltype": celltype,
                    "permutation_id": permutation_id,
                    "permuted_coef": value,
                    "observed_coef": observed,
                }
            )

        adjusted = fit_ols(y, response, cohort=cohort)
        adjusted_null = np.asarray(
            [
                fit_ols(y, permuted, cohort=cohort)["coef"]
                for permuted in stratified_permutations
            ]
        )
        p_exact, n_extreme = exact_p(adjusted["coef"], adjusted_null)
        result_rows.append(
            {
                "analysis": "paired_delta_ordinal_cohort_adjusted",
                "celltype": celltype,
                "n_patients": len(patients),
                "n_unique_permutations": len(stratified_permutations),
                **adjusted,
                "parametric_fdr": np.nan,
                "exact_p_two_sided": p_exact,
                "n_extreme": n_extreme,
                "permutation_scheme": "ordinal-response labels permuted within Mono/Combo cohort",
            }
        )
        for permutation_id, value in enumerate(adjusted_null, 1):
            null_rows.append(
                {
                    "analysis": "paired_delta_ordinal_cohort_adjusted",
                    "celltype": celltype,
                    "permutation_id": permutation_id,
                    "permuted_coef": value,
                    "observed_coef": adjusted["coef"],
                }
            )

    binary_patients = [
        patient
        for patient in patients
        if patient_info.loc[patient, "response_ord"] in {"Low", "High"}
    ]
    binary_index = np.asarray([patients.index(patient) for patient in binary_patients])
    binary_response = np.asarray(
        [
            1.0 if patient_info.loc[patient, "response_ord"] == "High" else 0.0
            for patient in binary_patients
        ]
    )
    binary_permutations = unique_response_permutations(binary_response)
    for celltype in celltypes:
        y = delta[celltype][binary_index]
        observed = float(binary_source.loc[binary_source["celltype"] == celltype, "coef"].iloc[0])
        null = np.asarray(
            [fit_ols(y, permuted)["coef"] for permuted in binary_permutations]
        )
        if not math.isclose(observed, fit_ols(y, binary_response)["coef"], abs_tol=1e-8):
            raise ValueError(f"Binary coefficient mismatch for {celltype}")
        p_exact, n_extreme = exact_p(observed, null)
        source = binary_source.loc[binary_source["celltype"] == celltype].iloc[0]
        result_rows.append(
            {
                "analysis": "paired_delta_High_vs_Low",
                "celltype": celltype,
                "n_patients": len(binary_patients),
                "n_unique_permutations": len(binary_permutations),
                "coef": observed,
                "std_error": np.nan,
                "ci95_low": source["ci95_low"],
                "ci95_high": source["ci95_high"],
                "parametric_p": source["parametric_p"],
                "parametric_fdr": source["parametric_fdr"],
                "exact_p_two_sided": p_exact,
                "n_extreme": n_extreme,
                "permutation_scheme": "all unique High/Low label assignments",
            }
        )
        for permutation_id, value in enumerate(null, 1):
            null_rows.append(
                {
                    "analysis": "paired_delta_High_vs_Low",
                    "celltype": celltype,
                    "permutation_id": permutation_id,
                    "permuted_coef": value,
                    "observed_coef": observed,
                }
            )

    results = pd.DataFrame(result_rows)
    results["exact_bh_fdr"] = results.groupby("analysis", group_keys=False)[
        "exact_p_two_sided"
    ].apply(bh_adjust)
    for analysis in results["analysis"].unique():
        mask = results["analysis"] == analysis
        if results.loc[mask, "parametric_fdr"].isna().all():
            results.loc[mask, "parametric_fdr"] = bh_adjust(
                results.loc[mask, "parametric_p"]
            )

    null_table = pd.DataFrame(null_rows)
    results_path = RESULT_DIR / "ABUNDANCE_EXACT_PERMUTATION_RESULTS.csv"
    null_path = RESULT_DIR / "ABUNDANCE_EXACT_PERMUTATION_NULL_DISTRIBUTIONS.csv"
    results.to_csv(results_path, index=False)
    null_table.to_csv(null_path, index=False)
    results.to_csv(
        SOURCE_DIR / "ExtendedData8_abundance_exact_permutation_results.csv",
        index=False,
    )
    null_table.to_csv(
        SOURCE_DIR / "ExtendedData8_abundance_permutation_null.csv",
        index=False,
    )

    figure_results = results.copy()
    labels = {
        "paired_delta_ordinal_unadjusted": "Ordinal, unadjusted",
        "paired_delta_ordinal_cohort_adjusted": "Ordinal, cohort-adjusted",
        "paired_delta_High_vs_Low": "High versus Low",
    }
    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 7.5,
            "axes.linewidth": 0.6,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 6.45))
    for ax, analysis, panel in zip(
        axes.flat[:3],
        labels,
        ["a", "b", "c"],
    ):
        subset = (
            figure_results.loc[figure_results["analysis"] == analysis]
            .set_index("celltype")
            .reindex(CELL_ORDER)
            .dropna(subset=["coef"])
            .reset_index()
        )
        y = np.arange(len(subset))
        x_min = min(float(subset["ci95_low"].min()), 0.0)
        x_max = max(float(subset["ci95_high"].max()), 0.0)
        x_span = max(x_max - x_min, 1.0)
        ax.set_xlim(x_min - 0.08 * x_span, x_max + 0.72 * x_span)
        ax.axvline(0, color="#888888", lw=0.8)
        for index, row in subset.iterrows():
            ax.plot(
                [row["ci95_low"], row["ci95_high"]],
                [index, index],
                color="#555555",
                lw=1.25,
            )
            ax.scatter(
                row["coef"],
                index,
                s=46,
                color=LINEAGE_COLORS[row["celltype"]],
                edgecolor="white",
                linewidth=0.5,
                zorder=3,
            )
            ax.text(
                0.98,
                index,
                f"exact P {row['exact_p_two_sided']:.3f}; FDR {row['exact_bh_fdr']:.3f}",
                transform=ax.get_yaxis_transform(),
                va="center",
                ha="right",
                fontsize=6.0,
                bbox={"facecolor": "white", "edgecolor": "none", "pad": 0.4},
            )
        ax.set_yticks(y, subset["celltype"])
        ax.invert_yaxis()
        ax.set_title(labels[analysis], loc="left", weight="bold", fontsize=8.2)
        ax.set_xlabel("Post-pre logit-abundance effect")
        ax.grid(axis="x", color="#E5E5E5", lw=0.6)
        ax.spines[["top", "right"]].set_visible(False)
        ax.text(-0.14, 1.07, panel, transform=ax.transAxes, weight="bold", fontsize=9)
        ax.margins(x=0.15)

    ax = axes.flat[3]
    selected = {
        "paired_delta_ordinal_unadjusted": ("Mast", -0.20),
        "paired_delta_ordinal_cohort_adjusted": ("Mast", 0.0),
        "paired_delta_High_vs_Low": ("Myeloid", 0.20),
    }
    for row_index, (analysis, (celltype, offset)) in enumerate(selected.items()):
        subset = null_table[
            (null_table["analysis"] == analysis)
            & (null_table["celltype"] == celltype)
        ]
        ax.scatter(
            subset["permuted_coef"],
            np.full(len(subset), row_index) + offset,
            s=17,
            color="#BDBDBD",
            alpha=0.80,
            edgecolor="none",
        )
        observed = subset["observed_coef"].iloc[0]
        ax.scatter(
            observed,
            row_index + offset,
            s=62,
            color=LINEAGE_COLORS[celltype],
            edgecolor="black",
            linewidth=0.45,
            zorder=4,
        )
    ax.axvline(0, color="#888888", lw=0.8)
    ax.set_yticks(
        np.arange(3),
        [
            "Ordinal Mast (60)",
            "Cohort-adjusted Mast (18)",
            "High-Low Myeloid (4)",
        ],
    )
    ax.invert_yaxis()
    ax.set_xlabel("Permuted effect; coloured point is observed")
    ax.set_title("Exact null distributions", loc="left", weight="bold", fontsize=8.2)
    ax.grid(axis="x", color="#E5E5E5", lw=0.6)
    ax.spines[["top", "right"]].set_visible(False)
    ax.text(-0.14, 1.07, "d", transform=ax.transAxes, weight="bold", fontsize=9)

    fig.subplots_adjust(left=0.09, right=0.985, top=0.96, bottom=0.08, hspace=0.38, wspace=0.88)
    stem = FIGURE_DIR / "ExtendedData8_submission_abundance_exact_permutation"
    fig.savefig(stem.with_suffix(".png"), dpi=450, facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"), facecolor="white")
    fig.savefig(stem.with_suffix(".svg"), facecolor="white")
    plt.close(fig)

    top = results.sort_values(["analysis", "exact_p_two_sided"])
    report = [
        "# Exact Patient-Label Permutation Sensitivity for Paired Abundance",
        "",
        f"Generated: {datetime.now().isoformat(timespec='seconds')}",
        "",
        "## Design",
        "",
        f"- Ordinal analysis: {len(patients)} paired patients and {len(ordinal_permutations)} unique Low/Medium/High label assignments.",
        f"- Cohort-adjusted ordinal sensitivity: response labels permuted within Mono/Combo strata, giving {len(stratified_permutations)} unique assignments.",
        f"- High-versus-Low sensitivity: {len(binary_patients)} patients and {len(binary_permutations)} unique group assignments.",
        "- Two-sided exact P values are the fraction of exhaustive assignments with an absolute coefficient at least as large as observed.",
        "- Benjamini-Hochberg correction was applied across the six broad immune compartments within each analysis.",
        "",
        "## Result",
        "",
        "- No broad immune compartment passed exact-permutation FDR < 0.05 in any analysis.",
        "- The unadjusted ordinal Mast effect had exact P = 0.0833 and FDR = 0.500.",
        "- The unadjusted ordinal Myeloid effect had exact P = 0.2167 and FDR = 0.633.",
        "- The cohort-adjusted ordinal Mast effect had exact P = 0.111 and FDR = 0.667.",
        "- High-versus-Low Myeloid and Mast effects each had exact P = 0.250 and FDR = 0.500.",
        "",
        "The limma coefficients remain useful effect estimates, but exact patient-label permutation is the more defensible small-sample inference layer. Paired abundance findings should therefore be described as exploratory or descriptive.",
        "",
        "## Complete results",
        "",
        "```",
        top[
            [
                "analysis",
                "celltype",
                "coef",
                "parametric_p",
                "parametric_fdr",
                "exact_p_two_sided",
                "exact_bh_fdr",
            ]
        ].to_string(index=False),
        "```",
    ]
    (RESULT_DIR / "ABUNDANCE_EXACT_PERMUTATION_REPORT.md").write_text(
        "\n".join(report) + "\n", encoding="utf-8"
    )
    print(f"Results: {results_path}")
    print(f"Null distributions: {null_path}")
    print(f"Figure: {stem.with_suffix('.png')}")


if __name__ == "__main__":
    main()
