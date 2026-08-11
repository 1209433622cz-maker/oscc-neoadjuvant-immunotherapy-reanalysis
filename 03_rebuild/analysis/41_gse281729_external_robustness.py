#!/usr/bin/env python
"""Patient-level robustness analyses for GSE281729 locked-module validation."""

from __future__ import annotations

import math
import os
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
from matplotlib.gridspec import GridSpec
from statsmodels.stats.multitest import multipletests


WORKSPACE = Path(os.environ.get("GSE200996_WORKSPACE", Path.cwd())).resolve()
VALIDATION_DIR = WORKSPACE / "03_rebuild" / "validation" / "GSE281729_bulk_module_validation"
SCORES_PATH = VALIDATION_DIR / "GSE281729_LOCKED_MODULE_SAMPLE_SCORES.csv"
PRIMARY_STATS_PATH = VALIDATION_DIR / "GSE281729_PAIRED_DELTA_MODULE_RESPONSE_MODELS.csv"
COVERAGE_PATH = VALIDATION_DIR / "GSE281729_LOCKED_MODULE_GENE_COVERAGE.csv"
FIGURE_DIR = WORKSPACE / "03_rebuild" / "figures" / "submission"
SOURCE_DIR = FIGURE_DIR / "source_data"
FIGURE_DIR.mkdir(parents=True, exist_ok=True)
SOURCE_DIR.mkdir(parents=True, exist_ok=True)

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

LINEAGE_COLORS = {"T_cell": "#3B6FB6", "Myeloid": "#B4584B"}


def bh(values: pd.Series) -> np.ndarray:
    out = np.full(len(values), np.nan)
    valid = values.notna().to_numpy()
    if valid.any():
        out[valid] = multipletests(values.loc[valid].astype(float), method="fdr_bh")[1]
    return out


def clean_label(signature: str, lineage: str) -> str:
    prefix = "T cell" if lineage == "T_cell" else "Myeloid"
    body = signature.replace("T_", "").replace("M_", "").replace("LE_", "")
    body = body.replace("INTERFERON_ALPHA_RESPONSE", "IFN-alpha")
    body = body.replace("INTERFERON_GAMMA_RESPONSE", "IFN-gamma")
    body = body.replace("MTORC1_SIGNALING", "mTORC1")
    body = body.replace("union_core", "Union core").replace("_", " ")
    return f"{prefix} | {body}"


def build_paired_table() -> pd.DataFrame:
    scores = pd.read_csv(SCORES_PATH)
    scores["timepoint"] = scores["timepoint"].astype(str).str.lower()
    wide = scores.pivot_table(
        index=["patient_id", "signature"],
        columns="timepoint",
        values="module_score",
        aggfunc="mean",
    ).reset_index()
    meta_cols = [
        "patient_id",
        "signature",
        "target_lineage",
        "response_harmonized_ordinal",
        "response_ord_num",
        "hpv",
        "second_drug",
        "n_genes_present",
    ]
    meta = scores[meta_cols].drop_duplicates(["patient_id", "signature"])
    paired = wide.merge(meta, on=["patient_id", "signature"], how="left")
    paired["post_minus_pre"] = paired["post"] - paired["pre"]
    paired = paired.dropna(subset=["post_minus_pre", "response_ord_num"]).copy()
    return paired


def design_matrix(data: pd.DataFrame, adjusted: bool) -> pd.DataFrame:
    parts = [data[["response_ord_num"]].astype(float)]
    if adjusted:
        for covariate in ["hpv", "second_drug"]:
            if data[covariate].nunique() > 1:
                parts.append(pd.get_dummies(data[covariate].astype(str), prefix=covariate, drop_first=True, dtype=float))
    return sm.add_constant(pd.concat(parts, axis=1), has_constant="add")


def fit_model(data: pd.DataFrame, adjusted: bool, covariance: str) -> tuple[dict[str, float | int | str], object]:
    keep = ["post_minus_pre", "response_ord_num", "hpv", "second_drug", "patient_id"]
    model_data = data[keep].replace([np.inf, -np.inf], np.nan).dropna().copy()
    X = design_matrix(model_data, adjusted)
    fit = sm.OLS(model_data["post_minus_pre"].astype(float), X.astype(float)).fit()
    if covariance == "HC3":
        fit = fit.get_robustcov_results(cov_type="HC3")
        names = list(X.columns)
        coef_index = names.index("response_ord_num")
        coef = float(fit.params[coef_index])
        se = float(fit.bse[coef_index])
        p_value = float(fit.pvalues[coef_index])
        ci = fit.conf_int(alpha=0.05)[coef_index]
    else:
        coef = float(fit.params["response_ord_num"])
        se = float(fit.bse["response_ord_num"])
        p_value = float(fit.pvalues["response_ord_num"])
        ci = fit.conf_int(alpha=0.05).loc["response_ord_num"].to_numpy(dtype=float)
    return (
        {
            "n": int(len(model_data)),
            "coef": coef,
            "std_error": se,
            "ci95_low": float(ci[0]),
            "ci95_high": float(ci[1]),
            "p_value": p_value,
        },
        fit,
    )


def run_robust_models(paired: pd.DataFrame) -> pd.DataFrame:
    rows = []
    model_specs = [
        ("unadjusted_classical", False, "classical"),
        ("unadjusted_HC3", False, "HC3"),
        ("adjusted_hpv_second_drug_classical", True, "classical"),
        ("adjusted_hpv_second_drug_HC3", True, "HC3"),
    ]
    for signature, sub in paired.groupby("signature"):
        lineage = str(sub["target_lineage"].iloc[0])
        for model, adjusted, covariance in model_specs:
            result, _ = fit_model(sub, adjusted=adjusted, covariance=covariance)
            rows.append(
                {
                    "signature": signature,
                    "target_lineage": lineage,
                    "model": model,
                    "adjusted": adjusted,
                    "covariance": covariance,
                    **result,
                }
            )
    out = pd.DataFrame(rows)
    out["fdr"] = np.nan
    for model in out["model"].unique():
        mask = out["model"] == model
        out.loc[mask, "fdr"] = bh(out.loc[mask, "p_value"])
    return out


def verify_against_primary(robust: pd.DataFrame) -> None:
    primary = pd.read_csv(PRIMARY_STATS_PATH)
    checks = [
        ("ordinal_unadjusted", "unadjusted_classical"),
        ("ordinal_adjusted_hpv_second_drug", "adjusted_hpv_second_drug_classical"),
    ]
    for primary_model, robust_model in checks:
        left = primary[
            (primary["analysis"] == "paired_post_minus_pre")
            & (primary["model"] == primary_model)
            & (primary["model_status"] == "ok")
        ][["signature", "coef"]]
        right = robust[robust["model"] == robust_model][["signature", "coef"]]
        merged = left.merge(right, on="signature", suffixes=("_primary", "_robust"))
        if len(merged) != len(left) or not np.allclose(merged["coef_primary"], merged["coef_robust"], atol=1e-10):
            raise RuntimeError(f"Robustness rerun does not reproduce {primary_model} coefficients")


def run_leave_one_out(paired: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    for signature in SELECTED:
        sub = paired[paired["signature"] == signature].copy()
        lineage = str(sub["target_lineage"].iloc[0])
        for adjusted, model in [(False, "unadjusted"), (True, "adjusted_hpv_second_drug")]:
            full, _ = fit_model(sub, adjusted=adjusted, covariance="classical")
            for patient in sorted(sub["patient_id"].unique()):
                result, _ = fit_model(sub[sub["patient_id"] != patient], adjusted=adjusted, covariance="classical")
                rows.append(
                    {
                        "signature": signature,
                        "target_lineage": lineage,
                        "model": model,
                        "left_out_patient": patient,
                        "full_coef": full["coef"],
                        **result,
                    }
                )
    loo = pd.DataFrame(rows)
    summary_rows = []
    for (signature, model), sub in loo.groupby(["signature", "model"]):
        summary_rows.append(
            {
                "signature": signature,
                "target_lineage": sub["target_lineage"].iloc[0],
                "model": model,
                "n_refits": len(sub),
                "full_coef": sub["full_coef"].iloc[0],
                "loo_coef_min": sub["coef"].min(),
                "loo_coef_q1": sub["coef"].quantile(0.25),
                "loo_coef_median": sub["coef"].median(),
                "loo_coef_q3": sub["coef"].quantile(0.75),
                "loo_coef_max": sub["coef"].max(),
                "negative_refit_fraction": float((sub["coef"] < 0).mean()),
            }
        )
    return loo, pd.DataFrame(summary_rows)


def run_influence(paired: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    for signature in SELECTED:
        sub = paired[paired["signature"] == signature].copy().reset_index(drop=True)
        lineage = str(sub["target_lineage"].iloc[0])
        for adjusted, model in [(False, "unadjusted"), (True, "adjusted_hpv_second_drug")]:
            X = design_matrix(sub, adjusted)
            fit = sm.OLS(sub["post_minus_pre"].astype(float), X.astype(float)).fit()
            cooks = fit.get_influence().cooks_distance[0]
            threshold = 4.0 / len(sub)
            for idx, cook in enumerate(cooks):
                rows.append(
                    {
                        "signature": signature,
                        "target_lineage": lineage,
                        "model": model,
                        "patient_id": sub.loc[idx, "patient_id"],
                        "cooks_distance": float(cook),
                        "threshold_4_over_n": threshold,
                        "above_threshold": bool(cook > threshold),
                    }
                )
    influence = pd.DataFrame(rows)
    summary_rows = []
    for (signature, model), sub in influence.groupby(["signature", "model"]):
        top = sub.loc[sub["cooks_distance"].idxmax()]
        summary_rows.append(
            {
                "signature": signature,
                "target_lineage": sub["target_lineage"].iloc[0],
                "model": model,
                "top_influence_patient": top["patient_id"],
                "max_cooks_distance": top["cooks_distance"],
                "threshold_4_over_n": top["threshold_4_over_n"],
                "n_above_threshold": int(sub["above_threshold"].sum()),
            }
        )
    return influence, pd.DataFrame(summary_rows)


def selected_coverage() -> pd.DataFrame:
    coverage = pd.read_csv(COVERAGE_PATH)
    coverage = coverage[coverage["signature"].isin(SELECTED)].copy()
    coverage["coverage_fraction"] = coverage["n_genes_present_in_GSE281729"] / coverage["n_genes_defined"]
    coverage["label"] = [clean_label(s, l) for s, l in zip(coverage["signature"], coverage["target_lineage"])]
    return coverage


def style_axis(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(labelsize=6.5, width=0.6, length=3)
    ax.grid(axis="x", color="#dedede", linewidth=0.45, zorder=0)


def make_figure(robust: pd.DataFrame, loo_summary: pd.DataFrame, coverage: pd.DataFrame) -> None:
    labels = {
        row["signature"]: clean_label(row["signature"], row["target_lineage"])
        for _, row in robust[robust["signature"].isin(SELECTED)].iterrows()
    }
    order = (
        robust[(robust["signature"].isin(SELECTED)) & (robust["model"] == "unadjusted_classical")]
        .sort_values("coef")["signature"]
        .tolist()
    )
    y = np.arange(len(order))

    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 7,
            "axes.linewidth": 0.6,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    fig = plt.figure(figsize=(7.2, 6.2))
    gs = GridSpec(2, 2, figure=fig, height_ratios=[1.05, 1.0], width_ratios=[1.12, 0.88], hspace=0.44, wspace=0.52)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])

    ax_a.axvline(0, color="#777777", linewidth=0.7, zorder=1)
    for offset, model, marker, name in [
        (-0.12, "unadjusted_classical", "o", "Unadjusted"),
        (0.12, "adjusted_hpv_second_drug_classical", "s", "HPV + second drug adjusted"),
    ]:
        sub = robust.set_index(["signature", "model"]).loc[[(sig, model) for sig in order]].reset_index()
        x = sub["coef"].to_numpy(float)
        low = sub["ci95_low"].to_numpy(float)
        high = sub["ci95_high"].to_numpy(float)
        ax_a.errorbar(x, y + offset, xerr=np.vstack([x - low, high - x]), fmt="none", ecolor="#4b4b4b", elinewidth=0.7, capsize=2, zorder=2)
        ax_a.scatter(x, y + offset, s=26, marker=marker, facecolor="white" if marker == "s" else "#333333", edgecolor="#333333", linewidth=0.8, label=name, zorder=3)
    ax_a.set_yticks(y)
    ax_a.set_yticklabels([labels[sig] for sig in order])
    ax_a.set_xlabel("Ordinal response slope (95% CI)")
    ax_a.set_title("Model-adjusted response associations", loc="left", fontsize=8, fontweight="bold", pad=13)
    ax_a.text(0.00, 1.01, "Filled circle, unadjusted; open square, HPV + second-drug adjusted", transform=ax_a.transAxes, fontsize=5.8, color="#444444", va="bottom")
    ax_a.text(-0.20, 1.08, "a", transform=ax_a.transAxes, fontsize=10, fontweight="bold", va="top")
    style_axis(ax_a)

    ax_b.axvline(0, color="#777777", linewidth=0.7, zorder=1)
    for offset, model, color, name in [
        (-0.12, "adjusted_hpv_second_drug_classical", "#9a9a9a", "Classical"),
        (0.12, "adjusted_hpv_second_drug_HC3", "#1f1f1f", "HC3 robust"),
    ]:
        sub = robust.set_index(["signature", "model"]).loc[[(sig, model) for sig in order]].reset_index()
        x = sub["coef"].to_numpy(float)
        low = sub["ci95_low"].to_numpy(float)
        high = sub["ci95_high"].to_numpy(float)
        ax_b.errorbar(x, y + offset, xerr=np.vstack([x - low, high - x]), fmt="none", ecolor=color, elinewidth=0.8, capsize=2, zorder=2)
        ax_b.scatter(x, y + offset, s=22, color=color, edgecolor="white", linewidth=0.4, label=name, zorder=3)
    ax_b.set_yticks(y)
    ax_b.set_yticklabels([])
    ax_b.set_xlabel("Adjusted response slope (95% CI)")
    ax_b.set_title("Heteroskedasticity-robust inference", loc="left", fontsize=8, fontweight="bold", pad=13)
    ax_b.text(0.00, 1.01, "Grey, classical CI; black, HC3 robust CI", transform=ax_b.transAxes, fontsize=5.8, color="#444444", va="bottom")
    ax_b.text(-0.13, 1.08, "b", transform=ax_b.transAxes, fontsize=10, fontweight="bold", va="top")
    style_axis(ax_b)

    ax_c.axvline(0, color="#777777", linewidth=0.7, zorder=1)
    loo_plot = loo_summary[loo_summary["model"] == "adjusted_hpv_second_drug"].set_index("signature").loc[order].reset_index()
    for yi, row in loo_plot.iterrows():
        color = LINEAGE_COLORS[row["target_lineage"]]
        ax_c.plot([row["loo_coef_min"], row["loo_coef_max"]], [yi, yi], color="#777777", linewidth=0.7, zorder=1)
        ax_c.plot([row["loo_coef_q1"], row["loo_coef_q3"]], [yi, yi], color=color, linewidth=4.0, solid_capstyle="butt", zorder=2)
        ax_c.scatter(row["full_coef"], yi, s=25, color=color, edgecolor="white", linewidth=0.5, zorder=3)
    ax_c.set_yticks(y)
    ax_c.set_yticklabels([labels[sig] for sig in order])
    ax_c.set_xlabel("Adjusted slope across leave-one-patient refits")
    ax_c.set_title("Patient influence on adjusted slopes", loc="left", fontsize=8, fontweight="bold", pad=5)
    ax_c.text(-0.20, 1.05, "c", transform=ax_c.transAxes, fontsize=10, fontweight="bold", va="top")
    style_axis(ax_c)

    cov = coverage.set_index("signature").loc[order].reset_index()
    colors = [LINEAGE_COLORS[lineage] for lineage in cov["target_lineage"]]
    ax_d.barh(y, cov["coverage_fraction"], color=colors, height=0.56, edgecolor="none")
    for yi, row in cov.iterrows():
        ax_d.text(min(row["coverage_fraction"] + 0.02, 1.01), yi, f"{int(row['n_genes_present_in_GSE281729'])}/{int(row['n_genes_defined'])}", va="center", fontsize=6.2)
    ax_d.set_yticks(y)
    ax_d.set_yticklabels([])
    ax_d.set_xlim(0, 1.16)
    ax_d.set_xlabel("Module gene coverage")
    ax_d.set_title("Locked-module coverage", loc="left", fontsize=8, fontweight="bold", pad=5)
    ax_d.text(-0.13, 1.05, "d", transform=ax_d.transAxes, fontsize=10, fontweight="bold", va="top")
    ax_d.set_xticks([0, 0.5, 1.0])
    ax_d.set_xticklabels(["0", "50%", "100%"])
    style_axis(ax_d)

    fig.subplots_adjust(left=0.25, right=0.97, top=0.96, bottom=0.08)
    stem = FIGURE_DIR / "ExtendedData7_submission_gse281729_robustness"
    fig.savefig(stem.with_suffix(".png"), dpi=450, bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
    plt.close(fig)


def write_report(
    paired: pd.DataFrame,
    robust: pd.DataFrame,
    loo_summary: pd.DataFrame,
    influence_summary: pd.DataFrame,
    coverage: pd.DataFrame,
) -> None:
    selected_robust = robust[robust["signature"].isin(SELECTED)].copy()
    hc3 = selected_robust[selected_robust["model"] == "adjusted_hpv_second_drug_HC3"]
    loo = loo_summary[loo_summary["model"] == "adjusted_hpv_second_drug"]
    influence = influence_summary[influence_summary["model"] == "adjusted_hpv_second_drug"]
    lines = [
        "# GSE281729 External Robustness Report",
        "",
        f"Generated: {datetime.now().astimezone().isoformat(timespec='seconds')}",
        "",
        "## Analysis Boundary",
        "",
        f"- Response-annotated paired patients: {paired['patient_id'].nunique()}",
        f"- Locked modules tested: {paired['signature'].nunique()}",
        f"- Leading modules displayed: {len(SELECTED)}",
        "- Primary covariate adjustment: HPV status and second-drug exposure.",
        "- Robust covariance: HC3.",
        "- Influence diagnostics: leave-one-patient refits and Cook's distance.",
        "",
        "## Main Findings",
        "",
        f"- Adjusted HC3 slopes negative among displayed modules: {int((hc3['coef'] < 0).sum())}/{len(hc3)}.",
        f"- Adjusted HC3 FDR < 0.05 among displayed modules: {int((hc3['fdr'] < 0.05).sum())}/{len(hc3)}.",
        f"- Modules with all adjusted leave-one-patient slopes negative: {int((loo['negative_refit_fraction'] == 1).sum())}/{len(loo)}.",
        f"- Modules with at least one adjusted Cook's-distance value above 4/n: {int((influence['n_above_threshold'] > 0).sum())}/{len(influence)}.",
        f"- Median displayed module gene coverage: {coverage['coverage_fraction'].median():.1%}.",
        "",
        "## Interpretation",
        "",
        "The disease-matched bulk response associations are not removed by HPV/second-drug adjustment or HC3 covariance.",
        "Leave-one-patient refits quantify sign stability but do not create an independent validation cohort.",
        "Cook's-distance flags identify observations requiring transparency; they do not by themselves invalidate a model.",
        "Because the assay is bulk RNA-seq, these analyses support pathway-level response association and cannot assign the signal to a specific cell lineage in this cohort.",
        "",
        "## Outputs",
        "",
        "- `GSE281729_PAIRED_MODULE_DELTAS_ROBUSTNESS.csv`",
        "- `GSE281729_ROBUST_RESPONSE_MODELS.csv`",
        "- `GSE281729_LEAVE_ONE_PATIENT_MODELS.csv`",
        "- `GSE281729_LEAVE_ONE_PATIENT_SUMMARY.csv`",
        "- `GSE281729_COOKS_DISTANCE.csv`",
        "- `GSE281729_COOKS_DISTANCE_SUMMARY.csv`",
        "- `ExtendedData7_submission_gse281729_robustness.png/pdf/svg`",
    ]
    (VALIDATION_DIR / "GSE281729_EXTERNAL_ROBUSTNESS_REPORT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    paired = build_paired_table()
    robust = run_robust_models(paired)
    verify_against_primary(robust)
    loo, loo_summary = run_leave_one_out(paired)
    influence, influence_summary = run_influence(paired)
    coverage = selected_coverage()

    paired.to_csv(VALIDATION_DIR / "GSE281729_PAIRED_MODULE_DELTAS_ROBUSTNESS.csv", index=False)
    robust.to_csv(VALIDATION_DIR / "GSE281729_ROBUST_RESPONSE_MODELS.csv", index=False)
    loo.to_csv(VALIDATION_DIR / "GSE281729_LEAVE_ONE_PATIENT_MODELS.csv", index=False)
    loo_summary.to_csv(VALIDATION_DIR / "GSE281729_LEAVE_ONE_PATIENT_SUMMARY.csv", index=False)
    influence.to_csv(VALIDATION_DIR / "GSE281729_COOKS_DISTANCE.csv", index=False)
    influence_summary.to_csv(VALIDATION_DIR / "GSE281729_COOKS_DISTANCE_SUMMARY.csv", index=False)

    robust[robust["signature"].isin(SELECTED)].to_csv(SOURCE_DIR / "ExtendedData7_model_comparison_source.csv", index=False)
    loo_summary[(loo_summary["signature"].isin(SELECTED)) & (loo_summary["model"] == "adjusted_hpv_second_drug")].to_csv(
        SOURCE_DIR / "ExtendedData7_leave_one_out_source.csv", index=False
    )
    influence_summary[(influence_summary["signature"].isin(SELECTED)) & (influence_summary["model"] == "adjusted_hpv_second_drug")].to_csv(
        SOURCE_DIR / "ExtendedData7_influence_source.csv", index=False
    )
    coverage.to_csv(SOURCE_DIR / "ExtendedData7_gene_coverage_source.csv", index=False)

    make_figure(robust, loo_summary, coverage)
    write_report(paired, robust, loo_summary, influence_summary, coverage)
    print(VALIDATION_DIR / "GSE281729_EXTERNAL_ROBUSTNESS_REPORT.md")
    print(FIGURE_DIR / "ExtendedData7_submission_gse281729_robustness.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
