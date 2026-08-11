from __future__ import annotations

import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec


WORKSPACE = Path(os.environ.get("GSE200996_WORKSPACE", Path.cwd())).resolve()
VALIDATION_DIR = (
    WORKSPACE
    / "03_rebuild"
    / "validation"
    / "GSE281729_bulk_module_validation"
)
FIGURE_DIR = WORKSPACE / "03_rebuild" / "figures" / "submission"
SOURCE_DIR = FIGURE_DIR / "source_data"
FIGURE_DIR.mkdir(parents=True, exist_ok=True)
SOURCE_DIR.mkdir(parents=True, exist_ok=True)

ROBUST_PATH = VALIDATION_DIR / "GSE281729_ROBUST_RESPONSE_MODELS.csv"
LOO_PATH = VALIDATION_DIR / "GSE281729_LEAVE_ONE_PATIENT_SUMMARY.csv"
COVERAGE_PATH = VALIDATION_DIR / "GSE281729_LOCKED_MODULE_GENE_COVERAGE.csv"
TIMING_PATH = VALIDATION_DIR / "GSE281729_RESPONSE_ADAPTIVE_TIMING_MODELS.csv"
STRATIFIED_PATH = VALIDATION_DIR / "GSE281729_RESPONSE_ADAPTIVE_TIMING_STRATIFIED.csv"

SELECTED = [
    "M_LE_INTERFERON_ALPHA_RESPONSE",
    "T_LE_INTERFERON_ALPHA_RESPONSE",
    "M_LE_INTERFERON_GAMMA_RESPONSE",
    "T_LE_INTERFERON_GAMMA_RESPONSE",
    "M_LE_MTORC1_SIGNALING",
    "T_LE_MTORC1_SIGNALING",
    "M_LE_union_core",
    "T_LE_union_core",
]
LINEAGE_COLORS = {"T_cell": "#3b6fb6", "Myeloid": "#b4584b"}


def clean_label(signature: str, lineage: str) -> str:
    prefix = "T cell" if lineage == "T_cell" else "Myeloid"
    body = signature.replace("T_", "").replace("M_", "").replace("LE_", "")
    body = body.replace("INTERFERON_ALPHA_RESPONSE", "IFN-alpha")
    body = body.replace("INTERFERON_GAMMA_RESPONSE", "IFN-gamma")
    body = body.replace("MTORC1_SIGNALING", "mTORC1")
    body = body.replace("union_core", "Union core")
    return f"{prefix} | {body.replace('_', ' ')}"


def style_axis(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(labelsize=6.2, width=0.6, length=3)
    ax.grid(axis="x", color="#dedede", linewidth=0.45, zorder=0)


def panel_label(ax: plt.Axes, label: str, x: float = -0.19) -> None:
    ax.text(x, 1.08, label, transform=ax.transAxes, fontsize=8, fontweight="bold", va="top")


def paired_forest(
    ax: plt.Axes,
    data: pd.DataFrame,
    order: list[str],
    labels: dict[str, str],
    models: list[tuple[str, float, str, str, str]],
    xlabel: str,
    title: str,
    subtitle: str,
    show_labels: bool,
) -> None:
    y = np.arange(len(order))
    ax.axvline(0, color="#777777", linewidth=0.7, zorder=1)
    indexed = data.set_index(["signature", "model"])
    for model, offset, marker, color, legend_label in models:
        sub = indexed.loc[[(sig, model) for sig in order]].reset_index()
        x = sub["coef"].to_numpy(float)
        low = sub["ci95_low"].to_numpy(float)
        high = sub["ci95_high"].to_numpy(float)
        ax.errorbar(
            x,
            y + offset,
            xerr=np.vstack([x - low, high - x]),
            fmt="none",
            ecolor=color,
            elinewidth=0.7,
            capsize=1.8,
            zorder=2,
        )
        ax.scatter(
            x,
            y + offset,
            s=22,
            marker=marker,
            facecolor="white" if marker == "s" else color,
            edgecolor=color,
            linewidth=0.8,
            label=legend_label,
            zorder=3,
        )
    ax.set_yticks(y)
    ax.set_yticklabels([labels[sig] for sig in order] if show_labels else [])
    ax.set_xlabel(xlabel)
    ax.set_title(title, loc="left", fontsize=8, fontweight="bold", pad=12)
    ax.text(0.00, 1.005, subtitle, transform=ax.transAxes, fontsize=5.6, color="#444444", va="bottom")
    style_axis(ax)


def main() -> None:
    robust = pd.read_csv(ROBUST_PATH)
    loo = pd.read_csv(LOO_PATH)
    coverage = pd.read_csv(COVERAGE_PATH)
    timing = pd.read_csv(TIMING_PATH)
    stratified = pd.read_csv(STRATIFIED_PATH)

    robust = robust[robust["signature"].isin(SELECTED)].copy()
    loo = loo[loo["signature"].isin(SELECTED)].copy()
    coverage = coverage[coverage["signature"].isin(SELECTED)].copy()
    timing = timing[timing["signature"].isin(SELECTED)].copy()
    stratified = stratified[stratified["signature"].isin(SELECTED)].copy()

    labels = {
        row["signature"]: clean_label(row["signature"], row["target_lineage"])
        for _, row in robust.drop_duplicates("signature").iterrows()
    }
    order = (
        robust[robust["model"].eq("unadjusted_classical")]
        .sort_values("coef")["signature"]
        .tolist()
    )
    y = np.arange(len(order))

    robust.to_csv(SOURCE_DIR / "ExtendedData7_robust_models_source.csv", index=False)
    loo.to_csv(SOURCE_DIR / "ExtendedData7_leave_one_out_source.csv", index=False)
    coverage.to_csv(SOURCE_DIR / "ExtendedData7_module_coverage_source.csv", index=False)
    timing.to_csv(SOURCE_DIR / "ExtendedData7_response_adaptive_timing_source.csv", index=False)
    stratified.to_csv(SOURCE_DIR / "ExtendedData7_timing_stratified_source.csv", index=False)

    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 7,
            "axes.linewidth": 0.6,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )
    fig = plt.figure(figsize=(7.2, 6.5))
    gs = GridSpec(3, 2, figure=fig, hspace=0.44, wspace=0.43)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])
    ax_e = fig.add_subplot(gs[2, 0])
    ax_f = fig.add_subplot(gs[2, 1])

    paired_forest(
        ax_a,
        robust,
        order,
        labels,
        [
            ("unadjusted_classical", -0.12, "o", "#333333", "Unadjusted"),
            ("adjusted_hpv_second_drug_classical", 0.12, "s", "#333333", "HPV + drug adjusted"),
        ],
        "Ordinal response slope (95% CI)",
        "Covariate adjustment",
        "Filled, unadjusted; open, HPV + second-drug adjusted",
        True,
    )
    panel_label(ax_a, "a")

    paired_forest(
        ax_b,
        robust,
        order,
        labels,
        [
            ("adjusted_hpv_second_drug_classical", -0.12, "o", "#999999", "Classical"),
            ("adjusted_hpv_second_drug_HC3", 0.12, "o", "#222222", "HC3"),
        ],
        "Adjusted response slope (95% CI)",
        "Robust covariance",
        "Grey, classical; black, HC3",
        False,
    )
    panel_label(ax_b, "b", -0.12)

    ax_c.axvline(0, color="#777777", linewidth=0.7, zorder=1)
    loo_plot = (
        loo[loo["model"].eq("adjusted_hpv_second_drug")]
        .set_index("signature")
        .loc[order]
        .reset_index()
    )
    for yi, row in loo_plot.iterrows():
        color = LINEAGE_COLORS[row["target_lineage"]]
        ax_c.plot([row["loo_coef_min"], row["loo_coef_max"]], [yi, yi], color="#777777", linewidth=0.7)
        ax_c.plot([row["loo_coef_q1"], row["loo_coef_q3"]], [yi, yi], color=color, linewidth=4.0)
        ax_c.scatter(row["full_coef"], yi, s=23, color=color, edgecolor="white", linewidth=0.4, zorder=3)
    ax_c.set_yticks(y)
    ax_c.set_yticklabels([labels[sig] for sig in order])
    ax_c.set_xlabel("Adjusted slope across leave-one-patient refits")
    ax_c.set_title("Patient influence", loc="left", fontsize=8, fontweight="bold", pad=5)
    style_axis(ax_c)
    panel_label(ax_c, "c")

    cov = coverage.set_index("signature").loc[order].reset_index()
    cov["coverage_fraction"] = cov["n_genes_present_in_GSE281729"] / cov["n_genes_defined"]
    colors = [LINEAGE_COLORS[lineage] for lineage in cov["target_lineage"]]
    ax_d.barh(y, cov["coverage_fraction"], color=colors, height=0.56, edgecolor="none")
    for yi, row in cov.iterrows():
        ax_d.text(
            min(row["coverage_fraction"] + 0.02, 1.01),
            yi,
            f"{int(row['n_genes_present_in_GSE281729'])}/{int(row['n_genes_defined'])}",
            va="center",
            fontsize=5.8,
        )
    ax_d.set_yticks(y)
    ax_d.set_yticklabels([])
    ax_d.set_xlim(0, 1.16)
    ax_d.set_xticks([0, 0.5, 1.0])
    ax_d.set_xticklabels(["0", "50%", "100%"])
    ax_d.set_xlabel("Module gene coverage")
    ax_d.set_title("Locked-module coverage", loc="left", fontsize=8, fontweight="bold", pad=5)
    style_axis(ax_d)
    panel_label(ax_d, "d", -0.12)

    timing_compare = pd.concat(
        [
            robust[robust["model"].eq("adjusted_hpv_second_drug_HC3")],
            timing[timing["model"].eq("ordinal_hpv_second_drug_doses_HC3")],
        ],
        ignore_index=True,
    )
    paired_forest(
        ax_e,
        timing_compare,
        order,
        labels,
        [
            ("adjusted_hpv_second_drug_HC3", -0.12, "o", "#333333", "HPV + drug"),
            ("ordinal_hpv_second_drug_doses_HC3", 0.12, "s", "#333333", "+ doses/timing"),
        ],
        "HC3-adjusted response slope (95% CI)",
        "Response-adaptive timing sensitivity",
        "Filled, HPV + drug; open, + doses/timing; minimum module-level FDR = 0.0537",
        True,
    )
    panel_label(ax_e, "e")

    paired_forest(
        ax_f,
        stratified.assign(model="dose_" + stratified["doses"].astype(str)),
        order,
        labels,
        [
            ("dose_1", -0.12, "s", "#777777", "One dose"),
            ("dose_2", 0.12, "o", "#222222", "Two doses"),
        ],
        "Within-stratum ordinal slope (95% CI)",
        "Timing-stratified associations",
        "Open, one dose (n = 16); filled, two doses (n = 14)",
        False,
    )
    panel_label(ax_f, "f", -0.12)

    fig.subplots_adjust(left=0.25, right=0.98, top=0.94, bottom=0.065)
    stem = FIGURE_DIR / "ExtendedData7_submission_gse281729_robustness"
    fig.savefig(stem.with_suffix(".png"), dpi=450, facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"), facecolor="white")
    fig.savefig(stem.with_suffix(".svg"), facecolor="white")
    plt.close(fig)

    print(stem.with_suffix(".png"))


if __name__ == "__main__":
    main()
