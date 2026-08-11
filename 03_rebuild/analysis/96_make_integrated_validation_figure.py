from __future__ import annotations

import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle


WORKSPACE = Path(os.environ.get("GSE200996_WORKSPACE", Path.cwd())).resolve()
GSE281729_STATS = (
    WORKSPACE
    / "03_rebuild"
    / "validation"
    / "GSE281729_bulk_module_validation"
    / "GSE281729_ROBUST_RESPONSE_MODELS.csv"
)
GSE179730_STATS = (
    WORKSPACE
    / "03_rebuild"
    / "validation"
    / "GSE179730_bulk_treatment_direction"
    / "GSE179730_LOCKED_MODULE_RESPONSE_EXACT.csv"
)
OUT_DIR = WORKSPACE / "03_rebuild" / "figures" / "external_validation"
SRC_DIR = OUT_DIR / "source_data"
OUT_DIR.mkdir(parents=True, exist_ok=True)
SRC_DIR.mkdir(parents=True, exist_ok=True)


def clean_label(signature: str, lineage: str) -> str:
    prefix = "T cell" if lineage == "T_cell" else "Myeloid"
    body = signature
    body = body.replace("T_", "").replace("M_", "")
    body = body.replace("LE_", "")
    body = body.replace("INTERFERON_ALPHA_RESPONSE", "IFN-alpha")
    body = body.replace("INTERFERON_GAMMA_RESPONSE", "IFN-gamma")
    body = body.replace("MTORC1_SIGNALING", "mTORC1")
    body = body.replace("TNFA_SIGNALING_VIA_NFKB", "TNFA/NF-kB")
    body = body.replace("union_core", "Union core")
    body = body.replace("_", " ")
    return f"{prefix} | {body}"


def load_gse281729() -> pd.DataFrame:
    df = pd.read_csv(GSE281729_STATS)
    df = df[
        df["model"].eq("adjusted_hpv_second_drug_HC3")
    ].copy()
    keep = [
        "M_LE_INTERFERON_ALPHA_RESPONSE",
        "T_LE_INTERFERON_ALPHA_RESPONSE",
        "M_LE_INTERFERON_GAMMA_RESPONSE",
        "T_LE_INTERFERON_GAMMA_RESPONSE",
        "M_LE_MTORC1_SIGNALING",
        "T_LE_MTORC1_SIGNALING",
        "M_LE_union_core",
        "T_LE_union_core",
    ]
    df = df[df["signature"].isin(keep)].copy()
    df["label"] = [clean_label(sig, lin) for sig, lin in zip(df["signature"], df["target_lineage"])]
    df = df.sort_values("coef")
    return df


def load_gse179730() -> pd.DataFrame:
    df = pd.read_csv(GSE179730_STATS)
    keep = [
        "T_LE_INTERFERON_ALPHA_RESPONSE",
        "M_LE_INTERFERON_ALPHA_RESPONSE",
        "T_LE_INTERFERON_GAMMA_RESPONSE",
        "M_LE_INTERFERON_GAMMA_RESPONSE",
        "T_LE_MTORC1_SIGNALING",
        "M_LE_MTORC1_SIGNALING",
        "T_LE_union_core",
        "M_LE_union_core",
    ]
    df = df[df["signature"].isin(keep)].copy()
    df["label"] = [clean_label(sig, lin) for sig, lin in zip(df["signature"], df["target_lineage"])]
    df = df.sort_values("responder_minus_non_responder", ascending=False)
    return df


def draw_cohort_matrix(ax: plt.Axes) -> None:
    rows = [
        ("GSE200996", "scRNA", "Discovery", "partial"),
        ("GSE281729", "Bulk", "Response/timing", "ready"),
        ("GSE179730", "Bulk", "Response check", "pending"),
        ("GSE301741", "scRNA/TCR", "Provisional/null", "pending"),
        ("GSE123813", "scRNA/TCR", "Boundary", "cross"),
    ]
    colors = {
        "ready": "#1f1f1f",
        "partial": "#555555",
        "pending": "#8f8f8f",
        "cross": "#c4c4c4",
    }
    x_positions = np.linspace(0.10, 0.90, len(rows))
    ax.plot([0.10, 0.90], [0.58, 0.58], color="#d0d0d0", lw=0.8, transform=ax.transAxes, zorder=1)
    for x, (cohort, modality, role, status) in zip(x_positions, rows):
        ax.scatter(x, 0.58, s=70, color=colors[status], edgecolor="white", linewidth=0.7, transform=ax.transAxes, zorder=2)
        ax.text(x, 0.35, cohort, ha="center", va="center", fontsize=7.5, fontweight="bold", transform=ax.transAxes)
        ax.text(x, 0.14, f"{modality} | {role}", ha="center", va="center", fontsize=6.1, color="#444444", transform=ax.transAxes)
    ax.text(-0.035, 1.02, "a", fontsize=10, fontweight="bold", va="top", transform=ax.transAxes)
    ax.set_title("External evidence hierarchy", loc="left", fontsize=8, fontweight="bold", pad=2)
    ax.axis("off")


def draw_gse281729(ax: plt.Axes, df: pd.DataFrame) -> None:
    y = np.arange(len(df))
    colors = np.where(df["target_lineage"].eq("T_cell"), "#3b6fb6", "#b4584b")
    ax.axvline(0, color="#888888", lw=0.7, zorder=1)
    x = df["coef"].to_numpy(dtype=float)
    xerr = np.vstack([x - df["ci95_low"].to_numpy(dtype=float), df["ci95_high"].to_numpy(dtype=float) - x])
    ax.errorbar(x, y, xerr=xerr, fmt="none", ecolor="#333333", elinewidth=0.7, capsize=2, zorder=2)
    ax.scatter(x, y, c=colors, s=np.where(df["fdr"] < 0.05, 42, 28), edgecolor="white", linewidth=0.5, zorder=3)
    for yi, (_, row) in zip(y, df.iterrows()):
        ax.text(0.045, yi, f"FDR {row['fdr']:.3g}", ha="left", va="center", fontsize=6, color="#444444")
    ax.set_yticks(y)
    ax.set_yticklabels(df["label"], fontsize=6.5)
    ax.set_xlabel("Adjusted ordinal response slope for paired post-pre delta (95% CI)")
    ax.text(-0.12, 1.04, "b", fontsize=10, fontweight="bold", va="top", transform=ax.transAxes)
    ax.set_title("GSE281729: HNSCC response association", loc="left", fontsize=8, fontweight="bold", pad=6)
    ax.text(
        0.00,
        1.005,
        "Module-level estimates; global-shift stress test in Extended Data Fig. 13",
        transform=ax.transAxes,
        fontsize=5.9,
        color="#444444",
        va="bottom",
    )
    ax.set_xlim(min(-0.70, float(df["ci95_low"].min()) - 0.04), 0.18)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="x", color="#dddddd", lw=0.4)


def draw_gse179730(ax: plt.Axes, df: pd.DataFrame) -> None:
    y = np.arange(len(df))
    colors = np.where(df["target_lineage"].eq("T_cell"), "#3b6fb6", "#b4584b")
    ax.axvline(0, color="#888888", lw=0.7, zorder=1)
    x = df["responder_minus_non_responder"].to_numpy()
    xerr = np.vstack([x - df["ci95_low"].to_numpy(), df["ci95_high"].to_numpy() - x])
    ax.errorbar(x, y, xerr=xerr, fmt="none", ecolor="#333333", elinewidth=0.7, capsize=2, zorder=2)
    ax.scatter(x, y, c=colors, s=30, edgecolor="white", linewidth=0.5, zorder=3)
    for yi, (_, row) in zip(y, df.iterrows()):
        ax.text(
            max(float(row["ci95_high"]) + 0.025, 0.34),
            yi,
            f"exact P {row['exact_p']:.2g}",
            ha="left",
            va="center",
            fontsize=5.7,
            color="#444444",
        )
    ax.set_yticks(y)
    ax.set_yticklabels(df["label"], fontsize=6.5)
    ax.set_xlabel("Responder - non-responder paired treatment delta")
    ax.text(-0.20, 1.04, "c", fontsize=10, fontweight="bold", va="top", transform=ax.transAxes)
    ax.set_title("GSE179730: OCSCC response check", loc="left", fontsize=8, fontweight="bold", pad=6)
    ax.text(
        0.00,
        1.005,
        "z-score modules; family composite is exact-null and rank-sensitive",
        transform=ax.transAxes,
        fontsize=5.9,
        color="#444444",
        va="bottom",
    )
    ax.set_xlim(min(-0.36, float(df["ci95_low"].min()) - 0.04), max(0.82, float(df["ci95_high"].max()) + 0.20))
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="x", color="#dddddd", lw=0.4)


def draw_interpretation(ax: plt.Axes) -> None:
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.text(-0.06, 1.04, "d", fontsize=10, fontweight="bold", va="top", transform=ax.transAxes)
    ax.text(0.02, 1.04, "Evidence roles and claim boundary", fontsize=8, fontweight="bold", va="top", transform=ax.transAxes)
    blocks = [
        (
            0.69,
            "Global-shift-sensitive response association",
            "GSE281729: the z-score family slope attenuates after response-blind global-PC adjustment; no single MCP-counter population explains the shift",
            "#1f1f1f",
        ),
        (0.40, "Independent response check, exact-null", "GSE179730: module direction under z scoring; family score is null and scale-sensitive", "#666666"),
        (0.11, "Generalization boundary", "GSE301741 and GSE123813: no validated clinical prediction claim", "#1f1f1f"),
    ]
    for y0, header, body, color in blocks:
        ax.add_patch(Rectangle((0.02, y0 - 0.01), 0.025, 0.11, color=color, transform=ax.transAxes, clip_on=False))
        ax.text(0.08, y0 + 0.075, header, fontsize=7.2, fontweight="bold", color="#222222", va="top")
        ax.text(0.08, y0 - 0.005, body, fontsize=6.2, color="#444444", va="top", wrap=True)
    ax.axis("off")


def main() -> None:
    gse281729 = load_gse281729()
    gse179730 = load_gse179730()
    gse281729.to_csv(SRC_DIR / "Figure3_gse281729_panel_source.csv", index=False)
    gse179730.to_csv(SRC_DIR / "Figure3_gse179730_panel_source.csv", index=False)

    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 7,
            "axes.linewidth": 0.6,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    fig = plt.figure(figsize=(7.2, 6.0))
    gs = GridSpec(3, 2, figure=fig, height_ratios=[0.52, 2.20, 1.72], width_ratios=[1.48, 1.0], hspace=0.50, wspace=0.38)
    ax_a = fig.add_subplot(gs[0, :])
    ax_b = fig.add_subplot(gs[1, :])
    ax_c = fig.add_subplot(gs[2, 0])
    ax_d = fig.add_subplot(gs[2, 1])
    draw_cohort_matrix(ax_a)
    draw_gse281729(ax_b, gse281729)
    draw_gse179730(ax_c, gse179730)
    draw_interpretation(ax_d)
    fig.subplots_adjust(top=0.97, bottom=0.08, left=0.22, right=0.98)
    fig.savefig(OUT_DIR / "Figure3_submission_external_validation.png", dpi=450, bbox_inches="tight")
    fig.savefig(OUT_DIR / "Figure3_submission_external_validation.pdf", bbox_inches="tight")
    fig.savefig(OUT_DIR / "Figure3_submission_external_validation.svg", bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
