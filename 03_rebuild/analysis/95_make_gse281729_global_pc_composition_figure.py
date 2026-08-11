from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


ROOT = Path(__file__).resolve().parents[2]
REBUILD = ROOT / "03_rebuild"
VALIDATION = REBUILD / "validation" / "gse281729_global_pc_composition"
OUT_DIR = REBUILD / "figures" / "submission"
SOURCE_DIR = OUT_DIR / "source_data"
OUT_DIR.mkdir(parents=True, exist_ok=True)
SOURCE_DIR.mkdir(parents=True, exist_ok=True)

PC_VARIANCE = VALIDATION / "GSE281729_GLOBAL_PC_VARIANCE.csv"
FAMILY_MODELS = VALIDATION / "GSE281729_FAMILY_GLOBAL_PC_COMPOSITION_MODELS.csv"
MCP_MODELS = VALIDATION / "GSE281729_MCP_COUNTER_RESPONSE_MODELS.csv"
CORRELATIONS = VALIDATION / "GSE281729_FAMILY_PC_COMPOSITION_CORRELATIONS.csv"
PC_QC = VALIDATION / "GSE281729_GLOBAL_PC_QC_CORRELATIONS.csv"
PAIRED_QC = VALIDATION / "GSE281729_PROCESSED_EXPRESSION_PAIRED_QC_DELTAS.csv"
PC_SCORES = VALIDATION / "GSE281729_GLOBAL_PC_PATIENT_SCORES.csv"
FAMILY_DELTAS = VALIDATION / "GSE281729_FAMILY_GLOBAL_PC_COMPOSITION_DELTAS.csv"

STEM = OUT_DIR / "ExtendedData13_submission_gse281729_global_pc_composition"

BLUE = "#3B6FB6"
RED = "#B95848"
GREEN = "#5B9E6F"
GREY = "#5F6368"
LIGHT_GREY = "#D9D9D9"


MODEL_LABELS = {
    "base": "Base",
    "global_pc1": "+ global PC1",
    "global_pc1_pc2": "+ global PCs 1-2",
    "composition_pc1": "+ composition PC1",
    "global_pc1_composition_pc1": "+ global PC1\n+ composition PC1",
}
MODEL_ORDER = list(MODEL_LABELS)


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.15,
        1.08,
        label,
        transform=ax.transAxes,
        fontsize=8,
        fontweight="bold",
        va="top",
    )


def forest_panel(
    ax: plt.Axes,
    frame: pd.DataFrame,
    method: str,
    color: str,
    title: str,
) -> None:
    current = (
        frame[
            frame["scoring_method"].eq(method)
            & frame["score_variant"].eq("module_mean_16")
        ]
        .set_index("model_spec")
        .loc[MODEL_ORDER]
        .reset_index()
    )
    y = np.arange(len(current))[::-1]
    ax.axvline(0, color="#777777", lw=0.7, zorder=0)
    ax.errorbar(
        current["effect"],
        y,
        xerr=np.vstack(
            [
                current["effect"] - current["ci95_low"],
                current["ci95_high"] - current["effect"],
            ]
        ),
        fmt="o",
        color=color,
        ecolor=color,
        markersize=4.5,
        capsize=2.5,
        lw=1.1,
    )
    ax.set_yticks(y)
    ax.set_yticklabels([MODEL_LABELS[value] for value in current["model_spec"]], fontsize=6.5)
    ax.set_xlabel("Adjusted ordinal response slope")
    ax.set_title(title, loc="left", fontsize=8, fontweight="bold", pad=5)
    ax.grid(axis="x", color=LIGHT_GREY, lw=0.45)
    ax.spines[["top", "right"]].set_visible(False)


def main() -> None:
    pc_variance = pd.read_csv(PC_VARIANCE)
    family_models = pd.read_csv(FAMILY_MODELS)
    mcp_models = pd.read_csv(MCP_MODELS)
    correlations = pd.read_csv(CORRELATIONS)
    pc_qc = pd.read_csv(PC_QC)
    paired_qc = pd.read_csv(PAIRED_QC)
    pc_scores = pd.read_csv(PC_SCORES)
    family_deltas = pd.read_csv(FAMILY_DELTAS)

    pc_variance.to_csv(SOURCE_DIR / "ExtendedData13_pc_variance_source.csv", index=False)
    family_models.to_csv(SOURCE_DIR / "ExtendedData13_family_models_source.csv", index=False)
    mcp_models.to_csv(SOURCE_DIR / "ExtendedData13_mcp_models_source.csv", index=False)
    correlations.to_csv(SOURCE_DIR / "ExtendedData13_family_correlations_source.csv", index=False)
    pc_qc.to_csv(SOURCE_DIR / "ExtendedData13_pc_qc_correlations_source.csv", index=False)

    primary_scores = pc_scores[
        pc_scores["gene_universe"].eq("primary_nonlocked_mean_ge_1")
    ][["patient_id", "global_pc1"]]
    scatter = (
        paired_qc[["patient_id", "delta_mean_logtpm"]]
        .merge(primary_scores, on="patient_id", validate="one_to_one")
        .merge(
            family_deltas[
                family_deltas["score_variant"].eq("module_mean_16")
                & family_deltas["scoring_method"].eq("z_score")
            ][["patient_id", "response_harmonized_ordinal"]].drop_duplicates(),
            on="patient_id",
            validate="one_to_one",
        )
    )
    scatter.to_csv(SOURCE_DIR / "ExtendedData13_pc1_qc_scatter_source.csv", index=False)

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
    fig, axes = plt.subplots(2, 3, figsize=(7.2, 5.7))
    ax_a, ax_b, ax_c, ax_d, ax_e, ax_f = axes.ravel()

    current_variance = pc_variance[pc_variance["pc"].isin(["global_pc1", "global_pc2", "global_pc3"])].copy()
    current_variance["label"] = current_variance["gene_universe"].map(
        {
            "primary_nonlocked_mean_ge_1": "Primary non-locked",
            "sensitivity_top5000_delta_variance": "Top-5,000 sensitivity",
        }
    )
    x = np.arange(3)
    width = 0.34
    for offset, (label, color) in enumerate(
        [("Primary non-locked", BLUE), ("Top-5,000 sensitivity", GREEN)]
    ):
        values = (
            current_variance[current_variance["label"].eq(label)]
            .set_index("pc")
            .loc[["global_pc1", "global_pc2", "global_pc3"], "explained_variance_ratio"]
            * 100
        )
        ax_a.bar(x + (offset - 0.5) * width, values, width=width, color=color, label=label)
    ax_a.set_xticks(x)
    ax_a.set_xticklabels(["PC1", "PC2", "PC3"])
    ax_a.set_ylabel("Explained paired-delta variance (%)")
    ax_a.set_title("Response-blind global PCs", loc="left", fontsize=8, fontweight="bold")
    ax_a.legend(frameon=False, fontsize=6, loc="upper right")
    ax_a.spines[["top", "right"]].set_visible(False)
    panel_label(ax_a, "a")

    forest_panel(ax_b, family_models, "z_score", BLUE, "Frozen family: z-score")
    panel_label(ax_b, "b")
    forest_panel(ax_c, family_models, "rank_mean", RED, "Frozen family: rank sensitivity")
    panel_label(ax_c, "c")

    mcp_order = (
        mcp_models.assign(abs_effect=lambda data: data["effect"].abs())
        .sort_values("abs_effect", ascending=True)["population"]
        .tolist()
    )
    current_mcp = mcp_models.set_index("population").loc[mcp_order].reset_index()
    y = np.arange(len(current_mcp))
    ax_d.axvline(0, color="#777777", lw=0.7)
    ax_d.errorbar(
        current_mcp["effect"],
        y,
        xerr=np.vstack(
            [
                current_mcp["effect"] - current_mcp["ci95_low"],
                current_mcp["ci95_high"] - current_mcp["effect"],
            ]
        ),
        fmt="o",
        color=GREY,
        ecolor=GREY,
        markersize=3.8,
        capsize=2.2,
        lw=1.0,
    )
    ax_d.set_yticks(y)
    population_labels = {
        "B lineage": "B lineage",
        "CD8 T cells": "CD8 T",
        "Fibroblasts": "Fibroblasts",
        "Myeloid dendritic cells": "Myeloid DC",
        "Cytotoxic lymphocytes": "Cytotoxic",
        "Monocytic lineage": "Monocytic",
        "T cells": "T cells",
        "NK cells": "NK",
        "Endothelial cells": "Endothelial",
        "Neutrophils": "Neutrophils",
    }
    ax_d.set_yticklabels(
        [population_labels.get(value, value) for value in current_mcp["population"]],
        fontsize=5.9,
    )
    ax_d.set_xlabel("Adjusted ordinal response slope")
    ax_d.set_title("MCP-counter composition deltas", loc="left", fontsize=8, fontweight="bold")
    ax_d.grid(axis="x", color=LIGHT_GREY, lw=0.45)
    ax_d.spines[["top", "right"]].set_visible(False)
    panel_label(ax_d, "d")

    selected_correlates = [
        "global_pc1",
        "global_pc2",
        "global_pc3",
        "composition_pc1",
        "Cytotoxic lymphocytes",
        "Monocytic lineage",
    ]
    heat = (
        correlations[correlations["correlate"].isin(selected_correlates)]
        .pivot(index="scoring_method", columns="correlate", values="pearson_correlation")
        .loc[["z_score", "rank_mean"], selected_correlates]
    )
    heat.index = ["z-score", "rank"]
    heat.columns = ["Global PC1", "Global PC2", "Global PC3", "Comp. PC1", "Cytotoxic", "Monocytic"]
    sns.heatmap(
        heat,
        ax=ax_e,
        cmap="coolwarm",
        vmin=-1,
        vmax=1,
        center=0,
        annot=True,
        fmt=".2f",
        annot_kws={"fontsize": 6.5},
        cbar_kws={"label": "Pearson r", "fraction": 0.04, "pad": 0.025},
        linewidths=0.5,
        linecolor="white",
    )
    ax_e.set_xlabel("")
    ax_e.set_ylabel("")
    ax_e.set_title("Family delta correlations", loc="left", fontsize=8, fontweight="bold")
    ax_e.tick_params(axis="x", labelsize=5.6, rotation=35)
    for label in ax_e.get_xticklabels():
        label.set_horizontalalignment("right")
    ax_e.tick_params(axis="y", labelsize=6.3, rotation=0)
    panel_label(ax_e, "e")

    response_colors = {"Low": "#777777", "Medium": "#D09B3D", "High": "#B95848"}
    for response, current in scatter.groupby("response_harmonized_ordinal"):
        ax_f.scatter(
            current["delta_mean_logtpm"],
            current["global_pc1"],
            s=22,
            color=response_colors[response],
            edgecolor="white",
            linewidth=0.4,
            label=response,
            zorder=3,
        )
    coefficients = np.polyfit(scatter["delta_mean_logtpm"], scatter["global_pc1"], 1)
    xx = np.linspace(scatter["delta_mean_logtpm"].min(), scatter["delta_mean_logtpm"].max(), 100)
    ax_f.plot(xx, coefficients[0] * xx + coefficients[1], color="#333333", lw=0.9)
    correlation = scatter["delta_mean_logtpm"].corr(scatter["global_pc1"])
    ax_f.text(
        0.04,
        0.96,
        f"Pearson r = {correlation:.2f}",
        transform=ax_f.transAxes,
        va="top",
        fontsize=6.5,
    )
    ax_f.set_xlabel("Post-pre mean logTPM")
    ax_f.set_ylabel("Global PC1 score")
    ax_f.set_title("Global PC1 marks logTPM shift", loc="left", fontsize=8, fontweight="bold")
    ax_f.legend(frameon=False, fontsize=6, loc="lower right", title="Response", title_fontsize=6)
    ax_f.spines[["top", "right"]].set_visible(False)
    panel_label(ax_f, "f")

    fig.subplots_adjust(
        left=0.115,
        right=0.985,
        top=0.94,
        bottom=0.09,
        wspace=0.58,
        hspace=0.48,
    )
    fig.savefig(STEM.with_suffix(".png"), dpi=450)
    fig.savefig(STEM.with_suffix(".pdf"))
    fig.savefig(STEM.with_suffix(".svg"))
    plt.close(fig)


if __name__ == "__main__":
    main()
