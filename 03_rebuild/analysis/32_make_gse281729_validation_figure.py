from __future__ import annotations

import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


WORKSPACE = Path(os.environ.get("GSE200996_WORKSPACE", Path.cwd())).resolve()
IN_DIR = WORKSPACE / "03_rebuild" / "validation" / "GSE281729_bulk_module_validation"
OUT_DIR = WORKSPACE / "03_rebuild" / "figures" / "external_validation"
SOURCE_DIR = OUT_DIR / "source_data"
OUT_DIR.mkdir(parents=True, exist_ok=True)
SOURCE_DIR.mkdir(parents=True, exist_ok=True)


def short_label(signature: str) -> str:
    labels = {
        "M_LE_INTERFERON_ALPHA_RESPONSE": "Myeloid IFN-alpha",
        "T_LE_INTERFERON_ALPHA_RESPONSE": "T cell IFN-alpha",
        "M_LE_INTERFERON_GAMMA_RESPONSE": "Myeloid IFN-gamma",
        "T_LE_INTERFERON_GAMMA_RESPONSE": "T cell IFN-gamma",
        "M_LE_MTORC1_SIGNALING": "Myeloid mTORC1",
        "T_LE_MTORC1_SIGNALING": "T cell mTORC1",
        "M_LE_union_core": "Myeloid union",
        "T_LE_union_core": "T cell union",
        "M_LE_TNFA_SIGNALING_VIA_NFKB": "Myeloid TNFA/NFKB",
        "T_LE_TNFA_SIGNALING_VIA_NFKB": "T cell TNFA/NFKB",
        "M_LE_INFLAMMATORY_RESPONSE": "Myeloid inflammatory",
        "M_LE_COMPLEMENT": "Myeloid complement",
    }
    return labels.get(signature, signature.replace("_", " "))


def main() -> int:
    stats = pd.read_csv(IN_DIR / "GSE281729_PAIRED_DELTA_MODULE_RESPONSE_MODELS.csv")
    scores = pd.read_csv(IN_DIR / "GSE281729_LOCKED_MODULE_SAMPLE_SCORES.csv")

    primary = stats[
        (stats["model"] == "ordinal_unadjusted")
        & (stats["model_status"] == "ok")
        & stats["signature"].str.contains("_LE_")
    ].copy()
    priority = [
        "M_LE_INTERFERON_ALPHA_RESPONSE",
        "T_LE_INTERFERON_ALPHA_RESPONSE",
        "M_LE_INTERFERON_GAMMA_RESPONSE",
        "T_LE_INTERFERON_GAMMA_RESPONSE",
        "M_LE_MTORC1_SIGNALING",
        "T_LE_MTORC1_SIGNALING",
        "M_LE_union_core",
        "T_LE_union_core",
        "M_LE_TNFA_SIGNALING_VIA_NFKB",
        "M_LE_INFLAMMATORY_RESPONSE",
        "M_LE_COMPLEMENT",
    ]
    primary = primary[primary["signature"].isin(priority)].copy()
    primary["label"] = primary["signature"].map(short_label)
    primary["order"] = primary["signature"].map({sig: i for i, sig in enumerate(priority)})
    primary = primary.sort_values("order", ascending=False)
    primary.to_csv(SOURCE_DIR / "GSE281729_validation_slope_source.csv", index=False)

    wide = scores.pivot_table(
        index=["patient_id", "signature"],
        columns="timepoint",
        values="module_score",
        aggfunc="mean",
    ).reset_index()
    meta = scores[
        [
            "patient_id",
            "signature",
            "response_harmonized_ordinal",
            "response_ord_num",
            "response_binary",
            "hpv",
            "second_drug",
        ]
    ].drop_duplicates(["patient_id", "signature"])
    wide = wide.merge(meta, on=["patient_id", "signature"], how="left")
    wide["post_minus_pre"] = wide["post"] - wide["pre"]
    selected = wide[
        wide["signature"].isin(["M_LE_INTERFERON_ALPHA_RESPONSE", "T_LE_INTERFERON_ALPHA_RESPONSE"])
        & wide["response_harmonized_ordinal"].isin(["Low", "Medium", "High"])
    ].copy()
    selected["label"] = selected["signature"].map(short_label)
    selected["response_harmonized_ordinal"] = pd.Categorical(
        selected["response_harmonized_ordinal"],
        categories=["Low", "Medium", "High"],
        ordered=True,
    )
    selected.to_csv(SOURCE_DIR / "GSE281729_validation_delta_distribution_source.csv", index=False)

    sns.set_theme(style="white", context="paper")
    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 8,
            "axes.linewidth": 0.8,
            "axes.labelsize": 8,
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
            "legend.fontsize": 7,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )

    fig = plt.figure(figsize=(8.2, 4.5))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.15, 1.0], wspace=0.48)
    ax0 = fig.add_subplot(gs[0, 0])
    ax1 = fig.add_subplot(gs[0, 1])

    colors = primary["target_lineage"].map({"T_cell": "#2F6C8E", "Myeloid": "#A95735"}).fillna("#606060")
    ax0.barh(primary["label"], primary["coef"], color=colors, height=0.72)
    ax0.axvline(0, color="#222222", linewidth=0.8)
    ax0.set_xlabel("Slope per response-depth step")
    ax0.set_ylabel("")
    ax0.set_title("Locked module association", loc="left", fontsize=9, fontweight="bold")
    ax0.set_xlim(min(-0.68, primary["coef"].min() - 0.08), 0.08)
    for spine in ["top", "right"]:
        ax0.spines[spine].set_visible(False)
    for i, (_, row) in enumerate(primary.iterrows()):
        if pd.notna(row["fdr"]) and row["fdr"] < 0.10:
            ax0.text(0.02, i, f"{row['fdr']:.2g}", ha="left", va="center", fontsize=6.5)

    palette = {"Low": "#606060", "Medium": "#9D8B3E", "High": "#B94F4F"}
    sns.boxplot(
        data=selected,
        x="response_harmonized_ordinal",
        y="post_minus_pre",
        hue="label",
        ax=ax1,
        width=0.62,
        fliersize=0,
        palette=["#A95735", "#2F6C8E"],
    )
    sns.stripplot(
        data=selected,
        x="response_harmonized_ordinal",
        y="post_minus_pre",
        hue="label",
        dodge=True,
        ax=ax1,
        size=3.2,
        linewidth=0.35,
        edgecolor="white",
        palette=["#A95735", "#2F6C8E"],
        alpha=0.88,
    )
    handles, labels = ax1.get_legend_handles_labels()
    ax1.legend(handles[:2], labels[:2], frameon=False, loc="upper right", ncol=1, title=None)
    ax1.axhline(0, color="#222222", linewidth=0.8)
    ax1.set_xlabel("Pathological response")
    ax1.set_ylabel("Post-pre module score")
    ax1.set_title("Representative IFN-alpha deltas", loc="left", fontsize=9, fontweight="bold")
    for spine in ["top", "right"]:
        ax1.spines[spine].set_visible(False)

    fig.text(0.012, 0.975, "a", fontsize=10, fontweight="bold", va="top")
    fig.text(0.555, 0.975, "b", fontsize=10, fontweight="bold", va="top")
    fig.savefig(OUT_DIR / "GSE281729_bulk_validation_locked_modules.png", dpi=450, bbox_inches="tight")
    fig.savefig(OUT_DIR / "GSE281729_bulk_validation_locked_modules.pdf", bbox_inches="tight")
    plt.close(fig)
    print(OUT_DIR / "GSE281729_bulk_validation_locked_modules.png")
    print(OUT_DIR / "GSE281729_bulk_validation_locked_modules.pdf")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
