from __future__ import annotations

import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


WORKSPACE = Path(os.environ.get("GSE200996_WORKSPACE", Path.cwd())).resolve()
VALIDATION = (
    WORKSPACE / "03_rebuild" / "validation" / "family_composite_stress_test"
)
TESTS_PATH = VALIDATION / "FAMILY_VARIANT_TESTS.csv"
WEIGHTS_PATH = VALIDATION / "FAMILY_VARIANT_WEIGHT_AUDIT.csv"
OUT_DIR = WORKSPACE / "03_rebuild" / "figures" / "submission"
SOURCE_DIR = OUT_DIR / "source_data"
OUT_DIR.mkdir(parents=True, exist_ok=True)
SOURCE_DIR.mkdir(parents=True, exist_ok=True)
OUT_STEM = OUT_DIR / "ExtendedData12_submission_family_composite_stress_test"


VARIANT_ORDER = [
    "module_mean_16",
    "unique_gene_equal",
    "inverse_membership_module_mean",
    "no_union_module_mean",
    "hallmark_only_module_mean",
    "dynamic_only_module_mean",
]
VARIANT_LABELS = {
    "module_mean_16": "Module mean (16)",
    "unique_gene_equal": "Unique-gene mean",
    "inverse_membership_module_mean": "Inverse membership",
    "no_union_module_mean": "No union modules",
    "hallmark_only_module_mean": "Hallmark only",
    "dynamic_only_module_mean": "Dynamic only",
}
COLORS = {
    False: "#3B6FB6",
    True: "#B4584B",
}


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.14,
        1.08,
        label,
        transform=ax.transAxes,
        fontsize=8,
        fontweight="bold",
        va="top",
    )


def clean_axis(ax: plt.Axes, grid_axis: str = "x") -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(width=0.6, length=3)
    ax.grid(axis=grid_axis, color="#E2E2E2", linewidth=0.45, zorder=0)


def draw_weight_audit(ax: plt.Axes, weights: pd.DataFrame) -> None:
    selected = (
        weights[
            (weights["cohort"] == "GSE281729")
            & weights["score_variant"].isin(VARIANT_ORDER[:3])
        ]
        .set_index("score_variant")
        .loc[VARIANT_ORDER[:3]]
    )
    y = np.arange(len(selected))
    bars = ax.barh(
        y,
        selected["effective_gene_number"],
        color=["#3B6FB6", "#579A68", "#B4584B"],
        height=0.58,
        zorder=2,
    )
    ax.set_yticks(y, [VARIANT_LABELS[name] for name in selected.index])
    ax.invert_yaxis()
    ax.set_xlabel("Effective number of genes")
    ax.set_title(
        "GSE281729: composite weight concentration",
        loc="left",
        fontsize=8,
        fontweight="bold",
        pad=5,
    )
    for bar, ratio in zip(bars, selected["max_to_min_weight_ratio"]):
        ax.text(
            bar.get_width() + 8,
            bar.get_y() + bar.get_height() / 2,
            f"max:min {ratio:.1f}",
            va="center",
            fontsize=6.2,
            color="#333333",
        )
    ax.set_xlim(0, 490)
    clean_axis(ax)


def draw_effects(
    ax: plt.Axes,
    tests: pd.DataFrame,
    method: str,
    title: str,
) -> None:
    selected = tests[
        (tests["cohort"] == "GSE281729")
        & (tests["scoring_method"] == method)
    ].copy()
    y = np.arange(len(VARIANT_ORDER))
    offsets = {False: -0.14, True: 0.14}
    for adjusted in [False, True]:
        current = (
            selected[selected["background_adjusted"] == adjusted]
            .set_index("score_variant")
            .loc[VARIANT_ORDER]
        )
        ax.errorbar(
            current["effect"],
            y + offsets[adjusted],
            xerr=np.vstack(
                [
                    current["effect"] - current["ci95_low"],
                    current["ci95_high"] - current["effect"],
                ]
            ),
            fmt="o",
            markersize=3.7,
            color=COLORS[adjusted],
            ecolor=COLORS[adjusted],
            elinewidth=0.8,
            capsize=2,
            label=(
                "Matched-background adjusted"
                if adjusted
                else "Unadjusted for background"
            ),
            zorder=3,
        )
    ax.axvline(0, color="#777777", linewidth=0.65, zorder=1)
    ax.set_yticks(y, [VARIANT_LABELS[name] for name in VARIANT_ORDER])
    ax.invert_yaxis()
    ax.set_xlabel("Adjusted ordinal response slope")
    ax.set_title(title, loc="left", fontsize=8, fontweight="bold", pad=5)
    clean_axis(ax)


def draw_correlations(ax: plt.Axes, tests: pd.DataFrame) -> None:
    correlations = (
        tests[tests["background_adjusted"] == False]  # noqa: E712
        .pivot_table(
            index="score_variant",
            columns=["cohort", "scoring_method"],
            values="family_background_delta_correlation",
            aggfunc="first",
        )
        .loc[VARIANT_ORDER]
    )
    column_order = [
        ("GSE179730", "z_score"),
        ("GSE179730", "rank_mean"),
        ("GSE281729", "z_score"),
        ("GSE281729", "rank_mean"),
    ]
    correlations = correlations.loc[:, column_order]
    image = ax.imshow(
        correlations.to_numpy(float),
        cmap="coolwarm",
        vmin=-1,
        vmax=1,
        aspect="auto",
    )
    ax.set_yticks(
        np.arange(len(VARIANT_ORDER)),
        [VARIANT_LABELS[name] for name in VARIANT_ORDER],
    )
    ax.set_xticks(
        np.arange(4),
        ["GSE179730\nz-score", "GSE179730\nrank", "GSE281729\nz-score", "GSE281729\nrank"],
    )
    ax.tick_params(length=0)
    for row in range(correlations.shape[0]):
        for column in range(correlations.shape[1]):
            value = correlations.iloc[row, column]
            ax.text(
                column,
                row,
                f"{value:.2f}",
                ha="center",
                va="center",
                fontsize=6.2,
                color="white" if abs(value) > 0.62 else "#222222",
            )
    ax.set_title(
        "Family versus matched-background delta",
        loc="left",
        fontsize=8,
        fontweight="bold",
        pad=5,
    )
    colorbar = plt.colorbar(image, ax=ax, fraction=0.032, pad=0.025, shrink=0.88)
    colorbar.set_label("Pearson r", fontsize=6.5)
    colorbar.ax.tick_params(labelsize=6, width=0.5, length=2)


def main() -> None:
    tests = pd.read_csv(TESTS_PATH)
    weights = pd.read_csv(WEIGHTS_PATH)
    if len(tests) != 48:
        raise ValueError(f"Expected 48 test rows, observed {len(tests)}")

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
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.7))
    draw_weight_audit(axes[0, 0], weights)
    draw_effects(
        axes[0, 1],
        tests,
        "z_score",
        "GSE281729: z-score family variants",
    )
    draw_effects(
        axes[1, 0],
        tests,
        "rank_mean",
        "GSE281729: rank-based family variants",
    )
    draw_correlations(axes[1, 1], tests)
    for ax, label in zip(axes.flat, "abcd"):
        panel_label(ax, label)
    handles, labels = axes[0, 1].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.995),
        frameon=False,
        ncol=2,
        fontsize=6.5,
        handletextpad=0.5,
        columnspacing=1.5,
    )
    fig.subplots_adjust(
        left=0.17,
        right=0.955,
        bottom=0.10,
        top=0.93,
        wspace=0.42,
        hspace=0.42,
    )
    fig.savefig(OUT_STEM.with_suffix(".png"), dpi=450)
    fig.savefig(OUT_STEM.with_suffix(".pdf"))
    fig.savefig(OUT_STEM.with_suffix(".svg"))
    plt.close(fig)

    tests.to_csv(
        SOURCE_DIR / "ExtendedData12_family_variant_tests_source.csv",
        index=False,
    )
    weights.to_csv(
        SOURCE_DIR / "ExtendedData12_family_weight_audit_source.csv",
        index=False,
    )


if __name__ == "__main__":
    main()
