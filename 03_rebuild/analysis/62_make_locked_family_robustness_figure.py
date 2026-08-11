from __future__ import annotations

import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


WORKSPACE = Path(os.environ.get("GSE200996_WORKSPACE", Path.cwd())).resolve()
VALIDATION_DIR = (
    WORKSPACE / "03_rebuild" / "validation" / "locked_family_robustness"
)
TEST_PATH = VALIDATION_DIR / "LOCKED_FAMILY_TESTS.csv"
RANDOM_PATH = VALIDATION_DIR / "LOCKED_FAMILY_MATCHED_RANDOM_EFFECTS.csv"
EXACT_PATH = VALIDATION_DIR / "GSE179730_LOCKED_FAMILY_EXACT_NULL.csv"
OUT_DIR = WORKSPACE / "03_rebuild" / "figures" / "submission"
OUT_DIR.mkdir(parents=True, exist_ok=True)
OUT_STEM = OUT_DIR / "ExtendedData11_submission_locked_family_robustness"


COLORS = {
    "null": "#b8b8b8",
    "null_edge": "#ffffff",
    "observed": "#111111",
    "z_score": "#3b6fb6",
    "rank_mean": "#b4584b",
}


def get_test(tests: pd.DataFrame, cohort: str, method: str) -> pd.Series:
    row = tests[
        tests["cohort"].eq(cohort) & tests["scoring_method"].eq(method)
    ]
    if len(row) != 1:
        raise ValueError(f"Expected one test row for {cohort}/{method}; found {len(row)}")
    return row.iloc[0]


def clean_axis(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(width=0.6, length=3)
    ax.grid(axis="y", color="#e2e2e2", linewidth=0.45, zorder=0)


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.14,
        1.08,
        label,
        transform=ax.transAxes,
        fontsize=10,
        fontweight="bold",
        va="top",
    )


def draw_exact_null(
    ax: plt.Axes,
    exact: pd.DataFrame,
    test: pd.Series,
    method: str,
    label: str,
) -> None:
    values = exact.loc[exact["scoring_method"].eq(method), "effect"].to_numpy(float)
    observed = float(test["effect"])
    ax.hist(
        values,
        bins=24,
        color=COLORS["null"],
        edgecolor=COLORS["null_edge"],
        linewidth=0.35,
        zorder=2,
    )
    ax.axvline(observed, color=COLORS[method], linewidth=1.5, zorder=4)
    ax.axvline(-observed, color=COLORS[method], linewidth=0.8, linestyle="--", zorder=3)
    ax.set_title(
        f"GSE179730: {label} exact null",
        loc="left",
        fontsize=8,
        fontweight="bold",
        pad=5,
    )
    ax.text(
        0.98,
        0.95,
        f"observed = {observed:.4f}\nexact P = {float(test['p_value']):.3f}",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=6.4,
        color="#333333",
    )
    ax.set_xlabel("Responder - non-responder family effect")
    ax.set_ylabel("Assignments")
    clean_axis(ax)


def draw_random_null(
    ax: plt.Axes,
    random_effects: pd.DataFrame,
    test: pd.Series,
    method: str,
    label: str,
) -> None:
    values = random_effects.loc[
        random_effects["cohort"].eq("GSE281729")
        & random_effects["scoring_method"].eq(method),
        "effect",
    ].to_numpy(float)
    observed = float(test["effect"])
    ax.hist(
        values,
        bins=30,
        color=COLORS["null"],
        edgecolor=COLORS["null_edge"],
        linewidth=0.35,
        zorder=2,
    )
    ax.axvline(observed, color=COLORS[method], linewidth=1.5, zorder=4)
    ax.axvspan(
        float(test["matched_random_q025"]),
        float(test["matched_random_q975"]),
        color="#777777",
        alpha=0.12,
        linewidth=0,
        zorder=1,
    )
    ax.set_title(
        f"GSE281729: {label} matched null",
        loc="left",
        fontsize=8,
        fontweight="bold",
        pad=5,
    )
    ax.text(
        0.98,
        0.95,
        (
            f"observed = {observed:.4f}\n"
            f"HC3 P = {float(test['p_value']):.3f}\n"
            f"empirical P = {float(test['empirical_specificity_p']):.4f}"
        ),
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=6.4,
        color="#333333",
    )
    ax.set_xlabel("Adjusted ordinal response slope")
    ax.set_ylabel("Matched random families")
    clean_axis(ax)


def main() -> None:
    tests = pd.read_csv(TEST_PATH)
    random_effects = pd.read_csv(RANDOM_PATH)
    exact = pd.read_csv(EXACT_PATH)

    expected_random = 2_000
    counts = random_effects.groupby(["cohort", "scoring_method"]).size()
    if not counts.eq(expected_random).all():
        raise ValueError(f"Matched-random iteration count mismatch: {counts.to_dict()}")
    exact_counts = exact.groupby("scoring_method").size()
    if not exact_counts.eq(462).all():
        raise ValueError(f"Exact-null assignment count mismatch: {exact_counts.to_dict()}")

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
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.5))
    draw_exact_null(
        axes[0, 0],
        exact,
        get_test(tests, "GSE179730", "z_score"),
        "z_score",
        "z-score",
    )
    draw_exact_null(
        axes[0, 1],
        exact,
        get_test(tests, "GSE179730", "rank_mean"),
        "rank_mean",
        "rank-mean",
    )
    draw_random_null(
        axes[1, 0],
        random_effects,
        get_test(tests, "GSE281729", "z_score"),
        "z_score",
        "z-score",
    )
    draw_random_null(
        axes[1, 1],
        random_effects,
        get_test(tests, "GSE281729", "rank_mean"),
        "rank_mean",
        "rank-mean",
    )
    for ax, label in zip(axes.flat, "abcd"):
        panel_label(ax, label)
    fig.suptitle(
        "Locked 16-module family: exact and matched-random robustness",
        x=0.08,
        y=0.995,
        ha="left",
        fontsize=9,
        fontweight="bold",
    )
    fig.subplots_adjust(left=0.10, right=0.985, bottom=0.10, top=0.91, wspace=0.30, hspace=0.42)
    fig.savefig(OUT_STEM.with_suffix(".png"), dpi=450)
    fig.savefig(OUT_STEM.with_suffix(".pdf"))
    fig.savefig(OUT_STEM.with_suffix(".svg"))
    plt.close(fig)


if __name__ == "__main__":
    main()
