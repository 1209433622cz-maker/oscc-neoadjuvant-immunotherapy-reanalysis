#!/usr/bin/env python
"""Orthogonal external TCR validation in GSE123813.

This analysis uses the BCC and SCC TCR/metadata components of GSE123813 to test
whether an independent anti-PD-1 setting shows treatment-associated T-cell state
and clonotype remodeling. It is intentionally not framed as OSCC response-depth
validation.
"""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

try:
    from scipy import stats as scipy_stats
except Exception:  # noqa: BLE001
    scipy_stats = None


SCRIPT_PATH = Path(__file__).resolve()
REBUILD_DIR = SCRIPT_PATH.parents[1]
ROOT = REBUILD_DIR.parent
RAW_DIR = ROOT / "00_raw_data" / "GSE123813_validation"
RESULTS_DIR = REBUILD_DIR / "results" / "external_tcr_validation"
FIGURE_DIR = REBUILD_DIR / "figures" / "submission"
SOURCE_DATA_DIR = FIGURE_DIR / "source_data"

FIG_STEM = "ExtendedData6_external_TCR_validation_GSE123813"
MIN_TCR_CELLS_PER_TIMEPOINT = 50

CD8_STATE_CLUSTERS = {"CD8_ex", "CD8_act", "CD8_ex_act", "CD8_eff"}
CD8_EX_CLUSTERS = {"CD8_ex", "CD8_ex_act"}
CLONOTYPE_COL = "cdr3s_aa"


@dataclass(frozen=True)
class DatasetSpec:
    disease: str
    metadata_file: str
    tcr_file: str


DATASETS = [
    DatasetSpec("BCC", "GSE123813_bcc_tcell_metadata.txt.gz", "GSE123813_bcc_tcr.txt.gz"),
    DatasetSpec("SCC", "GSE123813_scc_metadata.txt.gz", "GSE123813_scc_tcr.txt.gz"),
]

DELTA_METRICS = [
    "cd8_state_fraction",
    "cd8_ex_fraction",
    "expanded_cell_fraction",
    "hyperexpanded_cell_fraction",
    "top1_fraction",
    "top10_fraction",
    "normalized_shannon",
    "simpson_concentration",
]


def ensure_dirs() -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    SOURCE_DATA_DIR.mkdir(parents=True, exist_ok=True)


def bh_adjust(pvalues: list[float]) -> list[float]:
    arr = np.asarray([1.0 if pd.isna(p) else float(p) for p in pvalues], dtype=float)
    n = len(arr)
    if n == 0:
        return []
    order = np.argsort(arr)
    ranks = np.empty(n, dtype=int)
    ranks[order] = np.arange(1, n + 1)
    adjusted = arr * n / ranks
    adjusted_ordered = adjusted[order]
    for i in range(n - 2, -1, -1):
        adjusted_ordered[i] = min(adjusted_ordered[i], adjusted_ordered[i + 1])
    out = np.empty(n, dtype=float)
    out[order] = np.minimum(adjusted_ordered, 1.0)
    return out.tolist()


def two_sided_sign_test(values: pd.Series | np.ndarray, null: float = 0.0) -> float:
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    signs = arr - null
    pos = int(np.sum(signs > 0))
    neg = int(np.sum(signs < 0))
    n = pos + neg
    if n == 0:
        return math.nan
    k = min(pos, neg)
    prob = sum(math.comb(n, i) for i in range(k + 1)) / (2**n)
    return min(1.0, 2.0 * prob)


def ttest_1samp(values: pd.Series | np.ndarray, null: float = 0.0) -> tuple[float, float]:
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    if len(arr) < 2 or np.allclose(arr, arr[0]):
        return math.nan, math.nan
    if scipy_stats is None:
        return math.nan, math.nan
    res = scipy_stats.ttest_1samp(arr, popmean=null)
    return float(res.statistic), float(res.pvalue)


def wilcoxon_1samp(values: pd.Series | np.ndarray, null: float = 0.0) -> tuple[float, float]:
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)] - null
    if len(arr) < 3 or np.allclose(arr, 0):
        return math.nan, math.nan
    if scipy_stats is None:
        return math.nan, math.nan
    try:
        res = scipy_stats.wilcoxon(arr, zero_method="wilcox", alternative="two-sided")
        return float(res.statistic), float(res.pvalue)
    except ValueError:
        return math.nan, math.nan


def load_dataset(spec: DatasetSpec) -> pd.DataFrame:
    meta = pd.read_csv(RAW_DIR / spec.metadata_file, sep="\t", compression="gzip")
    tcr = (
        pd.read_csv(RAW_DIR / spec.tcr_file, sep="\t", compression="gzip", index_col=0)
        .reset_index()
        .rename(columns={"index": "cell.id"})
    )
    required = {"cell.id", "patient", "treatment", "cluster"}
    missing = required - set(meta.columns)
    if missing:
        raise ValueError(f"Missing metadata columns for {spec.disease}: {missing}")
    if CLONOTYPE_COL not in tcr.columns:
        raise ValueError(f"Missing {CLONOTYPE_COL} in {spec.tcr_file}")
    merged = meta.merge(tcr[["cell.id", CLONOTYPE_COL]], on="cell.id", how="inner")
    merged = merged[merged[CLONOTYPE_COL].notna()].copy()
    merged[CLONOTYPE_COL] = merged[CLONOTYPE_COL].astype(str)
    merged = merged[merged[CLONOTYPE_COL].str.len() > 0].copy()
    merged["disease"] = spec.disease
    merged["cd8_state"] = merged["cluster"].isin(CD8_STATE_CLUSTERS)
    merged["cd8_exhausted"] = merged["cluster"].isin(CD8_EX_CLUSTERS)
    return merged


def sample_metrics(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (disease, patient, treatment), sub in df.groupby(["disease", "patient", "treatment"], sort=True):
        counts = sub[CLONOTYPE_COL].value_counts()
        n = int(counts.sum())
        n_clones = int(len(counts))
        pvec = counts.to_numpy(dtype=float) / float(n)
        entropy = float(-np.sum(pvec * np.log(pvec))) if n_clones else math.nan
        expanded = counts[counts >= 2].index
        hyperexpanded = counts[counts >= 5].index
        rows.append(
            {
                "disease": disease,
                "patient": patient,
                "treatment": treatment,
                "n_tcr_cells": n,
                "n_clonotypes": n_clones,
                "top1_fraction": float(counts.iloc[0] / n) if n else math.nan,
                "top10_fraction": float(counts.iloc[:10].sum() / n) if n else math.nan,
                "expanded_cell_fraction": float(sub[CLONOTYPE_COL].isin(expanded).mean()) if n else math.nan,
                "hyperexpanded_cell_fraction": float(sub[CLONOTYPE_COL].isin(hyperexpanded).mean()) if n else math.nan,
                "normalized_shannon": float(entropy / math.log(n_clones)) if n_clones > 1 else math.nan,
                "simpson_concentration": float(np.sum(pvec**2)) if n else math.nan,
                "cd8_state_fraction": float(sub["cd8_state"].mean()) if n else math.nan,
                "cd8_ex_fraction": float(sub["cd8_exhausted"].mean()) if n else math.nan,
                "cd8_act_fraction": float((sub["cluster"] == "CD8_act").mean()) if n else math.nan,
            }
        )
    return pd.DataFrame(rows).sort_values(["disease", "patient", "treatment"]).reset_index(drop=True)


def paired_delta_stats(metrics: pd.DataFrame) -> pd.DataFrame:
    rows = []
    groups = [("ALL", metrics)]
    groups.extend((disease, sub) for disease, sub in metrics.groupby("disease", sort=True))
    for group_name, sub in groups:
        for metric in DELTA_METRICS:
            wide = sub.pivot_table(index=["disease", "patient"], columns="treatment", values=metric, aggfunc="mean")
            nwide = sub.pivot_table(index=["disease", "patient"], columns="treatment", values="n_tcr_cells", aggfunc="mean")
            if not {"pre", "post"}.issubset(set(wide.columns)):
                continue
            keep = (nwide["pre"] >= MIN_TCR_CELLS_PER_TIMEPOINT) & (nwide["post"] >= MIN_TCR_CELLS_PER_TIMEPOINT)
            deltas = (wide.loc[keep, "post"] - wide.loc[keep, "pre"]).dropna()
            t_stat, t_p = ttest_1samp(deltas, 0.0)
            w_stat, w_p = wilcoxon_1samp(deltas, 0.0)
            rows.append(
                {
                    "disease_group": group_name,
                    "metric": metric,
                    "n_paired_patients": int(len(deltas)),
                    "min_tcr_cells_per_timepoint": MIN_TCR_CELLS_PER_TIMEPOINT,
                    "mean_post_minus_pre": float(deltas.mean()) if len(deltas) else math.nan,
                    "median_post_minus_pre": float(deltas.median()) if len(deltas) else math.nan,
                    "sd_post_minus_pre": float(deltas.std(ddof=1)) if len(deltas) > 1 else math.nan,
                    "positive_delta_patients": int((deltas > 0).sum()) if len(deltas) else 0,
                    "positive_delta_fraction": float((deltas > 0).sum() / len(deltas)) if len(deltas) else math.nan,
                    "paired_t_stat": t_stat,
                    "paired_t_pvalue": t_p,
                    "wilcoxon_stat": w_stat,
                    "wilcoxon_pvalue": w_p,
                    "sign_test_pvalue": two_sided_sign_test(deltas, 0.0),
                    "interpretation_boundary": "Orthogonal TCR/cluster treatment-direction validation only; no OSCC response-depth labels.",
                }
            )
    out = pd.DataFrame(rows)
    for pcol in ["paired_t_pvalue", "wilcoxon_pvalue", "sign_test_pvalue"]:
        out[pcol.replace("_pvalue", "_fdr")] = bh_adjust(out[pcol].tolist())
    return out.sort_values(["disease_group", "metric"]).reset_index(drop=True)


def pair_turnover(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (disease, patient), sub in df.groupby(["disease", "patient"], sort=True):
        if not {"pre", "post"}.issubset(set(sub["treatment"])):
            continue
        pre = sub.loc[sub["treatment"] == "pre", CLONOTYPE_COL].value_counts()
        post = sub.loc[sub["treatment"] == "post", CLONOTYPE_COL].value_counts()
        if pre.sum() < MIN_TCR_CELLS_PER_TIMEPOINT or post.sum() < MIN_TCR_CELLS_PER_TIMEPOINT:
            continue
        all_clones = pre.index.union(post.index)
        pre_vec = pre.reindex(all_clones, fill_value=0).to_numpy(dtype=float)
        post_vec = post.reindex(all_clones, fill_value=0).to_numpy(dtype=float)
        pre_norm = pre_vec / pre_vec.sum()
        post_norm = post_vec / post_vec.sum()
        denom = float(np.sum(pre_norm**2) + np.sum(post_norm**2))
        morisita_horn = float(2 * np.sum(pre_norm * post_norm) / denom) if denom > 0 else math.nan
        shared = pre.index.intersection(post.index)
        top10_post = post.iloc[: min(10, len(post))]
        rows.append(
            {
                "disease": disease,
                "patient": patient,
                "pre_n_tcr_cells": int(pre.sum()),
                "post_n_tcr_cells": int(post.sum()),
                "pre_n_clonotypes": int(len(pre)),
                "post_n_clonotypes": int(len(post)),
                "shared_clonotypes": int(len(shared)),
                "jaccard_clonotype_overlap": float(len(shared) / len(all_clones)) if len(all_clones) else math.nan,
                "morisita_horn_overlap": morisita_horn,
                "post_new_cell_fraction": float(post[~post.index.isin(pre.index)].sum() / post.sum()),
                "pre_lost_cell_fraction": float(pre[~pre.index.isin(post.index)].sum() / pre.sum()),
                "replacement_asymmetry_post_new_minus_pre_lost": float(
                    post[~post.index.isin(pre.index)].sum() / post.sum()
                    - pre[~pre.index.isin(post.index)].sum() / pre.sum()
                ),
                "top10_post_new_fraction": float((~top10_post.index.isin(pre.index)).mean()) if len(top10_post) else math.nan,
            }
        )
    return pd.DataFrame(rows).sort_values(["disease", "patient"]).reset_index(drop=True)


def turnover_summary(turnover: pd.DataFrame) -> pd.DataFrame:
    rows = []
    metrics_and_nulls = {
        "post_new_cell_fraction": 0.5,
        "pre_lost_cell_fraction": 0.5,
        "replacement_asymmetry_post_new_minus_pre_lost": 0.0,
        "morisita_horn_overlap": math.nan,
        "top10_post_new_fraction": 0.5,
    }
    groups = [("ALL", turnover)]
    groups.extend((disease, sub) for disease, sub in turnover.groupby("disease", sort=True))
    for group_name, sub in groups:
        for metric, null in metrics_and_nulls.items():
            values = pd.to_numeric(sub[metric], errors="coerce").dropna()
            if math.isnan(null):
                t_stat = t_p = w_stat = w_p = sign_p = math.nan
            else:
                t_stat, t_p = ttest_1samp(values, null)
                w_stat, w_p = wilcoxon_1samp(values, null)
                sign_p = two_sided_sign_test(values, null)
            rows.append(
                {
                    "disease_group": group_name,
                    "metric": metric,
                    "n_pairs": int(len(values)),
                    "mean": float(values.mean()) if len(values) else math.nan,
                    "median": float(values.median()) if len(values) else math.nan,
                    "min": float(values.min()) if len(values) else math.nan,
                    "max": float(values.max()) if len(values) else math.nan,
                    "null_value_for_test": null,
                    "one_sample_t_stat": t_stat,
                    "one_sample_t_pvalue": t_p,
                    "wilcoxon_stat": w_stat,
                    "wilcoxon_pvalue": w_p,
                    "sign_test_pvalue": sign_p,
                }
            )
    out = pd.DataFrame(rows)
    for pcol in ["one_sample_t_pvalue", "wilcoxon_pvalue", "sign_test_pvalue"]:
        out[pcol.replace("_pvalue", "_fdr")] = bh_adjust(out[pcol].tolist())
    return out.sort_values(["disease_group", "metric"]).reset_index(drop=True)


def write_metadata_summary(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for disease, sub in df.groupby("disease", sort=True):
        rows.append(
            {
                "disease": disease,
                "n_tcr_matched_cells": int(len(sub)),
                "n_patients": int(sub["patient"].nunique()),
                "n_paired_patients": int(
                    sum({"pre", "post"}.issubset(set(x["treatment"])) for _, x in sub.groupby("patient"))
                ),
                "pre_cells": int((sub["treatment"] == "pre").sum()),
                "post_cells": int((sub["treatment"] == "post").sum()),
                "clusters": ";".join(sorted(sub["cluster"].dropna().unique())),
            }
        )
    summary = pd.DataFrame(rows)
    summary.to_csv(RESULTS_DIR / "GSE123813_TCR_metadata_summary.csv", index=False)
    return summary


def make_plot(metrics: pd.DataFrame, delta_stats: pd.DataFrame, turnover: pd.DataFrame) -> list[Path]:
    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 10,
            "axes.linewidth": 0.6,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    fig = plt.figure(figsize=(11.5, 8.4), constrained_layout=True)
    grid = fig.add_gridspec(2, 2)
    ax_state = fig.add_subplot(grid[0, 0])
    ax_expand = fig.add_subplot(grid[0, 1])
    ax_turnover = fig.add_subplot(grid[1, 0])
    ax_overlap = fig.add_subplot(grid[1, 1])

    colors = {"BCC": "#0072B2", "SCC": "#009E73"}

    def paired_panel(ax, metric: str, title: str, ylabel: str):
        offset = {"BCC": -0.08, "SCC": 0.08}
        for disease, sub in metrics.groupby("disease", sort=True):
            wide = sub.pivot(index="patient", columns="treatment", values=metric)
            nwide = sub.pivot(index="patient", columns="treatment", values="n_tcr_cells")
            if not {"pre", "post"}.issubset(wide.columns):
                continue
            keep = (nwide["pre"] >= MIN_TCR_CELLS_PER_TIMEPOINT) & (nwide["post"] >= MIN_TCR_CELLS_PER_TIMEPOINT)
            wide = wide.loc[keep]
            for _, row in wide.iterrows():
                ax.plot(
                    [0 + offset[disease], 1 + offset[disease]],
                    [row["pre"], row["post"]],
                    color=colors[disease],
                    alpha=0.35,
                    lw=1.0,
                )
                ax.scatter(
                    [0 + offset[disease], 1 + offset[disease]],
                    [row["pre"], row["post"]],
                    color=colors[disease],
                    s=24,
                    alpha=0.8,
                )
            if len(wide):
                ax.plot(
                    [0 + offset[disease], 1 + offset[disease]],
                    [wide["pre"].mean(), wide["post"].mean()],
                    color=colors[disease],
                    lw=2.6,
                    marker="o",
                    label=disease,
                )
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["pre", "post"])
        ax.set_title(title, fontsize=9.5, fontweight="semibold")
        ax.set_ylabel(ylabel, fontsize=9)
        ax.tick_params(labelsize=8)
        ax.spines[["top", "right"]].set_visible(False)

    paired_panel(ax_state, "cd8_ex_fraction", "TCR-matched CD8 exhausted-state fraction", "fraction of TCR cells")
    paired_panel(ax_expand, "expanded_cell_fraction", "Expanded-clonotype cell fraction", "fraction of TCR cells")
    ax_state.legend(frameon=False, fontsize=8)

    for disease, sub in turnover.groupby("disease", sort=True):
        ax_turnover.scatter(
            sub["pre_lost_cell_fraction"],
            sub["post_new_cell_fraction"],
            s=50,
            color=colors[disease],
            label=disease,
            edgecolor="black",
            linewidth=0.4,
        )
        for row in sub.itertuples(index=False):
            ax_turnover.text(
                row.pre_lost_cell_fraction + 0.006,
                row.post_new_cell_fraction + 0.006,
                row.patient,
                fontsize=6.5,
                color="#333333",
            )
    ax_turnover.plot([0, 1], [0, 1], color="#777777", lw=1, ls="--")
    ax_turnover.set_xlim(0.35, 1.02)
    ax_turnover.set_ylim(0.45, 1.02)
    ax_turnover.set_xlabel("pre-lost clonotype cell fraction", fontsize=9)
    ax_turnover.set_ylabel("post-new clonotype cell fraction", fontsize=9)
    ax_turnover.set_title("Pre/post clonotype turnover", fontsize=9.5, fontweight="semibold")
    ax_turnover.tick_params(labelsize=8)
    ax_turnover.spines[["top", "right"]].set_visible(False)

    y = np.arange(len(turnover))
    sorted_turnover = turnover.sort_values(["disease", "morisita_horn_overlap"]).reset_index(drop=True)
    ax_overlap.scatter(
        sorted_turnover["morisita_horn_overlap"],
        y,
        c=[colors[d] for d in sorted_turnover["disease"]],
        s=48,
        edgecolor="black",
        linewidth=0.4,
    )
    ax_overlap.set_yticks(y)
    ax_overlap.set_yticklabels(
        [f"{r.disease} {r.patient}" for r in sorted_turnover.itertuples(index=False)],
        fontsize=7,
    )
    ax_overlap.set_xlabel("Morisita-Horn repertoire overlap", fontsize=9)
    ax_overlap.set_title("Lower overlap indicates stronger clonal turnover", fontsize=9.5, fontweight="semibold")
    ax_overlap.tick_params(axis="x", labelsize=8)
    ax_overlap.spines[["top", "right"]].set_visible(False)

    for label, ax in zip(["a", "b", "c", "d"], [ax_state, ax_expand, ax_turnover, ax_overlap], strict=True):
        ax.text(
            -0.10,
            1.08,
            label,
            transform=ax.transAxes,
            fontsize=12,
            fontweight="bold",
            va="bottom",
            ha="left",
        )
    paths = []
    for ext in ["png", "pdf", "svg"]:
        out = FIGURE_DIR / f"{FIG_STEM}.{ext}"
        fig.savefig(out, dpi=300)
        paths.append(out)
    plt.close(fig)
    return paths


def update_extended_data_manifest() -> None:
    manifest_path = FIGURE_DIR / "EXTENDED_DATA_FIGURE_MANIFEST.csv"
    rows: list[dict[str, str]] = []
    fieldnames = ["generated_at", "figure", "stem", "primary_message"]
    if manifest_path.exists():
        with manifest_path.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle)
            if reader.fieldnames:
                fieldnames = reader.fieldnames
            rows = [row for row in reader if row.get("figure") != "ExtendedData6"]
    rows.append(
        {
            "generated_at": datetime.now().strftime("%Y-%m-%d %H:%M:%S HKT"),
            "figure": "ExtendedData6",
            "stem": FIG_STEM,
            "primary_message": "External GSE123813 BCC/SCC TCR data support treatment-associated T-cell state and repertoire remodeling, not OSCC response prediction.",
        }
    )
    rows = sorted(rows, key=lambda r: r.get("figure", ""))
    with manifest_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_source_csvs(
    metadata_summary: pd.DataFrame,
    metrics: pd.DataFrame,
    delta_stats: pd.DataFrame,
    turnover: pd.DataFrame,
    turnover_stats: pd.DataFrame,
) -> list[Path]:
    specs = [
        ("ExtendedData6_TCR_metadata_summary.csv", metadata_summary),
        ("ExtendedData6_TCR_sample_metrics.csv", metrics),
        ("ExtendedData6_TCR_paired_delta_stats.csv", delta_stats),
        ("ExtendedData6_TCR_pair_turnover.csv", turnover),
        ("ExtendedData6_TCR_turnover_summary.csv", turnover_stats),
    ]
    paths = []
    for name, df in specs:
        path = SOURCE_DATA_DIR / name
        df.to_csv(path, index=False)
        paths.append(path)
    return paths


def write_interpretation(
    metadata_summary: pd.DataFrame,
    delta_stats: pd.DataFrame,
    turnover: pd.DataFrame,
    turnover_stats: pd.DataFrame,
    figure_paths: list[Path],
) -> None:
    def stat_row(df: pd.DataFrame, **kwargs) -> pd.Series:
        sub = df.copy()
        for key, value in kwargs.items():
            sub = sub[sub[key] == value]
        return sub.iloc[0]

    all_cd8ex = stat_row(delta_stats, disease_group="ALL", metric="cd8_ex_fraction")
    bcc_postnew = stat_row(turnover_stats, disease_group="BCC", metric="post_new_cell_fraction")
    all_postnew = stat_row(turnover_stats, disease_group="ALL", metric="post_new_cell_fraction")
    all_mh = stat_row(turnover_stats, disease_group="ALL", metric="morisita_horn_overlap")

    lines = [
        "# Orthogonal external TCR validation: GSE123813 BCC/SCC anti-PD-1",
        "",
        f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S HKT')}",
        "",
        "## Dataset boundary",
        "",
        "- Source: GSE123813 BCC and SCC TCR plus T-cell metadata files already downloaded under `00_raw_data/GSE123813_validation/`.",
        "- Use: orthogonal external validation that paired anti-PD-1 data contain T-cell state/repertoire remodeling.",
        "- Boundary: this is not OSCC response-depth validation and should not be described as validating a fixed OSCC biomarker signature.",
        "",
        "## Local data summary",
        "",
        metadata_summary.to_string(index=False),
        "",
        "## Key statistics",
        "",
        "- Combined BCC/SCC TCR-matched CD8 exhausted-state fraction increased in {pos}/{n} paired patients; mean post-pre delta={mean:.4f}, Wilcoxon P={wp:.3g}, Wilcoxon FDR={wq:.3g}.".format(
            pos=int(all_cd8ex["positive_delta_patients"]),
            n=int(all_cd8ex["n_paired_patients"]),
            mean=float(all_cd8ex["mean_post_minus_pre"]),
            wp=float(all_cd8ex["wilcoxon_pvalue"]) if not pd.isna(all_cd8ex["wilcoxon_pvalue"]) else math.nan,
            wq=float(all_cd8ex["wilcoxon_fdr"]) if not pd.isna(all_cd8ex["wilcoxon_fdr"]) else math.nan,
        ),
        "- BCC post-treatment novel-clonotype cell fraction had median={median:.3f}, mean={mean:.3f} across {n} paired patients.".format(
            median=float(bcc_postnew["median"]),
            mean=float(bcc_postnew["mean"]),
            n=int(bcc_postnew["n_pairs"]),
        ),
        "- Across all BCC/SCC pairs, post-new clonotype cell fraction median={median:.3f}, mean={mean:.3f}; Morisita-Horn repertoire overlap median={mh_median:.3f}, mean={mh_mean:.3f}.".format(
            median=float(all_postnew["median"]),
            mean=float(all_postnew["mean"]),
            mh_median=float(all_mh["median"]),
            mh_mean=float(all_mh["mean"]),
        ),
        "",
        "## Interpretation",
        "",
        "This TCR analysis supports a broad anti-PD-1 treatment-remodeling context: TCR-matched cells show a directionally positive exhausted CD8 state shift and substantial pre/post repertoire turnover. However, clonal expansion metrics were mixed, SCC sample size was small, and both post-new and pre-lost clonotype fractions were high. The result should be used to calibrate the manuscript toward pathway/state/repertoire remodeling, not toward a fixed predictive signature.",
        "",
        "## Outputs",
        "",
    ]
    for path in figure_paths:
        lines.append(f"- Figure: `{path.relative_to(ROOT).as_posix()}`")
    lines.extend(
        [
            f"- Sample metrics: `{(RESULTS_DIR / 'GSE123813_TCR_sample_metrics.csv').relative_to(ROOT).as_posix()}`",
            f"- Paired delta stats: `{(RESULTS_DIR / 'GSE123813_TCR_paired_delta_stats.csv').relative_to(ROOT).as_posix()}`",
            f"- Pair turnover: `{(RESULTS_DIR / 'GSE123813_TCR_pair_turnover.csv').relative_to(ROOT).as_posix()}`",
            "",
        ]
    )
    (RESULTS_DIR / "GSE123813_TCR_validation_interpretation.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    ensure_dirs()
    print(f"Workspace: {ROOT}")
    print(f"Raw dir:   {RAW_DIR}")
    merged = pd.concat([load_dataset(spec) for spec in DATASETS], ignore_index=True)
    metadata_summary = write_metadata_summary(merged)

    metrics = sample_metrics(merged)
    delta_stats = paired_delta_stats(metrics)
    turnover = pair_turnover(merged)
    turnover_stats = turnover_summary(turnover)

    metrics.to_csv(RESULTS_DIR / "GSE123813_TCR_sample_metrics.csv", index=False)
    delta_stats.to_csv(RESULTS_DIR / "GSE123813_TCR_paired_delta_stats.csv", index=False)
    turnover.to_csv(RESULTS_DIR / "GSE123813_TCR_pair_turnover.csv", index=False)
    turnover_stats.to_csv(RESULTS_DIR / "GSE123813_TCR_turnover_summary.csv", index=False)

    write_source_csvs(metadata_summary, metrics, delta_stats, turnover, turnover_stats)
    figure_paths = make_plot(metrics, delta_stats, turnover)
    update_extended_data_manifest()
    write_interpretation(metadata_summary, delta_stats, turnover, turnover_stats, figure_paths)

    print("External TCR validation complete.")
    print(delta_stats[["disease_group", "metric", "n_paired_patients", "mean_post_minus_pre", "wilcoxon_pvalue", "wilcoxon_fdr"]].to_string(index=False))
    print(turnover_stats[["disease_group", "metric", "n_pairs", "mean", "median"]].to_string(index=False))


if __name__ == "__main__":
    main()
