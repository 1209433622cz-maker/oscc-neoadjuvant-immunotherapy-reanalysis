from __future__ import annotations

import gzip
import math
import os
import re
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests


WORKSPACE = Path(os.environ.get("GSE200996_WORKSPACE", Path.cwd())).resolve()
RAW_PATH = WORKSPACE / "00_raw_data" / "external_validation" / "GSE179730" / "GSE179730_RNAseq-combinedCPM.txt.gz"
MODULE_PATH = WORKSPACE / "03_rebuild" / "results" / "external_validation" / "GSE123813_gene_set_manifest.csv"
OUT_DIR = WORKSPACE / "03_rebuild" / "validation" / "GSE179730_bulk_treatment_direction"
FIG_DIR = WORKSPACE / "03_rebuild" / "figures" / "external_validation"
SRC_DIR = FIG_DIR / "source_data"
for path in [OUT_DIR, FIG_DIR, SRC_DIR]:
    path.mkdir(parents=True, exist_ok=True)


def bh(values: list[float]) -> list[float]:
    valid = [i for i, value in enumerate(values) if not math.isnan(value)]
    out = [math.nan] * len(values)
    if valid:
        adjusted = multipletests([values[i] for i in valid], method="fdr_bh")[1]
        for idx, value in zip(valid, adjusted):
            out[idx] = float(value)
    return out


def load_expression() -> pd.DataFrame:
    with gzip.open(RAW_PATH, "rt", encoding="utf-8", errors="replace") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        sample_ids = header[1:]
        records = []
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) != len(header):
                continue
            gene = fields[0].strip()
            if not gene:
                continue
            values = pd.to_numeric(pd.Series(fields[1:]), errors="coerce").to_numpy(dtype=float)
            records.append((gene, values))
    expr = pd.DataFrame(np.vstack([record[1] for record in records]), index=[record[0] for record in records], columns=sample_ids)
    return expr.groupby(expr.index).mean()


def sample_annotation(sample_ids: list[str]) -> pd.DataFrame:
    rows = []
    for sample_id in sample_ids:
        match = re.match(r"^(HN\d+)\.(Pre|Post|Recur)$", sample_id, flags=re.IGNORECASE)
        patient_id = match.group(1) if match else ""
        timepoint = match.group(2).lower() if match else ""
        rows.append(
            {
                "sample_id": sample_id,
                "patient_id": patient_id,
                "timepoint": timepoint,
                "dataset_role": "external_bulk_treatment_direction_response_pending",
                "response_status": "pending_patient_level_table",
            }
        )
    return pd.DataFrame(rows)


def load_modules(expr_index: pd.Index) -> pd.DataFrame:
    modules = pd.read_csv(MODULE_PATH)
    present = set(expr_index)
    rows = []
    for _, row in modules.iterrows():
        genes = [gene.strip() for gene in str(row["genes_defined"]).split(";") if gene.strip()]
        genes_present = [gene for gene in genes if gene in present]
        rows.append(
            {
                "signature": row["signature"],
                "target_lineage": row["target_lineage"],
                "source": row["source"],
                "n_genes_defined": len(genes),
                "n_genes_present_in_GSE179730": len(genes_present),
                "coverage_fraction": len(genes_present) / len(genes) if genes else math.nan,
                "genes_present": ";".join(genes_present),
                "genes_missing": ";".join([gene for gene in genes if gene not in present]),
            }
        )
    return pd.DataFrame(rows)


def score_modules(expr: pd.DataFrame, modules: pd.DataFrame, paired_sample_ids: list[str]) -> pd.DataFrame:
    log_expr = np.log2(expr[paired_sample_ids] + 1.0)
    z = log_expr.sub(log_expr.mean(axis=1), axis=0).div(log_expr.std(axis=1, ddof=0).replace(0, np.nan), axis=0)
    score_rows = []
    for _, module in modules.iterrows():
        genes = [gene for gene in str(module["genes_present"]).split(";") if gene]
        values = z.loc[genes].mean(axis=0, skipna=True) if genes else pd.Series(np.nan, index=paired_sample_ids)
        for sample_id, score in values.items():
            score_rows.append(
                {
                    "sample_id": sample_id,
                    "signature": module["signature"],
                    "target_lineage": module["target_lineage"],
                    "module_score": float(score) if pd.notna(score) else math.nan,
                    "n_genes_present": int(module["n_genes_present_in_GSE179730"]),
                }
            )
    return pd.DataFrame(score_rows)


def paired_delta_stats(scored: pd.DataFrame, samples: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    merged = scored.merge(samples, on="sample_id", how="left")
    pre = merged[merged["timepoint"] == "pre"]
    post = merged[merged["timepoint"] == "post"]
    pairs = pre.merge(
        post,
        on=["patient_id", "signature", "target_lineage"],
        suffixes=("_pre", "_post"),
        how="inner",
    )
    pairs["delta_post_minus_pre"] = pairs["module_score_post"] - pairs["module_score_pre"]

    rows = []
    for signature, sub in pairs.groupby("signature"):
        values = sub["delta_post_minus_pre"].replace([np.inf, -np.inf], np.nan).dropna().astype(float)
        info = sub.iloc[0]
        if len(values) >= 2:
            t_p = stats.ttest_1samp(values, 0.0, nan_policy="omit").pvalue
            try:
                w_p = stats.wilcoxon(values).pvalue
            except ValueError:
                w_p = math.nan
            sem = values.std(ddof=1) / math.sqrt(len(values))
            ci = stats.t.ppf(0.975, df=len(values) - 1) * sem if len(values) > 1 else math.nan
        else:
            t_p = w_p = ci = math.nan
        rows.append(
            {
                "signature": signature,
                "target_lineage": info["target_lineage"],
                "n_pairs": int(len(values)),
                "mean_delta_post_minus_pre": float(values.mean()) if len(values) else math.nan,
                "median_delta_post_minus_pre": float(values.median()) if len(values) else math.nan,
                "sd_delta": float(values.std(ddof=1)) if len(values) > 1 else math.nan,
                "ci95_low": float(values.mean() - ci) if len(values) > 1 else math.nan,
                "ci95_high": float(values.mean() + ci) if len(values) > 1 else math.nan,
                "t_test_p": float(t_p) if not math.isnan(t_p) else math.nan,
                "wilcoxon_p": float(w_p) if not math.isnan(w_p) else math.nan,
            }
        )
    stats_df = pd.DataFrame(rows)
    stats_df["t_test_fdr"] = bh(stats_df["t_test_p"].tolist())
    stats_df["wilcoxon_fdr"] = bh(stats_df["wilcoxon_p"].tolist())
    stats_df = stats_df.sort_values(["t_test_fdr", "t_test_p", "signature"], na_position="last")
    return pairs, stats_df


def make_figure(stats_df: pd.DataFrame, pairs: pd.DataFrame) -> None:
    focus = stats_df[
        stats_df["signature"].str.contains("INTERFERON|MTORC1|TNFA|P53|union_core", regex=True)
    ].copy()
    focus = focus.sort_values("mean_delta_post_minus_pre")
    clean_name = focus["signature"].str.replace(r"^[TM]_", "", regex=True).str.replace("_", " ")
    lineage_label = focus["target_lineage"].replace({"T_cell": "T cell", "Myeloid": "Myeloid"})
    focus["label"] = lineage_label + " | " + clean_name

    pairs_focus = pairs[pairs["signature"].isin(focus["signature"])].copy()
    pairs_focus["label"] = pairs_focus["signature"].map(dict(zip(focus["signature"], focus["label"])))
    pairs_focus["label"] = pd.Categorical(pairs_focus["label"], categories=focus["label"], ordered=True)
    pairs_focus.to_csv(SRC_DIR / "GSE179730_treatment_delta_distribution_source.csv", index=False)
    focus.to_csv(SRC_DIR / "GSE179730_treatment_delta_summary_source.csv", index=False)

    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 7,
            "axes.linewidth": 0.6,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    y = np.arange(len(focus))
    colors = np.where(focus["target_lineage"].eq("T_cell"), "#3b6fb6", "#b4584b")
    ax.axvline(0, color="#888888", lw=0.7, zorder=1)
    ax.errorbar(
        focus["mean_delta_post_minus_pre"],
        y,
        xerr=[
            focus["mean_delta_post_minus_pre"] - focus["ci95_low"],
            focus["ci95_high"] - focus["mean_delta_post_minus_pre"],
        ],
        fmt="none",
        ecolor="#333333",
        elinewidth=0.7,
        capsize=2,
        zorder=2,
    )
    ax.scatter(focus["mean_delta_post_minus_pre"], y, c=colors, s=28, edgecolor="white", linewidth=0.4, zorder=3)
    for i, (_, row) in enumerate(focus.iterrows()):
        vals = pairs_focus[pairs_focus["signature"] == row["signature"]]["delta_post_minus_pre"].dropna().to_numpy()
        if len(vals):
            jitter = np.linspace(-0.18, 0.18, len(vals))
            ax.scatter(vals, np.full(len(vals), i) + jitter, c="#111111", alpha=0.28, s=8, linewidth=0, zorder=2)
    ax.set_yticks(y)
    ax.set_yticklabels(focus["label"])
    ax.set_xlabel("Module score change after nivolumab (post - pre)")
    ax.set_title("GSE179730 paired bulk RNA-seq: locked module treatment direction", loc="left", fontsize=9, pad=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="x", color="#dddddd", lw=0.4)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "GSE179730_bulk_treatment_direction_locked_modules.png", dpi=450)
    fig.savefig(FIG_DIR / "GSE179730_bulk_treatment_direction_locked_modules.pdf")
    plt.close(fig)


def write_report(expr: pd.DataFrame, samples: pd.DataFrame, modules: pd.DataFrame, stats_df: pd.DataFrame) -> None:
    paired_n = samples[samples["timepoint"].isin(["pre", "post"])].groupby("patient_id")["timepoint"].nunique().eq(2).sum()
    top = stats_df.head(10)
    report = [
        "# GSE179730 locked-module treatment-direction validation",
        "",
        f"Generated: {datetime.now().isoformat(timespec='seconds')}",
        "",
        "## Inputs",
        "",
        f"- Expression: `{RAW_PATH}`",
        f"- Locked module manifest: `{MODULE_PATH}`",
        f"- Expression matrix: {expr.shape[0]} genes x {expr.shape[1]} samples.",
        f"- Paired pre/post patients used for treatment-direction testing: {paired_n}.",
        "",
        "## Response label gate",
        "",
        "The local GEO/SRA/expression files do not contain patient-level response labels. The associated open full text confirms the response definition and overall responder/non-responder counts, but the patient-level Table S1/S2 mapping is still required before using this cohort as response validation.",
        "",
        "Therefore, this analysis is intentionally restricted to paired treatment-direction validation of the locked modules.",
        "",
        "## Top paired post-minus-pre module changes",
        "",
        "| signature | lineage | n pairs | mean delta | t-test P | t-test FDR |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for _, row in top.iterrows():
        report.append(
            f"| {row['signature']} | {row['target_lineage']} | {int(row['n_pairs'])} | "
            f"{row['mean_delta_post_minus_pre']:.4f} | {row['t_test_p']:.3g} | {row['t_test_fdr']:.3g} |"
        )
    report.extend(
        [
            "",
            "## Interpretation",
            "",
            "Use this result as a response-pending, cross-cohort treatment perturbation check. It can support whether the locked module families are biologically mobile after neoadjuvant PD-1 blockade, but it must not be written as responder/non-responder validation until the patient-level response table is recovered.",
            "",
            "## Outputs",
            "",
            f"- Sample annotation: `{OUT_DIR / 'GSE179730_SAMPLE_ANNOTATION.csv'}`",
            f"- Module coverage: `{OUT_DIR / 'GSE179730_LOCKED_MODULE_GENE_COVERAGE.csv'}`",
            f"- Sample scores: `{OUT_DIR / 'GSE179730_LOCKED_MODULE_SAMPLE_SCORES.csv'}`",
            f"- Paired deltas: `{OUT_DIR / 'GSE179730_LOCKED_MODULE_PAIRED_DELTAS.csv'}`",
            f"- Paired statistics: `{OUT_DIR / 'GSE179730_PAIRED_TREATMENT_DELTA_STATS.csv'}`",
            f"- Figure PNG/PDF: `{FIG_DIR / 'GSE179730_bulk_treatment_direction_locked_modules.png'}`",
        ]
    )
    (OUT_DIR / "GSE179730_TREATMENT_DIRECTION_REPORT.md").write_text("\n".join(report) + "\n", encoding="utf-8")


def main() -> None:
    expr = load_expression()
    samples = sample_annotation(expr.columns.tolist())
    paired_sample_ids = samples[samples["timepoint"].isin(["pre", "post"])]["sample_id"].tolist()
    modules = load_modules(expr.index)
    score_values = score_modules(expr, modules, paired_sample_ids)
    scored = score_values.merge(samples, on="sample_id", how="left")
    pairs, stats_df = paired_delta_stats(score_values, samples)

    samples.to_csv(OUT_DIR / "GSE179730_SAMPLE_ANNOTATION.csv", index=False)
    modules.to_csv(OUT_DIR / "GSE179730_LOCKED_MODULE_GENE_COVERAGE.csv", index=False)
    scored.to_csv(OUT_DIR / "GSE179730_LOCKED_MODULE_SAMPLE_SCORES.csv", index=False)
    pairs.to_csv(OUT_DIR / "GSE179730_LOCKED_MODULE_PAIRED_DELTAS.csv", index=False)
    stats_df.to_csv(OUT_DIR / "GSE179730_PAIRED_TREATMENT_DELTA_STATS.csv", index=False)
    make_figure(stats_df, pairs)
    write_report(expr, samples, modules, stats_df)


if __name__ == "__main__":
    main()
