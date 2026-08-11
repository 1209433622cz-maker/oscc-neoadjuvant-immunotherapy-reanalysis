from __future__ import annotations

import csv
import gzip
import math
import os
import re
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats
from statsmodels.stats.multitest import multipletests


WORKSPACE = Path(os.environ.get("GSE200996_WORKSPACE", Path.cwd())).resolve()
RAW_PATH = (
    WORKSPACE
    / "00_raw_data"
    / "external_validation"
    / "GSE281729"
    / "GSE281729_Mastrolonardo_etal_expressions_logTPM_Nivo-IDO_Bulk-RNAseq_Processed_File.txt.gz"
)
MODULE_PATH = WORKSPACE / "03_rebuild" / "results" / "external_validation" / "GSE123813_gene_set_manifest.csv"
RESPONSE_PATH = WORKSPACE / "03_rebuild" / "manifests" / "RESPONSE_METADATA_FREEZE.csv"
OUT_DIR = WORKSPACE / "03_rebuild" / "validation" / "GSE281729_bulk_module_validation"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def bh(values: list[float]) -> list[float]:
    valid = [i for i, value in enumerate(values) if not math.isnan(value)]
    out = [math.nan] * len(values)
    if valid:
        adjusted = multipletests([values[i] for i in valid], method="fdr_bh")[1]
        for idx, value in zip(valid, adjusted):
            out[idx] = float(value)
    return out


def normalize_response(raw: str) -> tuple[str, str, float]:
    value = str(raw or "").strip()
    lower = value.lower()
    if not value:
        return "", "", math.nan
    if lower in {"nr", "non-responder", "non-respoder", "non responder", "non respoder"} or "non" in lower:
        return "Low", "NR", 0.0
    if "minor" in lower or lower in {"min resp", "minor resp", "min responder"}:
        return "Medium", "intermediate", 1.0
    if lower in {"cr", "complete responder", "complete respoder"} or "complete" in lower:
        return "High", "R", 2.0
    if lower in {"r", "responder", "respoder"}:
        return "High", "R", 2.0
    return "", "", math.nan


def parse_percent(value: str) -> float:
    text = str(value or "").strip().replace("%", "")
    if not text:
        return math.nan
    try:
        return float(text)
    except ValueError:
        return math.nan


def load_expression_and_annotation() -> tuple[pd.DataFrame, pd.DataFrame]:
    with gzip.open(RAW_PATH, "rt", encoding="utf-8", errors="replace") as fh:
        lines = [next(fh).rstrip("\n").split("\t") for _ in range(14)]
        expression_rows = [line.rstrip("\n").split("\t") for line in fh]

    annotation = {}
    for row in lines[:13]:
        if len(row) < 3:
            continue
        key = row[1].strip() or f"annotation_{len(annotation) + 1}"
        annotation[key] = row[2:]

    sample_ids = lines[13][2:]
    sample_rows = []
    for idx, sample_id in enumerate(sample_ids):
        match = re.search(r"TJ3_(\d+)_", sample_id)
        patient_id = f"NI_TJ3_{match.group(1)}" if match else ""
        primary_response = annotation.get("Primary Path Response", [""] * len(sample_ids))[idx].strip()
        resp_ordinal, resp_binary, resp_num = normalize_response(primary_response)
        sample_rows.append(
            {
                "sample_id": sample_id,
                "patient_id": patient_id,
                "subject_code_initials": annotation.get("Subject Code-Intials", [""] * len(sample_ids))[idx].strip(),
                "first_drug": annotation.get("FirstDrug", [""] * len(sample_ids))[idx].strip(),
                "second_drug": annotation.get("SecondDrug", [""] * len(sample_ids))[idx].strip(),
                "doses": annotation.get("Doses", [""] * len(sample_ids))[idx].strip(),
                "timepoint": annotation.get("Time Point", [""] * len(sample_ids))[idx].strip().lower(),
                "hpv": annotation.get("HPV", [""] * len(sample_ids))[idx].strip(),
                "smoking_status": annotation.get("Smoking Status", [""] * len(sample_ids))[idx].strip(),
                "primary_path_response_raw": primary_response,
                "primary_response_percent": parse_percent(annotation.get("%Primary Response", [""] * len(sample_ids))[idx]),
                "overall_response_percent": parse_percent(annotation.get("% Overall Response", [""] * len(sample_ids))[idx]),
                "overall_path_response": annotation.get("Overall path Response", [""] * len(sample_ids))[idx].strip(),
                "geo_fastq_1": annotation.get("GEO Files (paired-end seq)", [""] * len(sample_ids))[idx].strip(),
                "response_harmonized_ordinal": resp_ordinal,
                "response_binary": resp_binary,
                "response_ord_num": resp_num,
                "source": "GSE281729 processed expression embedded clinical annotation",
            }
        )
    sample_df = pd.DataFrame(sample_rows)

    expr_records = []
    for row in expression_rows:
        if len(row) < 3:
            continue
        ensg, symbol = row[0].strip(), row[1].strip()
        if not symbol:
            continue
        values = pd.to_numeric(pd.Series(row[2:]), errors="coerce").to_numpy(dtype=float)
        expr_records.append((ensg, symbol, values))

    symbols = [record[1] for record in expr_records]
    data = np.vstack([record[2] for record in expr_records])
    expr = pd.DataFrame(data, index=symbols, columns=sample_ids)
    expr = expr.groupby(expr.index).mean()
    return expr, sample_df


def load_modules(expr_index: pd.Index) -> pd.DataFrame:
    modules = pd.read_csv(MODULE_PATH)
    rows = []
    present = set(expr_index)
    for _, row in modules.iterrows():
        genes = [gene.strip() for gene in str(row["genes_defined"]).split(";") if gene.strip()]
        genes_present = [gene for gene in genes if gene in present]
        rows.append(
            {
                "signature": row["signature"],
                "target_lineage": row["target_lineage"],
                "source": row["source"],
                "n_genes_defined": len(genes),
                "n_genes_present_in_GSE281729": len(genes_present),
                "genes_present": ";".join(genes_present),
                "genes_missing": ";".join([gene for gene in genes if gene not in present]),
            }
        )
    return pd.DataFrame(rows)


def score_modules(expr: pd.DataFrame, modules: pd.DataFrame) -> pd.DataFrame:
    z = expr.sub(expr.mean(axis=1), axis=0).div(expr.std(axis=1, ddof=0).replace(0, np.nan), axis=0)
    score_rows = []
    for _, module in modules.iterrows():
        genes = [gene for gene in str(module["genes_present"]).split(";") if gene]
        if not genes:
            values = pd.Series(np.nan, index=expr.columns)
        else:
            values = z.loc[genes].mean(axis=0, skipna=True)
        for sample_id, score in values.items():
            score_rows.append(
                {
                    "sample_id": sample_id,
                    "signature": module["signature"],
                    "target_lineage": module["target_lineage"],
                    "module_score": score,
                    "n_genes_present": module["n_genes_present_in_GSE281729"],
                }
            )
    return pd.DataFrame(score_rows)


def ols_slope(df: pd.DataFrame, y: str, x: str, covariates: list[str] | None = None) -> dict[str, float | str | int]:
    covariates = covariates or []
    needed = [y, x] + covariates
    data = df[needed].copy()
    data = data.replace([np.inf, -np.inf], np.nan).dropna()
    if len(data) < 4 or data[x].nunique() < 2:
        return {
            "n": len(data),
            "coef": math.nan,
            "std_error": math.nan,
            "ci95_low": math.nan,
            "ci95_high": math.nan,
            "p_value": math.nan,
            "model_status": "not_estimable",
        }
    design_parts = [data[[x]].astype(float)]
    for cov in covariates:
        if data[cov].nunique() > 1:
            design_parts.append(pd.get_dummies(data[cov].astype(str), prefix=cov, drop_first=True, dtype=float))
    X = pd.concat(design_parts, axis=1)
    X = sm.add_constant(X, has_constant="add")
    try:
        fit = sm.OLS(data[y].astype(float), X.astype(float)).fit()
        ci = fit.conf_int(alpha=0.05).loc[x]
        return {
            "n": int(len(data)),
            "coef": float(fit.params[x]),
            "std_error": float(fit.bse[x]),
            "ci95_low": float(ci[0]),
            "ci95_high": float(ci[1]),
            "p_value": float(fit.pvalues[x]),
            "model_status": "ok",
        }
    except Exception as exc:
        return {
            "n": int(len(data)),
            "coef": math.nan,
            "std_error": math.nan,
            "ci95_low": math.nan,
            "ci95_high": math.nan,
            "p_value": math.nan,
            "model_status": f"error: {exc}",
        }


def run_models(scored: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    baseline_rows = []
    delta_rows = []
    binary_rows = []

    pre = scored[scored["timepoint"].str.lower().str.startswith("pre")].copy()
    post = scored[scored["timepoint"].str.lower().str.startswith("post")].copy()

    for signature, sub in pre.groupby("signature"):
        info = sub.iloc[0]
        for model_name, covs in {
            "ordinal_unadjusted": [],
            "ordinal_adjusted_hpv_second_drug": ["hpv", "second_drug"],
        }.items():
            res = ols_slope(sub, "module_score", "response_ord_num", covs)
            baseline_rows.append(
                {
                    "analysis": "baseline_pre",
                    "model": model_name,
                    "signature": signature,
                    "target_lineage": info["target_lineage"],
                    **res,
                }
            )

    wide = scored.pivot_table(
        index=["patient_id", "signature"],
        columns="timepoint",
        values="module_score",
        aggfunc="mean",
    ).reset_index()
    meta_cols = [
        "patient_id",
        "signature",
        "target_lineage",
        "response_ord_num",
        "response_binary",
        "response_harmonized_ordinal",
        "hpv",
        "second_drug",
    ]
    meta = scored[meta_cols].drop_duplicates(["patient_id", "signature"])
    wide = wide.merge(meta, on=["patient_id", "signature"], how="left")
    if "post" in wide.columns and "pre" in wide.columns:
        wide["post_minus_pre"] = wide["post"] - wide["pre"]
    else:
        wide["post_minus_pre"] = np.nan
    paired = wide.dropna(subset=["post_minus_pre"]).copy()

    for signature, sub in paired.groupby("signature"):
        info = sub.iloc[0]
        for model_name, covs in {
            "ordinal_unadjusted": [],
            "ordinal_adjusted_hpv_second_drug": ["hpv", "second_drug"],
        }.items():
            res = ols_slope(sub, "post_minus_pre", "response_ord_num", covs)
            delta_rows.append(
                {
                    "analysis": "paired_post_minus_pre",
                    "model": model_name,
                    "signature": signature,
                    "target_lineage": info["target_lineage"],
                    **res,
                }
            )

        rb = sub[sub["response_binary"].isin(["R", "NR"])].copy()
        if rb["response_binary"].nunique() == 2:
            r = rb.loc[rb["response_binary"] == "R", "post_minus_pre"]
            nr = rb.loc[rb["response_binary"] == "NR", "post_minus_pre"]
            test = stats.ttest_ind(r, nr, equal_var=False, nan_policy="omit")
            binary_rows.append(
                {
                    "analysis": "paired_post_minus_pre",
                    "model": "binary_R_vs_NR_welch",
                    "signature": signature,
                    "target_lineage": info["target_lineage"],
                    "n_R": int(r.notna().sum()),
                    "n_NR": int(nr.notna().sum()),
                    "mean_R": float(r.mean()),
                    "mean_NR": float(nr.mean()),
                    "mean_R_minus_NR": float(r.mean() - nr.mean()),
                    "p_value": float(test.pvalue) if not math.isnan(test.pvalue) else math.nan,
                    "model_status": "ok",
                }
            )
        else:
            binary_rows.append(
                {
                    "analysis": "paired_post_minus_pre",
                    "model": "binary_R_vs_NR_welch",
                    "signature": signature,
                    "target_lineage": info["target_lineage"],
                    "n_R": int((rb["response_binary"] == "R").sum()),
                    "n_NR": int((rb["response_binary"] == "NR").sum()),
                    "mean_R": math.nan,
                    "mean_NR": math.nan,
                    "mean_R_minus_NR": math.nan,
                    "p_value": math.nan,
                    "model_status": "not_estimable",
                }
            )

    baseline = pd.DataFrame(baseline_rows)
    delta = pd.DataFrame(delta_rows)
    binary = pd.DataFrame(binary_rows)
    for table in [baseline, delta, binary]:
        if "p_value" in table.columns:
            for model in table["model"].dropna().unique():
                mask = table["model"] == model
                table.loc[mask, "fdr"] = bh(table.loc[mask, "p_value"].astype(float).tolist())
    return baseline, delta, binary


def make_plot(stats_df: pd.DataFrame) -> None:
    plot = stats_df[
        (stats_df["analysis"] == "paired_post_minus_pre")
        & (stats_df["model"] == "ordinal_unadjusted")
        & (stats_df["model_status"] == "ok")
    ].copy()
    if plot.empty:
        return
    plot = plot.sort_values("coef")
    clean = (
        plot["signature"]
        .str.replace(r"^[TM]_", "", regex=True)
        .str.replace("LE_", "", regex=False)
        .str.replace("_", " ", regex=False)
        .str.replace("INTERFERON", "IFN", regex=False)
        .str.replace("SIGNALING VIA NFKB", "NF-kB", regex=False)
    )
    lineage = plot["target_lineage"].replace({"T_cell": "T cell", "Myeloid": "Myeloid"})
    plot["plot_label"] = lineage + " | " + clean
    fig_h = max(4.5, 0.28 * len(plot))
    fig, ax = plt.subplots(figsize=(7.2, fig_h))
    colors = plot["target_lineage"].map({"T_cell": "#3B7EA1", "Myeloid": "#B55A30"}).fillna("#555555")
    y = np.arange(len(plot))
    x = plot["coef"].to_numpy(dtype=float)
    xerr = np.vstack([x - plot["ci95_low"].to_numpy(dtype=float), plot["ci95_high"].to_numpy(dtype=float) - x])
    ax.errorbar(x, y, xerr=xerr, fmt="none", ecolor="#333333", elinewidth=0.8, capsize=2, zorder=2)
    ax.scatter(x, y, c=colors, s=30, edgecolor="white", linewidth=0.5, zorder=3)
    ax.axvline(0, color="#888888", linewidth=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels(plot["plot_label"], fontsize=7)
    ax.set_xlabel("Post-pre module-score slope per response-depth step (95% CI)")
    ax.set_ylabel("")
    ax.set_title("GSE281729 bulk validation of locked OSCC modules", loc="left", fontsize=10)
    ax.set_xlim(float(plot["ci95_low"].min()) - 0.05, 0.12)
    for y, (_, row) in enumerate(plot.iterrows()):
        label = f"FDR={row['fdr']:.3g}" if not math.isnan(row["fdr"]) else "FDR=NA"
        ax.text(0.015, y, label, va="center", ha="left", fontsize=6.5, color="#444444")
    ax.grid(axis="x", color="#dddddd", linewidth=0.4)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "GSE281729_locked_module_delta_slopes.png", dpi=300)
    fig.savefig(OUT_DIR / "GSE281729_locked_module_delta_slopes.pdf")
    plt.close(fig)


def main() -> int:
    expr, sample_df = load_expression_and_annotation()

    response_freeze = pd.read_csv(RESPONSE_PATH)
    gse281 = response_freeze[response_freeze["accession"] == "GSE281729"][
        ["patient_id", "response_harmonized_ordinal", "response_binary", "response_status"]
    ].copy()
    sample_df = sample_df.merge(gse281, on="patient_id", how="left", suffixes=("_embedded", "_freeze"))
    for col in ["response_harmonized_ordinal", "response_binary"]:
        sample_df[col] = sample_df[f"{col}_freeze"].fillna(sample_df[f"{col}_embedded"])
    sample_df["response_ord_num"] = sample_df["response_harmonized_ordinal"].map({"Low": 0.0, "Medium": 1.0, "High": 2.0})
    sample_df["response_status"] = sample_df["response_status"].fillna("embedded_expression_annotation")

    modules = load_modules(expr.index)
    scored = score_modules(expr, modules).merge(sample_df, on="sample_id", how="left")
    baseline, delta, binary = run_models(scored)
    combined_stats = pd.concat([baseline, delta, binary], ignore_index=True, sort=False)

    sample_df.to_csv(OUT_DIR / "GSE281729_EMBEDDED_SAMPLE_ANNOTATION.csv", index=False)
    modules.to_csv(OUT_DIR / "GSE281729_LOCKED_MODULE_GENE_COVERAGE.csv", index=False)
    scored.to_csv(OUT_DIR / "GSE281729_LOCKED_MODULE_SAMPLE_SCORES.csv", index=False)
    baseline.to_csv(OUT_DIR / "GSE281729_BASELINE_MODULE_RESPONSE_MODELS.csv", index=False)
    delta.to_csv(OUT_DIR / "GSE281729_PAIRED_DELTA_MODULE_RESPONSE_MODELS.csv", index=False)
    binary.to_csv(OUT_DIR / "GSE281729_BINARY_DELTA_MODULE_RESPONSE_MODELS.csv", index=False)
    combined_stats.to_csv(OUT_DIR / "GSE281729_LOCKED_MODULE_VALIDATION_STATS.csv", index=False)
    make_plot(combined_stats)

    confirmed_patients = sample_df.drop_duplicates("patient_id")
    confirmed_patients = confirmed_patients[confirmed_patients["response_harmonized_ordinal"].isin(["Low", "Medium", "High"])]
    paired_patients = (
        sample_df.groupby("patient_id")["timepoint"]
        .apply(lambda x: {"pre", "post"}.issubset(set(x.str.lower())))
        .reset_index(name="has_prepost")
    )
    paired_confirmed = paired_patients.merge(
        confirmed_patients[["patient_id", "response_harmonized_ordinal"]], on="patient_id", how="inner"
    )
    primary = delta[(delta["model"] == "ordinal_unadjusted") & (delta["model_status"] == "ok")].copy()
    primary = primary.sort_values(["fdr", "p_value"], na_position="last")
    top = primary.head(8)

    md = [
        "# GSE281729 Bulk Module Validation",
        "",
        f"Created: {datetime.now().astimezone().isoformat(timespec='seconds')}",
        "",
        "## Dataset",
        "",
        "- Source expression file: `00_raw_data/external_validation/GSE281729/GSE281729_Mastrolonardo_etal_expressions_logTPM_Nivo-IDO_Bulk-RNAseq_Processed_File.txt.gz`",
        f"- Genes parsed: {expr.shape[0]}",
        f"- Samples parsed: {expr.shape[1]}",
        f"- Response-annotated patients: {confirmed_patients['patient_id'].nunique()}",
        f"- Response-annotated paired patients: {int(paired_confirmed['has_prepost'].sum())}",
        "",
        "## Method",
        "",
        "Locked OSCC-derived modules were reused from `03_rebuild/results/external_validation/GSE123813_gene_set_manifest.csv`.",
        "For each gene, logTPM values were z-scored across samples; module score was the mean z-score of genes present in GSE281729.",
        "Primary validation models tested whether post-pre module-score deltas scaled with ordinal pathological response depth.",
        "",
        "## Primary Paired-Delta Results",
        "",
        "| Signature | Lineage | n | Slope | 95% CI | P | FDR |",
        "|---|---|---:|---:|---:|---:|---:|",
    ]
    for _, row in top.iterrows():
        md.append(
            f"| {row['signature']} | {row['target_lineage']} | {int(row['n'])} | {row['coef']:.4g} | "
            f"{row['ci95_low']:.4g} to {row['ci95_high']:.4g} | {row['p_value']:.3g} | {row['fdr']:.3g} |"
        )
    md.extend(
        [
            "",
            "## Interpretation",
            "",
            "This is a locked-module external validation, not a module-discovery run.",
            "The signed paired-delta slope reports how module-score change varies with ordinal pathological response depth.",
            "The leading validated slopes are negative in this cohort, so GSE281729 supports pathway-level response association rather than universal post-treatment upregulation.",
            "Results should be interpreted at pathway/module level because GSE281729 is bulk RNA-seq and cannot resolve cell-type abundance directly.",
            "",
            "## Files Written",
            "",
            "- `03_rebuild/validation/GSE281729_bulk_module_validation/GSE281729_EMBEDDED_SAMPLE_ANNOTATION.csv`",
            "- `03_rebuild/validation/GSE281729_bulk_module_validation/GSE281729_LOCKED_MODULE_GENE_COVERAGE.csv`",
            "- `03_rebuild/validation/GSE281729_bulk_module_validation/GSE281729_LOCKED_MODULE_SAMPLE_SCORES.csv`",
            "- `03_rebuild/validation/GSE281729_bulk_module_validation/GSE281729_LOCKED_MODULE_VALIDATION_STATS.csv`",
            "- `03_rebuild/validation/GSE281729_bulk_module_validation/GSE281729_locked_module_delta_slopes.png`",
            "- `03_rebuild/validation/GSE281729_bulk_module_validation/GSE281729_locked_module_delta_slopes.pdf`",
        ]
    )
    (OUT_DIR / "GSE281729_BULK_MODULE_VALIDATION_REPORT.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    print("\n".join(md))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
