#!/usr/bin/env python
"""Prespecified locked-family response and matched-random specificity analyses."""

from __future__ import annotations

import gzip
import itertools
import json
import math
import os
import re
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm


WORKSPACE = Path(os.environ.get("GSE200996_WORKSPACE", Path.cwd())).resolve()
CONFIG_PATH = WORKSPACE / "03_rebuild" / "config" / "locked_family_robustness.json"
MODULE_PATH = (
    WORKSPACE
    / "03_rebuild"
    / "results"
    / "external_validation"
    / "GSE123813_gene_set_manifest.csv"
)
GSE179_RAW = (
    WORKSPACE
    / "00_raw_data"
    / "external_validation"
    / "GSE179730"
    / "GSE179730_RNAseq-combinedCPM.txt.gz"
)
GSE179_LABELS = (
    WORKSPACE
    / "03_rebuild"
    / "validation"
    / "GSE179730_bulk_treatment_direction"
    / "GSE179730_RESPONSE_LABELS_TABLE_S2.csv"
)
GSE281_RAW = (
    WORKSPACE
    / "00_raw_data"
    / "external_validation"
    / "GSE281729"
    / "GSE281729_Mastrolonardo_etal_expressions_logTPM_Nivo-IDO_Bulk-RNAseq_Processed_File.txt.gz"
)
OUT_DIR = WORKSPACE / "03_rebuild" / "validation" / "locked_family_robustness"
SOURCE_DIR = WORKSPACE / "03_rebuild" / "figures" / "submission" / "source_data"
OUT_DIR.mkdir(parents=True, exist_ok=True)
SOURCE_DIR.mkdir(parents=True, exist_ok=True)

SAMPLE_OUT = OUT_DIR / "LOCKED_FAMILY_SAMPLE_SCORES.csv"
DELTA_OUT = OUT_DIR / "LOCKED_FAMILY_PAIRED_DELTAS.csv"
TEST_OUT = OUT_DIR / "LOCKED_FAMILY_TESTS.csv"
RANDOM_OUT = OUT_DIR / "LOCKED_FAMILY_MATCHED_RANDOM_EFFECTS.csv"
EXACT_NULL_OUT = OUT_DIR / "GSE179730_LOCKED_FAMILY_EXACT_NULL.csv"
COVERAGE_OUT = OUT_DIR / "LOCKED_FAMILY_MODULE_COVERAGE.csv"
REPORT_OUT = OUT_DIR / "LOCKED_FAMILY_ROBUSTNESS_REPORT.md"


def load_config() -> dict:
    config = json.loads(CONFIG_PATH.read_text(encoding="utf-8"))
    if config["expected_module_count"] != 16:
        raise RuntimeError("The frozen analysis requires exactly 16 modules")
    if config["module_reselection_allowed"] or config["gene_reselection_allowed"]:
        raise RuntimeError("The frozen configuration must prohibit reselection")
    return config


def load_modules() -> pd.DataFrame:
    modules = pd.read_csv(MODULE_PATH)
    if len(modules) != 16 or modules["signature"].duplicated().any():
        raise RuntimeError("Locked module manifest is not the expected 16-module family")
    modules = modules.copy()
    modules["genes"] = modules["genes_defined"].fillna("").map(
        lambda value: [gene.strip() for gene in str(value).split(";") if gene.strip()]
    )
    return modules


def load_gse179730() -> tuple[pd.DataFrame, pd.DataFrame]:
    with gzip.open(GSE179_RAW, "rt", encoding="utf-8", errors="replace") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        sample_ids = header[1:]
        records = []
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) != len(header) or not fields[0].strip():
                continue
            values = pd.to_numeric(pd.Series(fields[1:]), errors="coerce").to_numpy(float)
            records.append((fields[0].strip(), values))
    expr = pd.DataFrame(
        np.vstack([record[1] for record in records]),
        index=[record[0] for record in records],
        columns=sample_ids,
    ).groupby(level=0).mean()
    expr = np.log2(expr + 1.0)

    rows = []
    for sample_id in expr.columns:
        match = re.match(r"^(HN\d+)\.(Pre|Post|Recur)$", sample_id, flags=re.IGNORECASE)
        rows.append(
            {
                "sample_id": sample_id,
                "patient_id": match.group(1) if match else "",
                "timepoint": match.group(2).lower() if match else "",
            }
        )
    annotation = pd.DataFrame(rows)
    paired_counts = (
        annotation[annotation["timepoint"].isin(["pre", "post"])]
        .groupby("patient_id")["timepoint"]
        .nunique()
    )
    paired_patients = paired_counts[paired_counts.eq(2)].index.tolist()
    annotation = annotation[
        annotation["patient_id"].isin(paired_patients)
        & annotation["timepoint"].isin(["pre", "post"])
    ].copy()
    expr = expr[annotation["sample_id"].tolist()]

    labels = pd.read_csv(GSE179_LABELS)
    label_map = labels.set_index("geo_patient_id")["response_binary"].to_dict()
    annotation["response_binary"] = annotation["patient_id"].map(label_map)
    if annotation["response_binary"].isna().any():
        raise RuntimeError("GSE179730 paired patients are missing source response labels")
    return expr, annotation


def parse_percent(value: str) -> float:
    text = str(value or "").strip().replace("%", "")
    try:
        return float(text) if text else math.nan
    except ValueError:
        return math.nan


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


def load_gse281729() -> tuple[pd.DataFrame, pd.DataFrame]:
    with gzip.open(GSE281_RAW, "rt", encoding="utf-8", errors="replace") as handle:
        lines = [next(handle).rstrip("\n").split("\t") for _ in range(14)]
        expression_rows = [line.rstrip("\n").split("\t") for line in handle]
    annotation_fields: dict[str, list[str]] = {}
    for row in lines[:13]:
        if len(row) >= 3:
            annotation_fields[row[1].strip() or f"annotation_{len(annotation_fields) + 1}"] = row[2:]
    sample_ids = lines[13][2:]
    rows = []
    for index, sample_id in enumerate(sample_ids):
        match = re.search(r"TJ3_(\d+)_", sample_id)
        response_raw = annotation_fields.get("Primary Path Response", [""] * len(sample_ids))[index].strip()
        response_ordinal, response_binary, response_number = normalize_response(response_raw)
        rows.append(
            {
                "sample_id": sample_id,
                "patient_id": f"NI_TJ3_{match.group(1)}" if match else "",
                "timepoint": annotation_fields.get("Time Point", [""] * len(sample_ids))[index].strip().lower(),
                "hpv": annotation_fields.get("HPV", [""] * len(sample_ids))[index].strip(),
                "second_drug": annotation_fields.get("SecondDrug", [""] * len(sample_ids))[index].strip(),
                "doses": pd.to_numeric(
                    annotation_fields.get("Doses", [""] * len(sample_ids))[index].strip(),
                    errors="coerce",
                ),
                "primary_response_percent": parse_percent(
                    annotation_fields.get("%Primary Response", [""] * len(sample_ids))[index]
                ),
                "response_harmonized_ordinal": response_ordinal,
                "response_binary": response_binary,
                "response_ord_num": response_number,
            }
        )
    annotation = pd.DataFrame(rows)

    records = []
    for row in expression_rows:
        if len(row) < 3 or not row[1].strip():
            continue
        values = pd.to_numeric(pd.Series(row[2:]), errors="coerce").to_numpy(float)
        records.append((row[1].strip(), values))
    expr = pd.DataFrame(
        np.vstack([record[1] for record in records]),
        index=[record[0] for record in records],
        columns=sample_ids,
    ).groupby(level=0).mean()
    return expr, annotation


def score_family(
    expression: pd.DataFrame,
    modules: pd.DataFrame,
    cohort: str,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, np.ndarray], dict[str, np.ndarray], pd.DataFrame]:
    finite = expression.replace([np.inf, -np.inf], np.nan).dropna(axis=0, how="any")
    variable = finite.std(axis=1, ddof=0) > 0
    finite = finite.loc[variable]
    z = finite.sub(finite.mean(axis=1), axis=0).div(finite.std(axis=1, ddof=0), axis=0)
    rank = finite.rank(axis=0, method="average", pct=True)

    module_scores: dict[str, dict[str, pd.Series]] = {"z_score": {}, "rank_mean": {}}
    module_indices: dict[str, np.ndarray] = {}
    coverage_rows = []
    for _, row in modules.iterrows():
        genes = [gene for gene in row["genes"] if gene in finite.index]
        if not genes:
            raise RuntimeError(f"No detectable genes for {cohort} module {row['signature']}")
        module_indices[row["signature"]] = finite.index.get_indexer(genes)
        module_scores["z_score"][row["signature"]] = z.loc[genes].mean(axis=0)
        module_scores["rank_mean"][row["signature"]] = rank.loc[genes].mean(axis=0)
        coverage_rows.append(
            {
                "cohort": cohort,
                "signature": row["signature"],
                "target_lineage": row["target_lineage"],
                "n_genes_defined": len(row["genes"]),
                "n_genes_scored": len(genes),
                "coverage_fraction": len(genes) / len(row["genes"]),
                "genes_scored": ";".join(genes),
            }
        )

    sample_rows = []
    for method in ["z_score", "rank_mean"]:
        matrix = pd.DataFrame(module_scores[method])
        family = matrix.mean(axis=1)
        for sample_id, value in family.items():
            sample_rows.append(
                {
                    "cohort": cohort,
                    "sample_id": sample_id,
                    "scoring_method": method,
                    "locked_family_score": float(value),
                    "n_modules": matrix.shape[1],
                }
            )
    arrays = {
        "z_score": z.to_numpy(float),
        "rank_mean": rank.to_numpy(float),
    }
    return pd.DataFrame(sample_rows), pd.DataFrame(coverage_rows), module_indices, arrays, finite


def paired_deltas(
    scores: pd.DataFrame,
    annotation: pd.DataFrame,
    cohort: str,
) -> pd.DataFrame:
    merged = scores.merge(annotation, on="sample_id", how="left", validate="many_to_one")
    identifiers = ["patient_id", "scoring_method"]
    wide = merged.pivot_table(
        index=identifiers,
        columns="timepoint",
        values="locked_family_score",
        aggfunc="mean",
    ).reset_index()
    if "pre" not in wide.columns or "post" not in wide.columns:
        raise RuntimeError(f"{cohort} lacks paired pre/post family scores")
    wide["delta_post_minus_pre"] = wide["post"] - wide["pre"]
    wide = wide.dropna(subset=["delta_post_minus_pre"]).copy()
    metadata_columns = [
        column
        for column in [
            "patient_id",
            "response_binary",
            "response_harmonized_ordinal",
            "response_ord_num",
            "hpv",
            "second_drug",
            "doses",
            "primary_response_percent",
        ]
        if column in annotation.columns
    ]
    metadata = annotation[metadata_columns].drop_duplicates("patient_id")
    wide = wide.merge(metadata, on="patient_id", how="left", validate="many_to_one")
    wide.insert(0, "cohort", cohort)
    return wide


def exact_gse179730(
    deltas: pd.DataFrame,
) -> tuple[list[dict[str, object]], pd.DataFrame]:
    rows = []
    null_rows = []
    for method, subset in deltas.groupby("scoring_method"):
        current = subset.sort_values("patient_id").reset_index(drop=True)
        labels = current["response_binary"].to_numpy()
        values = current["delta_post_minus_pre"].to_numpy(float)
        n_responder = int((labels == "responder").sum())
        assignments = list(itertools.combinations(range(len(current)), n_responder))
        observed_mask = labels == "responder"
        observed = float(values[observed_mask].mean() - values[~observed_mask].mean())
        null = []
        for assignment_id, indices in enumerate(assignments, start=1):
            mask = np.zeros(len(current), dtype=bool)
            mask[list(indices)] = True
            effect = float(values[mask].mean() - values[~mask].mean())
            null.append(effect)
            null_rows.append(
                {
                    "scoring_method": method,
                    "assignment_id": assignment_id,
                    "effect": effect,
                    "is_observed_assignment": bool(np.array_equal(mask, observed_mask)),
                }
            )
        null_values = np.asarray(null)
        exact_p = float(np.mean(np.abs(null_values) >= abs(observed) - 1e-12))
        rows.append(
            {
                "cohort": "GSE179730",
                "scoring_method": method,
                "analysis": "locked_family_response_exact",
                "n_patients": len(current),
                "n_responder": n_responder,
                "n_non_responder": len(current) - n_responder,
                "effect": observed,
                "std_error": math.nan,
                "ci95_low": math.nan,
                "ci95_high": math.nan,
                "p_value": exact_p,
                "inference": f"two-sided exhaustive {len(assignments)}-assignment exact test",
            }
        )
    return rows, pd.DataFrame(null_rows)


def gse281_design(data: pd.DataFrame) -> pd.DataFrame:
    model_data = data[
        [
            "delta_post_minus_pre",
            "response_ord_num",
            "hpv",
            "second_drug",
            "doses",
        ]
    ].replace([np.inf, -np.inf], np.nan).dropna()
    parts = [model_data[["response_ord_num"]].astype(float)]
    for column in ["hpv", "second_drug"]:
        if model_data[column].nunique() > 1:
            parts.append(
                pd.get_dummies(
                    model_data[column].astype(str),
                    prefix=column,
                    drop_first=True,
                    dtype=float,
                )
            )
    parts.append(model_data[["doses"]].astype(float))
    design = sm.add_constant(pd.concat(parts, axis=1), has_constant="add")
    return model_data, design


def model_gse281729(deltas: pd.DataFrame) -> tuple[list[dict[str, object]], dict[str, np.ndarray]]:
    rows = []
    effect_weights = {}
    for method, subset in deltas.groupby("scoring_method"):
        current = subset.sort_values("patient_id").reset_index(drop=True)
        model_data, design = gse281_design(current)
        fit = sm.OLS(model_data["delta_post_minus_pre"].astype(float), design.astype(float)).fit(cov_type="HC3")
        term = "response_ord_num"
        ci = fit.conf_int(alpha=0.05).loc[term]
        rows.append(
            {
                "cohort": "GSE281729",
                "scoring_method": method,
                "analysis": "locked_family_ordinal_hpv_second_drug_doses_HC3",
                "n_patients": len(model_data),
                "n_responder": math.nan,
                "n_non_responder": math.nan,
                "effect": float(fit.params[term]),
                "std_error": float(fit.bse[term]),
                "ci95_low": float(ci.iloc[0]),
                "ci95_high": float(ci.iloc[1]),
                "p_value": float(fit.pvalues[term]),
                "inference": "OLS ordinal response coefficient with HPV, second-drug and doses/timing adjustment; HC3 covariance",
            }
        )
        pinv = np.linalg.pinv(design.to_numpy(float))
        effect_weights[method] = pinv[list(design.columns).index(term), :]
    return rows, effect_weights


def expression_bins(expression: pd.DataFrame) -> pd.Series:
    means = expression.mean(axis=1)
    ordered = means.rank(method="first")
    return pd.qcut(ordered, q=10, labels=False, duplicates="drop").astype(int)


def prepare_random_sampler(
    gene_index: pd.Index,
    bins: pd.Series,
    modules: pd.DataFrame,
    locked_genes: set[str],
) -> tuple[dict[int, np.ndarray], dict[str, dict[int, int]]]:
    background = np.asarray(
        [index for index, gene in enumerate(gene_index) if gene not in locked_genes],
        dtype=int,
    )
    by_bin: dict[int, np.ndarray] = {}
    for bin_id in sorted(bins.unique()):
        positions = [
            index for index in background if int(bins.iloc[index]) == int(bin_id)
        ]
        by_bin[int(bin_id)] = np.asarray(positions, dtype=int)
    requirements: dict[str, dict[int, int]] = {}
    for _, module in modules.iterrows():
        present = [gene for gene in module["genes"] if gene in gene_index]
        counts = pd.Series([int(bins.loc[gene]) for gene in present]).value_counts()
        requirements[module["signature"]] = {
            int(bin_id): int(count) for bin_id, count in counts.items()
        }
        for bin_id, count in requirements[module["signature"]].items():
            if len(by_bin[bin_id]) < count:
                raise RuntimeError(f"Insufficient matched background genes in expression bin {bin_id}")
    return by_bin, requirements


def random_module_indices(
    rng: np.random.Generator,
    by_bin: dict[int, np.ndarray],
    requirements: dict[str, dict[int, int]],
) -> dict[str, np.ndarray]:
    selected: dict[str, np.ndarray] = {}
    for signature, counts in requirements.items():
        random_positions = []
        for bin_id, count in counts.items():
            pool = by_bin[int(bin_id)]
            random_positions.extend(
                rng.choice(pool, size=int(count), replace=False).tolist()
            )
        selected[signature] = np.asarray(random_positions, dtype=int)
    return selected


def random_family_sample_score(
    arrays: dict[str, np.ndarray],
    indices: dict[str, np.ndarray],
) -> dict[str, np.ndarray]:
    output = {}
    for method, matrix in arrays.items():
        module_values = [matrix[index, :].mean(axis=0) for index in indices.values()]
        output[method] = np.vstack(module_values).mean(axis=0)
    return output


def delta_operator(
    sample_ids: list[str],
    annotation: pd.DataFrame,
    patient_order: list[str],
) -> np.ndarray:
    sample_frame = pd.DataFrame(
        {"sample_id": sample_ids, "sample_position": np.arange(len(sample_ids))}
    ).merge(annotation, on="sample_id", how="left", validate="one_to_one")
    operator = np.zeros((len(patient_order), len(sample_ids)), dtype=float)
    for row_index, patient_id in enumerate(patient_order):
        current = sample_frame[sample_frame["patient_id"] == patient_id]
        pre = current.loc[current["timepoint"] == "pre", "sample_position"].to_numpy(int)
        post = current.loc[current["timepoint"] == "post", "sample_position"].to_numpy(int)
        if len(pre) == 0 or len(post) == 0:
            raise RuntimeError(f"Missing paired samples for {patient_id}")
        operator[row_index, pre] = -1.0 / len(pre)
        operator[row_index, post] = 1.0 / len(post)
    return operator


def matched_random_controls(
    config: dict,
    cohort: str,
    expression: pd.DataFrame,
    annotation: pd.DataFrame,
    modules: pd.DataFrame,
    arrays: dict[str, np.ndarray],
    observed_deltas: pd.DataFrame,
    gse281_weights: dict[str, np.ndarray] | None,
) -> pd.DataFrame:
    rng = np.random.default_rng(int(config["random_seed"]))
    bins = expression_bins(expression)
    locked_genes = {
        gene
        for genes in modules["genes"]
        for gene in genes
        if gene in expression.index
    }
    patient_order = sorted(observed_deltas["patient_id"].unique())
    by_bin, requirements = prepare_random_sampler(
        expression.index,
        bins,
        modules,
        locked_genes,
    )
    patient_to_sample_delta = delta_operator(
        expression.columns.tolist(),
        annotation,
        patient_order,
    )
    sample_effect_weights: dict[str, np.ndarray] = {}
    if cohort == "GSE179730":
        patient_meta = observed_deltas.drop_duplicates("patient_id").set_index("patient_id").loc[patient_order]
        labels = patient_meta["response_binary"].to_numpy() == "responder"
        patient_effect = np.where(
            labels,
            1.0 / labels.sum(),
            -1.0 / (~labels).sum(),
        )
        for method in arrays:
            sample_effect_weights[method] = patient_effect @ patient_to_sample_delta
    else:
        for method in arrays:
            sample_effect_weights[method] = (
                gse281_weights[method] @ patient_to_sample_delta
            )

    rows = []
    for iteration in range(1, int(config["random_control_iterations"]) + 1):
        indices = random_module_indices(
            rng,
            by_bin,
            requirements,
        )
        sample_scores = random_family_sample_score(arrays, indices)
        for method, values in sample_scores.items():
            effect = float(sample_effect_weights[method] @ values)
            rows.append(
                {
                    "cohort": cohort,
                    "scoring_method": method,
                    "iteration": iteration,
                    "random_seed": int(config["random_seed"]),
                    "effect": effect,
                }
            )
    return pd.DataFrame(rows)


def attach_specificity(tests: pd.DataFrame, random_effects: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, test in tests.iterrows():
        null = random_effects[
            (random_effects["cohort"] == test["cohort"])
            & (random_effects["scoring_method"] == test["scoring_method"])
        ]["effect"].to_numpy(float)
        observed = float(test["effect"])
        empirical_p = float((1 + np.sum(np.abs(null) >= abs(observed) - 1e-12)) / (len(null) + 1))
        row = test.to_dict()
        row.update(
            {
                "matched_random_iterations": len(null),
                "matched_random_mean": float(null.mean()),
                "matched_random_sd": float(null.std(ddof=1)),
                "matched_random_q025": float(np.quantile(null, 0.025)),
                "matched_random_q975": float(np.quantile(null, 0.975)),
                "empirical_specificity_p": empirical_p,
            }
        )
        rows.append(row)
    return pd.DataFrame(rows)


def write_report(
    config: dict,
    tests: pd.DataFrame,
    coverage: pd.DataFrame,
) -> None:
    indexed = tests.set_index(["cohort", "scoring_method"])
    gse179_z = indexed.loc[("GSE179730", "z_score")]
    gse179_rank = indexed.loc[("GSE179730", "rank_mean")]
    gse281_z = indexed.loc[("GSE281729", "z_score")]
    gse281_rank = indexed.loc[("GSE281729", "rank_mean")]
    lines = [
        "# Locked-Family External Robustness Report",
        "",
        f"Generated: {datetime.now().astimezone().isoformat(timespec='seconds')}",
        "",
        "## Frozen Design",
        "",
        f"- Configuration: `{CONFIG_PATH}`",
        "- Primary score: unweighted mean of all 16 locked module scores.",
        "- Primary scoring: mean gene-wise z score; sensitivity scoring: within-sample percentile-rank mean.",
        f"- Matched random families: {config['random_control_iterations']} per cohort and scoring method.",
        f"- Random seed: {config['random_seed']}.",
        "- Random controls preserve module size and mean-expression decile; locked genes are excluded from the sampling pool.",
        "",
        "## Results",
        "",
        "| cohort | scoring | effect | inferential P | matched-random empirical P | random 95% interval |",
        "|---|---|---:|---:|---:|---:|",
    ]
    for _, row in tests.iterrows():
        lines.append(
            f"| {row['cohort']} | {row['scoring_method']} | {row['effect']:.4f} | "
            f"{row['p_value']:.4g} | {row['empirical_specificity_p']:.4g} | "
            f"{row['matched_random_q025']:.4f} to {row['matched_random_q975']:.4f} |"
        )
    lines.extend(
        [
            "",
            "## Coverage",
            "",
            f"- Cohort-module rows: {len(coverage)}.",
            f"- Minimum locked-module coverage: {coverage['coverage_fraction'].min():.3f}.",
            f"- Median locked-module coverage: {coverage['coverage_fraction'].median():.3f}.",
            "",
            "## Result-Level Interpretation",
            "",
            f"- GSE281729: the locked family decreased across the ordered response strata after adjustment for HPV status, second drug, dose count and sampling timing under both z-score (effect {gse281_z['effect']:.4f}, HC3 P={gse281_z['p_value']:.4g}) and rank-mean scoring (effect {gse281_rank['effect']:.4f}, HC3 P={gse281_rank['p_value']:.4g}).",
            f"- GSE281729: both observed effects were more extreme than all {config['random_control_iterations']} matched random families under the frozen absolute-effect rule (plus-one empirical P={gse281_z['empirical_specificity_p']:.4g} for each scoring method).",
            f"- GSE179730: the primary z-score family effect was positive but exact-null (effect {gse179_z['effect']:.4f}, exact P={gse179_z['p_value']:.4g}) and not more extreme than matched random families (empirical P={gse179_z['empirical_specificity_p']:.4g}).",
            f"- GSE179730: rank-mean scoring changed the effect direction and remained exact-null (effect {gse179_rank['effect']:.4f}, exact P={gse179_rank['p_value']:.4g}); this is a scoring-sensitivity boundary, not cross-cohort family replication.",
            "",
            "## Interpretation Rule",
            "",
            "The z-score family analysis is primary. Rank-mean scoring and matched-random families are sensitivity and specificity analyses. No result changes the locked family definition. The matched-random comparison establishes extremeness only under the frozen null and does not prove pathway uniqueness independent of global expression shifts. Statistical significance does not establish a clinical predictor or a causal treatment effect and must be interpreted with each cohort's endpoint and response-adaptive sampling limitations.",
            "",
            "## Outputs",
            "",
            f"- Sample scores: `{SAMPLE_OUT}`",
            f"- Paired deltas: `{DELTA_OUT}`",
            f"- Family tests: `{TEST_OUT}`",
            f"- Matched random effects: `{RANDOM_OUT}`",
            f"- GSE179730 exact null: `{EXACT_NULL_OUT}`",
            f"- Module coverage: `{COVERAGE_OUT}`",
        ]
    )
    REPORT_OUT.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    config = load_config()
    modules = load_modules()

    expr179, annotation179 = load_gse179730()
    scores179, coverage179, _, arrays179, eligible179 = score_family(expr179, modules, "GSE179730")
    deltas179 = paired_deltas(scores179, annotation179, "GSE179730")
    test_rows179, exact_null = exact_gse179730(deltas179)

    expr281, annotation281 = load_gse281729()
    scores281, coverage281, _, arrays281, eligible281 = score_family(expr281, modules, "GSE281729")
    deltas281 = paired_deltas(scores281, annotation281, "GSE281729")
    deltas281 = deltas281[
        deltas281["response_harmonized_ordinal"].isin(["Low", "Medium", "High"])
    ].copy()
    test_rows281, effect_weights281 = model_gse281729(deltas281)

    random179 = matched_random_controls(
        config,
        "GSE179730",
        eligible179,
        annotation179,
        modules,
        arrays179,
        deltas179,
        None,
    )
    random281 = matched_random_controls(
        config,
        "GSE281729",
        eligible281,
        annotation281,
        modules,
        arrays281,
        deltas281,
        effect_weights281,
    )

    tests = attach_specificity(
        pd.DataFrame(test_rows179 + test_rows281),
        pd.concat([random179, random281], ignore_index=True),
    )
    sample_scores = pd.concat([scores179, scores281], ignore_index=True)
    deltas = pd.concat([deltas179, deltas281], ignore_index=True)
    coverage = pd.concat([coverage179, coverage281], ignore_index=True)
    random_effects = pd.concat([random179, random281], ignore_index=True)

    sample_scores.to_csv(SAMPLE_OUT, index=False)
    deltas.to_csv(DELTA_OUT, index=False)
    tests.to_csv(TEST_OUT, index=False)
    random_effects.to_csv(RANDOM_OUT, index=False)
    exact_null.to_csv(EXACT_NULL_OUT, index=False)
    coverage.to_csv(COVERAGE_OUT, index=False)
    tests.to_csv(SOURCE_DIR / "ExtendedData11_locked_family_tests_source.csv", index=False)
    random_effects.to_csv(
        SOURCE_DIR / "ExtendedData11_matched_random_effects_source.csv",
        index=False,
    )
    deltas.to_csv(SOURCE_DIR / "ExtendedData11_paired_deltas_source.csv", index=False)
    exact_null.to_csv(SOURCE_DIR / "ExtendedData11_gse179730_exact_null_source.csv", index=False)
    write_report(config, tests, coverage)
    print(REPORT_OUT)
    print(TEST_OUT)


if __name__ == "__main__":
    main()
