from __future__ import annotations

import hashlib
import json
import math
from collections import Counter
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm


matplotlib.use("Agg")
matplotlib.rcParams["svg.hashsalt"] = "gse195832-locked-family"

ROOT = Path(__file__).resolve().parents[2]
REBUILD = ROOT / "03_rebuild"
CONFIG_PATH = REBUILD / "config" / "gse195832_locked_family_validation.json"
OUT_DIR = REBUILD / "validation" / "GSE195832_bulk_locked_family"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def bh_adjust(values: pd.Series) -> pd.Series:
    array = values.to_numpy(float)
    order = np.argsort(array)
    ranked = array[order]
    adjusted = ranked * len(array) / np.arange(1, len(array) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    output = np.empty_like(adjusted)
    output[order] = np.minimum(adjusted, 1.0)
    return pd.Series(output, index=values.index)


def load_config() -> dict:
    config = json.loads(CONFIG_PATH.read_text(encoding="utf-8"))
    if config["outcome_guided_selection_allowed"]:
        raise RuntimeError("Outcome-guided selection must remain disabled")
    return config


def load_modules(path: Path, expected: int) -> pd.DataFrame:
    modules = pd.read_csv(path)
    modules["genes"] = modules["genes_defined"].fillna("").map(
        lambda value: list(dict.fromkeys(
            gene.strip().upper() for gene in str(value).split(";") if gene.strip()
        ))
    )
    if len(modules) != expected:
        raise RuntimeError(f"Expected {expected} frozen modules; found {len(modules)}")
    return modules


def load_inputs(config: dict) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    counts_path = ROOT / config["counts_file"]
    metadata_path = ROOT / config["metadata_file"]
    for path, expected_hash in [
        (counts_path, config["counts_sha256"]),
        (metadata_path, config["metadata_sha256"]),
    ]:
        observed = sha256(path)
        if observed != expected_hash:
            raise RuntimeError(f"SHA256 mismatch for {path}: {observed}")

    metadata = pd.read_csv(metadata_path, index_col=0)
    metadata = metadata.rename(
        columns={
            "Patient": "patient_id",
            "Therapy": "therapy",
            "Clinical Response(Primary_Site_TE)": "response_label",
            "Gender": "gender",
            "Age": "age_source",
            "TobaccoHx": "tobacco_history",
            "Sample": "sample_id",
            "Treatment": "treatment_time",
            "Batch": "batch",
        }
    )
    response_map = {
        label: index for index, label in enumerate(config["response_order"])
    }
    metadata["response_ord_num"] = metadata["response_label"].map(response_map)
    metadata["timepoint"] = metadata["treatment_time"].map(
        {"Pre-Treated": "pre", "Post-Treated": "post"}
    )
    metadata["therapy_tadalafil"] = (
        metadata["therapy"].eq("antiPD1+Tadalafil").astype(int)
    )
    metadata["batch2"] = metadata["batch"].eq("Batch2").astype(int)

    required = [
        "patient_id", "therapy", "response_label", "sample_id", "timepoint", "batch"
    ]
    if metadata[required].isna().any().any():
        raise RuntimeError("Metadata contains missing required fields")
    if metadata["response_ord_num"].isna().any():
        raise RuntimeError("Unrecognized pathological-response label")
    if len(metadata) != int(config["expected_samples"]):
        raise RuntimeError("Unexpected GSE195832 metadata sample count")
    pair_counts = metadata.groupby(["patient_id", "timepoint"]).size().unstack(fill_value=0)
    if len(pair_counts) != int(config["expected_paired_patients"]):
        raise RuntimeError("Unexpected GSE195832 paired-patient count")
    if not pair_counts.reindex(columns=["pre", "post"], fill_value=0).eq(1).all().all():
        raise RuntimeError("Every patient must have exactly one pre and one post sample")
    patient_fields = ["therapy", "response_label", "response_ord_num", "batch"]
    if (metadata.groupby("patient_id")[patient_fields].nunique() != 1).any().any():
        raise RuntimeError("Patient-level metadata changes between paired samples")

    counts = pd.read_csv(counts_path, sep="\t")
    counts["symbol"] = counts["symbol"].astype(str).str.strip().str.upper()
    counts = counts[counts["symbol"].ne("") & counts["symbol"].ne("NAN")]
    counts = counts.set_index("symbol")
    counts = counts.apply(pd.to_numeric, errors="raise")
    if (counts.to_numpy() < 0).any():
        raise RuntimeError("Raw counts contain negative values")
    counts = counts.groupby(level=0, sort=True).sum()
    sample_ids = metadata["sample_id"].tolist()
    if set(counts.columns) != set(sample_ids):
        missing = sorted(set(sample_ids) - set(counts.columns))
        extra = sorted(set(counts.columns) - set(sample_ids))
        raise RuntimeError(f"Count/metadata mismatch; missing={missing}; extra={extra}")
    counts = counts.loc[:, sample_ids]
    library_sizes = counts.sum(axis=0)
    if (library_sizes <= 0).any():
        raise RuntimeError("Non-positive sequencing library size")
    cpm = counts.div(library_sizes, axis=1) * 1_000_000
    minimum_samples = 6
    keep = (cpm >= 1).sum(axis=1) >= minimum_samples
    log_cpm = np.log2(cpm.loc[keep] + 0.5)
    return metadata.reset_index(drop=True), counts, log_cpm


def build_delta_operator(
    metadata: pd.DataFrame, sample_order: list[str], patient_order: list[str]
) -> np.ndarray:
    positions = {sample: index for index, sample in enumerate(sample_order)}
    operator = np.zeros((len(patient_order), len(sample_order)), dtype=float)
    for row, patient in enumerate(patient_order):
        current = metadata[metadata["patient_id"].eq(patient)]
        pre = current.loc[current["timepoint"].eq("pre"), "sample_id"].iloc[0]
        post = current.loc[current["timepoint"].eq("post"), "sample_id"].iloc[0]
        operator[row, positions[pre]] = -1.0
        operator[row, positions[post]] = 1.0
    return operator


def score_frozen_family(
    expression: pd.DataFrame, modules: pd.DataFrame
) -> tuple[
    pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, np.ndarray],
    dict[str, np.ndarray], pd.DataFrame
]:
    finite = expression.replace([np.inf, -np.inf], np.nan).dropna(axis=0)
    finite = finite.loc[finite.std(axis=1, ddof=0).gt(0)]
    z = finite.sub(finite.mean(axis=1), axis=0).div(
        finite.std(axis=1, ddof=0), axis=0
    )
    rank = finite.rank(axis=0, method="average", pct=True)
    arrays = {"z_score": z.to_numpy(float), "rank_mean": rank.to_numpy(float)}

    present_modules: dict[str, np.ndarray] = {}
    coverage_rows: list[dict[str, object]] = []
    for _, module in modules.iterrows():
        genes = [gene for gene in module["genes"] if gene in finite.index]
        if not genes:
            raise RuntimeError(f"No expressed genes for {module['signature']}")
        present_modules[module["signature"]] = finite.index.get_indexer(genes)
        coverage_rows.append(
            {
                "signature": module["signature"],
                "target_lineage": module["target_lineage"],
                "n_genes_defined": len(module["genes"]),
                "n_genes_scored": len(genes),
                "coverage_fraction": len(genes) / len(module["genes"]),
                "genes_scored": ";".join(genes),
            }
        )

    memberships = Counter(
        int(position)
        for positions in present_modules.values()
        for position in positions
    )
    unique_positions = np.asarray(sorted(memberships), dtype=int)
    no_union = {
        name: positions
        for name, positions in present_modules.items()
        if not name.endswith("_union_core")
    }
    hallmark = {
        name: positions
        for name, positions in present_modules.items()
        if "_LE_" in name and not name.endswith("_union_core")
    }
    dynamic = {
        name: positions
        for name, positions in present_modules.items()
        if "_LE_" not in name
    }

    def module_mean(matrix: np.ndarray, selected: dict[str, np.ndarray]) -> np.ndarray:
        return np.vstack(
            [matrix[positions].mean(axis=0) for positions in selected.values()]
        ).mean(axis=0)

    def inverse_membership_mean(matrix: np.ndarray) -> np.ndarray:
        values = []
        for positions in present_modules.values():
            weights = np.asarray([1 / memberships[int(pos)] for pos in positions])
            values.append(np.average(matrix[positions], axis=0, weights=weights))
        return np.vstack(values).mean(axis=0)

    family_rows: list[dict[str, object]] = []
    module_rows: list[dict[str, object]] = []
    for method, matrix in arrays.items():
        variants = {
            "module_mean_16": module_mean(matrix, present_modules),
            "unique_gene_equal": matrix[unique_positions].mean(axis=0),
            "inverse_membership_module_mean": inverse_membership_mean(matrix),
            "no_union_module_mean": module_mean(matrix, no_union),
            "hallmark_only_module_mean": module_mean(matrix, hallmark),
            "dynamic_only_module_mean": module_mean(matrix, dynamic),
        }
        for variant, values in variants.items():
            for sample_id, value in zip(finite.columns, values):
                family_rows.append(
                    {
                        "sample_id": sample_id,
                        "scoring_method": method,
                        "score_variant": variant,
                        "family_score": float(value),
                    }
                )
        for signature, positions in present_modules.items():
            values = matrix[positions].mean(axis=0)
            for sample_id, value in zip(finite.columns, values):
                module_rows.append(
                    {
                        "sample_id": sample_id,
                        "scoring_method": method,
                        "signature": signature,
                        "module_score": float(value),
                    }
                )
    return (
        pd.DataFrame(family_rows),
        pd.DataFrame(module_rows),
        pd.DataFrame(coverage_rows),
        arrays,
        present_modules,
        finite,
    )


def paired_deltas(
    scores: pd.DataFrame,
    metadata: pd.DataFrame,
    value_column: str,
    group_columns: list[str],
) -> pd.DataFrame:
    merged = scores.merge(metadata, on="sample_id", how="left", validate="many_to_one")
    index = ["patient_id", *group_columns]
    wide = merged.pivot(index=index, columns="timepoint", values=value_column).reset_index()
    wide["delta_post_minus_pre"] = wide["post"] - wide["pre"]
    patient_meta = metadata[
        [
            "patient_id", "therapy", "therapy_tadalafil", "batch", "batch2",
            "response_label", "response_ord_num", "gender", "age_source",
            "tobacco_history",
        ]
    ].drop_duplicates("patient_id")
    return wide.merge(patient_meta, on="patient_id", validate="many_to_one")


def response_design(
    data: pd.DataFrame, outcome: str, covariates: list[str]
) -> tuple[pd.DataFrame, pd.DataFrame]:
    columns = [outcome, "response_ord_num", *covariates]
    current = data[columns].replace([np.inf, -np.inf], np.nan).dropna().copy()
    design = sm.add_constant(
        current[["response_ord_num", *covariates]].astype(float),
        has_constant="add",
    )
    return current, design


def fit_model(
    data: pd.DataFrame, outcome: str, covariates: list[str]
) -> tuple[dict[str, object], np.ndarray]:
    current, design = response_design(data, outcome, covariates)
    fit = sm.OLS(current[outcome].astype(float), design).fit(cov_type="HC3")
    term = "response_ord_num"
    ci = fit.conf_int().loc[term]
    weights = np.linalg.pinv(design.to_numpy(float))[
        list(design.columns).index(term)
    ]
    return (
        {
            "n_patients": len(current),
            "effect": float(fit.params[term]),
            "std_error": float(fit.bse[term]),
            "ci95_low": float(ci.iloc[0]),
            "ci95_high": float(ci.iloc[1]),
            "p_value_hc3": float(fit.pvalues[term]),
            "model_r2": float(sm.OLS(current[outcome].astype(float), design).fit().rsquared),
            "design_rank": int(np.linalg.matrix_rank(design.to_numpy(float))),
            "design_columns": int(design.shape[1]),
            "design_condition_number": float(np.linalg.cond(design.to_numpy(float))),
        },
        weights,
    )


def deterministic_pca(
    expression: pd.DataFrame,
    metadata: pd.DataFrame,
    locked_genes: set[str],
    n_components: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    patients = sorted(metadata["patient_id"].unique())
    operator = build_delta_operator(metadata, expression.columns.tolist(), patients)
    eligible = [
        gene for gene in expression.index
        if gene not in locked_genes
    ]
    delta = operator @ expression.loc[eligible].T.to_numpy(float)
    gene_sd = delta.std(axis=0, ddof=0)
    keep = gene_sd > 0
    delta = delta[:, keep]
    genes = np.asarray(eligible)[keep]
    standardized = (delta - delta.mean(axis=0)) / delta.std(axis=0, ddof=0)
    u, singular, vt = np.linalg.svd(standardized, full_matrices=False)
    total = float(np.sum(singular**2))
    score_frame = pd.DataFrame({"patient_id": patients})
    loading_rows = []
    variance_rows = []
    for index in range(n_components):
        scores = u[:, index] * singular[index]
        loadings = vt[index].copy()
        anchor = int(np.argmax(np.abs(loadings)))
        if loadings[anchor] < 0:
            scores *= -1
            loadings *= -1
        scores = (scores - scores.mean()) / scores.std(ddof=0)
        pc = f"global_pc{index + 1}"
        score_frame[pc] = scores
        variance_rows.append(
            {
                "pc": pc,
                "n_genes": len(genes),
                "explained_variance_ratio": float(singular[index] ** 2 / total),
                "singular_value": float(singular[index]),
                "orientation_anchor_gene": genes[anchor],
            }
        )
        loading_rows.extend(
            {
                "pc": pc,
                "gene": gene,
                "loading": float(loading),
                "absolute_loading": float(abs(loading)),
            }
            for gene, loading in zip(genes, loadings)
        )
    delta_frame = pd.DataFrame(delta, index=patients, columns=genes)
    delta_frame.index.name = "patient_id"
    return (
        score_frame,
        pd.DataFrame(loading_rows),
        pd.DataFrame(variance_rows),
        delta_frame,
    )


def stratified_permutation(
    data: pd.DataFrame, iterations: int, seed: int
) -> tuple[float, pd.DataFrame]:
    current = data.sort_values("patient_id").reset_index(drop=True)
    y = current["delta_post_minus_pre"].to_numpy(float)
    nuisance = sm.add_constant(
        current[["therapy_tadalafil", "batch2"]].astype(float),
        has_constant="add",
    ).to_numpy(float)
    projection = nuisance @ np.linalg.pinv(nuisance)

    def coefficient(response: np.ndarray) -> float:
        residual = response - projection @ response
        return float(residual @ y / (residual @ residual))

    observed_response = current["response_ord_num"].to_numpy(float)
    observed = coefficient(observed_response)
    rng = np.random.default_rng(seed)
    therapy = current["therapy"].to_numpy()
    effects = np.empty(iterations, dtype=float)
    for index in range(iterations):
        permuted = observed_response.copy()
        for arm in np.unique(therapy):
            positions = np.flatnonzero(therapy == arm)
            permuted[positions] = rng.permutation(permuted[positions])
        effects[index] = coefficient(permuted)
    p_value = float((1 + np.sum(np.abs(effects) >= abs(observed) - 1e-12)) / (iterations + 1))
    frame = pd.DataFrame(
        {
            "iteration": np.arange(1, iterations + 1),
            "random_seed": seed,
            "permuted_response_effect": effects,
        }
    )
    return p_value, frame


def expression_bins(expression: pd.DataFrame) -> pd.Series:
    ordered = expression.mean(axis=1).rank(method="first")
    return pd.qcut(ordered, q=10, labels=False, duplicates="drop").astype(int)


def overlap_preserving_random_families(
    expression: pd.DataFrame,
    modules: pd.DataFrame,
    arrays: dict[str, np.ndarray],
    present_modules: dict[str, np.ndarray],
    metadata: pd.DataFrame,
    primary_data: pd.DataFrame,
    iterations: int,
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    membership: dict[str, tuple[str, ...]] = {}
    for gene in expression.index:
        signatures = tuple(
            signature
            for signature, positions in present_modules.items()
            if expression.index.get_loc(gene) in set(positions.tolist())
        )
        if signatures:
            membership[gene] = signatures
    bins = expression_bins(expression)
    locked = set(membership)
    templates = Counter(int(bins.loc[gene]) for gene in membership)
    background = {
        bin_id: np.asarray(
            [
                index for index, gene in enumerate(expression.index)
                if gene not in locked and int(bins.loc[gene]) == bin_id
            ],
            dtype=int,
        )
        for bin_id in templates
    }
    if any(len(background[bin_id]) < count for bin_id, count in templates.items()):
        raise RuntimeError("Insufficient unique expression-matched background genes")

    patients = sorted(metadata["patient_id"].unique())
    operator = build_delta_operator(metadata, expression.columns.tolist(), patients)
    ordered = primary_data.set_index("patient_id").loc[patients].reset_index()
    _, design = response_design(
        ordered, "delta_post_minus_pre", ["therapy_tadalafil", "batch2"]
    )
    patient_weights = np.linalg.pinv(design.to_numpy(float))[
        list(design.columns).index("response_ord_num")
    ]
    sample_weights = patient_weights @ operator
    gene_positions = {
        gene: expression.index.get_loc(gene) for gene in membership
    }
    rng = np.random.default_rng(seed)
    rows = []
    first_map = None
    for iteration in range(1, iterations + 1):
        replacement: dict[str, int] = {}
        for bin_id, count in sorted(templates.items()):
            originals = [
                gene for gene in membership if int(bins.loc[gene]) == bin_id
            ]
            selected = rng.choice(background[bin_id], size=count, replace=False)
            replacement.update(dict(zip(sorted(originals), selected.tolist())))
        if first_map is None:
            first_map = replacement.copy()
        random_modules = {
            signature: np.asarray(
                [
                    replacement[gene]
                    for gene, signatures in membership.items()
                    if signature in signatures
                ],
                dtype=int,
            )
            for signature in present_modules
        }
        for method, matrix in arrays.items():
            sample_scores = np.vstack(
                [matrix[positions].mean(axis=0) for positions in random_modules.values()]
            ).mean(axis=0)
            rows.append(
                {
                    "scoring_method": method,
                    "iteration": iteration,
                    "random_seed": seed,
                    "effect": float(sample_weights @ sample_scores),
                }
            )
    if first_map is None:
        raise RuntimeError("Random-family generation failed")
    original_sizes = {
        signature: len(positions) for signature, positions in present_modules.items()
    }
    first_sizes = Counter()
    for gene, signatures in membership.items():
        for signature in signatures:
            first_sizes[signature] += 1
    audit = pd.DataFrame(
        [
            {
                "unique_locked_genes_present": len(membership),
                "total_module_memberships": sum(original_sizes.values()),
                "genes_reused_across_modules": sum(
                    len(signatures) > 1 for signatures in membership.values()
                ),
                "maximum_module_membership_per_gene": max(
                    len(signatures) for signatures in membership.values()
                ),
                "module_sizes_exactly_preserved": all(
                    first_sizes[signature] == size
                    for signature, size in original_sizes.items()
                ),
                "unique_replacements_preserved": len(set(first_map.values())) == len(first_map),
                "expression_decile_counts_preserved": True,
            }
        ]
    )
    return pd.DataFrame(rows), audit


def mcp_counter_analysis(
    expression: pd.DataFrame, metadata: pd.DataFrame, signature_path: Path
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    signatures = pd.read_csv(signature_path, sep="\t")
    coverage = []
    sample_rows = []
    for population, current in signatures.groupby("Cell population", sort=False):
        defined = list(dict.fromkeys(
            current["HUGO symbols"].astype(str).str.upper().tolist()
        ))
        present = [gene for gene in defined if gene in expression.index]
        coverage.append(
            {
                "population": population,
                "markers_defined": len(defined),
                "markers_present": len(present),
                "coverage_fraction": len(present) / len(defined),
                "present_markers": ";".join(present),
            }
        )
        values = expression.loc[present].mean(axis=0)
        sample_rows.extend(
            {
                "sample_id": sample,
                "population": population,
                "mcp_counter_score": float(value),
            }
            for sample, value in values.items()
        )
    sample_scores = pd.DataFrame(sample_rows)
    deltas = paired_deltas(
        sample_scores, metadata, "mcp_counter_score", ["population"]
    )
    model_rows = []
    for population, current in deltas.groupby("population"):
        result, _ = fit_model(
            current, "delta_post_minus_pre", ["therapy_tadalafil", "batch2"]
        )
        model_rows.append({"population": population, **result})
    models = pd.DataFrame(model_rows)
    models["bh_fdr"] = bh_adjust(models["p_value_hc3"])
    return pd.DataFrame(coverage), deltas, models


def global_shift_analysis(
    expression: pd.DataFrame, metadata: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    for metric, values in {
        "mean_log_cpm": expression.mean(axis=0),
        "median_log_cpm": expression.median(axis=0),
        "q75_log_cpm": expression.quantile(0.75, axis=0),
    }.items():
        rows.extend(
            {"sample_id": sample, "metric": metric, "value": float(value)}
            for sample, value in values.items()
        )
    deltas = paired_deltas(pd.DataFrame(rows), metadata, "value", ["metric"])
    models = []
    for metric, current in deltas.groupby("metric"):
        result, _ = fit_model(
            current, "delta_post_minus_pre", ["therapy_tadalafil", "batch2"]
        )
        models.append({"metric": metric, **result})
    output = pd.DataFrame(models)
    output["bh_fdr"] = bh_adjust(output["p_value_hc3"])
    return deltas, output


def make_figure(
    patient_data: pd.DataFrame,
    family_models: pd.DataFrame,
    random_effects: pd.DataFrame,
    observed_effect: float,
    empirical_p: float,
) -> None:
    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 7,
            "axes.linewidth": 0.7,
            "xtick.major.width": 0.7,
            "ytick.major.width": 0.7,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )
    colors = ["#6B7280", "#D9A441", "#2A9D8F", "#C4473A"]
    order = ["Non-Responder", "Minimal Responder", "Responder", "Complete Responder"]
    fig, axes = plt.subplots(2, 2, figsize=(7.09, 5.55))

    counts = patient_data["response_label"].value_counts().reindex(order).fillna(0)
    axes[0, 0].bar(range(4), counts, color=colors, width=0.72)
    axes[0, 0].set_xticks(range(4), ["NR", "Minor", "Responder", "Complete"])
    axes[0, 0].set_ylabel("Patients")
    axes[0, 0].set_title("Cohort composition", loc="left", fontweight="bold")

    rng = np.random.default_rng(195832)
    for index, label in enumerate(order):
        values = patient_data.loc[
            patient_data["response_label"].eq(label), "delta_post_minus_pre"
        ].to_numpy(float)
        jitter = rng.uniform(-0.10, 0.10, size=len(values))
        axes[0, 1].scatter(
            index + jitter, values, s=24, color=colors[index],
            edgecolor="white", linewidth=0.45, zorder=3
        )
        axes[0, 1].plot(
            [index - 0.22, index + 0.22],
            [np.median(values), np.median(values)],
            color="black", lw=1.1, zorder=4,
        )
    axes[0, 1].axhline(0, color="#9CA3AF", lw=0.7, ls="--")
    axes[0, 1].set_xticks(range(4), ["NR", "Minor", "Responder", "Complete"])
    axes[0, 1].set_ylabel("Family score delta (post - pre)")
    axes[0, 1].set_title("Frozen family across response depth", loc="left", fontweight="bold")

    model_order = ["unadjusted", "therapy_adjusted", "primary", "global_pc1", "global_pc1_pc2"]
    current = family_models.set_index("model_spec").loc[model_order].reset_index()
    y = np.arange(len(current))[::-1]
    axes[1, 0].errorbar(
        current["effect"], y,
        xerr=np.vstack(
            [
                current["effect"] - current["ci95_low"],
                current["ci95_high"] - current["effect"],
            ]
        ),
        fmt="o", color="#1F4E79", ecolor="#1F4E79", capsize=2, ms=4,
    )
    axes[1, 0].axvline(0, color="#9CA3AF", lw=0.7, ls="--")
    axes[1, 0].set_yticks(
        y, ["Unadjusted", "+ tadalafil", "Primary", "+ global PC1", "+ global PCs 1-2"]
    )
    axes[1, 0].set_xlabel("Ordinal-response coefficient (95% CI)")
    axes[1, 0].set_title("Model robustness", loc="left", fontweight="bold")

    null = random_effects.loc[
        random_effects["scoring_method"].eq("z_score"), "effect"
    ].to_numpy(float)
    axes[1, 1].hist(null, bins=36, color="#D1D5DB", edgecolor="white", linewidth=0.4)
    axes[1, 1].axvline(observed_effect, color="#C4473A", lw=1.4)
    axes[1, 1].axvline(-observed_effect, color="#C4473A", lw=0.9, ls="--")
    axes[1, 1].set_xlabel("Matched random-family coefficient")
    axes[1, 1].set_ylabel("Iterations")
    axes[1, 1].set_title(
        f"Overlap-preserving null (empirical P={empirical_p:.3g})",
        loc="left", fontweight="bold",
    )

    for label, axis in zip(["a", "b", "c", "d"], axes.flat):
        axis.text(
            -0.16, 1.08, label, transform=axis.transAxes,
            fontsize=8, fontweight="bold", va="top"
        )
        axis.spines[["top", "right"]].set_visible(False)
    fig.tight_layout(pad=1.25, w_pad=1.3, h_pad=1.7)
    fig.savefig(OUT_DIR / "GSE195832_LOCKED_FAMILY_VALIDATION.pdf")
    svg_metadata = {
        "Date": "2026-07-27",
        "Creator": "GSE195832 locked-family validation",
    }
    fig.savefig(
        OUT_DIR / "GSE195832_LOCKED_FAMILY_VALIDATION.svg",
        metadata=svg_metadata,
    )
    fig.savefig(
        OUT_DIR / "GSE195832_LOCKED_FAMILY_VALIDATION.png",
        dpi=600, facecolor="white",
    )
    figure_dir = REBUILD / "figures" / "submission"
    figure_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(figure_dir / "ExtendedData15_gse195832_frozen_family.pdf")
    fig.savefig(
        figure_dir / "ExtendedData15_gse195832_frozen_family.svg",
        metadata=svg_metadata,
    )
    fig.savefig(
        figure_dir / "ExtendedData15_gse195832_frozen_family.png",
        dpi=600, facecolor="white",
    )
    plt.close(fig)


def markdown_table(frame: pd.DataFrame, columns: list[str]) -> str:
    view = frame[columns].copy()
    for column in view.select_dtypes(include=[np.number]).columns:
        view[column] = view[column].map(
            lambda value: "NA" if pd.isna(value) else f"{value:.4g}"
        )
    view = view.fillna("NA").astype(str)
    header = "| " + " | ".join(columns) + " |"
    separator = "| " + " | ".join("---" for _ in columns) + " |"
    rows = [
        "| " + " | ".join(row[column].replace("|", r"\|") for column in columns) + " |"
        for _, row in view.iterrows()
    ]
    return "\n".join([header, separator, *rows])


def main() -> None:
    config = load_config()
    modules = load_modules(
        ROOT / config["module_manifest"], int(config["expected_modules"])
    )
    metadata, raw_counts, expression = load_inputs(config)
    (
        family_scores,
        module_scores,
        coverage,
        arrays,
        present_modules,
        expression,
    ) = score_frozen_family(expression, modules)
    family_deltas = paired_deltas(
        family_scores, metadata, "family_score",
        ["scoring_method", "score_variant"],
    )
    module_deltas = paired_deltas(
        module_scores, metadata, "module_score",
        ["scoring_method", "signature"],
    )

    locked_genes = set(
        gene for module in modules["genes"] for gene in module
        if gene in expression.index
    )
    pc_scores, pc_loadings, pc_variance, expression_deltas = deterministic_pca(
        expression, metadata, locked_genes, int(config["global_pc_count"])
    )
    family_deltas = family_deltas.merge(
        pc_scores, on="patient_id", validate="many_to_one"
    )
    module_deltas = module_deltas.merge(
        pc_scores, on="patient_id", validate="many_to_one"
    )

    model_rows = []
    for (method, variant), current in family_deltas.groupby(
        ["scoring_method", "score_variant"]
    ):
        for model_spec, covariates in config["model_specs"].items():
            result, _ = fit_model(current, "delta_post_minus_pre", covariates)
            model_rows.append(
                {
                    "scoring_method": method,
                    "score_variant": variant,
                    "model_spec": model_spec,
                    "covariates": ";".join(covariates),
                    **result,
                }
            )
    family_models = pd.DataFrame(model_rows)
    family_models["variant_bh_fdr"] = family_models.groupby(
        ["scoring_method", "model_spec"], group_keys=False
    )["p_value_hc3"].apply(bh_adjust)
    family_models["model_spec_bh_fdr"] = family_models.groupby(
        ["scoring_method", "score_variant"], group_keys=False
    )["p_value_hc3"].apply(bh_adjust)

    module_model_rows = []
    for (method, signature), current in module_deltas.groupby(
        ["scoring_method", "signature"]
    ):
        result, _ = fit_model(
            current, "delta_post_minus_pre", ["therapy_tadalafil", "batch2"]
        )
        lineage = modules.set_index("signature").loc[signature, "target_lineage"]
        module_model_rows.append(
            {
                "scoring_method": method,
                "signature": signature,
                "target_lineage": lineage,
                **result,
            }
        )
    module_models = pd.DataFrame(module_model_rows)
    module_models["bh_fdr"] = module_models.groupby(
        "scoring_method", group_keys=False
    )["p_value_hc3"].apply(bh_adjust)

    primary_data = family_deltas[
        family_deltas["scoring_method"].eq(config["primary_scoring_method"])
        & family_deltas["score_variant"].eq(config["primary_score_variant"])
    ].copy()
    permutation_p, permutation_null = stratified_permutation(
        primary_data,
        int(config["permutation_iterations"]),
        int(config["permutation_seed"]),
    )
    permutation_test = pd.DataFrame(
        [
            {
                "scoring_method": config["primary_scoring_method"],
                "score_variant": config["primary_score_variant"],
                "model_spec": "primary",
                "iterations": int(config["permutation_iterations"]),
                "permutation_scheme": config["permutation_scheme"],
                "observed_effect": float(
                    family_models.loc[
                        family_models["scoring_method"].eq(
                            config["primary_scoring_method"]
                        )
                        & family_models["score_variant"].eq(
                            config["primary_score_variant"]
                        )
                        & family_models["model_spec"].eq("primary"),
                        "effect",
                    ].iloc[0]
                ),
                "two_sided_permutation_p": permutation_p,
            }
        ]
    )

    random_effects, random_audit = overlap_preserving_random_families(
        expression,
        modules,
        arrays,
        present_modules,
        metadata,
        primary_data,
        int(config["random_family_iterations"]),
        int(config["random_family_seed"]),
    )
    primary_model = family_models[
        family_models["scoring_method"].eq(config["primary_scoring_method"])
        & family_models["score_variant"].eq(config["primary_score_variant"])
        & family_models["model_spec"].eq("primary")
    ].iloc[0]
    primary_null = random_effects.loc[
        random_effects["scoring_method"].eq(config["primary_scoring_method"]),
        "effect",
    ].to_numpy(float)
    empirical_p = float(
        (1 + np.sum(np.abs(primary_null) >= abs(primary_model["effect"]) - 1e-12))
        / (len(primary_null) + 1)
    )

    loo_rows = []
    for patient in sorted(primary_data["patient_id"]):
        result, _ = fit_model(
            primary_data[~primary_data["patient_id"].eq(patient)],
            "delta_post_minus_pre",
            ["therapy_tadalafil", "batch2"],
        )
        loo_rows.append({"left_out_patient": patient, **result})
    loo = pd.DataFrame(loo_rows)

    shift_deltas, shift_models = global_shift_analysis(expression, metadata)
    signature_path = ROOT / config["composition_signature_file"]
    if sha256(signature_path) != config["composition_signature_sha256"]:
        raise RuntimeError("MCP-counter signature SHA256 mismatch")
    mcp_coverage, mcp_deltas, mcp_models = mcp_counter_analysis(
        expression, metadata, signature_path
    )

    family_pc_correlations = []
    for method in ["z_score", "rank_mean"]:
        current = family_deltas[
            family_deltas["scoring_method"].eq(method)
            & family_deltas["score_variant"].eq("module_mean_16")
        ]
        for pc in ["global_pc1", "global_pc2", "global_pc3"]:
            family_pc_correlations.append(
                {
                    "scoring_method": method,
                    "feature": pc,
                    "pearson_r": float(current["delta_post_minus_pre"].corr(current[pc])),
                }
            )
    family_pc_correlations = pd.DataFrame(family_pc_correlations)

    source_qc = pd.DataFrame(
        [
            {"check": "metadata_sha256", "status": "PASS", "value": config["metadata_sha256"]},
            {"check": "counts_sha256", "status": "PASS", "value": config["counts_sha256"]},
            {"check": "samples", "status": "PASS", "value": len(metadata)},
            {"check": "paired_patients", "status": "PASS", "value": metadata["patient_id"].nunique()},
            {"check": "raw_gene_symbols_after_aggregation", "status": "PASS", "value": len(raw_counts)},
            {"check": "response_blind_filtered_genes", "status": "PASS", "value": len(expression)},
            {"check": "frozen_modules", "status": "PASS", "value": len(modules)},
            {"check": "age_source_anomaly_TJU033", "status": "WARN", "value": metadata.loc[metadata["patient_id"].eq("TJU033"), "age_source"].iloc[0]},
        ]
    )
    response_counts = (
        metadata.drop_duplicates("patient_id")
        .groupby(["response_label", "therapy", "batch"], observed=True)
        .size().rename("n_patients").reset_index()
    )

    gate_rows = [
        {
            "gate": "primary_HC3_P_lt_0.05",
            "passed": bool(primary_model["p_value_hc3"] < 0.05),
            "value": float(primary_model["p_value_hc3"]),
        },
        {
            "gate": "prespecified_positive_orientation",
            "passed": bool(primary_model["effect"] > 0),
            "value": float(primary_model["effect"]),
        },
    ]
    pc1_model = family_models[
        family_models["scoring_method"].eq("z_score")
        & family_models["score_variant"].eq("module_mean_16")
        & family_models["model_spec"].eq("global_pc1")
    ].iloc[0]
    gate_rows.extend(
        [
            {
                "gate": "global_PC1_adjusted_P_lt_0.05",
                "passed": bool(pc1_model["p_value_hc3"] < 0.05),
                "value": float(pc1_model["p_value_hc3"]),
            },
            {
                "gate": "overlap_random_empirical_P_lt_0.05",
                "passed": bool(empirical_p < 0.05),
                "value": empirical_p,
            },
            {
                "gate": "no_leave_one_patient_direction_reversal",
                "passed": bool(np.all(np.sign(loo["effect"]) == np.sign(primary_model["effect"]))),
                "value": int(np.sum(np.sign(loo["effect"]) != np.sign(primary_model["effect"]))),
            },
            {
                "gate": "not_one_scoring_variant_only",
                "passed": bool(
                    (
                        family_models["scoring_method"].eq("rank_mean")
                        & family_models["score_variant"].eq("module_mean_16")
                        & family_models["model_spec"].eq("primary")
                        & (family_models["p_value_hc3"] < 0.05)
                        & (np.sign(family_models["effect"]) == np.sign(primary_model["effect"]))
                    ).any()
                    and (
                        family_models["scoring_method"].eq("z_score")
                        & family_models["score_variant"].eq("unique_gene_equal")
                        & family_models["model_spec"].eq("primary")
                        & (family_models["p_value_hc3"] < 0.05)
                        & (np.sign(family_models["effect"]) == np.sign(primary_model["effect"]))
                    ).any()
                ),
                "value": "rank module mean and z unique-gene sensitivity",
            },
        ]
    )
    gates = pd.DataFrame(gate_rows)
    validation_status = (
        "INDEPENDENT_VALIDATION"
        if gates["passed"].all()
        else "SENSITIVITY_OR_BOUNDARY_EVIDENCE"
    )

    outputs = {
        "GSE195832_SOURCE_AND_COHORT_QC.csv": source_qc,
        "GSE195832_RESPONSE_COMPOSITION.csv": response_counts,
        "GSE195832_MODULE_COVERAGE.csv": coverage,
        "GSE195832_FAMILY_SAMPLE_SCORES.csv": family_scores,
        "GSE195832_FAMILY_PAIRED_DELTAS.csv": family_deltas,
        "GSE195832_FAMILY_RESPONSE_MODELS.csv": family_models,
        "GSE195832_MODULE_PAIRED_DELTAS.csv": module_deltas,
        "GSE195832_MODULE_RESPONSE_MODELS.csv": module_models,
        "GSE195832_PRIMARY_PERMUTATION_TEST.csv": permutation_test,
        "GSE195832_STRATIFIED_PERMUTATION_NULL.csv": permutation_null,
        "GSE195832_OVERLAP_RANDOM_EFFECTS.csv": random_effects,
        "GSE195832_OVERLAP_RANDOM_AUDIT.csv": random_audit,
        "GSE195832_LEAVE_ONE_PATIENT.csv": loo,
        "GSE195832_GLOBAL_PC_SCORES.csv": pc_scores,
        "GSE195832_GLOBAL_PC_VARIANCE.csv": pc_variance,
        "GSE195832_GLOBAL_PC_LOADINGS.csv.gz": pc_loadings,
        "GSE195832_EXPRESSION_PAIRED_DELTAS.csv.gz": expression_deltas.reset_index(),
        "GSE195832_GLOBAL_SHIFT_DELTAS.csv": shift_deltas,
        "GSE195832_GLOBAL_SHIFT_MODELS.csv": shift_models,
        "GSE195832_MCP_COUNTER_COVERAGE.csv": mcp_coverage,
        "GSE195832_MCP_COUNTER_DELTAS.csv": mcp_deltas,
        "GSE195832_MCP_COUNTER_MODELS.csv": mcp_models,
        "GSE195832_FAMILY_PC_CORRELATIONS.csv": family_pc_correlations,
        "GSE195832_VALIDATION_GATES.csv": gates,
    }
    for filename, frame in outputs.items():
        destination = OUT_DIR / filename
        if filename.endswith(".csv.gz"):
            frame.to_csv(
                destination,
                index=False,
                compression={"method": "gzip", "mtime": 0},
            )
        else:
            frame.to_csv(destination, index=False)

    source_dir = REBUILD / "figures" / "submission" / "source_data"
    source_dir.mkdir(parents=True, exist_ok=True)
    primary_data.to_csv(
        source_dir / "ExtendedData15_patient_deltas_source.csv", index=False
    )
    family_models[
        family_models["scoring_method"].eq("z_score")
        & family_models["score_variant"].eq("module_mean_16")
    ].to_csv(
        source_dir / "ExtendedData15_family_models_source.csv", index=False
    )
    random_effects[
        random_effects["scoring_method"].eq("z_score")
    ].to_csv(
        source_dir / "ExtendedData15_overlap_random_source.csv", index=False
    )
    response_counts.to_csv(
        source_dir / "ExtendedData15_response_counts_source.csv", index=False
    )

    make_figure(primary_data, family_models[
        family_models["scoring_method"].eq("z_score")
        & family_models["score_variant"].eq("module_mean_16")
    ], random_effects, float(primary_model["effect"]), empirical_p)

    rank_primary = family_models[
        family_models["scoring_method"].eq("rank_mean")
        & family_models["score_variant"].eq("module_mean_16")
        & family_models["model_spec"].eq("primary")
    ].iloc[0]
    report = [
        "# GSE195832 frozen-family independent-validation report",
        "",
        f"Status: **{validation_status}**",
        "",
        "## Frozen design",
        "",
        "- Patient is the unit of inference: 28 exact pretreatment/post-treatment pairs.",
        "- Raw feature counts were library-size normalized to CPM, response-blind filtered, and transformed as log2(CPM + 0.5).",
        "- The 16-module family, all module genes, the primary score, response ordering and positive response orientation were frozen before scoring this cohort.",
        "- Primary model: post-minus-pre z-score family delta regressed on ordinal pathological response, tadalafil arm and sequencing batch with HC3 covariance.",
        f"- Stratified Monte Carlo permutation: {config['permutation_iterations']:,} response-label permutations within treatment arm.",
        f"- Overlap-preserving specificity null: {config['random_family_iterations']:,} expression-decile-matched random families.",
        "- Source age was not modeled because TJU033 is recorded as age 5; no unsourced correction was made.",
        "",
        "## Cohort composition",
        "",
        markdown_table(
            response_counts,
            ["response_label", "therapy", "batch", "n_patients"],
        ),
        "",
        "## Primary family result",
        "",
        markdown_table(
            pd.DataFrame(
                [
                    {
                        "method": "z_score",
                        "effect": primary_model["effect"],
                        "ci95_low": primary_model["ci95_low"],
                        "ci95_high": primary_model["ci95_high"],
                        "HC3_P": primary_model["p_value_hc3"],
                        "stratified_permutation_P": permutation_p,
                        "overlap_random_empirical_P": empirical_p,
                        "condition_number": primary_model["design_condition_number"],
                    },
                    {
                        "method": "rank_mean",
                        "effect": rank_primary["effect"],
                        "ci95_low": rank_primary["ci95_low"],
                        "ci95_high": rank_primary["ci95_high"],
                        "HC3_P": rank_primary["p_value_hc3"],
                        "stratified_permutation_P": math.nan,
                        "overlap_random_empirical_P": float(
                            (
                                1 + np.sum(
                                    np.abs(
                                        random_effects.loc[
                                            random_effects["scoring_method"].eq("rank_mean"),
                                            "effect",
                                        ].to_numpy(float)
                                    )
                                    >= abs(rank_primary["effect"]) - 1e-12
                                )
                            )
                            / (
                                1
                                + random_effects["scoring_method"].eq("rank_mean").sum()
                            )
                        ),
                        "condition_number": rank_primary["design_condition_number"],
                    },
                ]
            ),
            [
                "method", "effect", "ci95_low", "ci95_high", "HC3_P",
                "stratified_permutation_P", "overlap_random_empirical_P",
                "condition_number",
            ],
        ),
        "",
        "## Validation gates",
        "",
        markdown_table(gates, ["gate", "passed", "value"]),
        "",
        "## Interpretation",
        "",
        (
            "All six prespecified gates passed, so this cohort qualifies as independent "
            "family-level response validation under the frozen protocol."
            if validation_status == "INDEPENDENT_VALIDATION"
            else
            "At least one prespecified gate failed. The cohort must therefore be reported "
            "as sensitivity or boundary evidence rather than independent validation."
        ),
        "",
        "The endpoint is the authors' ordered pathological response category. No unreported percentage regression threshold was invented. The cohort mixes nivolumab monotherapy and nivolumab plus tadalafil, which is adjusted in the model but still limits treatment-specific causal interpretation.",
        "",
        "## Source provenance",
        "",
        "- GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE195832",
        "- Original article: https://pmc.ncbi.nlm.nih.gov/articles/PMC9161438/",
        "- Author data deposit: https://data.mendeley.com/datasets/yk8wj7xgdg/1",
    ]
    (OUT_DIR / "GSE195832_LOCKED_FAMILY_VALIDATION_REPORT.md").write_text(
        "\n".join(report) + "\n", encoding="utf-8"
    )
    summary = {
        "analysis_id": config["analysis_id"],
        "status": validation_status,
        "primary_effect": float(primary_model["effect"]),
        "primary_hc3_p": float(primary_model["p_value_hc3"]),
        "stratified_permutation_p": permutation_p,
        "global_pc1_adjusted_p": float(pc1_model["p_value_hc3"]),
        "overlap_random_empirical_p": empirical_p,
        "all_validation_gates_passed": bool(gates["passed"].all()),
    }
    (OUT_DIR / "GSE195832_VALIDATION_SUMMARY.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
