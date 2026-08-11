#!/usr/bin/env python3
"""Frozen patient-level validation of the 16-module family in GSE232240.

The expression matrix is streamed by gene. Cells are aggregated immediately to
patient/timepoint/lineage pseudobulks, avoiding a dense cell-by-gene object.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import itertools
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import rankdata


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workspace", type=Path, required=True)
    parser.add_argument("--metadata-only", action="store_true")
    return parser.parse_args()


def read_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def exact_binary_permutation(values: np.ndarray, labels: np.ndarray) -> tuple[float, float, int]:
    labels = labels.astype(int)
    n = len(values)
    n_pos = int(labels.sum())
    observed = float(values[labels == 1].mean() - values[labels == 0].mean())
    extreme = 0
    denominator = 0
    all_indices = np.arange(n)
    for positive_indices in itertools.combinations(range(n), n_pos):
        mask = np.zeros(n, dtype=bool)
        mask[list(positive_indices)] = True
        effect = float(values[mask].mean() - values[~mask].mean())
        denominator += 1
        if abs(effect) >= abs(observed) - 1e-12:
            extreme += 1
    return observed, extreme / denominator, denominator


def exact_binary_distribution(
    values: np.ndarray,
    labels: np.ndarray,
    patient_ids: list[str],
) -> pd.DataFrame:
    labels = labels.astype(bool)
    observed_positive = frozenset(np.flatnonzero(labels))
    observed_effect = float(values[labels].mean() - values[~labels].mean())
    rows = []
    for assignment_id, positive_indices in enumerate(
        itertools.combinations(range(len(values)), int(labels.sum())), start=1
    ):
        positive_set = frozenset(positive_indices)
        mask = np.zeros(len(values), dtype=bool)
        mask[list(positive_indices)] = True
        effect = float(values[mask].mean() - values[~mask].mean())
        rows.append(
            {
                "assignment_id": assignment_id,
                "RE_patients": ";".join(patient_ids[index] for index in positive_indices),
                "effect_RE_minus_NR": effect,
                "absolute_effect_at_least_observed": abs(effect) >= abs(observed_effect) - 1e-12,
                "is_observed_assignment": positive_set == observed_positive,
            }
        )
    return pd.DataFrame(rows)


def bh_adjust(p_values: np.ndarray) -> np.ndarray:
    order = np.argsort(p_values)
    ranked = p_values[order]
    adjusted = np.minimum.accumulate((ranked * len(ranked) / np.arange(1, len(ranked) + 1))[::-1])[::-1]
    output = np.empty_like(adjusted)
    output[order] = np.minimum(adjusted, 1.0)
    return output


def audit_metadata(metadata_path: Path, config: dict, output_dir: Path) -> tuple[pd.DataFrame, list[str]]:
    metadata = pd.read_csv(metadata_path, sep="\t", compression="gzip", dtype=str)
    missing = sorted(set(config["required_metadata_columns"]) - set(metadata.columns))
    if missing:
        raise ValueError(f"Missing required metadata columns: {missing}")
    if metadata["cell_id"].duplicated().any():
        raise ValueError("Metadata cell_id values are not unique")

    patient_response_counts = metadata.groupby("patient")["response"].nunique()
    if patient_response_counts.max() != 1:
        raise ValueError("At least one patient has inconsistent response labels")

    lineage_map = config["lineage_map"]
    metadata["lineage"] = metadata["cell_type"].map(lineage_map)
    metadata["excluded_reason"] = metadata["patient"].map(config["excluded_patients"])
    eligible_cells = metadata[metadata["excluded_reason"].isna() & metadata["lineage"].notna()].copy()

    counts = (
        eligible_cells.groupby(["patient", "response", "lineage", "timepoint"])
        .size()
        .unstack(fill_value=0)
        .reset_index()
    )
    for timepoint in config["timepoint_levels"]:
        if timepoint not in counts:
            counts[timepoint] = 0
    minimum = int(config["minimum_cells_per_patient_timepoint_lineage"])
    counts["eligible_lineage_pair"] = counts["pre"].ge(minimum) & counts["post"].ge(minimum)

    eligible_sets = {
        lineage: set(counts.loc[(counts["lineage"] == lineage) & counts["eligible_lineage_pair"], "patient"])
        for lineage in sorted(set(lineage_map.values()))
    }
    primary_patients = sorted(set.intersection(*eligible_sets.values()))
    counts["eligible_primary_both_lineages"] = counts["patient"].isin(primary_patients)
    counts.to_csv(output_dir / "GSE232240_metadata_eligibility.csv", index=False)

    mapping = metadata[["patient", "response"]].drop_duplicates()
    primary_mapping = mapping[mapping["patient"].isin(primary_patients)].sort_values("patient")
    response_counts = primary_mapping["response"].value_counts().to_dict()
    denominator = math.comb(len(primary_mapping), int(response_counts.get("RE", 0)))
    summary = {
        "metadata_cells": int(len(metadata)),
        "metadata_patients": int(metadata["patient"].nunique()),
        "excluded_patients": config["excluded_patients"],
        "eligible_T_cell_patients": len(eligible_sets.get("T_cell", set())),
        "eligible_Myeloid_patients": len(eligible_sets.get("Myeloid", set())),
        "primary_both_lineages_patients": len(primary_patients),
        "primary_response_counts": response_counts,
        "exact_permutation_denominator": denominator,
        "primary_patients": primary_patients,
        "status": "metadata_gate_passed",
    }
    with (output_dir / "GSE232240_metadata_audit.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)
    return metadata, primary_patients


def stream_pseudobulk(
    count_path: Path,
    metadata: pd.DataFrame,
    primary_patients: list[str],
) -> tuple[list[str], np.ndarray, list[tuple[str, str, str]], np.ndarray]:
    group_keys = [
        (patient, timepoint, lineage)
        for patient in primary_patients
        for timepoint in ("pre", "post")
        for lineage in ("T_cell", "Myeloid")
    ]
    group_lookup = {key: index for index, key in enumerate(group_keys)}
    metadata_index = metadata.set_index("cell_id", drop=False)

    with gzip.open(count_path, "rt", encoding="utf-8-sig", errors="strict") as handle:
        header = handle.readline().rstrip("\r\n").split("\t")
        if len(header) < 2:
            raise ValueError("Count matrix header has fewer than two columns")
        # This GEO file starts directly with cell IDs, whereas some exported
        # matrices include a leading gene-column label. Support both layouts.
        if header[0] in metadata_index.index:
            cell_ids = header
        elif len(header) > 1 and header[1] in metadata_index.index:
            cell_ids = header[1:]
        else:
            raise ValueError("Count matrix header does not match metadata cell identifiers")
        if len(cell_ids) != len(set(cell_ids)):
            raise ValueError("Count matrix cell identifiers are not unique")

        group_index = np.full(len(cell_ids), -1, dtype=np.int32)
        missing_metadata = []
        for column_index, cell_id in enumerate(cell_ids):
            if cell_id not in metadata_index.index:
                missing_metadata.append(cell_id)
                continue
            row = metadata_index.loc[cell_id]
            lineage = row["lineage"]
            key = (row["patient"], row["timepoint"], lineage)
            if key in group_lookup:
                group_index[column_index] = group_lookup[key]
        if missing_metadata:
            raise ValueError(f"{len(missing_metadata)} count columns lack metadata; first={missing_metadata[0]}")
        included_columns = np.flatnonzero(group_index >= 0)
        included_groups = group_index[included_columns]
        if len(included_columns) == 0:
            raise ValueError("No count columns passed the frozen primary eligibility gate")

        gene_sums: dict[str, np.ndarray] = {}
        library_sizes = np.zeros(len(group_keys), dtype=np.float64)
        expected_values = len(cell_ids)
        for line_number, line in enumerate(handle, start=2):
            gene, separator, value_text = line.rstrip("\r\n").partition("\t")
            if not separator:
                raise ValueError(f"Malformed count row at line {line_number}")
            values = np.fromstring(value_text, sep="\t", dtype=np.float64)
            if len(values) != expected_values:
                raise ValueError(
                    f"Count row {line_number} has {len(values)} values; expected {expected_values}"
                )
            sums = np.bincount(
                included_groups,
                weights=values[included_columns],
                minlength=len(group_keys),
            ).astype(np.float64)
            gene_key = gene.strip().upper()
            if not gene_key:
                continue
            if gene_key in gene_sums:
                gene_sums[gene_key] += sums
            else:
                gene_sums[gene_key] = sums
            library_sizes += sums

    if np.any(library_sizes <= 0):
        bad = [group_keys[index] for index in np.flatnonzero(library_sizes <= 0)]
        raise ValueError(f"Zero pseudobulk library sizes: {bad}")
    genes = sorted(gene_sums)
    matrix = np.vstack([gene_sums[gene] for gene in genes])
    return genes, matrix, group_keys, library_sizes


def standardize_by_lineage(log_cpm: np.ndarray, group_keys: list[tuple[str, str, str]]) -> np.ndarray:
    z = np.full_like(log_cpm, np.nan, dtype=np.float64)
    for lineage in ("T_cell", "Myeloid"):
        columns = np.array([key[2] == lineage for key in group_keys])
        values = log_cpm[:, columns]
        means = values.mean(axis=1, keepdims=True)
        standard_deviations = values.std(axis=1, ddof=1, keepdims=True)
        standard_deviations[standard_deviations == 0] = np.nan
        z[:, columns] = (values - means) / standard_deviations
    return z


def module_scores(
    score_matrix: np.ndarray,
    genes: list[str],
    group_keys: list[tuple[str, str, str]],
    modules: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    gene_lookup = {gene: index for index, gene in enumerate(genes)}
    rows = []
    coverage = []
    for module in modules.itertuples(index=False):
        defined = [gene.strip().upper() for gene in str(module.genes_defined).split(";") if gene.strip()]
        present = [gene for gene in defined if gene in gene_lookup]
        if not present:
            raise ValueError(f"Frozen module has zero measured genes: {module.signature}")
        lineage_columns = [index for index, key in enumerate(group_keys) if key[2] == module.target_lineage]
        indices = [gene_lookup[gene] for gene in present]
        values = np.nanmean(score_matrix[np.ix_(indices, lineage_columns)], axis=0)
        for column, value in zip(lineage_columns, values):
            patient, timepoint, lineage = group_keys[column]
            rows.append(
                {
                    "patient": patient,
                    "timepoint": timepoint,
                    "lineage": lineage,
                    "signature": module.signature,
                    "score": float(value),
                }
            )
        coverage.append(
            {
                "signature": module.signature,
                "target_lineage": module.target_lineage,
                "genes_defined": len(defined),
                "genes_measured": len(present),
                "coverage_fraction": len(present) / len(defined),
                "measured_gene_symbols": ";".join(present),
            }
        )
    return pd.DataFrame(rows), pd.DataFrame(coverage)


def family_scores(module_table: pd.DataFrame, response_map: dict[str, str], score_name: str) -> pd.DataFrame:
    wide = (
        module_table.groupby(["patient", "timepoint", "signature"], as_index=False)["score"]
        .mean()
        .pivot(index=["patient", "timepoint"], columns="signature", values="score")
    )
    if wide.isna().any().any():
        raise ValueError("At least one primary patient/timepoint lacks a frozen module score")
    family = wide.mean(axis=1).rename(score_name).reset_index()
    family["response"] = family["patient"].map(response_map)
    return family


def patient_delta(table: pd.DataFrame, value_column: str) -> pd.DataFrame:
    wide = table.pivot(index=["patient", "response"], columns="timepoint", values=value_column).reset_index()
    if not {"pre", "post"}.issubset(wide.columns):
        raise ValueError("Missing pre or post score in primary family table")
    wide[f"{value_column}_delta"] = wide["post"] - wide["pre"]
    return wide


def lineage_family_results(
    module_table: pd.DataFrame,
    response_map: dict[str, str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    lineage_scores = (
        module_table.groupby(["patient", "timepoint", "lineage"], as_index=False)["score"]
        .mean()
        .rename(columns={"score": "lineage_family_score"})
    )
    lineage_scores["response"] = lineage_scores["patient"].map(response_map)
    delta_rows = []
    result_rows = []
    for lineage, subset in lineage_scores.groupby("lineage", sort=True):
        delta = patient_delta(subset, "lineage_family_score")
        delta["lineage"] = lineage
        labels = delta["response"].eq("RE").to_numpy(dtype=int)
        effect, exact_p, denominator = exact_binary_permutation(
            delta["lineage_family_score_delta"].to_numpy(dtype=float), labels
        )
        delta_rows.append(delta)
        result_rows.append(
            {
                "lineage": lineage,
                "n_patients": len(delta),
                "effect_RE_minus_NR": effect,
                "exact_p": exact_p,
                "exact_denominator": denominator,
            }
        )
    return pd.concat(delta_rows, ignore_index=True), pd.DataFrame(result_rows)


def global_pc1_delta(
    z_matrix: np.ndarray,
    group_keys: list[tuple[str, str, str]],
    primary_patients: list[str],
) -> pd.Series:
    lineage_delta = {}
    for lineage in ("T_cell", "Myeloid"):
        columns = [index for index, key in enumerate(group_keys) if key[2] == lineage]
        sample_by_gene = z_matrix[:, columns].T
        sample_by_gene = np.nan_to_num(sample_by_gene, nan=0.0)
        u, singular_values, _ = np.linalg.svd(sample_by_gene, full_matrices=False)
        pc1 = u[:, 0] * singular_values[0]
        keyed = {group_keys[column][:2]: pc1[index] for index, column in enumerate(columns)}
        for patient in primary_patients:
            lineage_delta[(patient, lineage)] = keyed[(patient, "post")] - keyed[(patient, "pre")]
    return pd.Series(
        {
            patient: np.mean([lineage_delta[(patient, "T_cell")], lineage_delta[(patient, "Myeloid")]])
            for patient in primary_patients
        },
        name="global_pc1_delta",
    )


def overlap_preserving_null(
    z_matrix: np.ndarray,
    log_cpm: np.ndarray,
    genes: list[str],
    group_keys: list[tuple[str, str, str]],
    modules: pd.DataFrame,
    response_map: dict[str, str],
    observed_effect: float,
    iterations: int,
    seed: int,
) -> tuple[float, pd.DataFrame]:
    gene_lookup = {gene: index for index, gene in enumerate(genes)}
    module_genes = {
        row.signature: [g.strip().upper() for g in str(row.genes_defined).split(";") if g.strip().upper() in gene_lookup]
        for row in modules.itertuples(index=False)
    }
    locked_unique = sorted({gene for values in module_genes.values() for gene in values})
    # Mean pseudobulk abundance is computed without response labels and is used
    # only to match the detectability of null genes to frozen genes.
    detectability = np.nanmean(log_cpm, axis=1)
    finite = np.isfinite(detectability)
    quantile_edges = np.unique(np.quantile(detectability[finite], np.linspace(0, 1, 11)))
    deciles = np.digitize(detectability, quantile_edges[1:-1], right=True)
    locked_set = set(locked_unique)
    candidate_by_decile = {
        decile: [genes[index] for index in np.flatnonzero((deciles == decile) & finite) if genes[index] not in locked_set]
        for decile in range(10)
    }
    source_by_decile = {
        decile: [gene for gene in locked_unique if deciles[gene_lookup[gene]] == decile]
        for decile in range(10)
    }
    for decile, sources in source_by_decile.items():
        if len(candidate_by_decile[decile]) < len(sources):
            raise ValueError(f"Insufficient matched background genes in detectability decile {decile}")

    rng = np.random.default_rng(seed)
    sample_keys = sorted({key[:2] for key in group_keys})
    sample_lookup = {key: index for index, key in enumerate(sample_keys)}
    group_sample_index = np.array([sample_lookup[key[:2]] for key in group_keys], dtype=np.int32)
    patients = sorted({key[0] for key in sample_keys})
    patient_sample_lookup = {
        patient: (sample_lookup[(patient, "pre")], sample_lookup[(patient, "post")]) for patient in patients
    }
    labels = np.array([response_map[patient] == "RE" for patient in patients], dtype=bool)
    effects = []
    for iteration in range(iterations):
        replacement = {}
        for decile, sources in source_by_decile.items():
            if not sources:
                continue
            selected = rng.choice(candidate_by_decile[decile], size=len(sources), replace=False)
            replacement.update(dict(zip(sources, selected)))
        sample_scores = np.zeros(len(sample_keys), dtype=np.float64)
        for module in modules.itertuples(index=False):
            indices = [gene_lookup[replacement[gene]] for gene in module_genes[module.signature]]
            columns = [index for index, key in enumerate(group_keys) if key[2] == module.target_lineage]
            values = np.nanmean(z_matrix[np.ix_(indices, columns)], axis=0)
            np.add.at(sample_scores, group_sample_index[columns], values)
        sample_scores /= len(modules)
        deltas = np.array(
            [sample_scores[post_index] - sample_scores[pre_index] for pre_index, post_index in patient_sample_lookup.values()]
        )
        effect = float(deltas[labels].mean() - deltas[~labels].mean())
        effects.append({"iteration": iteration + 1, "effect": effect})
    null_table = pd.DataFrame(effects)
    empirical_p = (1 + np.sum(np.abs(null_table["effect"]) >= abs(observed_effect))) / (iterations + 1)
    return float(empirical_p), null_table


def main() -> None:
    args = parse_args()
    workspace = args.workspace.resolve()
    config_path = workspace / "03_rebuild/config/gse232240_validation.json"
    config = read_json(config_path)
    output_dir = workspace / "03_rebuild/results/external_validation/GSE232240"
    output_dir.mkdir(parents=True, exist_ok=True)
    metadata_path = workspace / config["metadata_file"]
    count_path = workspace / config["count_file"]
    module_path = workspace / config["module_manifest"]

    if sha256(metadata_path) != config["metadata_sha256"]:
        raise ValueError("Metadata SHA256 does not match the frozen configuration")
    if sha256(module_path) != config["module_manifest_sha256"]:
        raise ValueError("Module-manifest SHA256 does not match the frozen configuration")

    metadata, primary_patients = audit_metadata(metadata_path, config, output_dir)
    print(f"Metadata gate passed: {len(primary_patients)} primary patients")
    if args.metadata_only:
        print("Metadata-only run complete; no expression result was estimated.")
        return
    if not count_path.exists() or count_path.stat().st_size == 0:
        raise FileNotFoundError(f"Count matrix is missing or empty: {count_path}")

    modules = pd.read_csv(module_path)
    if len(modules) != int(config["expected_module_count"]):
        raise ValueError(f"Expected {config['expected_module_count']} modules; found {len(modules)}")
    genes, counts, group_keys, library_sizes = stream_pseudobulk(count_path, metadata, primary_patients)
    eligible_lineage_cells = metadata[
        metadata["patient"].isin(primary_patients) & metadata["lineage"].isin(["T_cell", "Myeloid"])
    ]
    cell_counts = eligible_lineage_cells.groupby(["patient", "timepoint", "lineage"]).size().to_dict()
    group_totals = pd.DataFrame(
        [
            {
                "patient": patient,
                "timepoint": timepoint,
                "lineage": lineage,
                "cell_count": int(cell_counts.get((patient, timepoint, lineage), 0)),
                "library_size": float(library_sizes[index]),
            }
            for index, (patient, timepoint, lineage) in enumerate(group_keys)
        ]
    )
    log_cpm = np.log2(counts / library_sizes[np.newaxis, :] * 1_000_000 + 0.5)
    z_matrix = standardize_by_lineage(log_cpm, group_keys)

    response_map = metadata[["patient", "response"]].drop_duplicates().set_index("patient")["response"].to_dict()
    z_modules, coverage = module_scores(z_matrix, genes, group_keys, modules)
    z_family = family_scores(z_modules, response_map, "family_score")
    z_delta = patient_delta(z_family, "family_score")
    lineage_deltas, lineage_results = lineage_family_results(z_modules, response_map)
    labels = z_delta["response"].eq("RE").to_numpy(dtype=int)
    delta_values = z_delta["family_score_delta"].to_numpy(dtype=float)
    effect, exact_p, exact_denominator = exact_binary_permutation(delta_values, labels)
    exact_null = exact_binary_distribution(delta_values, labels, z_delta["patient"].tolist())

    ranks = np.column_stack([rankdata(log_cpm[:, column], method="average") / len(genes) for column in range(log_cpm.shape[1])])
    rank_modules, _ = module_scores(ranks, genes, group_keys, modules)
    rank_family = family_scores(rank_modules, response_map, "rank_family_score")
    rank_delta = patient_delta(rank_family, "rank_family_score")
    rank_labels = rank_delta["response"].eq("RE").to_numpy(dtype=int)
    rank_effect, rank_p, _ = exact_binary_permutation(
        rank_delta["rank_family_score_delta"].to_numpy(dtype=float), rank_labels
    )

    pc1 = global_pc1_delta(z_matrix, group_keys, primary_patients)
    z_delta["global_pc1_delta"] = z_delta["patient"].map(pc1)
    design = sm.add_constant(
        pd.DataFrame({"response_RE": labels, "global_pc1_delta": z_delta["global_pc1_delta"].to_numpy()})
    )
    hc3 = sm.OLS(delta_values, design).fit(cov_type="HC3")
    adjusted_effect = float(hc3.params["response_RE"])
    adjusted_p = float(hc3.pvalues["response_RE"])

    loo_effects = []
    for patient in z_delta["patient"]:
        keep = z_delta["patient"].ne(patient).to_numpy()
        loo_labels = labels[keep]
        loo_values = delta_values[keep]
        loo_effect = float(loo_values[loo_labels == 1].mean() - loo_values[loo_labels == 0].mean())
        loo_effects.append({"excluded_patient": patient, "effect": loo_effect})
    loo = pd.DataFrame(loo_effects)
    no_loo_reversal = bool((np.sign(loo["effect"]) == np.sign(effect)).all())

    module_delta = z_modules.pivot_table(
        index=["patient", "signature"], columns="timepoint", values="score"
    ).reset_index()
    module_delta["delta"] = module_delta["post"] - module_delta["pre"]
    module_results = []
    for signature, subset in module_delta.groupby("signature", sort=False):
        subset = subset.copy()
        subset["response"] = subset["patient"].map(response_map)
        module_labels = subset["response"].eq("RE").to_numpy(dtype=int)
        module_effect, module_p, module_denominator = exact_binary_permutation(
            subset["delta"].to_numpy(dtype=float), module_labels
        )
        module_results.append(
            {
                "signature": signature,
                "effect_RE_minus_NR": module_effect,
                "exact_p": module_p,
                "exact_denominator": module_denominator,
            }
        )
    module_results = pd.DataFrame(module_results)
    module_results["BH_FDR"] = bh_adjust(module_results["exact_p"].to_numpy())

    empirical_p, null_table = overlap_preserving_null(
        z_matrix,
        log_cpm,
        genes,
        group_keys,
        modules,
        response_map,
        observed_effect=effect,
        iterations=int(config.get("random_control_iterations", 2000)),
        seed=int(config.get("random_seed", 20260809)),
    )

    gates = {
        "exact_permutation_p_lt_0_05": exact_p < 0.05,
        "prespecified_positive_orientation": effect > 0,
        "global_pc1_adjusted_hc3_p_lt_0_05": adjusted_p < 0.05,
        "overlap_preserving_empirical_p_lt_0_05": empirical_p < 0.05,
        "no_leave_one_patient_out_sign_reversal": no_loo_reversal,
        "rank_score_direction_concordant": rank_effect > 0,
    }
    primary_status = "VALIDATION" if all(gates.values()) else "SENSITIVITY_OR_BOUNDARY"
    null_quantiles = null_table["effect"].quantile([0.025, 0.5, 0.975])
    result = {
        "analysis_id": config["analysis_id"],
        "n_patients": int(len(z_delta)),
        "n_RE": int(labels.sum()),
        "n_NR": int((labels == 0).sum()),
        "effect_RE_minus_NR": effect,
        "exact_permutation_p": exact_p,
        "exact_permutation_denominator": exact_denominator,
        "global_pc1_adjusted_response_effect": adjusted_effect,
        "global_pc1_adjusted_HC3_p": adjusted_p,
        "rank_effect_RE_minus_NR": rank_effect,
        "rank_exact_p": rank_p,
        "overlap_preserving_empirical_p": empirical_p,
        "overlap_null_effect_q025": float(null_quantiles.loc[0.025]),
        "overlap_null_effect_median": float(null_quantiles.loc[0.5]),
        "overlap_null_effect_q975": float(null_quantiles.loc[0.975]),
        "gates_passed": int(sum(gates.values())),
        "gates_total": int(len(gates)),
        "status": primary_status,
        **gates,
    }

    z_delta.to_csv(output_dir / "GSE232240_patient_family_deltas.csv", index=False)
    group_totals.to_csv(output_dir / "GSE232240_pseudobulk_group_totals.csv", index=False)
    exact_null.to_csv(output_dir / "GSE232240_exact_permutation_null.csv", index=False)
    rank_delta.to_csv(output_dir / "GSE232240_patient_rank_family_deltas.csv", index=False)
    lineage_deltas.to_csv(output_dir / "GSE232240_lineage_family_deltas.csv", index=False)
    lineage_results.to_csv(output_dir / "GSE232240_lineage_family_results.csv", index=False)
    z_modules.to_csv(output_dir / "GSE232240_module_scores.csv", index=False)
    coverage.to_csv(output_dir / "GSE232240_frozen_gene_coverage.csv", index=False)
    module_results.to_csv(output_dir / "GSE232240_module_results.csv", index=False)
    loo.to_csv(output_dir / "GSE232240_leave_one_patient_out.csv", index=False)
    null_table.to_csv(output_dir / "GSE232240_overlap_preserving_null.csv.gz", index=False, compression="gzip")
    pd.DataFrame([result]).to_csv(output_dir / "GSE232240_primary_result.csv", index=False)
    with (output_dir / "GSE232240_primary_result.json").open("w", encoding="utf-8") as handle:
        json.dump(result, handle, indent=2)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
