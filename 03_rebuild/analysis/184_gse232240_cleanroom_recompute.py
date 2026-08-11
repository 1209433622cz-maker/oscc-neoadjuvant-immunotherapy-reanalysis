#!/usr/bin/env python3
"""Independent clean-room recomputation of the frozen GSE232240 analysis.

This implementation reads only staged raw inputs and the frozen configuration.
It deliberately does not import or call the primary implementation and does not
read any baseline result until the separate comparator is executed.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import itertools
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from scipy.stats import norm, rankdata


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--independent-seed", type=int, default=20260810)
    return parser.parse_args()


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_csv(path: Path, rows: list[dict], fieldnames: list[str] | None = None) -> None:
    if not rows and fieldnames is None:
        raise ValueError(f"Cannot infer columns for empty table: {path}")
    columns = fieldnames or list(rows[0])
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def read_inputs(input_dir: Path) -> tuple[dict, list[dict], list[dict]]:
    config_path = input_dir / "gse232240_validation.json"
    module_path = input_dir / "locked_family_gene_set_manifest.csv"
    metadata_path = input_dir / "GSM7324295_Meta_data_IMCISION.txt.gz"
    config = json.loads(config_path.read_text(encoding="utf-8"))
    if file_sha256(metadata_path) != config["metadata_sha256"]:
        raise ValueError("Staged metadata hash differs from the frozen configuration")
    if file_sha256(module_path) != config["module_manifest_sha256"]:
        raise ValueError("Staged module-manifest hash differs from the frozen configuration")
    with module_path.open(newline="", encoding="utf-8-sig") as handle:
        modules = list(csv.DictReader(handle))
    if len(modules) != int(config["expected_module_count"]):
        raise ValueError(f"Expected {config['expected_module_count']} modules; found {len(modules)}")
    with gzip.open(metadata_path, "rt", encoding="utf-8-sig", newline="") as handle:
        metadata = list(csv.DictReader(handle, delimiter="\t"))
    return config, modules, metadata


def audit_metadata(metadata: list[dict], config: dict) -> tuple[dict, list[str], list[dict], dict]:
    required = set(config["required_metadata_columns"])
    if not metadata or not required.issubset(metadata[0]):
        raise ValueError(f"Missing metadata columns: {sorted(required - set(metadata[0] if metadata else []))}")
    cell_ids = [row["cell_id"] for row in metadata]
    if len(cell_ids) != len(set(cell_ids)):
        raise ValueError("Metadata cell IDs are not unique")
    responses: dict[str, set[str]] = defaultdict(set)
    cell_meta: dict[str, tuple[str, str, str | None]] = {}
    cell_counts: Counter = Counter()
    lineage_map = config["lineage_map"]
    excluded = set(config["excluded_patients"])
    for row in metadata:
        patient = row["patient"]
        responses[patient].add(row["response"])
        lineage = lineage_map.get(row["cell_type"])
        cell_meta[row["cell_id"]] = (patient, row["timepoint"], lineage)
        if patient not in excluded and lineage is not None:
            cell_counts[(patient, row["response"], lineage, row["timepoint"])] += 1
    inconsistent = [patient for patient, values in responses.items() if len(values) != 1]
    if inconsistent:
        raise ValueError(f"Inconsistent response labels: {inconsistent}")

    minimum = int(config["minimum_cells_per_patient_timepoint_lineage"])
    lineages = sorted(set(lineage_map.values()))
    eligible_by_lineage: dict[str, set[str]] = {lineage: set() for lineage in lineages}
    eligibility_rows = []
    combinations_present = sorted({key[:3] for key in cell_counts})
    for patient, response, lineage in combinations_present:
        pre = int(cell_counts.get((patient, response, lineage, "pre"), 0))
        post = int(cell_counts.get((patient, response, lineage, "post"), 0))
        eligible = pre >= minimum and post >= minimum
        if eligible:
            eligible_by_lineage[lineage].add(patient)
        eligibility_rows.append(
            {
                "patient": patient,
                "response": response,
                "lineage": lineage,
                "post": post,
                "pre": pre,
                "eligible_lineage_pair": eligible,
                "eligible_primary_both_lineages": False,
            }
        )
    primary_patients = sorted(set.intersection(*(eligible_by_lineage[lineage] for lineage in lineages)))
    for row in eligibility_rows:
        row["eligible_primary_both_lineages"] = row["patient"] in primary_patients
    response_map = {patient: next(iter(responses[patient])) for patient in primary_patients}
    response_counts = Counter(response_map.values())
    audit = {
        "metadata_cells": len(metadata),
        "metadata_patients": len(responses),
        "excluded_patients": config["excluded_patients"],
        "eligible_T_cell_patients": len(eligible_by_lineage.get("T_cell", set())),
        "eligible_Myeloid_patients": len(eligible_by_lineage.get("Myeloid", set())),
        "primary_both_lineages_patients": len(primary_patients),
        "primary_response_counts": dict(response_counts),
        "exact_permutation_denominator": math.comb(len(primary_patients), response_counts.get("RE", 0)),
        "primary_patients": primary_patients,
        "status": "cleanroom_metadata_gate_passed",
    }
    return cell_meta, primary_patients, eligibility_rows, {"audit": audit, "response_map": response_map}


def aggregate_counts(
    count_path: Path,
    cell_meta: dict,
    primary_patients: list[str],
) -> tuple[list[str], np.ndarray, list[tuple[str, str, str]], np.ndarray, np.ndarray]:
    group_keys = [
        (patient, timepoint, lineage)
        for patient in primary_patients
        for timepoint in ("pre", "post")
        for lineage in ("T_cell", "Myeloid")
    ]
    group_lookup = {key: index for index, key in enumerate(group_keys)}
    with gzip.open(count_path, "rt", encoding="utf-8-sig", errors="strict") as handle:
        header = handle.readline().rstrip("\r\n").split("\t")
        if header[0] in cell_meta:
            matrix_cells = header
        elif len(header) > 1 and header[1] in cell_meta:
            matrix_cells = header[1:]
        else:
            raise ValueError("Count header cannot be mapped to staged metadata")
        if len(matrix_cells) != len(set(matrix_cells)):
            raise ValueError("Count header contains duplicate cell IDs")
        missing = [cell for cell in matrix_cells if cell not in cell_meta]
        if missing:
            raise ValueError(f"Count columns without metadata: {len(missing)}; first={missing[0]}")
        column_groups = np.full(len(matrix_cells), -1, dtype=np.int32)
        for column, cell in enumerate(matrix_cells):
            key = cell_meta[cell]
            if key in group_lookup:
                column_groups[column] = group_lookup[key]
        selected_columns = np.flatnonzero(column_groups >= 0)
        selected_groups = column_groups[selected_columns]
        sort_order = np.argsort(selected_groups, kind="stable")
        ordered_columns = selected_columns[sort_order]
        ordered_groups = selected_groups[sort_order]
        segment_starts = np.r_[0, np.flatnonzero(np.diff(ordered_groups)) + 1]
        segment_groups = ordered_groups[segment_starts]
        group_cell_counts = np.bincount(selected_groups, minlength=len(group_keys)).astype(np.int64)

        by_gene: dict[str, np.ndarray] = {}
        library_sizes = np.zeros(len(group_keys), dtype=np.float64)
        for line_number, line in enumerate(handle, start=2):
            gene, separator, value_text = line.rstrip("\r\n").partition("\t")
            if not separator:
                raise ValueError(f"Malformed count row at line {line_number}")
            values = np.fromstring(value_text, sep="\t", dtype=np.float64)
            if values.size != len(matrix_cells):
                raise ValueError(f"Count row {line_number}: {values.size} values; expected {len(matrix_cells)}")
            sums = np.zeros(len(group_keys), dtype=np.float64)
            sums[segment_groups] = np.add.reduceat(values[ordered_columns], segment_starts)
            symbol = gene.strip().upper()
            if symbol:
                if symbol in by_gene:
                    by_gene[symbol] += sums
                else:
                    by_gene[symbol] = sums
            library_sizes += sums
    if np.any(library_sizes <= 0):
        raise ValueError("At least one clean-room pseudobulk has a zero library size")
    genes = sorted(by_gene)
    matrix = np.vstack([by_gene[gene] for gene in genes])
    return genes, matrix, group_keys, library_sizes, group_cell_counts


def lineage_z_scores(log_cpm: np.ndarray, group_keys: list[tuple[str, str, str]]) -> np.ndarray:
    output = np.full_like(log_cpm, np.nan)
    for lineage in ("T_cell", "Myeloid"):
        columns = np.array([key[2] == lineage for key in group_keys])
        block = log_cpm[:, columns]
        means = np.sum(block, axis=1, keepdims=True) / block.shape[1]
        centered = block - means
        standard_deviations = np.sqrt(np.sum(centered * centered, axis=1, keepdims=True) / (block.shape[1] - 1))
        standard_deviations[standard_deviations == 0] = np.nan
        output[:, columns] = centered / standard_deviations
    return output


def score_modules(
    score_matrix: np.ndarray,
    genes: list[str],
    group_keys: list[tuple[str, str, str]],
    modules: list[dict],
) -> tuple[list[dict], list[dict], dict[str, list[str]]]:
    gene_lookup = {gene: index for index, gene in enumerate(genes)}
    scores = []
    coverage = []
    measured_by_module: dict[str, list[str]] = {}
    for module in modules:
        defined = [value.strip().upper() for value in module["genes_defined"].split(";") if value.strip()]
        measured = [gene for gene in defined if gene in gene_lookup]
        if not measured:
            raise ValueError(f"No measured genes for {module['signature']}")
        measured_by_module[module["signature"]] = measured
        columns = [index for index, key in enumerate(group_keys) if key[2] == module["target_lineage"]]
        indices = [gene_lookup[gene] for gene in measured]
        values = np.nanmean(score_matrix[np.ix_(indices, columns)], axis=0)
        for column, value in zip(columns, values):
            patient, timepoint, lineage = group_keys[column]
            scores.append(
                {
                    "patient": patient,
                    "timepoint": timepoint,
                    "lineage": lineage,
                    "signature": module["signature"],
                    "score": float(value),
                }
            )
        coverage.append(
            {
                "signature": module["signature"],
                "target_lineage": module["target_lineage"],
                "genes_defined": len(defined),
                "genes_measured": len(measured),
                "coverage_fraction": len(measured) / len(defined),
                "measured_gene_symbols": ";".join(measured),
            }
        )
    return scores, coverage, measured_by_module


def family_tables(module_scores: list[dict], response_map: dict[str, str], value_name: str) -> tuple[list[dict], list[dict]]:
    collected: dict[tuple[str, str], list[float]] = defaultdict(list)
    for row in module_scores:
        collected[(row["patient"], row["timepoint"])].append(float(row["score"]))
    scores = []
    for patient in sorted(response_map):
        for timepoint in ("pre", "post"):
            values = collected[(patient, timepoint)]
            if len(values) != 16:
                raise ValueError(f"Expected 16 module scores for {patient}/{timepoint}; found {len(values)}")
            scores.append({"patient": patient, "timepoint": timepoint, value_name: float(np.mean(values))})
    indexed = {(row["patient"], row["timepoint"]): row[value_name] for row in scores}
    delta_name = f"{value_name}_delta"
    deltas = [
        {
            "patient": patient,
            "response": response_map[patient],
            "post": indexed[(patient, "post")],
            "pre": indexed[(patient, "pre")],
            delta_name: indexed[(patient, "post")] - indexed[(patient, "pre")],
        }
        for patient in sorted(response_map)
    ]
    return scores, deltas


def exact_distribution(values: np.ndarray, labels: np.ndarray, patients: list[str]) -> tuple[float, float, list[dict]]:
    labels = labels.astype(bool)
    observed = float(values[labels].mean() - values[~labels].mean())
    observed_indices = tuple(np.flatnonzero(labels))
    rows = []
    extreme = 0
    for assignment_id, positive_indices in enumerate(itertools.combinations(range(len(values)), int(labels.sum())), start=1):
        mask = np.zeros(len(values), dtype=bool)
        mask[list(positive_indices)] = True
        effect = float(values[mask].mean() - values[~mask].mean())
        is_extreme = abs(effect) >= abs(observed) - 1e-12
        extreme += int(is_extreme)
        rows.append(
            {
                "assignment_id": assignment_id,
                "RE_patients": ";".join(patients[index] for index in positive_indices),
                "effect_RE_minus_NR": effect,
                "absolute_effect_at_least_observed": is_extreme,
                "is_observed_assignment": tuple(positive_indices) == observed_indices,
            }
        )
    return observed, extreme / len(rows), rows


def exact_test(values: np.ndarray, labels: np.ndarray) -> tuple[float, float, int]:
    patients = [str(index) for index in range(len(values))]
    effect, p_value, rows = exact_distribution(values, labels, patients)
    return effect, p_value, len(rows)


def global_pc_delta(z_matrix: np.ndarray, group_keys: list[tuple[str, str, str]], patients: list[str]) -> dict[str, float]:
    lineage_deltas: dict[tuple[str, str], float] = {}
    for lineage in ("T_cell", "Myeloid"):
        columns = [index for index, key in enumerate(group_keys) if key[2] == lineage]
        sample_by_gene = np.nan_to_num(z_matrix[:, columns].T, nan=0.0)
        gram = sample_by_gene @ sample_by_gene.T
        eigenvalues, eigenvectors = np.linalg.eigh(gram)
        pc1 = eigenvectors[:, -1] * math.sqrt(max(float(eigenvalues[-1]), 0.0))
        keyed = {group_keys[column][:2]: pc1[index] for index, column in enumerate(columns)}
        for patient in patients:
            lineage_deltas[(patient, lineage)] = keyed[(patient, "post")] - keyed[(patient, "pre")]
    return {
        patient: float(np.mean([lineage_deltas[(patient, "T_cell")], lineage_deltas[(patient, "Myeloid")]]))
        for patient in patients
    }


def hc3_response(values: np.ndarray, labels: np.ndarray, covariate: np.ndarray) -> tuple[float, float]:
    design = np.column_stack([np.ones(len(values)), labels.astype(float), covariate.astype(float)])
    inverse = np.linalg.inv(design.T @ design)
    beta = inverse @ design.T @ values
    residuals = values - design @ beta
    leverage = np.einsum("ij,jk,ik->i", design, inverse, design)
    adjusted_squared = (residuals / (1.0 - leverage)) ** 2
    meat = design.T @ (design * adjusted_squared[:, None])
    covariance = inverse @ meat @ inverse
    standard_error = math.sqrt(float(covariance[1, 1]))
    z_value = float(beta[1] / standard_error)
    return float(beta[1]), float(2.0 * norm.sf(abs(z_value)))


def bh_adjust(p_values: list[float]) -> np.ndarray:
    values = np.asarray(p_values, dtype=float)
    order = np.argsort(values)
    ranked = values[order]
    adjusted_ranked = np.minimum.accumulate((ranked * len(ranked) / np.arange(1, len(ranked) + 1))[::-1])[::-1]
    adjusted = np.empty_like(adjusted_ranked)
    adjusted[order] = np.minimum(adjusted_ranked, 1.0)
    return adjusted


def lineage_and_module_results(
    module_scores: list[dict], response_map: dict[str, str]
) -> tuple[list[dict], list[dict], list[dict]]:
    lineage_values: dict[tuple[str, str, str], list[float]] = defaultdict(list)
    module_values: dict[tuple[str, str, str], float] = {}
    for row in module_scores:
        lineage_values[(row["patient"], row["timepoint"], row["lineage"])].append(float(row["score"]))
        module_values[(row["patient"], row["timepoint"], row["signature"])] = float(row["score"])
    lineage_deltas = []
    lineage_results = []
    for lineage in ("Myeloid", "T_cell"):
        rows = []
        for patient in sorted(response_map):
            pre = float(np.mean(lineage_values[(patient, "pre", lineage)]))
            post = float(np.mean(lineage_values[(patient, "post", lineage)]))
            row = {
                "patient": patient,
                "response": response_map[patient],
                "post": post,
                "pre": pre,
                "lineage_family_score_delta": post - pre,
                "lineage": lineage,
            }
            rows.append(row)
            lineage_deltas.append(row)
        values = np.array([row["lineage_family_score_delta"] for row in rows])
        labels = np.array([row["response"] == "RE" for row in rows])
        effect, p_value, denominator = exact_test(values, labels)
        lineage_results.append(
            {"lineage": lineage, "n_patients": len(rows), "effect_RE_minus_NR": effect, "exact_p": p_value, "exact_denominator": denominator}
        )

    signatures = []
    for row in module_scores:
        if row["signature"] not in signatures:
            signatures.append(row["signature"])
    module_rows = []
    for signature in signatures:
        values = []
        labels = []
        for patient in sorted(response_map):
            values.append(module_values[(patient, "post", signature)] - module_values[(patient, "pre", signature)])
            labels.append(response_map[patient] == "RE")
        effect, p_value, denominator = exact_test(np.asarray(values), np.asarray(labels))
        module_rows.append(
            {"signature": signature, "effect_RE_minus_NR": effect, "exact_p": p_value, "exact_denominator": denominator}
        )
    adjusted = bh_adjust([row["exact_p"] for row in module_rows])
    for row, fdr in zip(module_rows, adjusted):
        row["BH_FDR"] = float(fdr)
    return lineage_deltas, lineage_results, module_rows


def matched_null(
    z_matrix: np.ndarray,
    log_cpm: np.ndarray,
    genes: list[str],
    group_keys: list[tuple[str, str, str]],
    modules: list[dict],
    measured_by_module: dict[str, list[str]],
    response_map: dict[str, str],
    observed_effect: float,
    iterations: int,
    seed: int,
) -> tuple[float, list[dict]]:
    gene_lookup = {gene: index for index, gene in enumerate(genes)}
    locked_unique = sorted({gene for values in measured_by_module.values() for gene in values})
    detectability = np.nanmean(log_cpm, axis=1)
    finite = np.isfinite(detectability)
    quantile_edges = np.unique(np.quantile(detectability[finite], np.linspace(0, 1, 11)))
    deciles = np.digitize(detectability, quantile_edges[1:-1], right=True)
    locked_set = set(locked_unique)
    candidates = {
        decile: [genes[index] for index in np.flatnonzero((deciles == decile) & finite) if genes[index] not in locked_set]
        for decile in range(10)
    }
    sources = {
        decile: [gene for gene in locked_unique if deciles[gene_lookup[gene]] == decile]
        for decile in range(10)
    }
    rng = np.random.default_rng(seed)
    sample_keys = sorted({key[:2] for key in group_keys})
    sample_lookup = {key: index for index, key in enumerate(sample_keys)}
    group_sample = np.array([sample_lookup[key[:2]] for key in group_keys], dtype=np.int32)
    patients = sorted(response_map)
    labels = np.array([response_map[patient] == "RE" for patient in patients])
    patient_pairs = [(sample_lookup[(patient, "pre")], sample_lookup[(patient, "post")]) for patient in patients]
    rows = []
    for iteration in range(1, iterations + 1):
        replacement: dict[str, str] = {}
        for decile in range(10):
            if sources[decile]:
                selected = rng.choice(candidates[decile], size=len(sources[decile]), replace=False)
                replacement.update(zip(sources[decile], selected))
        sample_scores = np.zeros(len(sample_keys), dtype=float)
        for module in modules:
            indices = [gene_lookup[replacement[gene]] for gene in measured_by_module[module["signature"]]]
            columns = [index for index, key in enumerate(group_keys) if key[2] == module["target_lineage"]]
            values = np.nanmean(z_matrix[np.ix_(indices, columns)], axis=0)
            for group_column, value in zip(columns, values):
                sample_scores[group_sample[group_column]] += value
        sample_scores /= len(modules)
        deltas = np.array([sample_scores[post] - sample_scores[pre] for pre, post in patient_pairs])
        effect = float(deltas[labels].mean() - deltas[~labels].mean())
        rows.append({"iteration": iteration, "effect": effect})
    empirical = (1 + sum(abs(row["effect"]) >= abs(observed_effect) for row in rows)) / (iterations + 1)
    return float(empirical), rows


def main() -> int:
    args = parse_args()
    input_dir = args.input_dir.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    config, modules, metadata = read_inputs(input_dir)
    cell_meta, patients, eligibility, metadata_state = audit_metadata(metadata, config)
    response_map = metadata_state["response_map"]
    count_path = input_dir / "GSM7324294_Count_data_IMCISION.txt.gz"
    genes, counts, group_keys, library_sizes, cell_counts = aggregate_counts(count_path, cell_meta, patients)
    group_totals = [
        {
            "patient": patient,
            "timepoint": timepoint,
            "lineage": lineage,
            "cell_count": int(cell_counts[index]),
            "library_size": float(library_sizes[index]),
        }
        for index, (patient, timepoint, lineage) in enumerate(group_keys)
    ]
    log_cpm = np.log2(counts / library_sizes[np.newaxis, :] * 1_000_000.0 + 0.5)
    z_matrix = lineage_z_scores(log_cpm, group_keys)
    module_scores, coverage, measured_by_module = score_modules(z_matrix, genes, group_keys, modules)
    _, patient_deltas = family_tables(module_scores, response_map, "family_score")
    labels = np.array([row["response"] == "RE" for row in patient_deltas])
    delta_values = np.array([row["family_score_delta"] for row in patient_deltas])
    observed_effect, exact_p, exact_rows = exact_distribution(delta_values, labels, patients)

    rank_matrix = np.column_stack(
        [rankdata(log_cpm[:, column], method="average") / len(genes) for column in range(log_cpm.shape[1])]
    )
    rank_scores, _, _ = score_modules(rank_matrix, genes, group_keys, modules)
    _, rank_deltas = family_tables(rank_scores, response_map, "rank_family_score")
    rank_values = np.array([row["rank_family_score_delta"] for row in rank_deltas])
    rank_effect, rank_p, _ = exact_test(rank_values, labels)

    pc_delta = global_pc_delta(z_matrix, group_keys, patients)
    for row in patient_deltas:
        row["global_pc1_delta"] = pc_delta[row["patient"]]
    adjusted_effect, adjusted_p = hc3_response(delta_values, labels, np.array([pc_delta[patient] for patient in patients]))

    lineage_deltas, lineage_results, module_results = lineage_and_module_results(module_scores, response_map)
    loo_rows = []
    for excluded_index, patient in enumerate(patients):
        keep = np.arange(len(patients)) != excluded_index
        effect = float(delta_values[keep][labels[keep]].mean() - delta_values[keep][~labels[keep]].mean())
        loo_rows.append({"excluded_patient": patient, "effect": effect})
    no_loo_reversal = all(np.sign(row["effect"]) == np.sign(observed_effect) for row in loo_rows)

    frozen_seed = int(config.get("random_seed", 20260809))
    iterations = int(config.get("random_control_iterations", 2000))
    matched_p, matched_rows = matched_null(
        z_matrix, log_cpm, genes, group_keys, modules, measured_by_module, response_map, observed_effect, iterations, frozen_seed
    )
    independent_p, independent_rows = matched_null(
        z_matrix,
        log_cpm,
        genes,
        group_keys,
        modules,
        measured_by_module,
        response_map,
        observed_effect,
        iterations,
        int(args.independent_seed),
    )

    gates = {
        "exact_permutation_p_lt_0_05": exact_p < 0.05,
        "prespecified_positive_orientation": observed_effect > 0,
        "global_pc1_adjusted_hc3_p_lt_0_05": adjusted_p < 0.05,
        "overlap_preserving_empirical_p_lt_0_05": matched_p < 0.05,
        "no_leave_one_patient_out_sign_reversal": no_loo_reversal,
        "rank_score_direction_concordant": rank_effect > 0,
    }
    matched_effects = np.array([row["effect"] for row in matched_rows])
    independent_effects = np.array([row["effect"] for row in independent_rows])
    result = {
        "analysis_id": config["analysis_id"],
        "implementation": "independent_cleanroom",
        "n_patients": len(patients),
        "n_RE": int(labels.sum()),
        "n_NR": int((~labels).sum()),
        "effect_RE_minus_NR": observed_effect,
        "exact_permutation_p": exact_p,
        "exact_permutation_denominator": len(exact_rows),
        "global_pc1_adjusted_response_effect": adjusted_effect,
        "global_pc1_adjusted_HC3_p": adjusted_p,
        "rank_effect_RE_minus_NR": rank_effect,
        "rank_exact_p": rank_p,
        "overlap_preserving_empirical_p": matched_p,
        "overlap_null_effect_q025": float(np.quantile(matched_effects, 0.025)),
        "overlap_null_effect_median": float(np.quantile(matched_effects, 0.5)),
        "overlap_null_effect_q975": float(np.quantile(matched_effects, 0.975)),
        "independent_seed": int(args.independent_seed),
        "independent_seed_empirical_p": independent_p,
        "independent_seed_q025": float(np.quantile(independent_effects, 0.025)),
        "independent_seed_median": float(np.quantile(independent_effects, 0.5)),
        "independent_seed_q975": float(np.quantile(independent_effects, 0.975)),
        "gates_passed": int(sum(gates.values())),
        "gates_total": len(gates),
        "status": "VALIDATION" if all(gates.values()) else "SENSITIVITY_OR_BOUNDARY",
        **gates,
    }

    write_csv(output_dir / "cleanroom_metadata_eligibility.csv", eligibility)
    write_csv(output_dir / "cleanroom_pseudobulk_group_totals.csv", group_totals)
    write_csv(output_dir / "cleanroom_patient_family_deltas.csv", patient_deltas)
    write_csv(output_dir / "cleanroom_patient_rank_family_deltas.csv", rank_deltas)
    write_csv(output_dir / "cleanroom_exact_permutation_null.csv", exact_rows)
    write_csv(output_dir / "cleanroom_lineage_family_deltas.csv", lineage_deltas)
    write_csv(output_dir / "cleanroom_lineage_family_results.csv", lineage_results)
    write_csv(output_dir / "cleanroom_module_scores.csv", module_scores)
    write_csv(output_dir / "cleanroom_frozen_gene_coverage.csv", coverage)
    write_csv(output_dir / "cleanroom_module_results.csv", module_results)
    write_csv(output_dir / "cleanroom_leave_one_patient_out.csv", loo_rows)
    with gzip.open(output_dir / "cleanroom_overlap_preserving_null.csv.gz", "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["iteration", "effect"])
        writer.writeheader()
        writer.writerows(matched_rows)
    with gzip.open(output_dir / "cleanroom_independent_seed_null.csv.gz", "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["iteration", "effect"])
        writer.writeheader()
        writer.writerows(independent_rows)
    write_csv(output_dir / "cleanroom_primary_result.csv", [result])
    (output_dir / "cleanroom_primary_result.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    (output_dir / "cleanroom_metadata_audit.json").write_text(json.dumps(metadata_state["audit"], indent=2), encoding="utf-8")
    provenance = {
        "implementation": "independent_cleanroom",
        "input_dir": str(input_dir),
        "input_hashes": {path.name: file_sha256(path) for path in sorted(input_dir.iterdir()) if path.is_file()},
        "script_sha256": file_sha256(Path(__file__).resolve()),
        "baseline_results_read": False,
    }
    (output_dir / "cleanroom_provenance.json").write_text(json.dumps(provenance, indent=2), encoding="utf-8")
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
