#!/usr/bin/env python3
"""Unblind and compare GSE232240 clean-room outputs with the frozen baseline."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--clean-dir", required=True, type=Path)
    parser.add_argument("--baseline-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--absolute-tolerance", type=float, default=1e-10)
    parser.add_argument("--relative-tolerance", type=float, default=1e-8)
    return parser.parse_args()


def canonical_hash(values: np.ndarray, decimals: int = 12) -> str:
    normalized = np.sort(np.round(np.asarray(values, dtype=float), decimals=decimals))
    text = "\n".join(f"{value:.{decimals}f}" for value in normalized)
    return hashlib.sha256(text.encode("ascii")).hexdigest()


def main() -> int:
    args = parse_args()
    clean_dir = args.clean_dir.resolve()
    baseline_dir = args.baseline_dir.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    checks: list[dict[str, str]] = []

    def check(check_id: str, condition: bool, evidence: object) -> None:
        checks.append({"check_id": check_id, "status": "PASS" if condition else "FAIL", "evidence": str(evidence)})

    def close(left: np.ndarray, right: np.ndarray, atol: float | None = None, rtol: float | None = None) -> bool:
        return bool(
            np.allclose(
                np.asarray(left, dtype=float),
                np.asarray(right, dtype=float),
                atol=args.absolute_tolerance if atol is None else atol,
                rtol=args.relative_tolerance if rtol is None else rtol,
                equal_nan=True,
            )
        )

    table_pairs = {
        "metadata_eligibility": ("cleanroom_metadata_eligibility.csv", "GSE232240_metadata_eligibility.csv"),
        "pseudobulk_group_totals": ("cleanroom_pseudobulk_group_totals.csv", "GSE232240_pseudobulk_group_totals.csv"),
        "patient_family_deltas": ("cleanroom_patient_family_deltas.csv", "GSE232240_patient_family_deltas.csv"),
        "patient_rank_family_deltas": ("cleanroom_patient_rank_family_deltas.csv", "GSE232240_patient_rank_family_deltas.csv"),
        "lineage_family_deltas": ("cleanroom_lineage_family_deltas.csv", "GSE232240_lineage_family_deltas.csv"),
        "lineage_family_results": ("cleanroom_lineage_family_results.csv", "GSE232240_lineage_family_results.csv"),
        "module_scores": ("cleanroom_module_scores.csv", "GSE232240_module_scores.csv"),
        "frozen_gene_coverage": ("cleanroom_frozen_gene_coverage.csv", "GSE232240_frozen_gene_coverage.csv"),
        "module_results": ("cleanroom_module_results.csv", "GSE232240_module_results.csv"),
        "leave_one_patient_out": ("cleanroom_leave_one_patient_out.csv", "GSE232240_leave_one_patient_out.csv"),
    }
    tables: dict[str, tuple[pd.DataFrame, pd.DataFrame]] = {}
    for name, (clean_name, baseline_name) in table_pairs.items():
        clean_path = clean_dir / clean_name
        baseline_path = baseline_dir / baseline_name
        check(f"{name}_files_exist", clean_path.is_file() and baseline_path.is_file(), f"{clean_path}; {baseline_path}")
        if clean_path.is_file() and baseline_path.is_file():
            clean = pd.read_csv(clean_path)
            baseline = pd.read_csv(baseline_path)
            tables[name] = (clean, baseline)
            check(f"{name}_shape", clean.shape == baseline.shape, f"clean={clean.shape}; baseline={baseline.shape}")
            check(f"{name}_columns", list(clean.columns) == list(baseline.columns), f"clean={list(clean.columns)}")

    key_map = {
        "metadata_eligibility": ["patient", "response", "lineage"],
        "pseudobulk_group_totals": ["patient", "timepoint", "lineage"],
        "patient_family_deltas": ["patient"],
        "patient_rank_family_deltas": ["patient"],
        "lineage_family_deltas": ["lineage", "patient"],
        "lineage_family_results": ["lineage"],
        "module_scores": ["signature", "patient", "timepoint"],
        "frozen_gene_coverage": ["signature"],
        "module_results": ["signature"],
        "leave_one_patient_out": ["excluded_patient"],
    }
    pc_sign = 1.0
    for name, (clean, baseline) in tables.items():
        keys = key_map[name]
        clean = clean.sort_values(keys).reset_index(drop=True)
        baseline = baseline.sort_values(keys).reset_index(drop=True)
        check(f"{name}_keys", clean[keys].astype(str).equals(baseline[keys].astype(str)), keys)
        for column in clean.columns:
            if column in keys:
                continue
            if name == "patient_family_deltas" and column == "global_pc1_delta":
                correlation = float(np.corrcoef(clean[column].astype(float), baseline[column].astype(float))[0, 1])
                pc_sign = 1.0 if correlation >= 0 else -1.0
                check("global_pc1_absolute_correlation", abs(correlation) >= 1 - 1e-10, f"correlation={correlation}")
                check("global_pc1_sign_aligned_values", close(clean[column].astype(float) * pc_sign, baseline[column].astype(float)), f"sign={pc_sign}")
            elif pd.api.types.is_numeric_dtype(clean[column]) and pd.api.types.is_numeric_dtype(baseline[column]):
                difference = float(np.nanmax(np.abs(clean[column].astype(float) - baseline[column].astype(float))))
                check(f"{name}_{column}", close(clean[column], baseline[column]), f"max_abs_diff={difference}")
            else:
                clean_values = clean[column].astype(str).str.lower().tolist()
                baseline_values = baseline[column].astype(str).str.lower().tolist()
                check(f"{name}_{column}", clean_values == baseline_values, f"rows={len(clean_values)}")

    clean_exact = pd.read_csv(clean_dir / "cleanroom_exact_permutation_null.csv")
    baseline_exact = pd.read_csv(baseline_dir / "GSE232240_exact_permutation_null.csv")
    check("exact_null_rows_3003", len(clean_exact) == len(baseline_exact) == 3003, f"clean={len(clean_exact)}; baseline={len(baseline_exact)}")
    check("exact_null_assignment_ids", clean_exact["assignment_id"].equals(baseline_exact["assignment_id"]), "ordered assignment IDs")
    check("exact_null_patient_assignments", clean_exact["RE_patients"].equals(baseline_exact["RE_patients"]), "ordered RE sets")
    check("exact_null_effects", close(clean_exact["effect_RE_minus_NR"], baseline_exact["effect_RE_minus_NR"]), f"max_abs_diff={np.max(np.abs(clean_exact['effect_RE_minus_NR'] - baseline_exact['effect_RE_minus_NR']))}")
    check("exact_null_extreme_flags", clean_exact["absolute_effect_at_least_observed"].equals(baseline_exact["absolute_effect_at_least_observed"]), "all flags")
    check("exact_null_observed_flag", clean_exact["is_observed_assignment"].equals(baseline_exact["is_observed_assignment"]), "one observed assignment")
    clean_exact_hash = canonical_hash(clean_exact["effect_RE_minus_NR"].to_numpy())
    baseline_exact_hash = canonical_hash(baseline_exact["effect_RE_minus_NR"].to_numpy())
    check("exact_null_sorted_12dp_canonical_hash", clean_exact_hash == baseline_exact_hash, f"clean={clean_exact_hash}; baseline={baseline_exact_hash}")

    clean_matched = pd.read_csv(clean_dir / "cleanroom_overlap_preserving_null.csv.gz")
    baseline_matched = pd.read_csv(baseline_dir / "GSE232240_overlap_preserving_null.csv.gz")
    check("matched_null_rows_2000", len(clean_matched) == len(baseline_matched) == 2000, f"clean={len(clean_matched)}; baseline={len(baseline_matched)}")
    check("matched_null_iterations", clean_matched["iteration"].equals(baseline_matched["iteration"]), "ordered iterations")
    check("matched_null_effects", close(clean_matched["effect"], baseline_matched["effect"]), f"max_abs_diff={np.max(np.abs(clean_matched['effect'] - baseline_matched['effect']))}")
    clean_matched_hash = canonical_hash(clean_matched["effect"].to_numpy())
    baseline_matched_hash = canonical_hash(baseline_matched["effect"].to_numpy())
    check("matched_null_sorted_12dp_canonical_hash", clean_matched_hash == baseline_matched_hash, f"clean={clean_matched_hash}; baseline={baseline_matched_hash}")

    clean_result = json.loads((clean_dir / "cleanroom_primary_result.json").read_text(encoding="utf-8"))
    baseline_result = json.loads((baseline_dir / "GSE232240_primary_result.json").read_text(encoding="utf-8"))
    result_keys = [
        "analysis_id", "n_patients", "n_RE", "n_NR", "effect_RE_minus_NR", "exact_permutation_p",
        "exact_permutation_denominator", "global_pc1_adjusted_response_effect", "global_pc1_adjusted_HC3_p",
        "rank_effect_RE_minus_NR", "rank_exact_p", "overlap_preserving_empirical_p", "overlap_null_effect_q025",
        "overlap_null_effect_median", "overlap_null_effect_q975", "gates_passed", "gates_total", "status",
        "exact_permutation_p_lt_0_05", "prespecified_positive_orientation", "global_pc1_adjusted_hc3_p_lt_0_05",
        "overlap_preserving_empirical_p_lt_0_05", "no_leave_one_patient_out_sign_reversal",
        "rank_score_direction_concordant",
    ]
    for key in result_keys:
        clean_value = clean_result.get(key)
        baseline_value = baseline_result.get(key)
        if isinstance(baseline_value, float):
            condition = close(np.array([clean_value]), np.array([baseline_value]))
        else:
            condition = clean_value == baseline_value
        check(f"primary_result_{key}", condition, f"clean={clean_value}; baseline={baseline_value}")

    independent_null = pd.read_csv(clean_dir / "cleanroom_independent_seed_null.csv.gz")
    independent_effects = independent_null["effect"].to_numpy(dtype=float)
    observed = abs(float(clean_result["effect_RE_minus_NR"]))
    independent_empirical = (1 + int(np.sum(np.abs(independent_effects) >= observed))) / (len(independent_effects) + 1)
    check("independent_seed_rows_2000", len(independent_null) == 2000, f"rows={len(independent_null)}")
    check("independent_seed_empirical_recomputed", abs(independent_empirical - float(clean_result["independent_seed_empirical_p"])) <= 1e-12, f"p={independent_empirical}")
    check("independent_seed_matched_null_extreme", independent_empirical < 0.05, f"p={independent_empirical}")
    baseline_quantiles = np.quantile(baseline_matched["effect"], [0.025, 0.5, 0.975])
    independent_quantiles = np.quantile(independent_effects, [0.025, 0.5, 0.975])
    check("independent_seed_null_center_stable", abs(independent_quantiles[1] - baseline_quantiles[1]) < 0.02, f"baseline={baseline_quantiles}; independent={independent_quantiles}")
    check("independent_seed_null_width_stable", abs((independent_quantiles[2] - independent_quantiles[0]) - (baseline_quantiles[2] - baseline_quantiles[0])) < 0.03, f"baseline={baseline_quantiles}; independent={independent_quantiles}")

    provenance = json.loads((clean_dir / "cleanroom_provenance.json").read_text(encoding="utf-8"))
    check("cleanroom_declares_no_baseline_read", provenance.get("baseline_results_read") is False, provenance.get("baseline_results_read"))
    check("cleanroom_input_file_count", len(provenance.get("input_hashes", {})) == 4, sorted(provenance.get("input_hashes", {})))

    output_csv = output_dir / "GSE232240_CLEANROOM_COMPARISON_AUDIT.csv"
    with output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["check_id", "status", "evidence"])
        writer.writeheader()
        writer.writerows(checks)
    failed = [row for row in checks if row["status"] == "FAIL"]
    output_md = output_dir / "GSE232240_CLEANROOM_COMPARISON_AUDIT.md"
    output_md.write_text(
        "\n".join(
            [
                "# GSE232240 independent clean-room comparison",
                "",
                f"- Checks: {len(checks)}",
                f"- Passed: {len(checks) - len(failed)}",
                f"- Failed: {len(failed)}",
                f"- Decision: {'PASS' if not failed else 'FAIL'}",
                f"- PC1 sign alignment applied: {pc_sign}",
                f"- Exact-null sorted 12-decimal canonical SHA256: `{clean_exact_hash}`",
                f"- Frozen-seed matched-null sorted 12-decimal canonical SHA256: `{clean_matched_hash}`",
                "",
                "## Failed checks",
                "",
                *([f"- `{row['check_id']}`: {row['evidence']}" for row in failed] or ["None."]),
                "",
                "## Interpretation",
                "",
                "A PASS establishes implementation-level reproducibility from isolated raw inputs. It does not change the frozen scientific classification or convert boundary evidence into validation.",
                "",
            ]
        ),
        encoding="utf-8",
    )
    print(f"Clean-room comparison: {len(checks)} checks; {len(failed)} failed")
    print(output_md)
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
