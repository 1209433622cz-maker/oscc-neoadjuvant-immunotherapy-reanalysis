from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


COMPARE_ROOTS = [
    Path("03_rebuild/results"),
    Path("03_rebuild/validation"),
    Path("03_rebuild/manifests"),
    Path("03_rebuild/figures/submission/source_data"),
    Path("03_rebuild/figures/external_validation/source_data"),
]

REQUIRED_OUTPUTS = [
    "03_rebuild/results/data_audit/paired_patient_response_summary.csv",
    "03_rebuild/results/pre_baseline/baseline_pre_composition_trend_respOrd_limma_logit.csv",
    "03_rebuild/results/pre_baseline/Baseline_T_cell_GSEA_Hallmark.csv",
    "03_rebuild/results/pre_baseline/Baseline_Myeloid_GSEA_Hallmark.csv",
    "03_rebuild/results/dynamic_paired/Fig4A_composition_delta_logit_limma_respOrd_trend.csv",
    "03_rebuild/results/dynamic_paired/Fig4B_T_cell_GSEA_Hallmark.csv",
    "03_rebuild/results/dynamic_paired/Fig4B_Myeloid_GSEA_Hallmark.csv",
    "03_rebuild/results/sensitivity_exact_permutation/ABUNDANCE_EXACT_PERMUTATION_RESULTS.csv",
    "03_rebuild/results/sensitivity_cohort_adjusted_pseudobulk/KEY_PATHWAY_MODEL_COMPARISON.csv",
    "03_rebuild/results/external_validation/GSE123813_paired_delta_stats.csv",
    "03_rebuild/results/external_tcr_validation/GSE123813_TCR_paired_delta_stats.csv",
    "03_rebuild/validation/GSE281729_bulk_module_validation/GSE281729_ROBUST_RESPONSE_MODELS.csv",
    "03_rebuild/validation/GSE179730_bulk_treatment_direction/GSE179730_LOCKED_MODULE_RESPONSE_EXACT.csv",
    "03_rebuild/validation/GSE301741_lineage_aware_validation/GSE301741_LINEAGE_AWARE_RESPONSE_TESTS.csv",
    "03_rebuild/validation/locked_family_robustness/LOCKED_FAMILY_TESTS.csv",
]

VOLATILE_COLUMN = re.compile(
    r"(generated|created|updated|timestamp|run_time|mtime|modified)",
    flags=re.IGNORECASE,
)

KEY_CANDIDATES = [
    "pathway",
    "gene",
    "gene_symbol",
    "celltype",
    "cell_type",
    "signature",
    "module",
    "model",
    "contrast",
    "pair_id",
    "patient_id",
    "sample_id",
    "geo_accession",
    "accession",
    "timepoint",
    "cluster",
    "lineage",
    "target_lineage",
    "metric",
    "test",
]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest().upper()


def relative_or_absolute(path: Path, root: Path) -> str:
    try:
        return path.relative_to(root).as_posix()
    except ValueError:
        return path.as_posix()


def provenance_files(workspace: Path) -> tuple[list[Path], list[Path]]:
    code_files: list[Path] = []
    allowed_suffixes = {".py", ".r", ".ps1", ".json", ".csv", ".md"}
    for relative in ("03_rebuild/analysis", "03_rebuild/config"):
        root = workspace / relative
        if root.exists():
            code_files.extend(
                path
                for path in root.rglob("*")
                if path.is_file()
                and "__pycache__" not in path.parts
                and path.suffix.lower() in allowed_suffixes
            )

    input_patterns = [
        "00_raw_data/GSE200996_RAW/*.h5",
        "00_raw_data/GSE200996_metadata/GSE200996_CD45.tumor.single.cell.meta.data.txt.gz",
        "00_raw_data/GSE123813_validation/*",
        "00_raw_data/external_validation/GSE281729/*",
        "00_raw_data/external_validation/GSE179730/*",
        "00_raw_data/external_validation/GSE301741/GSE301741_RAW.tar",
        "00_raw_data/external_validation/GSE301741/RAW_extracted/**/*.h5",
        "00_raw_data/external_validation/GSE301741/*series_matrix*",
        "00_raw_data/external_validation/GSE301741/*family.soft*",
        "02_references/external_supplements/GSE301741/*",
    ]
    input_files: list[Path] = []
    for pattern in input_patterns:
        input_files.extend(
            path for path in workspace.glob(pattern) if path.is_file()
        )
    return sorted(set(code_files)), sorted(set(input_files))


def build_manifest(
    workspace: Path,
    output: Path,
    deep_input_hash: bool,
) -> None:
    code_files, input_files = provenance_files(workspace)
    rows: list[dict[str, object]] = []
    for category, files in (("code_or_config", code_files), ("input", input_files)):
        for path in files:
            stat = path.stat()
            do_hash = category == "code_or_config" or deep_input_hash
            rows.append(
                {
                    "category": category,
                    "relative_path": relative_or_absolute(path, workspace),
                    "size_bytes": stat.st_size,
                    "modified_utc": datetime.fromtimestamp(
                        stat.st_mtime, tz=timezone.utc
                    ).isoformat(),
                    "sha256": sha256_file(path) if do_hash else "",
                    "hash_status": "computed" if do_hash else "deferred",
                }
            )
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(
        f"Provenance manifest: {output} "
        f"({len(code_files)} code/config; {len(input_files)} inputs)"
    )


def normalize_text(value: object, baseline: Path, clean: Path) -> str:
    if pd.isna(value):
        return "<NA>"
    text = str(value).replace("\\", "/")
    for root in sorted((baseline, clean), key=lambda item: len(item.as_posix()), reverse=True):
        root_text = root.as_posix()
        text = re.sub(re.escape(root_text), "<WORKSPACE>", text, flags=re.I)
    return text


def numeric_series(series: pd.Series) -> tuple[bool, np.ndarray]:
    text = series.astype("string").fillna("")
    nonempty = text.str.strip() != ""
    if not nonempty.any():
        return False, np.array([])
    numeric = pd.to_numeric(text.where(nonempty), errors="coerce")
    if numeric[nonempty].notna().all():
        return True, numeric.to_numpy(dtype=float)
    return False, np.array([])


def align_rows(
    left: pd.DataFrame,
    right: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, str]:
    candidates = [column for column in KEY_CANDIDATES if column in left.columns]
    for width in range(1, len(candidates) + 1):
        keys = candidates[:width]
        left_keys = left[keys].astype(str).agg("\x1f".join, axis=1)
        right_keys = right[keys].astype(str).agg("\x1f".join, axis=1)
        if (
            left_keys.is_unique
            and right_keys.is_unique
            and Counter(left_keys) == Counter(right_keys)
        ):
            left = left.assign(__clean_key=left_keys).sort_values(
                "__clean_key", kind="stable"
            )
            right = right.assign(__clean_key=right_keys).sort_values(
                "__clean_key", kind="stable"
            )
            return (
                left.drop(columns="__clean_key").reset_index(drop=True),
                right.drop(columns="__clean_key").reset_index(drop=True),
                ",".join(keys),
            )
    return left.reset_index(drop=True), right.reset_index(drop=True), "row_order"


def compare_csv(
    baseline_path: Path,
    clean_path: Path,
    baseline: Path,
    clean: Path,
) -> dict[str, object]:
    relative = relative_or_absolute(clean_path, clean)
    if not baseline_path.exists():
        return {
            "relative_path": relative,
            "status": "MISSING_BASELINE",
            "baseline_rows": "",
            "clean_rows": "",
            "max_abs_numeric_diff": "",
            "detail": "No corresponding frozen baseline table.",
        }

    try:
        left = pd.read_csv(baseline_path, dtype=str, keep_default_na=False)
        right = pd.read_csv(clean_path, dtype=str, keep_default_na=False)
    except Exception as exc:
        return {
            "relative_path": relative,
            "status": "READ_ERROR",
            "baseline_rows": "",
            "clean_rows": "",
            "max_abs_numeric_diff": "",
            "detail": str(exc),
        }

    if list(left.columns) != list(right.columns):
        return {
            "relative_path": relative,
            "status": "SCHEMA_FAIL",
            "baseline_rows": len(left),
            "clean_rows": len(right),
            "max_abs_numeric_diff": "",
            "detail": "Column names or order differ.",
        }
    if len(left) != len(right):
        return {
            "relative_path": relative,
            "status": "ROW_COUNT_FAIL",
            "baseline_rows": len(left),
            "clean_rows": len(right),
            "max_abs_numeric_diff": "",
            "detail": "Row counts differ.",
        }

    left, right, alignment = align_rows(left, right)
    failures: list[str] = []
    max_abs_diff = 0.0
    for column in left.columns:
        if VOLATILE_COLUMN.search(column):
            continue
        left_numeric, left_values = numeric_series(left[column])
        right_numeric, right_values = numeric_series(right[column])
        if left_numeric and right_numeric:
            equal = np.isclose(
                left_values,
                right_values,
                rtol=1e-8,
                atol=1e-10,
                equal_nan=True,
            )
            if not equal.all():
                finite = np.isfinite(left_values) & np.isfinite(right_values)
                if finite.any():
                    max_abs_diff = max(
                        max_abs_diff,
                        float(
                            np.max(
                                np.abs(left_values[finite] - right_values[finite])
                            )
                        ),
                    )
                failures.append(f"{column}:{int((~equal).sum())}")
        else:
            left_text = [
                normalize_text(value, baseline, clean)
                for value in left[column]
            ]
            right_text = [
                normalize_text(value, baseline, clean)
                for value in right[column]
            ]
            mismatch = sum(a != b for a, b in zip(left_text, right_text))
            if mismatch:
                failures.append(f"{column}:{mismatch}")

    return {
        "relative_path": relative,
        "status": "PASS" if not failures else "CONTENT_FAIL",
        "baseline_rows": len(left),
        "clean_rows": len(right),
        "max_abs_numeric_diff": f"{max_abs_diff:.12g}",
        "detail": (
            f"All non-volatile values agree after alignment by {alignment}."
            if not failures
            else (
                f"Alignment by {alignment}; mismatched cells by column: "
                + ", ".join(failures[:12])
            )
        ),
    }


def compare_outputs(
    baseline: Path,
    clean: Path,
    output_csv: Path,
    output_md: Path,
) -> None:
    clean_files: list[Path] = []
    for relative_root in COMPARE_ROOTS:
        root = clean / relative_root
        if root.exists():
            clean_files.extend(root.rglob("*.csv"))

    rows = [
        compare_csv(
            baseline / path.relative_to(clean),
            path,
            baseline,
            clean,
        )
        for path in sorted(set(clean_files))
    ]

    present = {row["relative_path"] for row in rows}
    for relative in REQUIRED_OUTPUTS:
        if relative not in present:
            rows.append(
                {
                    "relative_path": relative,
                    "status": "MISSING_CLEAN_REQUIRED",
                    "baseline_rows": "",
                    "clean_rows": "",
                    "max_abs_numeric_diff": "",
                    "detail": "Required clean-room output was not generated.",
                }
            )
    rows.sort(key=lambda row: str(row["relative_path"]))

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    with output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    counts = Counter(str(row["status"]) for row in rows)
    failure_count = sum(
        count for status, count in counts.items() if status != "PASS"
    )
    lines = [
        "# Clean-room result comparison",
        "",
        f"- Baseline: `{baseline}`",
        f"- Clean workspace: `{clean}`",
        f"- Tables assessed: {len(rows)}",
        f"- PASS: {counts.get('PASS', 0)}",
        f"- Non-PASS: {failure_count}",
        "",
        "## Status counts",
        "",
    ]
    lines.extend(f"- {status}: {count}" for status, count in sorted(counts.items()))
    lines.extend(["", "## Non-PASS tables", ""])
    non_pass = [row for row in rows if row["status"] != "PASS"]
    if non_pass:
        lines.extend(
            f"- `{row['relative_path']}`: {row['status']} - {row['detail']}"
            for row in non_pass
        )
    else:
        lines.append("- None.")
    output_md.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"Comparison CSV: {output_csv}")
    print(f"Comparison Markdown: {output_md}")
    print(json.dumps(dict(counts), ensure_ascii=True, sort_keys=True))
    if failure_count:
        raise SystemExit(2)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)

    manifest = subparsers.add_parser("manifest")
    manifest.add_argument("--workspace", type=Path, required=True)
    manifest.add_argument("--output", type=Path, required=True)
    manifest.add_argument("--deep-input-hash", action="store_true")

    compare = subparsers.add_parser("compare")
    compare.add_argument("--baseline", type=Path, required=True)
    compare.add_argument("--clean", type=Path, required=True)
    compare.add_argument("--output-csv", type=Path, required=True)
    compare.add_argument("--output-md", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.command == "manifest":
        build_manifest(
            args.workspace.resolve(),
            args.output.resolve(),
            args.deep_input_hash,
        )
        return 0
    compare_outputs(
        args.baseline.resolve(),
        args.clean.resolve(),
        args.output_csv.resolve(),
        args.output_md.resolve(),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
