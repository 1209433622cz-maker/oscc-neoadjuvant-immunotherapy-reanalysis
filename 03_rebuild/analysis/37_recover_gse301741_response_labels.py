#!/usr/bin/env python
"""Recover and audit patient-level GSE301741 response labels.

This script is deliberately conservative. It only promotes a response label when
it is present in an explicit source table, either extracted from the GSE301741
Seurat metadata or manually curated from the publication supplementary files.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import re
from collections import Counter, defaultdict
from pathlib import Path


PATIENT_RE = re.compile(r"\b(?:OCSCC|LSCC|OPSCC)\d+\b", re.IGNORECASE)
RESPONSE_COL_RE = re.compile(
    r"response|responder|outcome|p[_ -]?tr|path|regression|tumou?r[_ -]?regression|"
    r"m[_ -]?pr|c[_ -]?pr|non[_ -]?response|nonresponder",
    re.IGNORECASE,
)
PATIENT_COL_RE = re.compile(r"patient|subject|sample|donor|library|orig\.ident", re.IGNORECASE)
EXCLUDE_EVIDENCE_COL_RE = re.compile(
    r"definition|source|file|detail|curator|curation|note|recommended|status|provenance",
    re.IGNORECASE,
)
MANUAL_LABEL_COL_RE = re.compile(r"^response_label$|^pTR_class$|^pTR_percent$", re.IGNORECASE)


def workspace_root() -> Path:
    return Path(__file__).resolve().parents[2]


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="", encoding="utf-8-sig", errors="replace") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def extract_patient(row: dict[str, str]) -> str:
    for key, value in row.items():
        if value is None:
            continue
        if PATIENT_COL_RE.search(key) or PATIENT_RE.search(str(value)):
            match = PATIENT_RE.search(str(value))
            if match:
                return match.group(0).upper()
    return ""


def normalize_response(value: str) -> tuple[str, str, str]:
    raw = (value or "").strip()
    lower = raw.lower()
    if not raw:
        return "", "", ""
    if "pending" in lower or "unknown" in lower or "tbd" in lower:
        return "", "", ""
    if ("responders are" in lower or "defined as" in lower) and (
        "non-responder" in lower or "non-response" in lower or "non responder" in lower
    ):
        return "", "", ""

    ptr_match = re.search(r"ptr\s*[-_ ]?\s*([012])", lower)
    if ptr_match:
        pclass = f"pTR-{ptr_match.group(1)}"
        if pclass in {"pTR-1", "pTR-2"}:
            return "responder", pclass, "pTR response definition"
        return "non_responder", pclass, "pTR response definition"

    percent_match = re.search(r"([-+]?\d+(?:\.\d+)?)\s*%", lower)
    if percent_match and ("regression" in lower or "response" in lower or "ptr" in lower):
        pct = float(percent_match.group(1))
        if pct >= 10:
            inferred_class = "pTR-2" if pct > 50 else "pTR-1"
            return "responder", inferred_class, "percent tumor regression"
        return "non_responder", "pTR-0", "percent tumor regression"

    non_tokens = [
        "non-response",
        "non response",
        "non_response",
        "nonresponder",
        "non-responder",
        "non_responder",
        "no response",
        "nr",
    ]
    resp_tokens = ["responder", "response", "responders", "r"]
    if any(token in lower for token in non_tokens):
        return "non_responder", "", "explicit non-response text"
    if lower in {"r", "pr", "mpr", "cpr"} or any(token in lower for token in resp_tokens):
        return "responder", "", "explicit response text"

    return "", "", ""


def candidate_rows_from_table(path: Path, source_label: str) -> list[dict[str, str]]:
    rows = read_csv_rows(path)
    if not rows:
        return []

    headers = rows[0].keys()
    if source_label in {"manual_publication_curated_table", "publication_supplement_figure"}:
        response_cols = [col for col in headers if MANUAL_LABEL_COL_RE.search(col)]
    else:
        response_cols = [
            col
            for col in headers
            if RESPONSE_COL_RE.search(col) and not EXCLUDE_EVIDENCE_COL_RE.search(col)
        ]
    possible_patient_cols = [col for col in headers if PATIENT_COL_RE.search(col)]
    candidates: list[dict[str, str]] = []

    for row in rows:
        patient = extract_patient(row)
        if not patient:
            continue
        for col in response_cols:
            value = (row.get(col) or "").strip()
            label, pclass, rule = normalize_response(value)
            if not value:
                continue
            candidates.append(
                {
                    "patient_id": patient,
                    "response_label": label,
                    "pTR_class": pclass,
                    "source_type": source_label,
                    "source_file": str(path),
                    "evidence_column": col,
                    "evidence_value": value,
                    "normalization_rule": rule,
                    "candidate_status": "usable" if label else "needs_manual_review",
                    "patient_columns_detected": ";".join(possible_patient_cols),
                }
            )
    return candidates


def collapse_candidates(candidates: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in candidates:
        grouped[row["patient_id"]].append(row)

    collapsed: dict[str, dict[str, str]] = {}
    for patient, rows in grouped.items():
        usable = [row for row in rows if row["candidate_status"] == "usable" and row["response_label"]]
        labels = Counter(row["response_label"] for row in usable)
        ptrs = Counter(row["pTR_class"] for row in usable if row["pTR_class"])
        if len(labels) == 1:
            chosen = usable[0]
            figure_derived = chosen["source_type"] == "publication_supplement_figure"
            collapsed[patient] = {
                "response_status": "response_label_recovered_figure_derived"
                if figure_derived
                else "response_label_recovered",
                "response_label": next(iter(labels)),
                "pTR_class": ptrs.most_common(1)[0][0] if ptrs else "",
                "response_source": chosen["source_type"],
                "source_file": chosen["source_file"],
                "evidence_column": chosen["evidence_column"],
                "evidence_value": chosen["evidence_value"],
                "provenance_status": "figure_derived_requires_rds_or_table_crosscheck"
                if figure_derived
                else "single_consistent_label",
            }
        elif len(labels) > 1:
            collapsed[patient] = {
                "response_status": "conflicting_response_candidates",
                "response_label": "",
                "pTR_class": "",
                "response_source": ";".join(sorted({row["source_type"] for row in rows})),
                "source_file": ";".join(sorted({row["source_file"] for row in rows})),
                "evidence_column": ";".join(sorted({row["evidence_column"] for row in rows})),
                "evidence_value": " | ".join(sorted({row["evidence_value"] for row in rows})),
                "provenance_status": "manual_resolution_required",
            }
        else:
            collapsed[patient] = {
                "response_status": "candidate_values_need_manual_review",
                "response_label": "",
                "pTR_class": "",
                "response_source": ";".join(sorted({row["source_type"] for row in rows})),
                "source_file": ";".join(sorted({row["source_file"] for row in rows})),
                "evidence_column": ";".join(sorted({row["evidence_column"] for row in rows})),
                "evidence_value": " | ".join(sorted({row["evidence_value"] for row in rows})),
                "provenance_status": "manual_review_required",
            }
    return collapsed


def markdown_report(path: Path, patient_rows: list[dict[str, object]], candidates: list[dict[str, str]]) -> None:
    counts = Counter(str(row.get("response_status", "")) for row in patient_rows)
    recovered = [
        row
        for row in patient_rows
        if str(row.get("response_status", "")).startswith("response_label_recovered")
    ]
    lines = [
        "# GSE301741 Response Label Recovery",
        "",
        "## Decision",
        "",
    ]
    if recovered:
        lines.append(
            "Patient-level response labels were recovered for at least one patient. "
            "Use RDS/table-confirmed labels for final response-stratified claims when available; "
            "figure-derived labels should be treated as provisional until cross-checked."
        )
    else:
        lines.append(
            "No patient-level response label has been recovered from an explicit local source yet. "
            "Do not run or report GSE301741 response-stratified validation until the publication supplement or RDS metadata yields labels."
        )
    lines.extend(["", "## Response Status Counts", "", "| Status | Patients |", "|---|---:|"])
    for status, n in sorted(counts.items()):
        lines.append(f"| {status} | {n} |")
    lines.extend(
        [
            "",
            "## Source Hierarchy",
            "",
            "1. Publication supplementary table/PDF with patient-level outcome or pTR annotations.",
            "2. GSE301741 Seurat RDS metadata extracted locally from the deposited object.",
            "3. Manual curation template populated from the supplement with exact source/provenance fields.",
            "4. Figure-only labels are acceptable only as a temporary audit note, not as a final analysis source.",
            "",
            "## Files Written",
            "",
            "- `03_rebuild/validation/GSE301741_response_recovery/GSE301741_RESPONSE_LABEL_CANDIDATES.csv`",
            "- `03_rebuild/validation/GSE301741_response_recovery/GSE301741_RESPONSE_LABEL_STATUS.csv`",
            "- `03_rebuild/validation/GSE301741_response_recovery/GSE301741_VALIDATION_STRATA_WITH_RESPONSE.csv`",
        ]
    )
    if candidates:
        lines.extend(["", "## Candidate Evidence", "", f"- Candidate evidence rows detected: {len(candidates)}"])
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def merge_strata_with_response(
    strata_rows: list[dict[str, str]],
    patient_rows: list[dict[str, object]],
) -> list[dict[str, object]]:
    response_by_patient = {
        str(row["patient_id"]): row for row in patient_rows
    }
    merged_rows: list[dict[str, object]] = []
    response_fields = [
        "response_label",
        "pTR_class",
        "source_file",
        "evidence_column",
        "evidence_value",
        "provenance_status",
    ]
    for strata in strata_rows:
        patient = strata["patient_id"]
        response = response_by_patient.get(patient, {})
        merged: dict[str, object] = dict(strata)
        merged["response_status"] = response.get(
            "response_status", strata.get("response_status", "")
        )
        merged["response_source"] = response.get(
            "response_source", strata.get("response_source", "")
        )
        for field in response_fields:
            merged[field] = response.get(field, "")

        if not merged.get("response_label"):
            readiness = "response_label_not_ready"
        elif merged.get("validation_tier") == "A_same_fraction_same_library":
            readiness = "tierA_provisional_response_ready_needs_crosscheck"
        elif merged.get("validation_tier") in {
            "B_same_fraction_mixed_library",
            "C_prepost_mixed_fraction",
        }:
            readiness = "sensitivity_only_due_fraction_or_library_mismatch"
        else:
            readiness = "cross_sectional_context_only"
        merged["response_validation_readiness"] = readiness
        merged_rows.append(merged)
    return merged_rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workspace", default=str(workspace_root()))
    args = parser.parse_args()

    workspace = Path(args.workspace).resolve()
    strata_path = workspace / "03_rebuild" / "validation" / "GSE301741_stratification" / "GSE301741_VALIDATION_STRATA.csv"
    manual_path = workspace / "03_rebuild" / "config" / "GSE301741_RESPONSE_LABELS_MANUAL.csv"
    figure_path = workspace / "03_rebuild" / "validation" / "GSE301741_response_recovery" / "GSE301741_RESPONSE_LABELS_FIGURE_DERIVED.csv"
    manifest_dir = workspace / "03_rebuild" / "manifests"
    rds_summary_path = manifest_dir / "GSE301741_RDS_METADATA_RESPONSE_FIELDS.csv"
    rds_cell_meta_path = manifest_dir / "GSE301741_RDS_CELL_METADATA.csv.gz"
    out_dir = workspace / "03_rebuild" / "validation" / "GSE301741_response_recovery"

    strata_rows = read_csv_rows(strata_path)
    if not strata_rows:
        raise SystemExit(f"Missing strata table: {strata_path}")

    candidates: list[dict[str, str]] = []
    candidates.extend(candidate_rows_from_table(manual_path, "manual_publication_curated_table"))
    candidates.extend(candidate_rows_from_table(figure_path, "publication_supplement_figure"))
    candidates.extend(candidate_rows_from_table(rds_summary_path, "rds_metadata_summary"))
    candidates.extend(candidate_rows_from_table(rds_cell_meta_path, "rds_cell_metadata"))
    collapsed = collapse_candidates(candidates)

    patient_rows: list[dict[str, object]] = []
    for row in strata_rows:
        patient = row["patient_id"]
        recovered = collapsed.get(patient)
        if recovered is None:
            recovered = {
                "response_status": "pending_response_source_required",
                "response_label": "",
                "pTR_class": "",
                "response_source": "pending: publication supplementary table or GSE301741 RDS metadata",
                "source_file": "",
                "evidence_column": "",
                "evidence_value": "",
                "provenance_status": "not_recovered",
            }
        patient_rows.append(
            {
                "patient_id": patient,
                "site": row.get("site", ""),
                "validation_tier": row.get("validation_tier", ""),
                "recommended_use": row.get("recommended_use", ""),
                **recovered,
            }
        )

    candidate_fields = [
        "patient_id",
        "response_label",
        "pTR_class",
        "source_type",
        "source_file",
        "evidence_column",
        "evidence_value",
        "normalization_rule",
        "candidate_status",
        "patient_columns_detected",
    ]
    patient_fields = [
        "patient_id",
        "site",
        "validation_tier",
        "recommended_use",
        "response_status",
        "response_label",
        "pTR_class",
        "response_source",
        "source_file",
        "evidence_column",
        "evidence_value",
        "provenance_status",
    ]
    write_csv(out_dir / "GSE301741_RESPONSE_LABEL_CANDIDATES.csv", candidates, candidate_fields)
    write_csv(out_dir / "GSE301741_RESPONSE_LABEL_STATUS.csv", patient_rows, patient_fields)
    merged_rows = merge_strata_with_response(strata_rows, patient_rows)
    merged_fields = list(strata_rows[0])
    merged_fields.extend(
        field
        for field in [
            "response_label",
            "pTR_class",
            "source_file",
            "evidence_column",
            "evidence_value",
            "provenance_status",
            "response_validation_readiness",
        ]
        if field not in merged_fields
    )
    write_csv(
        out_dir / "GSE301741_VALIDATION_STRATA_WITH_RESPONSE.csv",
        merged_rows,
        merged_fields,
    )
    markdown_report(out_dir / "GSE301741_RESPONSE_LABEL_RECOVERY_REPORT.md", patient_rows, candidates)

    recovered_n = sum(
        1 for row in patient_rows if str(row["response_status"]).startswith("response_label_recovered")
    )
    print(f"Patients audited: {len(patient_rows)}")
    print(f"Response labels recovered: {recovered_n}")
    print(f"Candidate evidence rows: {len(candidates)}")
    print(f"Report: {out_dir / 'GSE301741_RESPONSE_LABEL_RECOVERY_REPORT.md'}")


if __name__ == "__main__":
    main()
