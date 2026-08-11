from __future__ import annotations

import csv
import math
import os
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path


WORKSPACE = Path(os.environ.get("GSE200996_WORKSPACE", Path.cwd())).resolve()
MANIFEST_DIR = WORKSPACE / "03_rebuild" / "manifests"
GSE281729_ANNOTATION = (
    WORKSPACE
    / "03_rebuild"
    / "validation"
    / "GSE281729_bulk_module_validation"
    / "GSE281729_EMBEDDED_SAMPLE_ANNOTATION.csv"
)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8-sig") as fh:
        return list(csv.DictReader(fh))


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    fieldnames: list[str] = []
    seen = set()
    for row in rows:
        for key in row:
            if key not in seen:
                seen.add(key)
                fieldnames.append(key)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def join_unique(values: list[str]) -> str:
    return ";".join(sorted({v for v in values if v and v != "NA" and str(v).lower() != "nan"}))


def norm_float(value: str) -> str:
    if value is None:
        return ""
    text = str(value).strip()
    if not text or text.lower() == "nan":
        return ""
    try:
        number = float(text)
    except ValueError:
        return ""
    if math.isnan(number):
        return ""
    if number.is_integer():
        return str(int(number))
    return f"{number:.4f}".rstrip("0").rstrip(".")


def main() -> int:
    v1_path = MANIFEST_DIR / "RESPONSE_METADATA_FREEZE.csv"
    if not v1_path.exists():
        raise FileNotFoundError(v1_path)
    if not GSE281729_ANNOTATION.exists():
        raise FileNotFoundError(GSE281729_ANNOTATION)

    rows = read_csv(v1_path)
    embedded = read_csv(GSE281729_ANNOTATION)
    by_patient: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in embedded:
        patient = row.get("patient_id", "")
        if patient:
            by_patient[patient].append(row)

    updates = {}
    for patient, patient_rows in by_patient.items():
        ordinals = join_unique([row.get("response_harmonized_ordinal", "") for row in patient_rows])
        binaries = join_unique([row.get("response_binary", "") for row in patient_rows])
        raw = join_unique([row.get("primary_path_response_raw", "") for row in patient_rows])
        primary_pct = join_unique([norm_float(row.get("primary_response_percent", "")) for row in patient_rows])
        overall_pct = join_unique([norm_float(row.get("overall_response_percent", "")) for row in patient_rows])
        overall_path = join_unique([row.get("overall_path_response", "") for row in patient_rows])
        hpv = join_unique([row.get("hpv", "") for row in patient_rows])
        second_drug = join_unique([row.get("second_drug", "") for row in patient_rows])
        timepoints = join_unique([row.get("timepoint", "") for row in patient_rows])
        if ordinals:
            status = "category_and_percent_from_processed_expression_annotation"
        else:
            status = "pending_response_source_required"
        updates[patient] = {
            "response_original": raw,
            "response_harmonized_ordinal": ordinals,
            "response_binary": binaries,
            "path_response_percent": primary_pct,
            "overall_response_percent": overall_pct,
            "overall_path_response": overall_path,
            "response_status": status,
            "response_source": "GSE281729 processed expression embedded annotation rows: Primary Path Response and %Primary Response",
            "response_notes": "Recovered from processed-expression embedded clinical annotation; response-adaptive timing cohort.",
            "hpv_status": hpv,
            "second_drug": second_drug,
            "timepoints": timepoints,
        }

    out = []
    for row in rows:
        new = dict(row)
        if new.get("accession") == "GSE281729" and new.get("patient_id") in updates:
            new.update(updates[new["patient_id"]])
        else:
            new.setdefault("overall_response_percent", "")
            new.setdefault("overall_path_response", "")
            new.setdefault("hpv_status", "")
            new.setdefault("second_drug", "")
        out.append(new)

    out_path = MANIFEST_DIR / "RESPONSE_METADATA_FREEZE.csv"
    write_csv(out_path, out)

    audit = []
    for accession in sorted({row.get("accession", "") for row in out}):
        subset = [row for row in out if row.get("accession") == accession]
        status_counts = Counter(row.get("response_status", "") for row in subset)
        ordinal_counts = Counter(row.get("response_harmonized_ordinal", "") for row in subset)
        for status, count in sorted(status_counts.items()):
            audit.append({"accession": accession, "field": "response_status", "value": status, "count": count})
        for value, count in sorted(ordinal_counts.items()):
            audit.append({"accession": accession, "field": "response_harmonized_ordinal", "value": value, "count": count})
    audit_path = MANIFEST_DIR / "RESPONSE_METADATA_FREEZE_AUDIT.csv"
    write_csv(audit_path, audit)

    md = [
        "# Response Metadata Freeze",
        "",
        f"Created: {datetime.now().astimezone().isoformat(timespec='seconds')}",
        "",
        "## Decision",
        "",
        "GSE281729 response metadata is upgraded using the processed expression file's embedded clinical annotation rows.",
        "GSE301741 and GSE179730 response labels remain pending until authoritative clinical metadata are recovered.",
        "",
        "## Cohort Status",
        "",
        "| Cohort | Patients | Confirmed response patients | Pending patients |",
        "|---|---:|---:|---:|",
    ]
    for accession in sorted({row.get("accession", "") for row in out}):
        subset = [row for row in out if row.get("accession") == accession]
        confirmed = [
            row
            for row in subset
            if row.get("response_harmonized_ordinal") in {"Low", "Medium", "High"}
            and "pending" not in row.get("response_status", "").lower()
        ]
        pending = [row for row in subset if "pending" in row.get("response_status", "").lower()]
        md.append(f"| {accession} | {len(subset)} | {len(confirmed)} | {len(pending)} |")
    md.extend(
        [
            "",
            "## GSE281729 Upgrade",
            "",
            "- Patients confirmed from GEO surgery titles: 31.",
            "- Patients confirmed from processed-expression embedded clinical annotation: 37.",
            "- Added `path_response_percent`, `overall_response_percent`, `overall_path_response`, `hpv_status` and `second_drug` where available.",
            "- This upgrades GSE281729 to the first response-ready external validation cohort.",
            "",
            "## Files Written",
            "",
            "- `03_rebuild/manifests/RESPONSE_METADATA_FREEZE.csv`",
            "- `03_rebuild/manifests/RESPONSE_METADATA_FREEZE_AUDIT.csv`",
            "- `03_rebuild/manifests/RESPONSE_METADATA_FREEZE_REPORT.md`",
        ]
    )
    report_path = MANIFEST_DIR / "RESPONSE_METADATA_FREEZE_REPORT.md"
    report_path.write_text("\n".join(md) + "\n", encoding="utf-8")
    print("\n".join(md))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

