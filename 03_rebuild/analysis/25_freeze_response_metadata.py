from __future__ import annotations

import csv
import json
import os
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path


WORKSPACE = Path(os.environ.get("GSE200996_WORKSPACE", Path.cwd())).resolve()
MANIFEST_DIR = WORKSPACE / "03_rebuild" / "manifests"
OUT_DIR = WORKSPACE / "03_rebuild" / "manifests"
CONFIG_DIR = WORKSPACE / "03_rebuild" / "config"
CONFIG_DIR.mkdir(parents=True, exist_ok=True)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8-sig") as fh:
        return list(csv.DictReader(fh))


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str] | None = None) -> None:
    if fieldnames is None:
        fieldnames = []
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


def gse200996_patient_response() -> list[dict[str, object]]:
    rows = read_csv(MANIFEST_DIR / "GSE200996_LOCAL_MULTIMODAL_MANIFEST.csv")
    by_patient: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        if row.get("patient_id"):
            by_patient[row["patient_id"]].append(row)

    out: list[dict[str, object]] = []
    for patient, patient_rows in sorted(by_patient.items()):
        categories = sorted(set(r.get("path_response_category", "") for r in patient_rows if r.get("path_response_category")))
        bins = sorted(set(r.get("path_response_binary", "") for r in patient_rows if r.get("path_response_binary")))
        ords = sorted(set(r.get("path_response_ord_num", "") for r in patient_rows if r.get("path_response_ord_num")))
        arms = sorted(set(r.get("treatment_arm", "") for r in patient_rows if r.get("treatment_arm")))
        out.append(
            {
                "accession": "GSE200996",
                "patient_id": patient,
                "treatment_arm": ";".join(arms),
                "response_original": ";".join(categories),
                "response_harmonized_ordinal": ";".join(categories),
                "response_ord_num": ";".join(ords),
                "response_binary": ";".join(bins),
                "path_response_percent": "",
                "response_status": "category_available" if categories else "pending_not_in_current_response_table",
                "response_source": "03_rebuild/tables/submission/csv/S1_patient_response.csv",
                "response_notes": "Continuous pathological response percentage not yet recovered from authoritative source.",
                "n_manifest_rows": len(patient_rows),
                "modalities": ";".join(sorted(set(r.get("modality", "") for r in patient_rows if r.get("modality")))),
                "timepoints": ";".join(sorted(set(r.get("timepoint_harmonized", "") for r in patient_rows if r.get("timepoint_harmonized")))),
            }
        )
    return out


def normalize_gse281729_response(raw: str) -> tuple[str, str, str]:
    value = raw.strip()
    lower = value.lower()
    if not value:
        return "", "", "pending_no_surgery_response_in_sample_title"
    if "non" in lower:
        return "Low", "NR", "category_from_surgery_sample_title"
    if "minor" in lower:
        return "Medium", "intermediate", "category_from_surgery_sample_title"
    if "complete" in lower:
        return "High", "R", "category_from_surgery_sample_title"
    if "responder" in lower:
        return "High", "R", "category_from_surgery_sample_title"
    return "", "", "pending_unmapped_response_label"


def external_patient_response() -> list[dict[str, object]]:
    rows = read_csv(MANIFEST_DIR / "EXTERNAL_COHORT_SAMPLE_MANIFEST.csv")
    by_accession_patient: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        patient = row.get("patient_id", "")
        accession = row.get("accession", "")
        if patient and accession:
            by_accession_patient[(accession, patient)].append(row)

    out: list[dict[str, object]] = []
    for (accession, patient), patient_rows in sorted(by_accession_patient.items()):
        response_values = sorted(set(r.get("response_raw", "") for r in patient_rows if r.get("response_raw", "")))
        treatment_values = sorted(set(r.get("treatment", "") for r in patient_rows if r.get("treatment", "")))
        site_values = sorted(set(r.get("site", "") for r in patient_rows if r.get("site", "")))

        response_original = ";".join(response_values)
        response_harmonized = ""
        response_binary = ""
        response_status = "pending_response_source_required"
        response_source = "not available in current GEO sample metadata"
        notes = ""

        if accession == "GSE281729":
            response_harmonized, response_binary, response_status = normalize_gse281729_response(response_original)
            response_source = "GEO series matrix sample titles; surgery-sample label propagated to patient"
            if not response_original:
                notes = "No surgery response label detected for this patient; exclude from response-linked validation unless recovered."
        elif accession == "GSE301741":
            response_source = "pending: GSE301741 RDS metadata or publication supplementary tables"
            notes = "Do not load full RDS on current occupied/34 GB RAM state; use idle/high-memory extraction script."
        elif accession == "GSE179730":
            response_source = "pending: publication/source supplementary clinical table"
            notes = "GEO matrix provides longitudinal sample titles but not patient-level response class."

        out.append(
            {
                "accession": accession,
                "patient_id": patient,
                "site": ";".join(site_values),
                "treatment": ";".join(treatment_values),
                "response_original": response_original,
                "response_harmonized_ordinal": response_harmonized,
                "response_binary": response_binary,
                "path_response_percent": "",
                "response_status": response_status,
                "response_source": response_source,
                "response_notes": notes,
                "n_manifest_rows": len(patient_rows),
                "modalities": ";".join(sorted(set(r.get("modality", "") for r in patient_rows if r.get("modality")))),
                "timepoints": ";".join(sorted(set(r.get("timepoint_harmonized", "") for r in patient_rows if r.get("timepoint_harmonized")))),
                "geo_accessions": ";".join(sorted(set(r.get("geo_accession", "") for r in patient_rows if r.get("geo_accession")))),
            }
        )
    return out


def main() -> int:
    gse200996 = gse200996_patient_response()
    external = external_patient_response()
    all_rows = gse200996 + external

    write_csv(OUT_DIR / "RESPONSE_METADATA_FREEZE.csv", all_rows)

    audit = []
    for accession in sorted(set(str(r["accession"]) for r in all_rows)):
        subset = [r for r in all_rows if r["accession"] == accession]
        status_counts = Counter(str(r["response_status"]) for r in subset)
        ordinal_counts = Counter(str(r["response_harmonized_ordinal"]) for r in subset)
        for status, count in sorted(status_counts.items()):
            audit.append({"accession": accession, "field": "response_status", "value": status, "count": count})
        for value, count in sorted(ordinal_counts.items()):
            audit.append({"accession": accession, "field": "response_harmonized_ordinal", "value": value, "count": count})
    write_csv(OUT_DIR / "RESPONSE_METADATA_FREEZE_AUDIT.csv", audit)

    config = {
        "created": datetime.now().astimezone().isoformat(timespec="seconds"),
        "status": "response metadata freeze; category-level only for confirmed cohorts",
        "primary_response_variable": {
            "preferred": "path_response_percent",
            "fallback": "response_harmonized_ordinal",
            "fallback_rule": "Use ordered categorical response only when authoritative continuous pathological response cannot be recovered.",
            "binary_rule": "Binary response is sensitivity only and must be cohort-specific, predeclared, and clinically justified.",
        },
        "cohort_response_rules": {
            "GSE200996": {
                "current_status": "Low/Medium/High categories available for 19 tumor-response patients; continuous percent pending.",
                "allowed_primary_before_percent_recovery": False,
            },
            "GSE301741": {
                "current_status": "sample manifest ready; patient response pending RDS metadata or supplementary table extraction.",
                "allowed_primary_before_response_recovery": False,
            },
            "GSE281729": {
                "current_status": "surgery sample titles provide Non-Responder/minor Responder/Responder/Complete Responder for many patients.",
                "harmonization": "Non-Responder=Low; minor Responder=Medium; Responder or Complete Responder=High for module-level direction checks.",
                "allowed_primary_before_percent_recovery": "secondary_validation_only",
            },
            "GSE179730": {
                "current_status": "sample manifest ready; patient response pending publication/source table.",
                "allowed_primary_before_response_recovery": False,
            },
        },
        "locked_outputs": [
            "03_rebuild/manifests/RESPONSE_METADATA_FREEZE.csv",
            "03_rebuild/manifests/RESPONSE_METADATA_FREEZE_AUDIT.csv",
        ],
    }
    (CONFIG_DIR / "response_harmonization_rules.json").write_text(json.dumps(config, indent=2), encoding="utf-8")

    md = [
        "# Response Metadata Freeze",
        "",
        f"Created: {config['created']}",
        "",
        "## Decision",
        "",
        "Response metadata is frozen only where it is available from local audited tables or GEO sample titles.",
        "Missing response labels are explicitly marked pending and must not be inferred from expression data.",
        "",
        "## Cohort Status",
        "",
        "| Cohort | Patients | Confirmed response patients | Pending patients |",
        "|---|---:|---:|---:|",
    ]
    for accession in sorted(set(str(r["accession"]) for r in all_rows)):
        subset = [r for r in all_rows if r["accession"] == accession]
        confirmed = [r for r in subset if "pending" not in str(r["response_status"]).lower() and r.get("response_harmonized_ordinal")]
        pending = [r for r in subset if "pending" in str(r["response_status"]).lower()]
        md.append(f"| {accession} | {len(subset)} | {len(confirmed)} | {len(pending)} |")
    md.extend(
        [
            "",
            "## Files Written",
            "",
            "- `03_rebuild/manifests/RESPONSE_METADATA_FREEZE.csv`",
            "- `03_rebuild/manifests/RESPONSE_METADATA_FREEZE_AUDIT.csv`",
            "- `03_rebuild/config/response_harmonization_rules.json`",
        ]
    )
    (OUT_DIR / "RESPONSE_METADATA_FREEZE_REPORT.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    print("\n".join(md))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
