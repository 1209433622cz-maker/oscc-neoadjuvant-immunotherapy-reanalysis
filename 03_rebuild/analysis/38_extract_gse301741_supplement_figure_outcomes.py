#!/usr/bin/env python
"""Extract provisional GSE301741 response labels from supplement figure outcome bars.

The deposited GEO/SRA metadata does not expose patient-level response labels. The
publication supplement contains heatmap annotation bars with explicit Patient,
Timepoint and Outcome rows. This script quantifies the Outcome row in two
independent panels from Supplementary Fig. 1 and writes a provisional,
figure-derived response table with audit counts.
"""

from __future__ import annotations

import csv
import os
import shutil
import subprocess
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path

from PIL import Image


PATIENTS = [
    "LSCC20",
    "LSCC21",
    "OCSCC1",
    "OCSCC12",
    "OCSCC13",
    "OCSCC19",
    "OCSCC3",
    "OCSCC30",
    "OCSCC31",
    "OCSCC33",
    "OCSCC34",
    "OCSCC35",
    "OCSCC36",
    "OCSCC37",
    "OCSCC38",
    "OPSCC41",
]

PATIENT_COLORS = dict(
    zip(
        PATIENTS,
        [
            (253, 215, 15),
            (174, 37, 196),
            (223, 157, 118),
            (86, 64, 138),
            (220, 97, 95),
            (141, 239, 143),
            (97, 67, 232),
            (70, 231, 232),
            (198, 207, 244),
            (113, 136, 84),
            (133, 197, 123),
            (222, 162, 229),
            (126, 233, 104),
            (219, 168, 70),
            (228, 65, 217),
            (230, 242, 150),
        ],
    )
)

OUTCOME_COLORS = {
    "non_responder": (38, 94, 15),
    "responder": (253, 101, 179),
}

TIMEPOINT_COLORS = {
    "post": (253, 23, 15),
    "pre": (125, 253, 15),
}


def workspace_root() -> Path:
    return Path(__file__).resolve().parents[2]


def locate_pdftoppm() -> Path:
    configured = os.environ.get("POPPLER_PDFTOPPM", "")
    candidates = [
        Path(configured) if configured else None,
        Path(shutil.which("pdftoppm")) if shutil.which("pdftoppm") else None,
        Path.home()
        / ".cache"
        / "codex-runtimes"
        / "codex-primary-runtime"
        / "dependencies"
        / "native"
        / "poppler"
        / "Library"
        / "bin"
        / "pdftoppm.exe",
    ]
    for candidate in candidates:
        if candidate and candidate.exists():
            return candidate.resolve()
    raise FileNotFoundError(
        "pdftoppm was not found. Set POPPLER_PDFTOPPM to the Poppler executable."
    )


def render_evidence_panels(
    workspace: Path,
    image_dir: Path,
    panel_paths: dict[str, Path],
) -> None:
    if all(path.exists() for path in panel_paths.values()):
        return

    source_pdf = (
        workspace
        / "02_references"
        / "external_supplements"
        / "GSE301741"
        / "mmc1.pdf"
    )
    if not source_pdf.exists():
        raise FileNotFoundError(source_pdf)

    image_dir.mkdir(parents=True, exist_ok=True)
    page_prefix = image_dir / "mmc1_page2_360dpi"
    page_image = page_prefix.with_suffix(".png")
    command = [
        str(locate_pdftoppm()),
        "-f",
        "2",
        "-l",
        "2",
        "-singlefile",
        "-r",
        "360",
        "-png",
        str(source_pdf),
        str(page_prefix),
    ]
    subprocess.run(command, check=True)
    if not page_image.exists():
        raise RuntimeError(f"Rendered PDF page was not created: {page_image}")

    page = Image.open(page_image).convert("RGB")
    if page.size != (3060, 3960):
        raise RuntimeError(
            f"Unexpected 360-dpi page size {page.size}; expected (3060, 3960)."
        )
    crops = {
        "supp_fig1_panelG_Th_cells": (200, 2180, 1080, 3890),
        "supp_fig1_panelI_Treg_cells": (2040, 2180, 2920, 3890),
    }
    for panel, bounds in crops.items():
        page.crop(bounds).save(panel_paths[panel])


def distance(a: tuple[int, int, int], b: tuple[int, int, int]) -> float:
    return sum((a[i] - b[i]) ** 2 for i in range(3)) ** 0.5


def mean_color(img: Image.Image, x: int, y0: int, y1: int) -> tuple[int, int, int]:
    vals = [img.getpixel((x, y)) for y in range(y0, y1)]
    return tuple(round(sum(v[i] for v in vals) / len(vals)) for i in range(3))


def nearest(color: tuple[int, int, int], palette: dict[str, tuple[int, int, int]]) -> tuple[str, float]:
    return min(((key, distance(color, value)) for key, value in palette.items()), key=lambda item: item[1])


def analyze_panel(path: Path) -> dict[str, dict[str, object]]:
    img = Image.open(path).convert("RGB")
    counts: dict[str, Counter[str]] = defaultdict(Counter)
    timepoints: dict[str, Counter[str]] = defaultdict(Counter)
    x_positions: dict[str, list[int]] = defaultdict(list)

    for x in range(100, 780):
        patient_color = mean_color(img, x, 225, 260)
        if max(patient_color) - min(patient_color) < 25 or sum(patient_color) > 720:
            continue
        patient, patient_distance = nearest(patient_color, PATIENT_COLORS)
        if patient_distance > 45:
            continue

        outcome_color = mean_color(img, x, 295, 330)
        outcome, outcome_distance = nearest(outcome_color, OUTCOME_COLORS)
        if outcome_distance > 55:
            continue

        timepoint_color = mean_color(img, x, 260, 295)
        timepoint, timepoint_distance = nearest(timepoint_color, TIMEPOINT_COLORS)
        if timepoint_distance > 60:
            timepoint = "unclassified"

        counts[patient][outcome] += 1
        timepoints[patient][timepoint] += 1
        x_positions[patient].append(x)

    panel_result: dict[str, dict[str, object]] = {}
    for patient in PATIENTS:
        total = sum(counts[patient].values())
        if total == 0:
            label = ""
            support_fraction = 0.0
        else:
            label, support_n = counts[patient].most_common(1)[0]
            support_fraction = support_n / total
        panel_result[patient] = {
            "label": label,
            "total_columns": total,
            "support_fraction": support_fraction,
            "outcome_counts": dict(counts[patient]),
            "timepoint_counts": dict(timepoints[patient]),
            "x_min": min(x_positions[patient]) if x_positions[patient] else "",
            "x_max": max(x_positions[patient]) if x_positions[patient] else "",
        }
    return panel_result


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def main() -> None:
    workspace = workspace_root()
    image_dir = workspace / "03_rebuild" / "tmp" / "pdf_page_checks" / "GSE301741"
    panel_paths = {
        "supp_fig1_panelG_Th_cells": image_dir / "panelG_full_360dpi.png",
        "supp_fig1_panelI_Treg_cells": image_dir / "panelI_full_360dpi.png",
    }
    render_evidence_panels(workspace, image_dir, panel_paths)
    missing = [str(path) for path in panel_paths.values() if not path.exists()]
    if missing:
        raise SystemExit("Missing rendered evidence images: " + "; ".join(missing))

    panel_results = {panel: analyze_panel(path) for panel, path in panel_paths.items()}
    rows: list[dict[str, object]] = []
    for patient in PATIENTS:
        labels = [panel_results[panel][patient]["label"] for panel in panel_paths]
        concordance = len(set(labels)) == 1 and labels[0] != ""
        response_label = labels[0] if concordance else ""
        rows.append(
            {
                "patient_id": patient,
                "response_label": response_label,
                "pTR_class": "",
                "pTR_percent": "",
                "panelG_label": panel_results["supp_fig1_panelG_Th_cells"][patient]["label"],
                "panelG_columns": panel_results["supp_fig1_panelG_Th_cells"][patient]["total_columns"],
                "panelG_support_fraction": panel_results["supp_fig1_panelG_Th_cells"][patient]["support_fraction"],
                "panelG_outcome_counts": panel_results["supp_fig1_panelG_Th_cells"][patient]["outcome_counts"],
                "panelI_label": panel_results["supp_fig1_panelI_Treg_cells"][patient]["label"],
                "panelI_columns": panel_results["supp_fig1_panelI_Treg_cells"][patient]["total_columns"],
                "panelI_support_fraction": panel_results["supp_fig1_panelI_Treg_cells"][patient]["support_fraction"],
                "panelI_outcome_counts": panel_results["supp_fig1_panelI_Treg_cells"][patient]["outcome_counts"],
                "concordance_status": "concordant" if concordance else "discordant_or_missing",
                "response_source": "publication_supplement_figure",
                "source_file": str(
                    workspace
                    / "02_references"
                    / "external_supplements"
                    / "GSE301741"
                    / "mmc1.pdf"
                ),
                "source_detail": "Supplementary Fig. 1, panels G and I, page 2; Outcome row color bar.",
                "method": "pixel_color_extraction_from_360dpi_rendered_supplement_figure",
                "provenance_status": "figure_derived_requires_rds_or_table_crosscheck",
            }
        )

    out_dir = workspace / "03_rebuild" / "validation" / "GSE301741_response_recovery"
    out_csv = out_dir / "GSE301741_RESPONSE_LABELS_FIGURE_DERIVED.csv"
    fields = [
        "patient_id",
        "response_label",
        "pTR_class",
        "pTR_percent",
        "panelG_label",
        "panelG_columns",
        "panelG_support_fraction",
        "panelG_outcome_counts",
        "panelI_label",
        "panelI_columns",
        "panelI_support_fraction",
        "panelI_outcome_counts",
        "concordance_status",
        "response_source",
        "source_file",
        "source_detail",
        "method",
        "provenance_status",
    ]
    write_csv(out_csv, rows, fields)

    counts = Counter(row["response_label"] for row in rows)
    md = out_dir / "GSE301741_RESPONSE_LABELS_FIGURE_DERIVED_REPORT.md"
    lines = [
        "# GSE301741 Supplement-Figure Response Label Extraction",
        "",
        f"Created: {datetime.now().isoformat(timespec='seconds')}",
        "",
        "## Decision",
        "",
        "Patient-level responder/non-responder labels were extracted from the explicit Outcome annotation bars in Supplementary Fig. 1 panels G and I. These labels are publication-supplement-derived but still require cross-checking against RDS metadata or a tabular supplement before they carry final manuscript response-validation claims.",
        "",
        "## Counts",
        "",
        f"- Responders: {counts.get('responder', 0)}",
        f"- Non-responders: {counts.get('non_responder', 0)}",
        f"- Missing/discordant: {counts.get('', 0)}",
        "",
        "## Extracted Labels",
        "",
        "| Patient | Label | Panel G | Panel I | Concordance |",
        "|---|---|---|---|---|",
    ]
    for row in rows:
        lines.append(
            f"| {row['patient_id']} | {row['response_label']} | {row['panelG_label']} | {row['panelI_label']} | {row['concordance_status']} |"
        )
    lines.extend(
        [
            "",
            "## Files",
            "",
            "- `03_rebuild/validation/GSE301741_response_recovery/GSE301741_RESPONSE_LABELS_FIGURE_DERIVED.csv`",
            "- Rendered evidence: `03_rebuild/tmp/pdf_page_checks/GSE301741/panelG_full_360dpi.png` and `panelI_full_360dpi.png`",
        ]
    )
    md.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"Wrote: {out_csv}")
    print(f"Wrote: {md}")
    print(f"Responders: {counts.get('responder', 0)}")
    print(f"Non-responders: {counts.get('non_responder', 0)}")
    print(f"Missing/discordant: {counts.get('', 0)}")


if __name__ == "__main__":
    main()
