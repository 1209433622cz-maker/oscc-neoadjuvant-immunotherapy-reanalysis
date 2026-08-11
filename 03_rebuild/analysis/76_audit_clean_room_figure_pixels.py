from __future__ import annotations

import csv
import hashlib
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image


BASELINE = Path(__file__).resolve().parents[2]
CLEAN = (
    BASELINE
    / "03_rebuild"
    / "clean_room"
    / "v11_clean_workspace"
)
QC_PATH = BASELINE / "03_rebuild" / "manuscript" / "NATURE_STYLE_FIGURE_QC.csv"
OUT_CSV = (
    BASELINE
    / "03_rebuild"
    / "manuscript"
    / "CLEAN_ROOM_FIGURE_PIXEL_AUDIT.csv"
)
OUT_MD = (
    BASELINE
    / "03_rebuild"
    / "manuscript"
    / "CLEAN_ROOM_FIGURE_PIXEL_AUDIT.md"
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest().upper()


def main() -> int:
    qc = pd.read_csv(QC_PATH)
    rows: list[dict[str, object]] = []
    for item in qc.to_dict("records"):
        relative = (
            Path("03_rebuild")
            / "figures"
            / str(item["directory"])
            / f"{item['stem']}.png"
        )
        baseline_path = BASELINE / relative
        clean_path = CLEAN / relative
        row: dict[str, object] = {
            "figure": item["figure"],
            "relative_path": relative.as_posix(),
            "baseline_exists": baseline_path.exists(),
            "clean_exists": clean_path.exists(),
            "baseline_size": "",
            "clean_size": "",
            "pixel_exact": False,
            "mean_absolute_pixel_difference": "",
            "max_absolute_pixel_difference": "",
            "baseline_sha256": "",
            "clean_sha256": "",
            "status": "FAIL",
        }
        if baseline_path.exists() and clean_path.exists():
            left = np.asarray(Image.open(baseline_path).convert("RGBA"))
            right = np.asarray(Image.open(clean_path).convert("RGBA"))
            row["baseline_size"] = f"{left.shape[1]}x{left.shape[0]}"
            row["clean_size"] = f"{right.shape[1]}x{right.shape[0]}"
            row["baseline_sha256"] = sha256(baseline_path)
            row["clean_sha256"] = sha256(clean_path)
            if left.shape == right.shape:
                difference = np.abs(left.astype(np.int16) - right.astype(np.int16))
                row["pixel_exact"] = bool(np.array_equal(left, right))
                row["mean_absolute_pixel_difference"] = float(difference.mean())
                row["max_absolute_pixel_difference"] = int(difference.max())
                row["status"] = "PASS" if row["pixel_exact"] else "CHECK"
            else:
                row["status"] = "DIMENSION_FAIL"
        rows.append(row)

    with OUT_CSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    pass_count = sum(row["status"] == "PASS" for row in rows)
    check_count = sum(row["status"] == "CHECK" for row in rows)
    fail_count = len(rows) - pass_count - check_count
    lines = [
        "# Clean-room Figure Pixel Audit",
        "",
        f"- Figures checked: {len(rows)}",
        f"- Pixel-exact PASS: {pass_count}",
        f"- Pixel-different CHECK: {check_count}",
        f"- Missing/dimension FAIL: {fail_count}",
        "",
        "## Results",
        "",
        "| Figure | Status | Mean absolute difference | Maximum difference |",
        "|---|---|---:|---:|",
    ]
    lines.extend(
        "| {figure} | {status} | {mean} | {maximum} |".format(
            figure=row["figure"],
            status=row["status"],
            mean=row["mean_absolute_pixel_difference"],
            maximum=row["max_absolute_pixel_difference"],
        )
        for row in rows
    )
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"PASS={pass_count} CHECK={check_count} FAIL={fail_count}")
    print(OUT_MD)
    return 0 if check_count == 0 and fail_count == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
