import importlib
import json
import os
import platform
import sys
from pathlib import Path


def module_status(name: str) -> dict:
    try:
        mod = importlib.import_module(name)
        version = getattr(mod, "__version__", "unknown")
        return {"package": name, "available": True, "version": str(version)}
    except Exception as exc:
        return {"package": name, "available": False, "version": "", "error": str(exc)}


def main() -> int:
    workspace = Path(os.environ.get("GSE200996_WORKSPACE", os.getcwd())).resolve()
    env_dir = workspace / "03_rebuild" / "env"
    log_dir = workspace / "03_rebuild" / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)

    req_file = env_dir / "requirements-python.txt"
    packages = []
    if req_file.exists():
        packages = [
            line.strip().replace("-", "_")
            for line in req_file.read_text(encoding="utf-8").splitlines()
            if line.strip() and not line.startswith("#")
        ]

    # Import names that differ from pip names.
    import_names = {
        "python_docx": "docx",
        "Pillow": "PIL",
        "biopython": "Bio",
    }
    statuses = []
    for pkg in packages:
        import_name = import_names.get(pkg, pkg)
        status = module_status(import_name)
        status["requirement"] = pkg
        statuses.append(status)

    payload = {
        "python": sys.executable,
        "version": sys.version,
        "platform": platform.platform(),
        "workspace": str(workspace),
        "packages": statuses,
    }

    json_path = log_dir / "python_environment_status_latest.json"
    json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    md_lines = [
        "# Python Environment Status",
        "",
        f"- Python: `{sys.executable}`",
        f"- Version: `{sys.version.splitlines()[0]}`",
        f"- Platform: `{platform.platform()}`",
        "",
        "| Requirement | Import | Available | Version |",
        "|---|---|---:|---:|",
    ]
    for item in statuses:
        md_lines.append(
            f"| {item['requirement']} | {item['package']} | {item['available']} | {item.get('version', '')} |"
        )

    md_path = log_dir / "PYTHON_ENVIRONMENT_STATUS_LATEST.md"
    md_path.write_text("\n".join(md_lines) + "\n", encoding="utf-8")
    print("\n".join(md_lines))

    missing = [item["requirement"] for item in statuses if not item["available"]]
    if missing:
        print("Missing Python packages:", ", ".join(missing), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
