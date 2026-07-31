#!/usr/bin/env python
"""Report docking and visualization tools visible in the current environment."""

from __future__ import annotations

import argparse
import hashlib
import importlib
import importlib.metadata
import json
import platform
from pathlib import Path
import shutil
import subprocess
import sys
from typing import Any


MODULES = ("vina", "meeko", "rdkit", "pymol", "numpy")
EXECUTABLES = ("vina", "mk_prepare_ligand.py", "mk_prepare_receptor.py", "pymol", "conda")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def module_version(name: str) -> dict[str, Any]:
    found = importlib.util.find_spec(name) is not None
    result: dict[str, Any] = {"available": found, "version": None}
    if not found:
        return result
    candidates = {
        "rdkit": ("rdkit",),
        "pymol": ("pymol-open-source", "pymol"),
    }.get(name, (name,))
    for distribution in candidates:
        try:
            result["version"] = importlib.metadata.version(distribution)
            break
        except importlib.metadata.PackageNotFoundError:
            continue
    if result["version"] is None:
        try:
            module = importlib.import_module(name)
            result["version"] = getattr(module, "__version__", None)
        except Exception as exc:  # importing optional scientific packages may fail
            result["import_error"] = f"{type(exc).__name__}: {exc}"
    return result


def conda_environments() -> list[dict[str, str]]:
    conda = shutil.which("conda")
    if not conda:
        return []
    try:
        completed = subprocess.run(
            [conda, "env", "list", "--json"],
            check=True,
            capture_output=True,
            text=True,
            timeout=30,
        )
        payload = json.loads(completed.stdout)
    except Exception:
        return []
    return [{"path": value, "name": Path(value).name} for value in payload.get("envs", [])]


def build_report(vina_path: Path | None = None) -> dict[str, Any]:
    modules = {name: module_version(name) for name in MODULES}
    executables = {name: shutil.which(name) for name in EXECUTABLES}
    external_vina = None
    if vina_path is not None:
        if not vina_path.is_file():
            raise FileNotFoundError(f"Vina executable does not exist: {vina_path}")
        completed = subprocess.run(
            [str(vina_path), "--version"],
            check=True,
            capture_output=True,
            text=True,
            timeout=30,
        )
        version_output = (completed.stdout or completed.stderr).strip()
        executables["vina"] = str(vina_path.resolve())
        external_vina = {
            "path": str(vina_path.resolve()),
            "sha256": sha256_file(vina_path),
            "version_output": version_output,
        }
    missing_core = [
        name
        for name in ("vina", "meeko")
        if not modules[name]["available"] and not executables.get(name)
    ]
    missing_visualization = (
        not modules["pymol"]["available"] and not executables.get("pymol")
    )
    return {
        "python": {
            "executable": sys.executable,
            "version": platform.python_version(),
            "platform": platform.platform(),
        },
        "modules": modules,
        "executables": executables,
        "external_vina": external_vina,
        "conda_environments": conda_environments(),
        "readiness": {
            "docking_engine_ready_in_current_environment": not missing_core,
            "pymol_ready_in_current_environment": not missing_visualization,
            "missing_docking_components": missing_core,
        },
        "note": (
            "Tools may exist in separate Conda environments. Run this script inside each "
            "intended environment and record both reports."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", type=Path, help="Optional JSON output path")
    parser.add_argument("--vina-path", type=Path, help="Explicit Vina executable to verify")
    args = parser.parse_args()
    report = build_report(args.vina_path)
    text = json.dumps(report, indent=2, ensure_ascii=False)
    print(text)
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(text + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
