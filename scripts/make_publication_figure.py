"""Run the verified PyMOL/RDKit publication-figure workflow across Conda environments."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
from typing import Any

from publication_common import (
    PUBLICATION_DPI,
    load_residue_map,
    parse_pdb_atoms,
    sha256_file,
    validate_report_atom_indices,
    validate_verification_report,
    write_json,
)


SCRIPT_DIR = Path(__file__).resolve().parent


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--complex", type=Path, required=True)
    parser.add_argument("--ligand-sdf", type=Path, required=True)
    parser.add_argument("--verification-report", type=Path, required=True)
    parser.add_argument("--residue-map", type=Path)
    parser.add_argument("--score-report", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--output-prefix", default="docking_publication")
    parser.add_argument("--protein-name", default="Protein")
    parser.add_argument("--ligand-name", default="Ligand")
    parser.add_argument("--title")
    parser.add_argument("--pymol-env", default="pymol")
    parser.add_argument("--chem-env", default="rdkit-vis")
    parser.add_argument("--conda-exe", default="conda")
    parser.add_argument("--overview-transparency", type=float, default=0.32)
    parser.add_argument("--pocket-transparency", type=float, default=0.55)
    parser.add_argument("--pocket-cutoff", type=float, default=4.0)
    parser.add_argument("--max-pocket-residues", type=int, default=10)
    parser.add_argument("--max-labels", type=int, default=6)
    parser.add_argument("--atom-map-tolerance-A", type=float, default=0.05)
    parser.add_argument(
        "--ray",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Use PyMOL ray tracing (default: enabled)",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Validate paths and print commands without creating outputs",
    )
    args = parser.parse_args()
    if not re.fullmatch(r"[A-Za-z0-9._-]+", args.output_prefix):
        parser.error("--output-prefix may contain only letters, digits, dot, underscore, hyphen")
    for name in ("overview_transparency", "pocket_transparency"):
        value = getattr(args, name)
        if not 0.0 <= value < 1.0:
            parser.error(f"--{name.replace('_', '-')} must be in [0, 1)")
    if min(
        args.pocket_cutoff,
        args.max_pocket_residues,
        args.max_labels,
        args.atom_map_tolerance_A,
    ) <= 0:
        parser.error("Cutoffs, tolerances, and item counts must be positive")
    return args


def require_file(path: Path, label: str) -> Path:
    resolved = path.resolve()
    if not resolved.is_file():
        raise FileNotFoundError(f"{label} does not exist: {resolved}")
    return resolved


def input_record(path: Path) -> dict[str, str]:
    return {"path": str(path), "sha256": sha256_file(path)}


def command_text(command: list[str]) -> str:
    return subprocess.list2cmdline(command) if sys.platform == "win32" else shlex.join(command)


def run_stage(
    name: str,
    command: list[str],
    logs_dir: Path,
) -> subprocess.CompletedProcess[str]:
    completed = subprocess.run(
        command,
        text=True,
        encoding="utf-8",
        errors="replace",
        capture_output=True,
        check=False,
    )
    log = (
        f"COMMAND\n{command_text(command)}\n\n"
        f"RETURN CODE\n{completed.returncode}\n\n"
        f"STDOUT\n{completed.stdout}\n\nSTDERR\n{completed.stderr}\n"
    )
    (logs_dir / f"{name}.log").write_text(log, encoding="utf-8")
    if completed.returncode != 0:
        raise RuntimeError(
            f"Stage '{name}' failed with exit code {completed.returncode}; "
            f"see {logs_dir / f'{name}.log'}"
        )
    return completed


def main() -> int:
    args = parse_args()
    complex_path = require_file(args.complex, "Complex")
    ligand_sdf = require_file(args.ligand_sdf, "Ligand SDF")
    verification_report = require_file(args.verification_report, "Verification report")
    residue_map = require_file(args.residue_map, "Residue map") if args.residue_map else None
    score_report = require_file(args.score_report, "Score report") if args.score_report else None
    verification = validate_verification_report(verification_report, complex_path)
    validate_report_atom_indices(verification, parse_pdb_atoms(complex_path))
    load_residue_map(residue_map)
    output_dir = args.output_dir.resolve()
    if output_dir.exists():
        if not output_dir.is_dir():
            raise NotADirectoryError(f"Output path is not a directory: {output_dir}")
        if any(output_dir.iterdir()):
            raise FileExistsError(
                f"Output directory is not empty: {output_dir}. Use a new run directory."
            )

    conda_executable = shutil.which(args.conda_exe)
    if conda_executable is None:
        candidate = Path(args.conda_exe)
        if not candidate.is_file():
            raise FileNotFoundError(f"Conda executable was not found: {args.conda_exe}")
        conda_executable = str(candidate.resolve())

    prefix = args.output_prefix
    common_render = [
        "--complex", str(complex_path),
        "--verification-report", str(verification_report),
        "--output-dir", str(output_dir),
        "--output-prefix", prefix,
    ]
    render_command = [
        conda_executable,
        "run", "-n", args.pymol_env, "python",
        str(SCRIPT_DIR / "pymol_publication_render.py"),
        *common_render,
        "--protein-name", args.protein_name,
        "--ligand-name", args.ligand_name,
        "--overview-transparency", str(args.overview_transparency),
        "--pocket-transparency", str(args.pocket_transparency),
        "--pocket-cutoff", str(args.pocket_cutoff),
        "--max-pocket-residues", str(args.max_pocket_residues),
        "--ray" if args.ray else "--no-ray",
    ]
    draw_command = [
        conda_executable,
        "run", "-n", args.chem_env, "python",
        str(SCRIPT_DIR / "draw_2d_interactions.py"),
        "--complex", str(complex_path),
        "--ligand-sdf", str(ligand_sdf),
        "--verification-report", str(verification_report),
        "--output-dir", str(output_dir),
        "--output-prefix", prefix,
        "--max-labels", str(args.max_labels),
        "--atom-map-tolerance-A", str(args.atom_map_tolerance_A),
    ]
    if residue_map is not None:
        render_command.extend(["--residue-map", str(residue_map)])
        draw_command.extend(["--residue-map", str(residue_map)])
    if score_report is not None:
        draw_command.extend(["--score-report", str(score_report)])

    compose_command = [
        conda_executable,
        "run", "-n", args.chem_env, "python",
        str(SCRIPT_DIR / "compose_publication_figure.py"),
        "--overview", str(output_dir / f"{prefix}_overall.png"),
        "--pocket", str(output_dir / f"{prefix}_pocket.png"),
        "--ligand-mask", str(output_dir / f"{prefix}_ligand_mask.png"),
        "--contacts-2d", str(output_dir / f"{prefix}_contacts_2d.png"),
        "--pymol-manifest", str(output_dir / f"{prefix}_pymol_manifest.json"),
        "--contact-report", str(output_dir / f"{prefix}_contacts.json"),
        "--output-dir", str(output_dir),
        "--output-prefix", prefix,
        "--protein-name", args.protein_name,
        "--ligand-name", args.ligand_name,
    ]
    if args.title:
        compose_command.extend(["--title", args.title])
    if score_report is not None:
        compose_command.extend(["--score-report", str(score_report)])

    commands = {
        "render_3d": render_command,
        "draw_2d": draw_command,
        "compose": compose_command,
    }
    if args.dry_run:
        print(
            json.dumps(
                {name: command_text(command) for name, command in commands.items()},
                indent=2,
                ensure_ascii=False,
            )
        )
        return 0

    output_dir.mkdir(parents=True, exist_ok=True)
    logs_dir = output_dir / "logs"
    logs_dir.mkdir()
    inputs: dict[str, Any] = {
        "complex": input_record(complex_path),
        "ligand_sdf": input_record(ligand_sdf),
        "verification_report": input_record(verification_report),
        "residue_map": input_record(residue_map) if residue_map else None,
        "score_report": input_record(score_report) if score_report else None,
    }
    manifest_path = output_dir / f"{prefix}_run_manifest.json"
    manifest: dict[str, Any] = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "running",
        "inputs": inputs,
        "parameters": {
            "protein_name": args.protein_name,
            "ligand_name": args.ligand_name,
            "title": args.title,
            "dpi_metadata": PUBLICATION_DPI,
            "overview_transparency": args.overview_transparency,
            "pocket_transparency": args.pocket_transparency,
            "pocket_cutoff_A": args.pocket_cutoff,
            "max_pocket_residues": args.max_pocket_residues,
            "max_labels": args.max_labels,
            "atom_map_tolerance_A": args.atom_map_tolerance_A,
            "ray": args.ray,
        },
        "environments": {"pymol": args.pymol_env, "chemistry": args.chem_env},
        "commands": {name: command_text(value) for name, value in commands.items()},
    }
    write_json(manifest_path, manifest)
    try:
        for stage, command in commands.items():
            print(f"Running {stage} ...", flush=True)
            run_stage(stage, command, logs_dir)
    except Exception as exc:
        manifest["status"] = "failed"
        manifest["finished_utc"] = datetime.now(timezone.utc).isoformat()
        manifest["error"] = str(exc)
        write_json(manifest_path, manifest)
        raise

    expected = [
        output_dir / f"{prefix}_overall.png",
        output_dir / f"{prefix}_pocket.png",
        output_dir / f"{prefix}_ligand_mask.png",
        output_dir / f"{prefix}.pse",
        output_dir / f"{prefix}_pymol_manifest.json",
        output_dir / f"{prefix}_contacts_2d.png",
        output_dir / f"{prefix}_contacts.json",
        output_dir / f"{prefix}_triptych.png",
        output_dir / f"{prefix}_figure_manifest.json",
    ]
    missing = [path for path in expected if not path.is_file()]
    if missing:
        error = f"Workflow completed but outputs are missing: {missing}"
        manifest["status"] = "failed"
        manifest["finished_utc"] = datetime.now(timezone.utc).isoformat()
        manifest["error"] = error
        write_json(manifest_path, manifest)
        raise RuntimeError(error)
    manifest["status"] = "complete"
    manifest["finished_utc"] = datetime.now(timezone.utc).isoformat()
    manifest["outputs"] = [
        {"path": str(path), "sha256": sha256_file(path)} for path in expected
    ]
    write_json(manifest_path, manifest)
    print(f"Wrote publication figure: {output_dir / f'{prefix}_triptych.png'}")
    print(f"Wrote editable PyMOL session: {output_dir / f'{prefix}.pse'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
