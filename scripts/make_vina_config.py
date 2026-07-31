#!/usr/bin/env python
"""Create an auditable AutoDock Vina config and parameter-plan JSON."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import re
from typing import Iterable


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atom_element_from_pdb(line: str) -> str:
    element = line[76:78].strip() if len(line) >= 78 else ""
    if element:
        return element.upper()
    name = line[12:16].strip().lstrip("0123456789")
    return (name[:2] if len(name) > 1 and name[:2].upper() in {"CL", "BR"} else name[:1]).upper()


def read_pdb_coordinates(path: Path, resn: str | None) -> list[tuple[float, float, float]]:
    coordinates: list[tuple[float, float, float]] = []
    saw_model = False
    in_first_model = True
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            record = line[:6].strip().upper()
            if record == "MODEL":
                if saw_model:
                    in_first_model = False
                saw_model = True
                continue
            if record == "ENDMDL" and saw_model:
                break
            if record not in {"ATOM", "HETATM"} or not in_first_model:
                continue
            if resn and line[17:20].strip().upper() != resn.upper():
                continue
            if atom_element_from_pdb(line) == "H":
                continue
            try:
                coordinates.append(
                    (float(line[30:38]), float(line[38:46]), float(line[46:54]))
                )
            except ValueError as exc:
                raise ValueError(f"Invalid coordinate record in {path}: {line.rstrip()}") from exc
    return coordinates


def read_sdf_coordinates(path: Path) -> list[tuple[float, float, float]]:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    if any("V3000" in line for line in lines[:10]):
        coordinates: list[tuple[float, float, float]] = []
        in_atoms = False
        for line in lines:
            if "M  V30 BEGIN ATOM" in line:
                in_atoms = True
                continue
            if "M  V30 END ATOM" in line:
                break
            if not in_atoms:
                continue
            tokens = line.split()
            if len(tokens) >= 7 and tokens[:2] == ["M", "V30"]:
                element = tokens[3].upper()
                if element != "H":
                    coordinates.append(tuple(float(value) for value in tokens[4:7]))
        return coordinates
    if len(lines) < 4:
        return []
    try:
        atom_count = int(lines[3][:3])
    except ValueError as exc:
        raise ValueError(f"Cannot parse V2000 atom count in {path}") from exc
    coordinates = []
    for line in lines[4 : 4 + atom_count]:
        if len(line) < 34:
            continue
        element = line[31:34].strip().upper()
        if element == "H":
            continue
        coordinates.append((float(line[0:10]), float(line[10:20]), float(line[20:30])))
    return coordinates


def read_mol2_coordinates(path: Path) -> list[tuple[float, float, float]]:
    coordinates: list[tuple[float, float, float]] = []
    in_atoms = False
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("@<TRIPOS>ATOM"):
            in_atoms = True
            continue
        if line.startswith("@<TRIPOS>") and in_atoms:
            break
        if not in_atoms or not line.strip():
            continue
        tokens = line.split()
        if len(tokens) < 6:
            continue
        element = re.split(r"[.0-9]", tokens[5], maxsplit=1)[0].upper()
        if element != "H":
            coordinates.append(tuple(float(value) for value in tokens[2:5]))
    return coordinates


def read_coordinates(path: Path, resn: str | None = None) -> list[tuple[float, float, float]]:
    suffix = path.suffix.lower()
    if suffix in {".pdb", ".pdbqt", ".ent"}:
        coordinates = read_pdb_coordinates(path, resn)
    elif suffix in {".sdf", ".mol"}:
        if resn:
            raise ValueError("--reference-resn is only supported for PDB/PDBQT references")
        coordinates = read_sdf_coordinates(path)
    elif suffix == ".mol2":
        if resn:
            raise ValueError("--reference-resn is only supported for PDB/PDBQT references")
        coordinates = read_mol2_coordinates(path)
    else:
        raise ValueError(f"Unsupported coordinate format: {path.suffix}")
    if not coordinates:
        qualifier = f" for residue {resn}" if resn else ""
        raise ValueError(f"No heavy-atom coordinates found in {path}{qualifier}")
    for xyz in coordinates:
        if not all(math.isfinite(value) for value in xyz):
            raise ValueError(f"Non-finite coordinate found in {path}")
    return coordinates


def extent_and_center(
    coordinates: Iterable[tuple[float, float, float]],
) -> tuple[list[float], list[float], list[float], list[float]]:
    values = list(coordinates)
    minima = [min(point[axis] for point in values) for axis in range(3)]
    maxima = [max(point[axis] for point in values) for axis in range(3)]
    extent = [maxima[axis] - minima[axis] for axis in range(3)]
    center = [(minima[axis] + maxima[axis]) / 2.0 for axis in range(3)]
    return minima, maxima, extent, center


def ensure_writeable(path: Path, force: bool) -> None:
    if path.exists() and not force:
        raise FileExistsError(f"Refusing to overwrite {path}; pass --force to replace it")
    path.parent.mkdir(parents=True, exist_ok=True)


def validate_docking_file(path: Path, role: str) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"{role} file does not exist: {path}")
    if path.suffix.lower() != ".pdbqt":
        raise ValueError(f"{role} must be docking-ready PDBQT, received: {path}")


def format_number(value: float) -> str:
    return f"{value:.4f}".rstrip("0").rstrip(".")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--receptor", type=Path, required=True)
    parser.add_argument("--ligand", type=Path, required=True)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--reference", type=Path, help="Reference ligand/site coordinate file")
    source.add_argument("--center", nargs=3, type=float, metavar=("X", "Y", "Z"))
    parser.add_argument("--size", nargs=3, type=float, metavar=("SX", "SY", "SZ"))
    parser.add_argument(
        "--box-rationale",
        help=(
            "Independent evidence for the selected site/box. Required with an explicit "
            "--center; do not cite the final docked pose itself."
        ),
    )
    parser.add_argument("--reference-resn", help="Residue name to select from a PDB reference")
    parser.add_argument(
        "--padding",
        type=float,
        default=5.0,
        help="Padding per side in Angstrom when --reference is used (default: 5)",
    )
    parser.add_argument("--scoring", choices=("vina", "vinardo", "ad4"), default="vina")
    parser.add_argument("--exhaustiveness", type=int, default=8)
    parser.add_argument("--num-modes", type=int, default=9)
    parser.add_argument("--energy-range", type=float, default=3.0)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--cpu", type=int, default=0)
    parser.add_argument("--poses", type=Path, default=Path("docked_poses.pdbqt"))
    parser.add_argument("--output", type=Path, required=True, help="Vina config output")
    parser.add_argument("--plan", type=Path, required=True, help="JSON parameter-plan output")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    validate_docking_file(args.receptor, "Receptor")
    validate_docking_file(args.ligand, "Ligand")
    if args.exhaustiveness < 1 or args.num_modes < 1:
        parser.error("--exhaustiveness and --num-modes must be positive")
    if args.energy_range <= 0 or args.padding < 0 or args.cpu < 0:
        parser.error("--energy-range must be positive; --padding and --cpu cannot be negative")

    reference_details = None
    if args.reference:
        if not args.reference.is_file():
            raise FileNotFoundError(f"Reference file does not exist: {args.reference}")
        coordinates = read_coordinates(args.reference, args.reference_resn)
        minima, maxima, reference_extent, center = extent_and_center(coordinates)
        size = [value + 2.0 * args.padding for value in reference_extent]
        reference_details = {
            "path": str(args.reference.resolve()),
            "sha256": sha256_file(args.reference),
            "residue_filter": args.reference_resn,
            "heavy_atom_count": len(coordinates),
            "minimum_A": minima,
            "maximum_A": maxima,
            "extent_A": reference_extent,
            "padding_per_side_A": args.padding,
            "formula": "size_axis = extent_axis + 2 * padding_per_side",
        }
        box_rationale = (
            args.box_rationale.strip()
            if args.box_rationale
            else "reference-coordinate extent plus per-side padding"
        )
    else:
        if args.size is None:
            parser.error("--size SX SY SZ is required with an explicit --center")
        if not args.box_rationale or not args.box_rationale.strip():
            parser.error(
                "--box-rationale is required with an explicit --center; state the cognate "
                "ligand, experimental residues, or independent pocket evidence"
            )
        center = list(args.center)
        size = list(args.size)
        box_rationale = args.box_rationale.strip()

    if any(value <= 0 or not math.isfinite(value) for value in size):
        parser.error("All box dimensions must be finite and positive")
    if any(not math.isfinite(value) for value in center):
        parser.error("All center coordinates must be finite")

    ligand_coordinates = read_coordinates(args.ligand)
    _, _, ligand_extent, _ = extent_and_center(ligand_coordinates)
    volume = size[0] * size[1] * size[2]
    warnings: list[str] = []
    if any(size[axis] < ligand_extent[axis] + 4.0 for axis in range(3)):
        warnings.append(
            "At least one box axis leaves <2 Å clearance per side relative to the input "
            "ligand extent; visually inspect and justify."
        )
    if volume > 27000:
        warnings.append(
            "Box volume exceeds 27,000 Å^3; search convergence may require substantially more "
            "sampling or a better-localized site."
        )
    if args.scoring == "ad4":
        warnings.append(
            "AD4 requires its compatible affinity-map workflow; this config alone may not be "
            "sufficient, and AD4 scores are not comparable to Vina/Vinardo scores."
        )
    if args.exhaustiveness == 8:
        warnings.append(
            "Exhaustiveness 8 is the official default, not evidence of convergence; run "
            "multi-seed/exhaustiveness sensitivity checks."
        )

    config_lines = [
        f"receptor = {args.receptor.resolve()}",
        f"ligand = {args.ligand.resolve()}",
        f"center_x = {format_number(center[0])}",
        f"center_y = {format_number(center[1])}",
        f"center_z = {format_number(center[2])}",
        f"size_x = {format_number(size[0])}",
        f"size_y = {format_number(size[1])}",
        f"size_z = {format_number(size[2])}",
        f"scoring = {args.scoring}",
        f"exhaustiveness = {args.exhaustiveness}",
        f"num_modes = {args.num_modes}",
        f"energy_range = {format_number(args.energy_range)}",
        f"seed = {args.seed}",
        f"cpu = {args.cpu}",
        f"out = {args.poses.resolve()}",
    ]

    ensure_writeable(args.output, args.force)
    ensure_writeable(args.plan, args.force)
    args.output.write_text("\n".join(config_lines) + "\n", encoding="utf-8")

    plan = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "receptor": {
                "path": str(args.receptor.resolve()),
                "sha256": sha256_file(args.receptor),
            },
            "ligand": {
                "path": str(args.ligand.resolve()),
                "sha256": sha256_file(args.ligand),
                "heavy_atom_extent_A": ligand_extent,
            },
            "reference": reference_details,
        },
        "box": {
            "center_A": center,
            "size_A": size,
            "volume_A3": volume,
            "rationale": box_rationale,
        },
        "parameters": {
            "scoring": args.scoring,
            "exhaustiveness": args.exhaustiveness,
            "num_modes": args.num_modes,
            "energy_range_kcal_mol": args.energy_range,
            "seed": args.seed,
            "cpu": args.cpu,
        },
        "outputs": {
            "config": str(args.output.resolve()),
            "poses": str(args.poses.resolve()),
        },
        "warnings": warnings,
        "interpretation": (
            "These parameters define one docking run. Validate search convergence across seeds "
            "and, when possible, by cognate-ligand redocking."
        ),
    }
    args.plan.write_text(
        json.dumps(plan, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(f"Wrote {args.output}")
    print(f"Wrote {args.plan}")
    for warning in warnings:
        print(f"WARNING: {warning}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
