"""Draw an audited RDKit 2D contact panel from a verified docking complex."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
from io import BytesIO
import json
import math
from pathlib import Path
import re
from typing import Any

from PIL import Image, ImageDraw, ImageFont
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D

from publication_common import (
    PUBLICATION_DPI,
    display_residue,
    ensure_targets,
    ligand_residue_key,
    load_json,
    load_residue_map,
    parse_pdb_atoms,
    residue_key,
    sha256_file,
    source_atom_label,
    validate_report_atom_indices,
    validate_verification_report,
    write_json,
)


CANVAS_SIZE = (1400, 900)
MOLECULE_SIZE = (860, 620)
MOLECULE_ORIGIN = (270, 140)


def resolve_font(bold: bool, size: int) -> ImageFont.FreeTypeFont:
    candidates = (
        [Path(r"C:\Windows\Fonts\arialbd.ttf"), Path("DejaVuSans-Bold.ttf")]
        if bold
        else [Path(r"C:\Windows\Fonts\arial.ttf"), Path("DejaVuSans.ttf")]
    )
    for candidate in candidates:
        try:
            return ImageFont.truetype(str(candidate), size)
        except OSError:
            continue
    raise OSError("No usable Arial or DejaVu Sans font was found")


def dashed_line(
    draw: ImageDraw.ImageDraw,
    start: tuple[float, float],
    end: tuple[float, float],
    *,
    fill: str,
    width: int,
    dash: int,
    gap: int,
) -> None:
    x1, y1 = start
    x2, y2 = end
    dx = x2 - x1
    dy = y2 - y1
    length = math.hypot(dx, dy)
    if length == 0:
        return
    ux, uy = dx / length, dy / length
    position = 0.0
    while position < length:
        segment_end = min(position + dash, length)
        draw.line(
            (
                x1 + ux * position,
                y1 + uy * position,
                x1 + ux * segment_end,
                y1 + uy * segment_end,
            ),
            fill=fill,
            width=width,
        )
        position += dash + gap


def draw_badge(
    draw: ImageDraw.ImageDraw,
    center: tuple[int, int],
    title: str,
    detail: str | None,
    *,
    candidate: bool,
) -> None:
    title_font = resolve_font(True, 28)
    detail_font = resolve_font(False, 21)
    title_box = draw.textbbox((0, 0), title, font=title_font)
    title_width = title_box[2] - title_box[0]
    width = max(170, title_width + 42)
    height = 82 if detail else 62
    left = center[0] - width // 2
    top = center[1] - height // 2
    right, bottom = left + width, top + height
    fill = "#FFF3C7" if candidate else "#F0F3F7"
    outline = "#D89B00" if candidate else "#A8B2BE"
    title_color = "#725000" if candidate else "#3E4B59"
    draw.rounded_rectangle(
        (left, top, right, bottom),
        radius=18,
        fill=fill,
        outline=outline,
        width=3,
    )
    title_y = top + 8 if detail else top + 13
    draw.text(
        (center[0] - title_width / 2, title_y),
        title,
        font=title_font,
        fill=title_color,
    )
    if detail:
        detail_box = draw.textbbox((0, 0), detail, font=detail_font)
        detail_width = detail_box[2] - detail_box[0]
        draw.text(
            (center[0] - detail_width / 2, top + 46),
            detail,
            font=detail_font,
            fill="#7A5A12",
        )


def residue_proximities(
    receptor_atoms: list[dict[str, Any]],
    ligand_atoms: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    residues: dict[tuple[str, str, str], list[dict[str, Any]]] = {}
    for atom in receptor_atoms:
        if not atom["is_hydrogen"]:
            residues.setdefault(residue_key(atom), []).append(atom)
    ligand_heavy = [atom for atom in ligand_atoms if not atom["is_hydrogen"]]
    output = []
    for key, atoms in residues.items():
        best = min(
            (
                (math.dist(rec["coord_A"], lig["coord_A"]), rec, lig)
                for rec in atoms
                for lig in ligand_heavy
            ),
            key=lambda item: item[0],
        )
        output.append(
            {
                "source_residue": key,
                "distance_A": best[0],
                "receptor_atom": best[1],
                "ligand_atom": best[2],
                "kind": "pocket_proximity",
            }
        )
    return sorted(output, key=lambda item: item["distance_A"])


def candidate_interactions(report: dict[str, Any]) -> list[dict[str, Any]]:
    shortest_by_residue: dict[tuple[str, str, str], dict[str, Any]] = {}
    for contact in report["geometry"]["candidate_polar_contacts"]:
        key = residue_key(contact["receptor_atom"])
        item = {
            "source_residue": key,
            "distance_A": float(contact["distance_A"]),
            "receptor_atom": contact["receptor_atom"],
            "ligand_atom": contact["ligand_atom"],
            "kind": "candidate_polar_contact",
        }
        prior = shortest_by_residue.get(key)
        if prior is None or item["distance_A"] < prior["distance_A"]:
            shortest_by_residue[key] = item
    return sorted(shortest_by_residue.values(), key=lambda item: item["distance_A"])


def coordinate_atom_map(
    molecule: Chem.Mol,
    pdb_heavy: list[dict[str, Any]],
    tolerance_A: float,
) -> tuple[dict[str, int], float]:
    if molecule.GetNumConformers() == 0:
        raise ValueError("Ligand SDF has no coordinates for atom mapping")
    conformer = molecule.GetConformer()
    available = set(range(molecule.GetNumAtoms()))
    mapping: dict[str, int] = {}
    distances: list[float] = []
    for pdb_atom in pdb_heavy:
        matches = []
        for index in available:
            rd_atom = molecule.GetAtomWithIdx(index)
            if rd_atom.GetSymbol().upper() != pdb_atom["element"].upper():
                continue
            position = conformer.GetAtomPosition(index)
            distance = math.dist(
                pdb_atom["coord_A"],
                (position.x, position.y, position.z),
            )
            if distance <= tolerance_A:
                matches.append((distance, index))
        if len(matches) != 1:
            raise ValueError(
                "Could not derive a unique coordinate-based SDF/PDB atom mapping for "
                f"{pdb_atom['name']}: found {len(matches)} matches within {tolerance_A} Å"
            )
        distance, index = matches[0]
        if pdb_atom["name"] in mapping:
            raise ValueError(f"Ligand atom name is not unique: {pdb_atom['name']}")
        mapping[pdb_atom["name"]] = index
        available.remove(index)
        distances.append(distance)
    if available:
        raise ValueError("SDF contains heavy atoms not mapped to the verified complex ligand")
    return mapping, max(distances, default=0.0)


def pdbqt_atom_element(line: str) -> str:
    atom_type = line.split()[-1].upper()
    mapping = {
        "A": "C",
        "C": "C",
        "N": "N",
        "NA": "N",
        "O": "O",
        "OA": "O",
        "S": "S",
        "SA": "S",
        "P": "P",
        "F": "F",
        "CL": "CL",
        "BR": "BR",
        "I": "I",
        "H": "H",
        "HD": "H",
        "HS": "H",
    }
    if atom_type not in mapping:
        raise ValueError(f"Unsupported PDBQT atom type in score pose: {atom_type}")
    return mapping[atom_type]


def first_pdbqt_pose_atoms(path: Path) -> list[dict[str, Any]]:
    atoms: list[dict[str, Any]] = []
    in_first_model = False
    saw_model = False
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        record = line[:6].strip().upper()
        if record == "MODEL":
            if saw_model:
                break
            saw_model = True
            in_first_model = True
            continue
        if record == "ENDMDL" and in_first_model:
            break
        if saw_model and not in_first_model:
            continue
        if record not in {"ATOM", "HETATM"}:
            continue
        try:
            coordinate = [
                float(line[30:38]),
                float(line[38:46]),
                float(line[46:54]),
            ]
        except ValueError as exc:
            raise ValueError(f"Invalid PDBQT coordinate line: {line}") from exc
        atoms.append(
            {
                "element": pdbqt_atom_element(line),
                "coord_A": coordinate,
            }
        )
    if not atoms:
        raise ValueError(f"No pose atoms were found in score output: {path}")
    return atoms


def verify_score_pose_association(
    score_report_path: Path,
    molecule: Chem.Mol,
    tolerance_A: float,
) -> dict[str, Any]:
    score_report = load_json(score_report_path)
    best_score = score_report.get("best_affinity_kcal_mol")
    poses = score_report.get("poses")
    if best_score is None or not isinstance(poses, list) or not poses:
        raise ValueError("Score report lacks best_affinity_kcal_mol or pose table")
    first_pose = poses[0]
    if int(first_pose.get("mode", -1)) != 1:
        raise ValueError("The first score-table entry is not Vina mode 1")
    if abs(float(first_pose["affinity_kcal_mol"]) - float(best_score)) > 1e-6:
        raise ValueError("Best score does not match the mode-1 score-table entry")

    pose_record = score_report.get("pose_output_file")
    if not isinstance(pose_record, dict) or not pose_record.get("path"):
        return {
            "status": "caller_asserted",
            "reason": "score report has no accessible pose_output_file record",
            "mode": 1,
            "affinity_kcal_mol": float(best_score),
        }
    pose_path = Path(str(pose_record["path"]))
    if not pose_path.is_file():
        return {
            "status": "caller_asserted",
            "reason": f"recorded pose output is unavailable: {pose_path}",
            "mode": 1,
            "affinity_kcal_mol": float(best_score),
        }
    recorded_hash = pose_record.get("sha256")
    actual_hash = sha256_file(pose_path)
    if not recorded_hash or str(recorded_hash).lower() != actual_hash.lower():
        raise ValueError("Score report pose-output SHA-256 does not match the current file")

    pose_atoms = [
        atom for atom in first_pdbqt_pose_atoms(pose_path) if atom["element"] != "H"
    ]
    if molecule.GetNumConformers() == 0:
        raise ValueError("Ligand SDF has no coordinates for score-pose verification")
    if len(pose_atoms) != molecule.GetNumAtoms():
        raise ValueError(
            "Vina mode 1 and ligand SDF have different heavy-atom counts: "
            f"{len(pose_atoms)} vs {molecule.GetNumAtoms()}"
        )
    conformer = molecule.GetConformer()
    available = set(range(molecule.GetNumAtoms()))
    distances: list[float] = []
    for pose_atom in pose_atoms:
        matches = []
        for index in available:
            rd_atom = molecule.GetAtomWithIdx(index)
            if rd_atom.GetSymbol().upper() != pose_atom["element"]:
                continue
            point = conformer.GetAtomPosition(index)
            distance = math.dist(
                pose_atom["coord_A"],
                (point.x, point.y, point.z),
            )
            if distance <= tolerance_A:
                matches.append((distance, index))
        if len(matches) != 1:
            raise ValueError(
                "Could not uniquely match a Vina mode-1 atom to the ligand SDF within "
                f"{tolerance_A} Å; matches={len(matches)}"
            )
        distance, index = matches[0]
        available.remove(index)
        distances.append(distance)
    return {
        "status": "coordinate_verified_mode_1",
        "mode": 1,
        "affinity_kcal_mol": float(best_score),
        "pose_output_file": {"path": str(pose_path.resolve()), "sha256": actual_hash},
        "atom_mapping_method": "unique element-and-3D-coordinate match",
        "tolerance_A": tolerance_A,
        "mapped_heavy_atoms": len(pose_atoms),
        "maximum_observed_distance_A": max(distances, default=0.0),
    }


def assign_label_centers(
    interactions: list[dict[str, Any]],
    target_points: dict[int, tuple[float, float]],
) -> dict[int, tuple[int, int]]:
    indexed = list(enumerate(interactions))
    indexed.sort(key=lambda item: target_points[item[0]][0])
    left_count = (len(indexed) + 1) // 2
    sides = (indexed[:left_count], indexed[left_count:])
    result: dict[int, tuple[int, int]] = {}
    for side_index, entries in enumerate(sides):
        ordered = sorted(entries, key=lambda item: target_points[item[0]][1])
        count = len(ordered)
        for position, (interaction_index, _) in enumerate(ordered):
            y = 160 if count == 1 else round(150 + position * 600 / (count - 1))
            result[interaction_index] = (120 if side_index == 0 else 1280, y)
    return result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--complex", type=Path, required=True)
    parser.add_argument("--ligand-sdf", type=Path, required=True)
    parser.add_argument("--verification-report", type=Path, required=True)
    parser.add_argument("--residue-map", type=Path)
    parser.add_argument("--score-report", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--output-prefix", default="docking_publication")
    parser.add_argument("--max-labels", type=int, default=6)
    parser.add_argument("--atom-map-tolerance-A", type=float, default=0.05)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    if not re.fullmatch(r"[A-Za-z0-9._-]+", args.output_prefix):
        parser.error("--output-prefix may contain only letters, digits, dot, underscore, hyphen")
    if args.max_labels < 1 or args.atom_map_tolerance_A <= 0:
        parser.error("--max-labels and --atom-map-tolerance-A must be positive")
    return args


def main() -> int:
    args = parse_args()
    if not args.ligand_sdf.is_file():
        raise FileNotFoundError(f"Ligand SDF does not exist: {args.ligand_sdf}")
    if args.score_report is not None and not args.score_report.is_file():
        raise FileNotFoundError(f"Score report does not exist: {args.score_report}")
    report = validate_verification_report(args.verification_report, args.complex)
    atoms = parse_pdb_atoms(args.complex)
    validate_report_atom_indices(report, atoms)
    residue_map, residue_map_metadata = load_residue_map(args.residue_map)
    ligand_key = ligand_residue_key(report)
    ligand_atoms = [atom for atom in atoms if residue_key(atom) == ligand_key]
    receptor_atoms = [atom for atom in atoms if residue_key(atom) != ligand_key]
    if not ligand_atoms or not receptor_atoms:
        raise ValueError("Could not split receptor and ligand from the verified complex")

    supplier = Chem.SDMolSupplier(str(args.ligand_sdf), removeHs=False)
    molecule_3d = supplier[0] if len(supplier) else None
    if molecule_3d is None:
        raise ValueError(f"RDKit could not parse {args.ligand_sdf}")
    molecule = Chem.RemoveHs(molecule_3d)
    pdb_heavy = [atom for atom in ligand_atoms if not atom["is_hydrogen"]]
    if molecule.GetNumAtoms() != len(pdb_heavy):
        raise ValueError(
            "Ligand SDF and verified complex have different heavy-atom counts: "
            f"{molecule.GetNumAtoms()} vs {len(pdb_heavy)}"
        )
    ligand_name_to_index, maximum_mapping_distance = coordinate_atom_map(
        molecule,
        pdb_heavy,
        args.atom_map_tolerance_A,
    )
    score_association = (
        verify_score_pose_association(
            args.score_report,
            molecule,
            args.atom_map_tolerance_A,
        )
        if args.score_report is not None
        else None
    )

    candidates = candidate_interactions(report)
    if len(candidates) > args.max_labels:
        raise ValueError(
            f"Verification report contains {len(candidates)} candidate-contact residues, "
            f"more than --max-labels={args.max_labels}"
        )
    candidate_keys = {item["source_residue"] for item in candidates}
    proximities = [
        item
        for item in residue_proximities(receptor_atoms, ligand_atoms)
        if item["source_residue"] not in candidate_keys
    ]
    interactions = candidates + proximities[: args.max_labels - len(candidates)]
    if not interactions:
        raise ValueError("No residues were available for the 2D interaction panel")
    if args.residue_map is not None:
        unmapped = sorted(
            {
                interaction["source_residue"]
                for interaction in interactions
                if interaction["source_residue"] not in residue_map
            }
        )
        if unmapped:
            labels = ", ".join(
                f"{chain or '-'}:{resn}{resi}" for chain, resn, resi in unmapped
            )
            raise ValueError(
                "The supplied residue map does not cover every displayed 2D residue: "
                f"{labels}. Add the mappings or reduce --max-labels; source and display "
                "numbering will not be mixed silently."
            )

    AllChem.Compute2DCoords(molecule)
    drawer = rdMolDraw2D.MolDraw2DCairo(*MOLECULE_SIZE)
    options = drawer.drawOptions()
    options.padding = 0.08
    options.bondLineWidth = 3.0
    options.minFontSize = 18
    options.maxFontSize = 28
    options.useMolBlockWedging = True
    drawer.DrawMolecule(molecule)
    drawer.FinishDrawing()
    molecule_image = Image.open(BytesIO(drawer.GetDrawingText())).convert("RGB")
    canvas = Image.new("RGB", CANVAS_SIZE, "white")
    canvas.paste(molecule_image, MOLECULE_ORIGIN)
    draw = ImageDraw.Draw(canvas)

    target_points: dict[int, tuple[float, float]] = {}
    for index, interaction in enumerate(interactions):
        ligand_atom_name = str(interaction["ligand_atom"]["name"])
        if ligand_atom_name not in ligand_name_to_index:
            raise ValueError(
                f"Interaction ligand atom is absent from the SDF/PDB atom map: "
                f"{ligand_atom_name}"
            )
        draw_point = drawer.GetDrawCoords(ligand_name_to_index[ligand_atom_name])
        target_points[index] = (
            MOLECULE_ORIGIN[0] + draw_point.x,
            MOLECULE_ORIGIN[1] + draw_point.y,
        )
    label_centers = assign_label_centers(interactions, target_points)

    reported_interactions = []
    for index, interaction in enumerate(interactions):
        candidate = interaction["kind"] == "candidate_polar_contact"
        line_color = "#D99500" if candidate else "#8A98A8"
        dashed_line(
            draw,
            label_centers[index],
            target_points[index],
            fill=line_color,
            width=4 if candidate else 3,
            dash=13 if candidate else 9,
            gap=8,
        )
        display = display_residue(interaction["source_residue"], residue_map)
        draw_badge(
            draw,
            label_centers[index],
            display["label"],
            f"{interaction['distance_A']:.1f} Å" if candidate else None,
            candidate=candidate,
        )
        reported_interactions.append(
            {
                "kind": interaction["kind"],
                "claim_status": "unverified" if candidate else "observed_proximity",
                "source_residue": {
                    "chain": interaction["source_residue"][0],
                    "resn": interaction["source_residue"][1],
                    "resi": interaction["source_residue"][2],
                },
                "display_residue": display,
                "distance_A": interaction["distance_A"],
                "ligand_atom": interaction["ligand_atom"],
                "receptor_atom": interaction["receptor_atom"],
                "source_atom_pair": (
                    f"{source_atom_label(interaction['ligand_atom'])} -- "
                    f"{source_atom_label(interaction['receptor_atom'])}"
                ),
            }
        )

    output_png = args.output_dir / f"{args.output_prefix}_contacts_2d.png"
    output_json = args.output_dir / f"{args.output_prefix}_contacts.json"
    ensure_targets([output_png, output_json], args.force)
    canvas.save(
        output_png,
        format="PNG",
        dpi=(PUBLICATION_DPI, PUBLICATION_DPI),
        optimize=True,
    )
    result = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "complex": {
                "path": str(args.complex.resolve()),
                "sha256": sha256_file(args.complex),
            },
            "ligand_sdf": {
                "path": str(args.ligand_sdf.resolve()),
                "sha256": sha256_file(args.ligand_sdf),
            },
            "verification_report": {
                "path": str(args.verification_report.resolve()),
                "sha256": sha256_file(args.verification_report),
            },
            "residue_map": residue_map_metadata,
            "score_report": (
                {
                    "path": str(args.score_report.resolve()),
                    "sha256": sha256_file(args.score_report),
                }
                if args.score_report is not None
                else None
            ),
        },
        "software": {"rdkit": rdkit.__version__},
        "output": {
            "path": str(output_png.resolve()),
            "pixels": list(CANVAS_SIZE),
            "mode": "RGB",
            "dpi_metadata": PUBLICATION_DPI,
        },
        "coordinate_handling": (
            "RDKit generated new 2D depiction coordinates on a molecule copy; source 3D "
            "docking coordinates were not modified. SDF/PDB atoms were matched uniquely by "
            "element and 3D coordinate before 2D depiction."
        ),
        "atom_mapping": {
            "method": "unique element-and-3D-coordinate match",
            "tolerance_A": args.atom_map_tolerance_A,
            "maximum_observed_distance_A": maximum_mapping_distance,
            "mapped_heavy_atoms": len(ligand_name_to_index),
        },
        "score_pose_association": score_association,
        "contact_policy": (
            "Gold is restricted to candidate_polar_contacts imported from the hash-matched "
            "PyMOL verification report. Gray lines are heavy-atom pocket proximities."
        ),
        "display_scope": {
            "kind": "curated_display_subset",
            "maximum_residue_labels": args.max_labels,
            "candidate_rule": "shortest reported candidate contact per receptor residue",
            "proximity_rule": (
                "nearest remaining receptor residues by minimum heavy-atom distance"
            ),
            "exhaustive_pocket_claim": False,
        },
        "interactions": reported_interactions,
        "output_png": str(output_png.resolve()),
    }
    write_json(output_json, result)
    print(f"RDKit {rdkit.__version__}")
    print(f"Wrote {output_png}")
    print(f"Wrote {output_json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
