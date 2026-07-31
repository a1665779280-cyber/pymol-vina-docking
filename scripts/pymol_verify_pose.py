#!/usr/bin/env python
"""Load a docked complex in PyMOL and emit coordinate-level verification artifacts."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
from pymol import cmd


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atom_record(atom: Any) -> dict[str, Any]:
    return {
        "object": atom.model,
        "index": atom.index,
        "chain": atom.chain,
        "resn": atom.resn,
        "resi": atom.resi,
        "name": atom.name,
        "element": atom.symbol,
        "coord_A": [float(value) for value in atom.coord],
    }


def atom_label(record: dict[str, Any]) -> str:
    chain = record["chain"] or "-"
    return (
        f"{record['object']}:{chain}/{record['resn']}{record['resi']}/{record['name']}"
    )


def create_box(center: list[float], size: list[float]) -> None:
    half = [value / 2.0 for value in size]
    x0, y0, z0 = [center[i] - half[i] for i in range(3)]
    x1, y1, z1 = [center[i] + half[i] for i in range(3)]
    vertices = [
        (x0, y0, z0),
        (x1, y0, z0),
        (x1, y1, z0),
        (x0, y1, z0),
        (x0, y0, z1),
        (x1, y0, z1),
        (x1, y1, z1),
        (x0, y1, z1),
    ]
    edges = (
        (0, 1),
        (1, 2),
        (2, 3),
        (3, 0),
        (4, 5),
        (5, 6),
        (6, 7),
        (7, 4),
        (0, 4),
        (1, 5),
        (2, 6),
        (3, 7),
    )
    cmd.delete("vina_search_box")
    for index, vertex in enumerate(vertices, start=1):
        cmd.pseudoatom(
            "vina_search_box",
            name=f"B{index}",
            elem="C",
            pos=vertex,
        )
    for start, end in edges:
        cmd.bond(
            f"vina_search_box and name B{start + 1}",
            f"vina_search_box and name B{end + 1}",
        )
    cmd.show("sticks", "vina_search_box")
    cmd.set("stick_radius", 0.075, "vina_search_box")
    cmd.color("marine", "vina_search_box")


def residue_sort_key(item: tuple[str, str, str]) -> tuple[str, int, str, str]:
    chain, resn, resi = item
    digits = "".join(character for character in resi if character.isdigit() or character == "-")
    suffix = resi[len(digits) :] if digits else resi
    try:
        number = int(digits)
    except ValueError:
        number = 10**9
    return chain, number, suffix, resn


def pair_metrics(
    receptor_atoms: list[Any],
    ligand_atoms: list[Any],
    clash_cutoff: float,
    contact_cutoff: float,
) -> dict[str, Any]:
    receptor_xyz = np.asarray([atom.coord for atom in receptor_atoms], dtype=float)
    ligand_xyz = np.asarray([atom.coord for atom in ligand_atoms], dtype=float)
    deltas = ligand_xyz[:, None, :] - receptor_xyz[None, :, :]
    distances = np.linalg.norm(deltas, axis=2)
    nearest_flat = int(np.argmin(distances))
    ligand_index, receptor_index = np.unravel_index(nearest_flat, distances.shape)
    nearest_distance = float(distances[ligand_index, receptor_index])

    clash_indices = np.argwhere(distances < clash_cutoff)
    severe_indices = np.argwhere(distances < 1.5)
    contact_indices = np.argwhere(
        (distances >= clash_cutoff) & (distances <= contact_cutoff)
    )
    clash_pairs = []
    sorted_clashes = sorted(
        (
            (float(distances[lig_i, rec_i]), int(lig_i), int(rec_i))
            for lig_i, rec_i in clash_indices
        ),
        key=lambda item: item[0],
    )
    for distance, lig_i, rec_i in sorted_clashes[:50]:
        ligand_record = atom_record(ligand_atoms[lig_i])
        receptor_record = atom_record(receptor_atoms[rec_i])
        clash_pairs.append(
            {
                "distance_A": distance,
                "ligand_atom": ligand_record,
                "receptor_atom": receptor_record,
                "label": f"{atom_label(ligand_record)} -- {atom_label(receptor_record)}",
            }
        )
    return {
        "nearest_heavy_atom_distance_A": nearest_distance,
        "nearest_pair": {
            "ligand_atom": atom_record(ligand_atoms[ligand_index]),
            "receptor_atom": atom_record(receptor_atoms[receptor_index]),
        },
        "heavy_atom_pairs_below_clash_cutoff": int(len(clash_indices)),
        "heavy_atom_pairs_below_1_5_A": int(len(severe_indices)),
        "heavy_atom_pairs_between_clash_and_contact_cutoffs": int(len(contact_indices)),
        "clash_pairs_first_50": clash_pairs,
    }


def candidate_polar_pairs(atom_lookup: dict[tuple[str, int], Any]) -> list[dict[str, Any]]:
    pairs: list[tuple[tuple[str, int], tuple[str, int], str, str]] = []
    for ligand_role, receptor_role in (("donor", "acceptor"), ("acceptor", "donor")):
        try:
            found = cmd.find_pairs(
                f"ligand and {ligand_role}",
                f"receptor and {receptor_role}",
                cutoff=3.6,
                mode=1,
                angle=55.0,
            )
            pairs.extend(
                (left, right, ligand_role, receptor_role)
                for left, right in found
            )
        except Exception:
            continue
    output = []
    seen: set[tuple[tuple[str, int], tuple[str, int]]] = set()
    for left, right, ligand_role, receptor_role in pairs:
        key = ((left[0], left[1]), (right[0], right[1]))
        if key in seen:
            continue
        seen.add(key)
        ligand_atom = atom_lookup.get((left[0], left[1]))
        receptor_atom = atom_lookup.get((right[0], right[1]))
        if ligand_atom is None or receptor_atom is None:
            continue
        distance = math.dist(ligand_atom.coord, receptor_atom.coord)
        ligand_record = atom_record(ligand_atom)
        receptor_record = atom_record(receptor_atom)
        output.append(
            {
                "distance_A": distance,
                "kind": "candidate_polar_contact",
                "claim_status": "unverified",
                "ligand_role": ligand_role,
                "receptor_role": receptor_role,
                "ligand_atom": ligand_record,
                "receptor_atom": receptor_record,
                "label": f"{atom_label(ligand_record)} -- {atom_label(receptor_record)}",
            }
        )
    return sorted(output, key=lambda item: item["distance_A"])


def configure_visuals() -> None:
    cmd.bg_color("white")
    cmd.set("orthoscopic", 1)
    cmd.set("depth_cue", 0)
    cmd.set("ray_trace_fog", 0)
    cmd.set("ray_shadows", 0)
    cmd.set("ray_opaque_background", 1)
    cmd.set("antialias", 2)
    cmd.set("ambient", 0.55)
    cmd.set("direct", 0.70)
    cmd.set("specular", 0.15)
    cmd.set("two_sided_lighting", 1)
    cmd.set("cartoon_fancy_helices", 1)
    cmd.set("cartoon_smooth_loops", 1)
    cmd.set("dash_radius", 0.05)
    cmd.set("dash_length", 0.22)
    cmd.set("dash_gap", 0.14)
    cmd.set_color("protein_blue", [0.48, 0.66, 0.88])
    cmd.set_color("ligand_orange", [0.98, 0.39, 0.02])
    cmd.set_color("pocket_gray", [0.68, 0.70, 0.74])
    cmd.set_color("contact_gold", [0.98, 0.70, 0.03])


def prepare_representation() -> None:
    cmd.hide("everything", "all")
    cmd.show("cartoon", "receptor and polymer.protein")
    cmd.color("protein_blue", "receptor and polymer.protein")
    cmd.set("cartoon_transparency", 0.42, "receptor and polymer.protein")
    cmd.select("pocket4", "byres (receptor within 4.0 of ligand)")
    cmd.show("sticks", "pocket4 and not elem H")
    cmd.set("stick_radius", 0.14, "pocket4")
    cmd.color("pocket_gray", "pocket4 and elem C")
    cmd.color("blue", "pocket4 and elem N")
    cmd.color("red", "pocket4 and elem O")
    cmd.color("yellow", "pocket4 and elem S")
    cmd.show("sticks", "ligand and not elem H")
    cmd.set("stick_radius", 0.24, "ligand")
    cmd.color("ligand_orange", "ligand and elem C")
    cmd.color("blue", "ligand and elem N")
    cmd.color("red", "ligand and elem O")
    cmd.color("yellow", "ligand and elem S")
    cmd.color("green", "ligand and elem F+CL+BR+I")
    for name, ligand_role, receptor_role in (
        ("ligand_donor_contacts", "donor", "acceptor"),
        ("ligand_acceptor_contacts", "acceptor", "donor"),
    ):
        try:
            cmd.distance(
                name,
                f"ligand and {ligand_role}",
                f"receptor and {receptor_role}",
                cutoff=3.6,
                mode=2,
            )
            cmd.hide("labels", name)
            cmd.set("dash_color", "contact_gold", name)
        except Exception:
            continue


def load_inputs(args: argparse.Namespace) -> tuple[dict[str, Any], str]:
    inputs: dict[str, Any] = {}
    if args.complex:
        cmd.load(str(args.complex), "docked_complex")
        cmd.select("ligand", f"docked_complex and ({args.ligand_selection})")
        cmd.select("receptor", "docked_complex and not ligand")
        inputs["complex"] = {
            "path": str(args.complex.resolve()),
            "sha256": sha256_file(args.complex),
        }
        source_mode = "complex"
    else:
        cmd.load(str(args.receptor), "receptor_object")
        cmd.load(str(args.ligand), "ligand_object")
        cmd.select("receptor", "receptor_object")
        cmd.select("ligand", "ligand_object")
        inputs["receptor"] = {
            "path": str(args.receptor.resolve()),
            "sha256": sha256_file(args.receptor),
        }
        inputs["ligand"] = {
            "path": str(args.ligand.resolve()),
            "sha256": sha256_file(args.ligand),
        }
        source_mode = "separate receptor and ligand"
    return inputs, source_mode


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--complex", type=Path)
    source.add_argument("--receptor", type=Path)
    parser.add_argument("--ligand", type=Path)
    parser.add_argument("--ligand-selection", default="resn LIG")
    parser.add_argument("--center", nargs=3, type=float, metavar=("X", "Y", "Z"))
    parser.add_argument("--size", nargs=3, type=float, metavar=("SX", "SY", "SZ"))
    parser.add_argument("--clash-cutoff", type=float, default=2.0)
    parser.add_argument("--contact-cutoff", type=float, default=4.0)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--ray", action="store_true", help="Ray trace PNGs")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    if args.complex:
        if args.ligand:
            parser.error("--ligand cannot be used with --complex")
        input_paths = [args.complex]
    else:
        if args.ligand is None:
            parser.error("--ligand is required with --receptor")
        input_paths = [args.receptor, args.ligand]
    for path in input_paths:
        if path is None or not path.is_file():
            parser.error(f"Input file does not exist: {path}")
    if (args.center is None) != (args.size is None):
        parser.error("--center and --size must be supplied together")
    if args.size and any(value <= 0 for value in args.size):
        parser.error("All box dimensions must be positive")
    if args.clash_cutoff <= 0 or args.contact_cutoff <= args.clash_cutoff:
        parser.error("--contact-cutoff must be greater than positive --clash-cutoff")
    if args.output_dir.exists() and any(args.output_dir.iterdir()) and not args.force:
        parser.error(f"Output directory is not empty: {args.output_dir}; pass --force")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    cmd.reinitialize()
    configure_visuals()
    inputs, source_mode = load_inputs(args)
    receptor_count = cmd.count_atoms("receptor")
    ligand_count = cmd.count_atoms("ligand")
    if receptor_count == 0:
        raise RuntimeError("PyMOL loaded no receptor atoms")
    if ligand_count == 0:
        raise RuntimeError(
            f"PyMOL loaded no ligand atoms for selection: {args.ligand_selection}"
        )
    ligand_states = cmd.count_states("ligand")
    receptor_protein_count = cmd.count_atoms("receptor and polymer.protein")
    receptor_context_count = cmd.count_atoms("receptor and not polymer.protein")
    receptor_context_resnames = sorted(
        {
            atom.resn
            for atom in cmd.get_model(
                "receptor and not polymer.protein and not elem H", state=1
            ).atom
        }
    )

    receptor_model = cmd.get_model("receptor and not elem H", state=1)
    ligand_model = cmd.get_model("ligand and not elem H", state=1)
    receptor_atoms = receptor_model.atom
    ligand_atoms = ligand_model.atom
    if not receptor_atoms or not ligand_atoms:
        raise RuntimeError("Heavy-atom selections are empty")
    lookup = {
        (atom.model, atom.index): atom
        for atom in cmd.get_model("receptor or ligand", state=1).atom
    }

    geometry = pair_metrics(
        receptor_atoms,
        ligand_atoms,
        args.clash_cutoff,
        args.contact_cutoff,
    )
    polar_pairs = candidate_polar_pairs(lookup)
    cmd.select("pocket4", "byres (receptor within 4.0 of ligand)")
    pocket_residues = sorted(
        {
            (atom.chain or "-", atom.resn, atom.resi)
            for atom in cmd.get_model("pocket4 and name CA", state=1).atom
        },
        key=residue_sort_key,
    )

    box_report = None
    if args.center and args.size:
        center = list(args.center)
        size = list(args.size)
        half = np.asarray(size, dtype=float) / 2.0
        ligand_xyz = np.asarray([atom.coord for atom in ligand_atoms], dtype=float)
        outside_mask = np.any(
            (ligand_xyz < np.asarray(center) - half)
            | (ligand_xyz > np.asarray(center) + half),
            axis=1,
        )
        create_box(center, size)
        box_report = {
            "center_A": center,
            "size_A": size,
            "ligand_heavy_atoms_outside": int(outside_mask.sum()),
            "all_ligand_heavy_atoms_inside": not bool(outside_mask.any()),
        }

    prepare_representation()
    if args.center and args.size:
        create_box(list(args.center), list(args.size))
    overview_path = args.output_dir / "overview.png"
    pocket_path = args.output_dir / "pocket.png"
    session_path = args.output_dir / "verification.pse"
    report_path = args.output_dir / "verification_report.json"

    cmd.orient("receptor")
    cmd.zoom("receptor", 4.0)
    cmd.scene("overview", "store")
    cmd.png(
        str(overview_path),
        width=1800,
        height=1400,
        dpi=300,
        ray=int(args.ray),
        quiet=1,
    )
    if args.center and args.size:
        cmd.disable("vina_search_box")
    cmd.orient("ligand")
    cmd.zoom("ligand or pocket4", 3.0)
    cmd.scene("pocket", "store")
    cmd.png(
        str(pocket_path),
        width=1800,
        height=1400,
        dpi=300,
        ray=int(args.ray),
        quiet=1,
    )
    if args.center and args.size:
        cmd.enable("vina_search_box")
    cmd.save(str(session_path))

    warnings: list[str] = []
    if ligand_states > 1:
        warnings.append(
            f"The ligand has {ligand_states} coordinate states; geometry, contacts, and images "
            "in this report analyze state 1 only. Inspect and cluster the remaining poses."
        )
    if receptor_context_count:
        warnings.append(
            f"The receptor includes {receptor_context_count} non-protein atoms "
            f"({', '.join(receptor_context_resnames) or 'unlabeled'}). They are included in "
            "distance checks; confirm whether waters, metals, or cofactors were intentionally "
            "retained."
        )
    if geometry["heavy_atom_pairs_below_1_5_A"]:
        warnings.append(
            "One or more receptor–ligand heavy-atom pairs are below 1.5 Å; inspect for severe "
            "steric overlap or coordinate/bonding artifacts."
        )
    elif geometry["heavy_atom_pairs_below_clash_cutoff"]:
        warnings.append(
            "One or more receptor–ligand heavy-atom pairs are below the selected clash cutoff; "
            "inspect atom types, protonation, and local geometry."
        )
    if box_report and not box_report["all_ligand_heavy_atoms_inside"]:
        warnings.append("At least one ligand heavy atom lies outside the supplied search box.")
    if not polar_pairs:
        warnings.append(
            "No candidate polar contacts were identified by PyMOL's donor/acceptor geometry "
            "heuristic; atom typing or missing hydrogens may affect this result."
        )

    report = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "software": {"pymol": cmd.get_version()[0], "numpy": np.__version__},
        "source_mode": source_mode,
        "inputs": inputs,
        "selection": {
            "ligand_selection": args.ligand_selection if args.complex else "ligand_object",
            "receptor_atoms": receptor_count,
            "receptor_heavy_atoms": len(receptor_atoms),
            "receptor_protein_atoms": receptor_protein_count,
            "receptor_nonprotein_atoms": receptor_context_count,
            "receptor_nonprotein_resnames": receptor_context_resnames,
            "ligand_atoms": ligand_count,
            "ligand_heavy_atoms": len(ligand_atoms),
            "ligand_states": ligand_states,
            "analyzed_ligand_state": 1,
        },
        "box": box_report,
        "geometry": {
            "clash_cutoff_A": args.clash_cutoff,
            "contact_cutoff_A": args.contact_cutoff,
            "candidate_polar_contact_method": {
                "engine": "PyMOL cmd.find_pairs",
                "distance_cutoff_A": 3.6,
                "mode": 1,
                "angle_parameter_deg": 55.0,
                "claim_policy": (
                    "Geometry-screened candidates only; protonation, bond order, and explicit "
                    "hydrogen geometry are not fully established."
                ),
            },
            **geometry,
            "candidate_polar_contacts": polar_pairs,
            "pocket_residues_within_4_A": [
                {"chain": chain, "resn": resn, "resi": resi}
                for chain, resn, resi in pocket_residues
            ],
        },
        "coordinate_load_verified": True,
        "gui_rotation_verified": False,
        "warnings": warnings,
        "outputs": {
            "overview_png": str(overview_path.resolve()),
            "pocket_png": str(pocket_path.resolve()),
            "pymol_session": str(session_path.resolve()),
        },
        "interpretation": (
            "This verifies coordinate loading and geometric plausibility only. It does not prove "
            "binding, affinity, efficacy, or experimental validity."
        ),
    }
    report_path.write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(f"Wrote {report_path}")
    print(f"Wrote {session_path}")
    print(f"Wrote {overview_path}")
    print(f"Wrote {pocket_path}")
    for warning in warnings:
        print(f"WARNING: {warning}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
