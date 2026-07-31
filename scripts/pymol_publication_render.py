"""Render transparent publication-style PyMOL views from a verified complex."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import re
from typing import Any

from pymol import cmd

from publication_common import (
    PUBLICATION_DPI,
    display_residue,
    ensure_targets,
    ligand_residue_key,
    load_residue_map,
    parse_pdb_atoms,
    residue_key,
    sha256_file,
    validate_report_atom_indices,
    validate_verification_report,
    write_json,
)


def configure_rendering() -> None:
    cmd.bg_color("white")
    cmd.set("orthoscopic", 1)
    cmd.set("depth_cue", 0)
    cmd.set("ray_trace_fog", 0)
    cmd.set("ray_shadows", 0)
    cmd.set("ray_opaque_background", 1)
    cmd.set("antialias", 2)
    cmd.set("ray_trace_mode", 0)
    cmd.set("ambient", 0.52)
    cmd.set("direct", 0.70)
    cmd.set("specular", 0.16)
    cmd.set("shininess", 16)
    cmd.set("reflect", 0.14)
    cmd.set("two_sided_lighting", 1)
    cmd.set("cartoon_fancy_helices", 1)
    cmd.set("cartoon_smooth_loops", 1)
    cmd.set("cartoon_flat_sheets", 1)
    cmd.set("cartoon_loop_radius", 0.20)
    cmd.set("cartoon_oval_length", 1.25)
    cmd.set("cartoon_oval_width", 0.26)
    cmd.set("dash_radius", 0.052)
    cmd.set("dash_length", 0.20)
    cmd.set("dash_gap", 0.13)

    cmd.set_color("ligand_orange", [0.98, 0.39, 0.02])
    cmd.set_color("pocket_gray", [0.67, 0.69, 0.73])
    cmd.set_color("pocket_light", [0.78, 0.79, 0.82])
    cmd.set_color("nitrogen_blue", [0.12, 0.30, 0.82])
    cmd.set_color("nitrogen_blue_soft", [0.42, 0.54, 0.84])
    cmd.set_color("oxygen_red", [0.88, 0.12, 0.12])
    cmd.set_color("oxygen_red_soft", [0.91, 0.42, 0.42])
    cmd.set_color("sulfur_gold", [0.95, 0.67, 0.05])
    cmd.set_color("halogen_green", [0.18, 0.58, 0.28])
    cmd.set_color("contact_gold", [0.90, 0.59, 0.00])


def color_ligand() -> None:
    cmd.color("ligand_orange", "ligand and elem C")
    cmd.color("nitrogen_blue", "ligand and elem N")
    cmd.color("oxygen_red", "ligand and elem O")
    cmd.color("sulfur_gold", "ligand and elem S")
    cmd.color("halogen_green", "ligand and elem F+CL+BR+I")


def color_pocket(selection: str, *, soft: bool = False) -> None:
    cmd.color("pocket_light" if soft else "pocket_gray", f"{selection} and elem C")
    cmd.color(
        "nitrogen_blue_soft" if soft else "nitrogen_blue",
        f"{selection} and elem N",
    )
    cmd.color(
        "oxygen_red_soft" if soft else "oxygen_red",
        f"{selection} and elem O",
    )
    cmd.color("sulfur_gold", f"{selection} and elem S")


def set_common_representation() -> None:
    cmd.hide("everything", "all")
    cmd.show("cartoon", "receptor and polymer.protein")
    cmd.spectrum("count", "rainbow", "receptor and polymer.protein")
    cmd.show("sticks", "ligand and not elem H")
    cmd.show("spheres", "ligand and not elem H")
    cmd.set("stick_radius", 0.26, "ligand")
    cmd.set("sphere_scale", 0.24, "ligand")
    color_ligand()


def residue_distances(
    receptor_atoms: list[Any],
    ligand_atoms: list[Any],
) -> tuple[dict[tuple[str, str, str], float], dict[tuple[str, str, str], int]]:
    minima: dict[tuple[str, str, str], float] = {}
    representative: dict[tuple[str, str, str], int] = {}
    for receptor_atom in receptor_atoms:
        if receptor_atom.symbol == "H":
            continue
        key = (receptor_atom.chain or "", receptor_atom.resn, receptor_atom.resi)
        representative.setdefault(key, receptor_atom.index)
        nearest = min(
            math.dist(receptor_atom.coord, ligand_atom.coord)
            for ligand_atom in ligand_atoms
            if ligand_atom.symbol != "H"
        )
        minima[key] = min(minima.get(key, math.inf), nearest)
    return minima, representative


def choose_display_residues(
    report: dict[str, Any],
    minima: dict[tuple[str, str, str], float],
    cutoff: float,
    maximum: int,
) -> list[tuple[str, str, str]]:
    ranked = [
        key
        for key, distance in sorted(minima.items(), key=lambda item: item[1])
        if distance <= cutoff
    ]
    selected = ranked[:maximum]
    candidate_keys = [
        residue_key(contact["receptor_atom"])
        for contact in report["geometry"]["candidate_polar_contacts"]
    ]
    for key in candidate_keys:
        if key not in selected:
            selected.append(key)
    return selected


def render_mask(path: Path, view: tuple[float, ...], width: int, height: int) -> None:
    cmd.hide("everything", "all")
    cmd.bg_color("black")
    cmd.show("spheres", "ligand and not elem H")
    cmd.color("white", "ligand and not elem H")
    cmd.set("sphere_scale", 0.32, "ligand")
    cmd.set_view(view)
    cmd.png(
        str(path),
        width=width,
        height=height,
        dpi=PUBLICATION_DPI,
        ray=0,
        quiet=1,
    )
    cmd.bg_color("white")


def contact_objects(report: dict[str, Any]) -> list[dict[str, Any]]:
    rendered: list[dict[str, Any]] = []
    for number, contact in enumerate(
        report["geometry"]["candidate_polar_contacts"], start=1
    ):
        ligand_atom = contact["ligand_atom"]
        receptor_atom = contact["receptor_atom"]
        name = f"candidate_polar_contact_{number}"
        cmd.distance(
            name,
            f"complex_object and index {int(ligand_atom['index'])}",
            f"complex_object and index {int(receptor_atom['index'])}",
        )
        cmd.color("contact_gold", name)
        cmd.set("dash_color", "contact_gold", name)
        cmd.set("dash_width", 3.0, name)
        cmd.hide("labels", name)
        rendered.append(
            {
                "object": name,
                "distance_A": float(contact["distance_A"]),
                "ligand_atom": ligand_atom,
                "receptor_atom": receptor_atom,
                "classification": (
                    "geometry-screened candidate polar contact; not a confirmed hydrogen bond"
                ),
            }
        )
    return rendered


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--complex", type=Path, required=True)
    parser.add_argument("--verification-report", type=Path, required=True)
    parser.add_argument("--residue-map", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--output-prefix", default="docking_publication")
    parser.add_argument("--protein-name", default="Protein")
    parser.add_argument("--ligand-name", default="Ligand")
    parser.add_argument("--overview-transparency", type=float, default=0.32)
    parser.add_argument("--pocket-transparency", type=float, default=0.55)
    parser.add_argument("--pocket-cutoff", type=float, default=4.0)
    parser.add_argument("--max-pocket-residues", type=int, default=10)
    parser.add_argument(
        "--overview-turn",
        nargs=3,
        type=float,
        default=(-8.0, 18.0, -6.0),
        metavar=("X", "Y", "Z"),
    )
    parser.add_argument(
        "--pocket-turn",
        nargs=3,
        type=float,
        default=(18.0, -15.0, -12.0),
        metavar=("X", "Y", "Z"),
    )
    parser.add_argument(
        "--ray",
        action=argparse.BooleanOptionalAction,
        default=True,
    )
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    if not re.fullmatch(r"[A-Za-z0-9._-]+", args.output_prefix):
        parser.error("--output-prefix may contain only letters, digits, dot, underscore, hyphen")
    for name in ("overview_transparency", "pocket_transparency"):
        value = getattr(args, name)
        if not 0.0 <= value < 1.0:
            parser.error(f"--{name.replace('_', '-')} must be in [0, 1)")
    if args.pocket_cutoff <= 0 or args.max_pocket_residues < 1:
        parser.error("Pocket cutoff and maximum residue count must be positive")
    return args


def main() -> int:
    args = parse_args()
    if not args.complex.is_file():
        raise FileNotFoundError(f"Complex does not exist: {args.complex}")
    report = validate_verification_report(args.verification_report, args.complex)
    source_atoms = parse_pdb_atoms(args.complex)
    validate_report_atom_indices(report, source_atoms)
    residue_map, residue_map_metadata = load_residue_map(args.residue_map)

    prefix = args.output_prefix
    overview_path = args.output_dir / f"{prefix}_overall.png"
    pocket_path = args.output_dir / f"{prefix}_pocket.png"
    mask_path = args.output_dir / f"{prefix}_ligand_mask.png"
    session_path = args.output_dir / f"{prefix}.pse"
    manifest_path = args.output_dir / f"{prefix}_pymol_manifest.json"
    ensure_targets(
        [overview_path, pocket_path, mask_path, session_path, manifest_path],
        args.force,
    )

    cmd.reinitialize()
    configure_rendering()
    cmd.load(str(args.complex), "complex_object")
    nearest_ligand_index = int(
        report["geometry"]["nearest_pair"]["ligand_atom"]["index"]
    )
    cmd.select(
        "ligand",
        f"byres (complex_object and index {nearest_ligand_index})",
    )
    cmd.select("receptor", "complex_object and not ligand")
    receptor_count = cmd.count_atoms("receptor")
    ligand_count = cmd.count_atoms("ligand")
    expected_receptor = int(report["selection"]["receptor_atoms"])
    expected_ligand = int(report["selection"]["ligand_atoms"])
    if (receptor_count, ligand_count) != (expected_receptor, expected_ligand):
        raise RuntimeError(
            "Loaded atom counts do not match the verification report: "
            f"observed receptor/ligand {receptor_count}/{ligand_count}, expected "
            f"{expected_receptor}/{expected_ligand}"
        )

    receptor_atoms = cmd.get_model("receptor and not elem H", state=1).atom
    ligand_atoms = cmd.get_model("ligand and not elem H", state=1).atom
    minima, representative = residue_distances(receptor_atoms, ligand_atoms)
    selected_keys = choose_display_residues(
        report,
        minima,
        args.pocket_cutoff,
        args.max_pocket_residues,
    )
    cmd.select("display_residues", "none")
    for key in selected_keys:
        index = representative.get(key)
        if index is None:
            raise RuntimeError(f"Selected residue is absent from PyMOL receptor: {key}")
        cmd.select(
            "display_residues",
            f"display_residues or byres (complex_object and index {index})",
        )
    cmd.select(
        "pocket",
        f"byres (receptor within {args.pocket_cutoff:.3f} of ligand)",
    )

    set_common_representation()
    cmd.set("cartoon_transparency", args.overview_transparency, "receptor")
    cmd.set("sphere_scale", 0.27, "ligand")
    cmd.show("sticks", "pocket and not elem H")
    cmd.set("stick_radius", 0.095, "pocket")
    color_pocket("pocket", soft=True)
    cmd.orient("receptor")
    for axis, angle in zip(("x", "y", "z"), args.overview_turn):
        cmd.turn(axis, angle)
    cmd.zoom("receptor", 3.5)
    overview_view = cmd.get_view()
    cmd.scene("overview", "store")
    cmd.png(
        str(overview_path),
        width=2000,
        height=2000,
        dpi=PUBLICATION_DPI,
        ray=int(args.ray),
        quiet=1,
    )
    render_mask(mask_path, overview_view, 2000, 2000)

    set_common_representation()
    cmd.set("cartoon_transparency", args.pocket_transparency, "receptor")
    cmd.show("sticks", "display_residues and not elem H")
    cmd.set("stick_radius", 0.17, "display_residues")
    color_pocket("display_residues")
    rendered_contacts = contact_objects(report)
    cmd.orient("ligand")
    for axis, angle in zip(("x", "y", "z"), args.pocket_turn):
        cmd.turn(axis, angle)
    cmd.zoom("ligand or display_residues", 2.8)
    cmd.clip("slab", 28)
    cmd.scene("binding_site", "store")
    cmd.png(
        str(pocket_path),
        width=1800,
        height=1200,
        dpi=PUBLICATION_DPI,
        ray=int(args.ray),
        quiet=1,
    )
    cmd.save(str(session_path))

    selected_report = []
    for key in selected_keys:
        display = display_residue(key, residue_map)
        selected_report.append(
            {
                "source": {"chain": key[0], "resn": key[1], "resi": key[2]},
                "display": display,
                "nearest_ligand_heavy_atom_distance_A": minima.get(key),
            }
        )
    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "complex": {
                "path": str(args.complex.resolve()),
                "sha256": sha256_file(args.complex),
            },
            "verification_report": {
                "path": str(args.verification_report.resolve()),
                "sha256": sha256_file(args.verification_report),
            },
            "residue_map": residue_map_metadata,
        },
        "software": {"pymol": cmd.get_version()[0]},
        "labels": {"protein": args.protein_name, "ligand": args.ligand_name},
        "atom_counts": {"receptor": receptor_count, "ligand": ligand_count},
        "render": {
            "dpi_metadata": PUBLICATION_DPI,
            "overview_cartoon_transparency": args.overview_transparency,
            "pocket_cartoon_transparency": args.pocket_transparency,
            "pocket_cutoff_A": args.pocket_cutoff,
            "ray_traced": bool(args.ray),
            "protein_color": (
                "N-to-C rainbow sequence gradient; aesthetic sequence cue, not a quantitative scale"
            ),
            "ligand_carbons": "orange",
            "candidate_polar_contacts": "gold dashed lines",
        },
        "selected_pocket_residues": selected_report,
        "candidate_polar_contacts": rendered_contacts,
        "warnings": [
            "Gold dashes are candidate polar contacts, not confirmed hydrogen bonds."
        ],
        "outputs": {
            "overview_png": str(overview_path.resolve()),
            "pocket_png": str(pocket_path.resolve()),
            "ligand_mask_png": str(mask_path.resolve()),
            "pymol_session": str(session_path.resolve()),
        },
    }
    write_json(manifest_path, manifest)
    print(f"Wrote {overview_path}")
    print(f"Wrote {pocket_path}")
    print(f"Wrote {session_path}")
    print(f"Wrote {manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
