"""Shared validation and provenance helpers for docking publication figures."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any, Iterable


PUBLICATION_DPI = 600


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise FileNotFoundError(f"JSON file does not exist: {path}")
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ValueError(f"Expected a JSON object in {path}")
    return data


def write_json(path: Path, value: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )


def ensure_targets(targets: Iterable[Path], force: bool) -> None:
    paths = list(targets)
    existing = [path for path in paths if path.exists()]
    if existing and not force:
        names = ", ".join(str(path) for path in existing)
        raise FileExistsError(f"Refusing to overwrite existing output(s): {names}; pass --force")
    for path in paths:
        path.parent.mkdir(parents=True, exist_ok=True)


def atom_element_from_pdb(line: str) -> str:
    element = line[76:78].strip() if len(line) >= 78 else ""
    if element:
        return element.upper()
    name = line[12:16].strip().lstrip("0123456789")
    two_letter = name[:2].upper()
    return two_letter if two_letter in {"BR", "CL"} else name[:1].upper()


def parse_pdb_atoms(path: Path) -> list[dict[str, Any]]:
    if not path.is_file():
        raise FileNotFoundError(f"Coordinate file does not exist: {path}")
    atoms: list[dict[str, Any]] = []
    saw_model = False
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            record = line[:6].strip().upper()
            if record == "MODEL":
                if saw_model:
                    break
                saw_model = True
                continue
            if record == "ENDMDL" and saw_model:
                break
            if record not in {"ATOM", "HETATM"} or len(line) < 54:
                continue
            try:
                coord = [
                    float(line[30:38]),
                    float(line[38:46]),
                    float(line[46:54]),
                ]
            except ValueError as exc:
                raise ValueError(f"Invalid PDB coordinate record: {line.rstrip()}") from exc
            element = atom_element_from_pdb(line)
            residue_number = line[22:26].strip()
            insertion_code = line[26:27].strip()
            atoms.append(
                {
                    "index": len(atoms) + 1,
                    "record": record,
                    "name": line[12:16].strip(),
                    "resn": line[17:20].strip(),
                    "chain": line[21:22].strip(),
                    "resi": f"{residue_number}{insertion_code}",
                    "insertion_code": insertion_code,
                    "altloc": line[16:17].strip(),
                    "element": element,
                    "is_hydrogen": element == "H",
                    "coord_A": coord,
                }
            )
    if not atoms:
        raise ValueError(f"No ATOM/HETATM coordinates found in {path}")
    return atoms


def residue_key(atom: dict[str, Any]) -> tuple[str, str, str]:
    return str(atom.get("chain", "")), str(atom["resn"]), str(atom["resi"])


def source_atom_label(atom: dict[str, Any]) -> str:
    chain = atom.get("chain") or "-"
    return f"{chain}/{atom['resn']}{atom['resi']}/{atom['name']}"


def validate_verification_report(
    report_path: Path,
    complex_path: Path,
) -> dict[str, Any]:
    report = load_json(report_path)
    if report.get("source_mode") != "complex":
        raise ValueError(
            "Publication rendering currently requires a complex-mode verification report"
        )
    if report.get("coordinate_load_verified") is not True:
        raise ValueError("Verification report does not assert successful coordinate loading")
    complex_record = report.get("inputs", {}).get("complex")
    if not isinstance(complex_record, dict) or not complex_record.get("sha256"):
        raise ValueError("Verification report does not contain inputs.complex.sha256")
    actual_hash = sha256_file(complex_path)
    if str(complex_record["sha256"]).lower() != actual_hash.lower():
        raise ValueError(
            "Verification report/input hash mismatch; regenerate verification for this exact "
            "complex before rendering"
        )
    geometry = report.get("geometry")
    if not isinstance(geometry, dict):
        raise ValueError("Verification report is missing geometry")
    candidates = geometry.get("candidate_polar_contacts")
    if not isinstance(candidates, list):
        raise ValueError("Verification report is missing candidate_polar_contacts")
    return report


def ligand_residue_key(report: dict[str, Any]) -> tuple[str, str, str]:
    nearest = report.get("geometry", {}).get("nearest_pair", {}).get("ligand_atom")
    if not isinstance(nearest, dict):
        raise ValueError("Verification report does not identify the ligand residue")
    return residue_key(nearest)


def validate_report_atom_indices(
    report: dict[str, Any],
    atoms: list[dict[str, Any]],
) -> None:
    """Resolve every reported contact atom uniquely by identity and coordinate.

    A PyMOL ``atom.index`` is an internal load index and is not necessarily the same as the
    atom's ordinal record in a PDB file. The exact-file hash plus identity/coordinate check is
    therefore the portable integrity test; the recorded PyMOL index remains available to
    reproduce the distance object after the same file is loaded in PyMOL.
    """
    for contact in report["geometry"]["candidate_polar_contacts"]:
        for role in ("ligand_atom", "receptor_atom"):
            recorded = contact.get(role)
            if not isinstance(recorded, dict) or "index" not in recorded:
                raise ValueError(f"Candidate contact is missing {role}.index")
            expected = (
                str(recorded.get("chain", "")),
                str(recorded.get("resn", "")),
                str(recorded.get("resi", "")),
                str(recorded.get("name", "")),
            )
            candidates = [
                atom
                for atom in atoms
                if (
                    atom["chain"],
                    atom["resn"],
                    atom["resi"],
                    atom["name"],
                )
                == expected
            ]
            coordinate = recorded.get("coord_A")
            if not isinstance(coordinate, list) or len(coordinate) != 3:
                raise ValueError(f"Candidate contact is missing {role}.coord_A")
            coordinate_matches = [
                atom
                for atom in candidates
                if math.dist(atom["coord_A"], [float(value) for value in coordinate])
                <= 0.01
            ]
            if len(coordinate_matches) != 1:
                raise ValueError(
                    "Reported atom identity/coordinate does not resolve uniquely in the exact "
                    f"PDB: {expected}; matches={len(coordinate_matches)}"
                )
            expected_element = str(recorded.get("element", "")).upper()
            if expected_element and coordinate_matches[0]["element"] != expected_element:
                raise ValueError(
                    f"Reported element mismatch for {expected}: expected {expected_element}, "
                    f"observed {coordinate_matches[0]['element']}"
                )


def load_residue_map(
    path: Path | None,
) -> tuple[dict[tuple[str, str, str], dict[str, str]], dict[str, Any]]:
    if path is None:
        return {}, {
            "source_description": "source complex numbering",
            "display_description": "source complex numbering",
            "path": None,
        }
    data = load_json(path)
    entries = data.get("residues")
    if not isinstance(entries, list):
        raise ValueError("Residue-map JSON requires a residues array")
    mapping: dict[tuple[str, str, str], dict[str, str]] = {}
    for entry in entries:
        if not isinstance(entry, dict):
            raise ValueError("Each residue-map entry must be a JSON object")
        required = ("source_resn", "source_resi", "display_resn", "display_resi")
        missing = [key for key in required if not str(entry.get(key, "")).strip()]
        if missing:
            raise ValueError(f"Residue-map entry is missing: {', '.join(missing)}")
        key = (
            str(entry.get("source_chain", "")),
            str(entry["source_resn"]),
            str(entry["source_resi"]),
        )
        if key in mapping:
            raise ValueError(f"Duplicate residue-map source key: {key}")
        mapping[key] = {
            "chain": str(entry.get("display_chain", "")),
            "resn": str(entry["display_resn"]),
            "resi": str(entry["display_resi"]),
            "label": str(
                entry.get(
                    "display_label",
                    f"{entry['display_resn']}{entry['display_resi']}",
                )
            ),
        }
    metadata = {
        "path": str(path.resolve()),
        "sha256": sha256_file(path),
        "source_description": str(
            data.get("source_description", "source complex numbering")
        ),
        "display_description": str(
            data.get("display_description", "mapped residue numbering")
        ),
        "entries": len(mapping),
    }
    return mapping, metadata


def display_residue(
    key: tuple[str, str, str],
    mapping: dict[tuple[str, str, str], dict[str, str]],
) -> dict[str, str]:
    if key in mapping:
        return mapping[key]
    chain, resn, resi = key
    return {
        "chain": chain,
        "resn": resn,
        "resi": resi,
        "label": f"{resn}{resi}",
    }
