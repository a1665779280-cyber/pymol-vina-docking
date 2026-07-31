#!/usr/bin/env python
"""Parse an AutoDock Vina log into a provenance-focused JSON record."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import shlex
from typing import Any


FLOAT = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_text_auto(path: Path) -> str:
    data = path.read_bytes()
    if data.startswith((b"\xff\xfe", b"\xfe\xff")):
        return data.decode("utf-16")
    if data.startswith(b"\xef\xbb\xbf"):
        return data.decode("utf-8-sig")
    return data.decode("utf-8", errors="replace")


def first_match(pattern: str, text: str, cast: Any = str) -> Any:
    match = re.search(pattern, text, flags=re.MULTILINE | re.IGNORECASE)
    if not match:
        return None
    value = match.group(1).strip()
    return cast(value)


def vector_match(label: str, text: str) -> list[float] | None:
    match = re.search(
        rf"^{re.escape(label)}\s*:\s*X\s*({FLOAT})\s*Y\s*({FLOAT})\s*Z\s*({FLOAT})",
        text,
        flags=re.MULTILINE | re.IGNORECASE,
    )
    return [float(match.group(index)) for index in range(1, 4)] if match else None


def parse_command(command: str | None) -> dict[str, Any]:
    if not command:
        return {}
    try:
        tokens = shlex.split(command, posix=True)
    except ValueError:
        tokens = command.split()
    parsed: dict[str, Any] = {"executable": tokens[0] if tokens else None, "options": {}}
    index = 1
    while index < len(tokens):
        token = tokens[index]
        if token.startswith("--"):
            key_value = token[2:].split("=", 1)
            key = key_value[0].replace("-", "_")
            if len(key_value) == 2:
                value: Any = key_value[1]
            elif index + 1 < len(tokens) and not tokens[index + 1].startswith("--"):
                index += 1
                value = tokens[index]
            else:
                value = True
            if isinstance(value, str):
                try:
                    value = float(value) if any(char in value for char in ".eE") else int(value)
                except ValueError:
                    pass
            parsed["options"][key] = value
        index += 1
    return parsed


def parse_scores(text: str) -> list[dict[str, float | int]]:
    scores: list[dict[str, float | int]] = []
    row_pattern = re.compile(
        rf"^\s*(\d+)\s+({FLOAT})\s+({FLOAT})\s+({FLOAT})\s*$",
        flags=re.MULTILINE,
    )
    for match in row_pattern.finditer(text):
        scores.append(
            {
                "mode": int(match.group(1)),
                "affinity_kcal_mol": float(match.group(2)),
                "rmsd_lower_bound_A": float(match.group(3)),
                "rmsd_upper_bound_A": float(match.group(4)),
            }
        )
    return scores


def parse_config(path: Path | None) -> dict[str, Any]:
    if path is None:
        return {}
    options: dict[str, Any] = {}
    for raw_line in read_text_auto(path).splitlines():
        line = raw_line.split("#", 1)[0].strip()
        if not line or "=" not in line:
            continue
        key, value = (part.strip() for part in line.split("=", 1))
        normalized_key = key.replace("-", "_")
        parsed_value: Any = value
        try:
            parsed_value = float(value) if any(char in value for char in ".eE") else int(value)
        except ValueError:
            pass
        options[normalized_key] = parsed_value
    return {
        "path": str(path.resolve()),
        "sha256": sha256_file(path),
        "options": options,
    }


def parse_log(
    path: Path,
    config_path: Path | None = None,
    poses_path: Path | None = None,
) -> dict[str, Any]:
    text = read_text_auto(path)
    command_match = re.search(r"^\$\s+(.+)$", text, flags=re.MULTILINE)
    command = command_match.group(1).strip() if command_match else None
    scores = parse_scores(text)
    config = parse_config(config_path)
    config_options = config.get("options", {})
    if poses_path is None and config_options.get("out"):
        configured_path = Path(str(config_options["out"]))
        if not configured_path.is_absolute() and config_path is not None:
            configured_path = config_path.parent / configured_path
        poses_path = configured_path
    pose_output = None
    if poses_path is not None:
        resolved_pose_path = poses_path.resolve()
        pose_output = {
            "path": str(resolved_pose_path),
            "exists": resolved_pose_path.is_file(),
        }
        if resolved_pose_path.is_file():
            pose_output.update(
                {
                    "sha256": sha256_file(resolved_pose_path),
                    "bytes": resolved_pose_path.stat().st_size,
                }
            )
    result: dict[str, Any] = {
        "parsed_utc": datetime.now(timezone.utc).isoformat(),
        "log": {
            "path": str(path.resolve()),
            "sha256": sha256_file(path),
            "bytes": path.stat().st_size,
        },
        "software_banner": first_match(r"^(AutoDock Vina[^\r\n]*)", text),
        "scoring_function": first_match(r"^Scoring function\s*:\s*(\S+)", text),
        "receptor": first_match(r"^(?:Rigid receptor|Receptor)\s*:\s*(.+)$", text),
        "ligand": first_match(r"^Ligand\s*:\s*(.+)$", text),
        "grid": {
            "center_A": vector_match("Grid center", text),
            "size_A": vector_match("Grid size", text),
            "spacing_A": first_match(r"^Grid space\s*:\s*(" + FLOAT + r")", text, float),
        },
        "parameters": {
            "exhaustiveness": first_match(r"^Exhaustiveness\s*:\s*(\d+)", text, int),
            "cpu": first_match(r"^CPU\s*:\s*(\d+)", text, int),
            "verbosity": first_match(r"^Verbosity\s*:\s*(\d+)", text, int),
            "seed": first_match(
                r"Performing docking\s*\(random seed:\s*([-+]?\d+)\)", text, int
            ),
        },
        "command": command,
        "command_parsed": parse_command(command),
        "config": config or None,
        "pose_output_file": pose_output,
        "poses": scores,
        "best_affinity_kcal_mol": scores[0]["affinity_kcal_mol"] if scores else None,
        "warnings": [],
        "interpretation": (
            "Vina affinity is a model score for ranking under this protocol, not an "
            "experimental binding free energy."
        ),
    }
    command_options = result["command_parsed"].get("options", {})
    for key in ("num_modes", "energy_range", "seed", "exhaustiveness", "cpu", "min_rmsd"):
        value = command_options.get(key, config_options.get(key))
        if value is None:
            continue
        if result["parameters"].get(key) is None or key not in result["parameters"]:
            result["parameters"][key] = value

    if result["scoring_function"] is None:
        result["scoring_function"] = command_options.get(
            "scoring", config_options.get("scoring")
        )

    required = {
        "scoring function": result["scoring_function"],
        "grid center": result["grid"]["center_A"],
        "grid size": result["grid"]["size_A"],
        "seed": result["parameters"]["seed"],
        "pose table": scores,
    }
    for label, value in required.items():
        if value is None or value == []:
            result["warnings"].append(f"Could not parse {label} from the supplied log.")
    if command is None:
        result["warnings"].append(
            "The log does not contain the exact invocation command; preserve the command or "
            "a hashed config alongside future runs."
        )
    requested_modes = result["parameters"].get("num_modes")
    if requested_modes is not None and scores and int(requested_modes) != len(scores):
        result["warnings"].append(
            f"Requested {int(requested_modes)} modes but the score table contains "
            f"{len(scores)}; Vina may return fewer distinct poses, so inspect the output."
        )
    if pose_output is not None and not pose_output["exists"]:
        result["warnings"].append(
            "The referenced pose-output file was not found when the log was parsed."
        )
    if result["scoring_function"] and result["scoring_function"].lower() != "vina":
        result["warnings"].append(
            "Do not compare this raw score numerically with Vina scores from another scoring "
            "function."
        )
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("log", type=Path)
    parser.add_argument("--config", type=Path, help="Optional Vina config used for the run")
    parser.add_argument(
        "--poses",
        type=Path,
        help="Optional docked-pose output to hash (otherwise inferred from config out=)",
    )
    parser.add_argument("--output", type=Path)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    if not args.log.is_file():
        parser.error(f"Log file does not exist: {args.log}")
    if args.config is not None and not args.config.is_file():
        parser.error(f"Config file does not exist: {args.config}")
    if args.poses is not None and not args.poses.is_file():
        parser.error(f"Pose output file does not exist: {args.poses}")
    result = parse_log(args.log, args.config, args.poses)
    text = json.dumps(result, indent=2, ensure_ascii=False) + "\n"
    print(text)
    if args.output:
        if args.output.exists() and not args.force:
            parser.error(f"Refusing to overwrite {args.output}; pass --force")
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")
        print(f"Wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
