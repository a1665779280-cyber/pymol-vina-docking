"""Compose verified PyMOL and RDKit panels into a publication-style triptych."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import re

from PIL import Image, ImageDraw, ImageFont, ImageOps, __version__ as pillow_version

from publication_common import (
    PUBLICATION_DPI,
    ensure_targets,
    load_json,
    sha256_file,
    write_json,
)


CANVAS_SIZE = (4000, 2600)
BACKGROUND = "#F4F7FB"
INK = "#172033"
MUTED = "#5D6878"
PANEL_BORDER = "#D8E0EA"
NAVY = "#253A68"
LEFT_PANEL = (90, 270, 2470, 2320)
RIGHT_TOP_PANEL = (2550, 270, 3910, 1240)
RIGHT_BOTTOM_PANEL = (2550, 1340, 3910, 2320)


def font(bold: bool, size: int) -> ImageFont.FreeTypeFont:
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
    dash: int = 16,
    gap: int = 11,
) -> None:
    x1, y1 = start
    x2, y2 = end
    dx, dy = x2 - x1, y2 - y1
    length = (dx * dx + dy * dy) ** 0.5
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


def dashed_rectangle(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    *,
    fill: str,
    width: int,
) -> None:
    left, top, right, bottom = box
    dashed_line(draw, (left, top), (right, top), fill=fill, width=width)
    dashed_line(draw, (right, top), (right, bottom), fill=fill, width=width)
    dashed_line(draw, (right, bottom), (left, bottom), fill=fill, width=width)
    dashed_line(draw, (left, bottom), (left, top), fill=fill, width=width)


def draw_panel(
    canvas: Image.Image,
    box: tuple[int, int, int, int],
    *,
    letter: str,
    title: str,
    subtitle: str,
) -> None:
    draw = ImageDraw.Draw(canvas)
    left, top, right, bottom = box
    draw.rounded_rectangle(
        (left + 12, top + 16, right + 12, bottom + 16),
        radius=28,
        fill="#DFE5ED",
    )
    draw.rounded_rectangle(
        box,
        radius=28,
        fill="white",
        outline=PANEL_BORDER,
        width=3,
    )
    draw.rounded_rectangle(
        (left + 30, top + 28, left + 104, top + 102),
        radius=18,
        fill=NAVY,
    )
    letter_font = font(True, 46)
    letter_box = draw.textbbox((0, 0), letter, font=letter_font)
    letter_width = letter_box[2] - letter_box[0]
    draw.text(
        (left + 67 - letter_width / 2, top + 34),
        letter,
        font=letter_font,
        fill="white",
    )
    draw.text((left + 132, top + 30), title, font=font(True, 39), fill=INK)
    draw.text((left + 132, top + 76), subtitle, font=font(False, 25), fill=MUTED)


def rounded_pill(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    text: str,
    *,
    fill: str,
    outline: str,
    text_color: str,
    font_size: int,
) -> None:
    draw.rounded_rectangle(
        box,
        radius=(box[3] - box[1]) // 2,
        fill=fill,
        outline=outline,
        width=2,
    )
    pill_font = font(True, font_size)
    text_box = draw.textbbox((0, 0), text, font=pill_font)
    text_width = text_box[2] - text_box[0]
    text_height = text_box[3] - text_box[1]
    draw.text(
        (
            (box[0] + box[2] - text_width) / 2,
            (box[1] + box[3] - text_height) / 2 - text_box[1],
        ),
        text,
        font=pill_font,
        fill=text_color,
    )


def draw_sequence_gradient(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int]) -> None:
    stops = [
        (0, 32, 255),
        (0, 204, 255),
        (0, 229, 120),
        (159, 240, 0),
        (255, 227, 0),
        (255, 117, 0),
        (247, 24, 24),
    ]
    left, top, right, bottom = box
    width = right - left
    segment_count = len(stops) - 1
    for x in range(width):
        position = x / max(width - 1, 1) * segment_count
        segment = min(int(position), segment_count - 1)
        fraction = position - segment
        c1, c2 = stops[segment], stops[segment + 1]
        color = tuple(
            round(c1[channel] + (c2[channel] - c1[channel]) * fraction)
            for channel in range(3)
        )
        draw.line((left + x, top, left + x, bottom), fill=color)
    draw.rounded_rectangle(box, radius=8, outline="#9AA5B1", width=2)


def trim_white_margin(
    image: Image.Image,
    padding: int,
) -> tuple[Image.Image, tuple[int, int, int, int]]:
    grayscale = ImageOps.grayscale(image.convert("RGB"))
    content_mask = grayscale.point(lambda value: 255 if value < 248 else 0)
    content_box = content_mask.getbbox()
    if content_box is None:
        return image, (0, 0, image.width, image.height)
    left = max(0, content_box[0] - padding)
    top = max(0, content_box[1] - padding)
    right = min(image.width, content_box[2] + padding)
    bottom = min(image.height, content_box[3] + padding)
    crop_box = (left, top, right, bottom)
    return image.crop(crop_box), crop_box


def mask_bbox(path: Path, padding: int = 65) -> tuple[int, int, int, int]:
    mask = Image.open(path).convert("L")
    binary = mask.point(lambda value: 255 if value > 32 else 0)
    box = binary.getbbox()
    if box is None:
        raise ValueError("Ligand mask contains no foreground pixels")
    return (
        max(0, box[0] - padding),
        max(0, box[1] - padding),
        min(mask.width, box[2] + padding),
        min(mask.height, box[3] + padding),
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--overview", type=Path, required=True)
    parser.add_argument("--pocket", type=Path, required=True)
    parser.add_argument("--ligand-mask", type=Path, required=True)
    parser.add_argument("--contacts-2d", type=Path, required=True)
    parser.add_argument("--pymol-manifest", type=Path, required=True)
    parser.add_argument("--contact-report", type=Path, required=True)
    parser.add_argument("--score-report", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--output-prefix", default="docking_publication")
    parser.add_argument("--protein-name", default="Protein")
    parser.add_argument("--ligand-name", default="Ligand")
    parser.add_argument("--title")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    if not re.fullmatch(r"[A-Za-z0-9._-]+", args.output_prefix):
        parser.error("--output-prefix may contain only letters, digits, dot, underscore, hyphen")
    for path in (
        args.overview,
        args.pocket,
        args.ligand_mask,
        args.contacts_2d,
        args.pymol_manifest,
        args.contact_report,
    ):
        if not path.is_file():
            parser.error(f"Input does not exist: {path}")
    if args.score_report is not None and not args.score_report.is_file():
        parser.error(f"Score report does not exist: {args.score_report}")
    return args


def main() -> int:
    args = parse_args()
    pymol_manifest = load_json(args.pymol_manifest)
    contact_report = load_json(args.contact_report)
    pymol_hash = pymol_manifest["inputs"]["complex"]["sha256"]
    contact_hash = contact_report["inputs"]["complex"]["sha256"]
    if str(pymol_hash).lower() != str(contact_hash).lower():
        raise ValueError("3D and 2D stage manifests refer to different complex hashes")
    pymol_verification_hash = pymol_manifest["inputs"]["verification_report"]["sha256"]
    contact_verification_hash = contact_report["inputs"]["verification_report"]["sha256"]
    if str(pymol_verification_hash).lower() != str(contact_verification_hash).lower():
        raise ValueError("3D and 2D stages used different verification reports")

    overview_original = Image.open(args.overview).convert("RGB")
    pocket_original = Image.open(args.pocket).convert("RGB")
    interaction_original = Image.open(args.contacts_2d).convert("RGB")
    overview, overview_crop = trim_white_margin(overview_original, 70)
    pocket, _ = trim_white_margin(pocket_original, 35)
    interaction_2d, _ = trim_white_margin(interaction_original, 25)

    canvas = Image.new("RGB", CANVAS_SIZE, BACKGROUND)
    draw = ImageDraw.Draw(canvas)
    title = args.title or f"{args.protein_name}–{args.ligand_name} docking visualization"
    draw.text((105, 57), title, font=font(True, 72), fill=INK)
    draw.text(
        (108, 148),
        "Ray-traced PyMOL views and an RDKit-derived 2D contact map from a verified pose",
        font=font(False, 32),
        fill=MUTED,
    )

    score_record = None
    if args.score_report is not None:
        score_data = load_json(args.score_report)
        score_value = score_data.get("best_affinity_kcal_mol")
        if score_value is None:
            raise ValueError("Score report has no best_affinity_kcal_mol")
        contact_score_input = contact_report.get("inputs", {}).get("score_report")
        if not isinstance(contact_score_input, dict):
            association = {
                "status": "caller_asserted",
                "reason": "2D stage did not receive the score report",
            }
        elif str(contact_score_input.get("sha256", "")).lower() != sha256_file(
            args.score_report
        ).lower():
            raise ValueError("2D stage and composer used different score reports")
        else:
            association = contact_report.get("score_pose_association") or {
                "status": "caller_asserted",
                "reason": "score-pose association was not recorded",
            }
        score_record = {
            "path": str(args.score_report.resolve()),
            "sha256": sha256_file(args.score_report),
            "best_affinity_kcal_mol": float(score_value),
            "pose_association": association,
        }
        score_prefix = (
            "Vina score"
            if association.get("status") == "coordinate_verified_mode_1"
            else "Reported Vina score"
        )
        rounded_pill(
            draw,
            (3140, 72, 3890, 154),
            f"{score_prefix}  {float(score_value):.3f} kcal/mol",
            fill="#EAF1FF",
            outline="#9EB5E4",
            text_color=NAVY,
            font_size=29,
        )

    draw_panel(
        canvas,
        LEFT_PANEL,
        letter="A",
        title="Overall binding pose",
        subtitle=(
            f"{args.protein_name} colored from N to C terminus • "
            f"{args.ligand_name} shown in orange"
        ),
    )
    draw_panel(
        canvas,
        RIGHT_TOP_PANEL,
        letter="B",
        title="Binding-pocket close-up",
        subtitle="Gray sticks: pocket residues • gold dashes: candidate polar contacts",
    )
    draw_panel(
        canvas,
        RIGHT_BOTTOM_PANEL,
        letter="C",
        title="2D interaction summary",
        subtitle="Gold: candidate polar contact • gray: nearest-residue proximity",
    )

    overview_display = ImageOps.contain(overview, (2160, 1880), Image.Resampling.LANCZOS)
    overview_origin = (
        LEFT_PANEL[0] + (LEFT_PANEL[2] - LEFT_PANEL[0] - overview_display.width) // 2,
        LEFT_PANEL[1] + 120,
    )
    canvas.paste(overview_display, overview_origin)
    pocket_display = ImageOps.contain(pocket, (1290, 830), Image.Resampling.LANCZOS)
    pocket_origin = (
        RIGHT_TOP_PANEL[0]
        + (RIGHT_TOP_PANEL[2] - RIGHT_TOP_PANEL[0] - pocket_display.width) // 2,
        RIGHT_TOP_PANEL[1] + 122,
    )
    canvas.paste(pocket_display, pocket_origin)
    interaction_display = ImageOps.contain(
        interaction_2d,
        (1300, 835),
        Image.Resampling.LANCZOS,
    )
    interaction_origin = (
        RIGHT_BOTTOM_PANEL[0]
        + (RIGHT_BOTTOM_PANEL[2] - RIGHT_BOTTOM_PANEL[0] - interaction_display.width) // 2,
        RIGHT_BOTTOM_PANEL[1] + 126,
    )
    canvas.paste(interaction_display, interaction_origin)

    source_ligand_box = mask_bbox(args.ligand_mask)
    overview_scale = overview_display.width / overview.width
    ligand_box = (
        round(overview_origin[0] + (source_ligand_box[0] - overview_crop[0]) * overview_scale),
        round(overview_origin[1] + (source_ligand_box[1] - overview_crop[1]) * overview_scale),
        round(overview_origin[0] + (source_ligand_box[2] - overview_crop[0]) * overview_scale),
        round(overview_origin[1] + (source_ligand_box[3] - overview_crop[1]) * overview_scale),
    )
    ligand_box = (
        max(LEFT_PANEL[0] + 20, ligand_box[0]),
        max(LEFT_PANEL[1] + 115, ligand_box[1]),
        min(LEFT_PANEL[2] - 20, ligand_box[2]),
        min(LEFT_PANEL[3] - 20, ligand_box[3]),
    )
    dashed_rectangle(draw, ligand_box, fill="#697586", width=4)
    rounded_pill(
        draw,
        (ligand_box[0], ligand_box[1] - 54, ligand_box[0] + 218, ligand_box[1] - 8),
        "binding site",
        fill="#FFFFFF",
        outline="#A7B0BC",
        text_color="#475467",
        font_size=22,
    )
    dashed_line(
        draw,
        (ligand_box[2], ligand_box[1]),
        (RIGHT_TOP_PANEL[0], RIGHT_TOP_PANEL[1] + 250),
        fill="#A0A9B5",
        width=3,
        dash=12,
        gap=10,
    )
    dashed_line(
        draw,
        (ligand_box[2], ligand_box[3]),
        (RIGHT_TOP_PANEL[0], RIGHT_TOP_PANEL[3] - 110),
        fill="#A0A9B5",
        width=3,
        dash=12,
        gap=10,
    )

    legend_y = 2392
    draw_sequence_gradient(draw, (120, legend_y + 4, 360, legend_y + 34))
    draw.text((382, legend_y - 2), "Sequence (N→C)", font=font(False, 27), fill=INK)
    draw.ellipse((870, legend_y, 904, legend_y + 34), fill="#FA6408")
    draw.text((922, legend_y - 2), args.ligand_name, font=font(False, 27), fill=INK)
    draw.line((1210, legend_y + 17, 1274, legend_y + 17), fill="#AEB3BB", width=9)
    draw.text((1290, legend_y - 2), "Pocket residues", font=font(False, 27), fill=INK)
    dashed_line(
        draw,
        (1740, legend_y + 17),
        (1812, legend_y + 17),
        fill="#D99500",
        width=6,
        dash=12,
        gap=8,
    )
    draw.text(
        (1830, legend_y - 2),
        "Candidate polar contact",
        font=font(False, 27),
        fill=INK,
    )
    dashed_line(
        draw,
        (2420, legend_y + 17),
        (2492, legend_y + 17),
        fill="#8A98A8",
        width=4,
        dash=9,
        gap=7,
    )
    draw.text(
        (2510, legend_y - 2),
        "Pocket proximity (2D)",
        font=font(False, 27),
        fill=INK,
    )

    numbering_description = contact_report["inputs"]["residue_map"].get(
        "display_description",
        "source complex numbering",
    )
    note = (
        "Distances come from the supplied docked coordinates. Gold marks a "
        "geometry-screened candidate polar contact, not a confirmed hydrogen bond. "
        f"Labels use {numbering_description}."
    )
    note_font = font(False, 24)
    note_box = draw.textbbox((0, 0), note, font=note_font)
    note_width = note_box[2] - note_box[0]
    draw.text(
        ((CANVAS_SIZE[0] - note_width) / 2, 2506),
        note,
        font=note_font,
        fill="#667085",
    )

    output_path = args.output_dir / f"{args.output_prefix}_triptych.png"
    manifest_path = args.output_dir / f"{args.output_prefix}_figure_manifest.json"
    ensure_targets([output_path, manifest_path], args.force)
    canvas.save(
        output_path,
        format="PNG",
        dpi=(PUBLICATION_DPI, PUBLICATION_DPI),
        optimize=True,
    )
    candidates = [
        item
        for item in contact_report.get("interactions", [])
        if item.get("kind") == "candidate_polar_contact"
    ]
    candidate_text = ", ".join(
        f"{item['display_residue']['label']} at {item['distance_A']:.1f} Å"
        for item in candidates
    ) or "no candidate polar contacts"
    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "overview": {"path": str(args.overview.resolve()), "sha256": sha256_file(args.overview)},
            "pocket": {"path": str(args.pocket.resolve()), "sha256": sha256_file(args.pocket)},
            "ligand_mask": {"path": str(args.ligand_mask.resolve()), "sha256": sha256_file(args.ligand_mask)},
            "contacts_2d": {"path": str(args.contacts_2d.resolve()), "sha256": sha256_file(args.contacts_2d)},
            "pymol_manifest": {"path": str(args.pymol_manifest.resolve()), "sha256": sha256_file(args.pymol_manifest)},
            "contact_report": {"path": str(args.contact_report.resolve()), "sha256": sha256_file(args.contact_report)},
            "score_report": score_record,
        },
        "complex_sha256": pymol_hash,
        "verification_report_sha256": pymol_verification_hash,
        "software": {"pillow": pillow_version},
        "output": {
            "path": str(output_path.resolve()),
            "pixels": list(CANVAS_SIZE),
            "mode": "RGB",
            "dpi_metadata": PUBLICATION_DPI,
        },
        "whole_image_operations": [
            "Lossless source PNGs were resized with Lanczos resampling for layout.",
            "Panel frames, titles, legend, and mask-derived zoom guides were added.",
            "Source molecular coordinates and measured distances were not modified.",
        ],
        "alt_text": (
            f"Three-panel {args.protein_name}-{args.ligand_name} docking figure. "
            f"Panel A shows a semi-transparent N-to-C sequence-colored protein with "
            f"{args.ligand_name} in orange. Panel B shows the 3D pocket and gold "
            f"candidate polar contacts. Panel C shows an RDKit 2D contact map with "
            f"{candidate_text}; other displayed residues are gray proximities."
        ),
        "limitations": (
            "This is a presentation of verified coordinates and geometric candidates; it "
            "does not establish binding, affinity, efficacy, or confirmed hydrogen bonds."
        ),
    }
    write_json(manifest_path, manifest)
    print(f"Wrote {output_path}")
    print(f"Wrote {manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
