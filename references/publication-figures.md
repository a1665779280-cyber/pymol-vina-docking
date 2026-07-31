# Publication-style docking figures

Use this workflow only after `pymol_verify_pose.py` has generated a verification report for the
exact combined complex. The publication renderer presents verified coordinates; it does not
validate biological binding or convert a Vina score into experimental affinity.

## Required inputs

| Input | Requirement |
|---|---|
| Combined complex PDB | The exact file audited in complex mode by `pymol_verify_pose.py`. |
| Verification report | Must contain `inputs.complex.sha256`, successful coordinate loading, atom counts, nearest-pair identity, and `candidate_polar_contacts`. |
| Docked ligand SDF | Same 3D pose and chemistry as the ligand in the complex. Atom order may differ; atoms are matched by element and 3D coordinate. |
| Residue map | Optional, but recommended when the complex uses local or renumbered residues. Every residue displayed in the 2D panel must be mapped if this file is supplied. |
| Parsed Vina score report | Optional. The score badge is omitted when no report is supplied. When the recorded multi-pose PDBQT is accessible, mode 1 is matched to the SDF by element and 3D coordinate. |

The renderer rejects a verification report whose stored complex SHA-256 does not match the
current PDB. It also requires every reported contact atom identity and 3D coordinate to resolve
uniquely in that PDB. PyMOL's recorded internal atom index is retained for session reproduction;
it is not assumed to equal the atom's ordinal PDB record number.

## Environments

The orchestrator itself uses only the Python standard library, but it deliberately calls two
separate environments:

```text
PyMOL rendering:  conda run -n pymol ...
RDKit + Pillow:   conda run -n rdkit-vis ...
```

Use different environment names with `--pymol-env` and `--chem-env`. Use `--conda-exe` when
Conda is not on `PATH`.

## Reference style and resolution

Use `assets/reference-style/publication_triptych_v1_reference.png` as the visual reference for
panel layout, hierarchy, colors, labels, and annotation density. It is an example artifact, not
scientific input, and its protein orientation must not be copied onto unrelated structures.

Publication PNGs use 600-DPI metadata. The final triptych remains 4000 x 2600 pixels, corresponding
to approximately 6.67 x 4.33 inches at 600 DPI. Increasing the metadata from 300 to 600 DPI does
not add pixels; use a larger canvas only when a larger physical print size is required.

## One-command workflow

Always choose a new, empty output directory:

```powershell
python scripts/make_publication_figure.py `
  --complex best_complex.pdb `
  --ligand-sdf best_pose.sdf `
  --verification-report pymol_verification\verification_report.json `
  --residue-map residue_map.json `
  --score-report vina_parsed.json `
  --protein-name EGFR `
  --ligand-name Erlotinib `
  --output-prefix EGFR_Erlotinib `
  --output-dir publication_run_001
```

Useful presentation parameters:

```text
--overview-transparency 0.32   protein-cartoon transparency in the overall view
--pocket-transparency 0.55     protein-cartoon transparency in the close-up
--pocket-cutoff 4.0            receptor residues considered around the ligand, Å
--max-pocket-residues 10       maximum residue sticks in the close-up
--max-labels 6                 maximum residue badges in the 2D panel
--atom-map-tolerance-A 0.05    SDF-to-complex coordinate matching tolerance, Å
--no-ray                       fast preview; keep ray tracing for final output
--dry-run                      validate paths and print resolved commands only
```

Direct worker scripts expose `--force`, but the orchestrator intentionally does not. A figure
run must not mix new intermediates with stale files.

## Residue-map schema

Keys use source chain, residue name, and full residue identifier. Include the insertion code in
`source_resi` or `display_resi`, for example `123A`.

```json
{
  "source_description": "chainless local numbering in the docking complex",
  "display_description": "native EGFR chain A numbering",
  "residues": [
    {
      "source_chain": "",
      "source_resn": "MET",
      "source_resi": "83",
      "display_chain": "A",
      "display_resn": "MET",
      "display_resi": "793",
      "display_label": "MET793"
    }
  ]
}
```

If no map is supplied, labels explicitly use the source-complex numbering. If a map is supplied,
the 2D stage fails rather than silently mixing source and mapped labels.

## Scientific and visual policy

- Gold dashed lines come only from the hash-matched verification report and are labeled
  **candidate polar contacts**, never confirmed hydrogen bonds.
- Gray dashed lines in the 2D panel are nearest-residue proximities, not bonds.
- The 2D panel is a curated display subset limited by `--max-labels`; it does not claim to list
  every pocket residue. If several candidates involve one residue, the shortest is displayed.
- A contact may be promoted to a hydrogen-bond claim only after donor/acceptor chemistry,
  protonation, bond order, donor-hydrogen geometry, and explicit distance/angle criteria are
  independently established. This workflow does not perform that promotion.
- The N-to-C rainbow is an aesthetic sequence-position cue, not a quantitative value scale.
- Ligand carbons are orange; heteroatom colors follow conventional element colors.
- “Transparent” refers to the semi-transparent protein cartoon. Final PNG files use an opaque
  white or very light background for predictable journal and slide rendering; do not describe
  them as transparent-background PNGs.
- RDKit creates a new 2D depiction on a molecule copy. The docked 3D coordinates are not changed.
- If the score report's hashed pose output is accessible, the workflow verifies that Vina mode 1
  and the displayed SDF have the same heavy-atom coordinates. Otherwise the badge says
  “Reported Vina score” and the manifest records the pose association as caller-asserted.

## Output contract

For prefix `docking_publication`, a successful run writes:

```text
output_dir/
  logs/
    render_3d.log
    draw_2d.log
    compose.log
  docking_publication_overall.png
  docking_publication_pocket.png
  docking_publication_ligand_mask.png
  docking_publication.pse
  docking_publication_pymol_manifest.json
  docking_publication_contacts_2d.png
  docking_publication_contacts.json
  docking_publication_triptych.png
  docking_publication_figure_manifest.json
  docking_publication_run_manifest.json
```

The manifests record input/output hashes, parameters, environment names, commands, residue and
atom mappings, candidate-contact policy, image metadata, alt text, and interpretation limits.
The final triptych must be RGB, 4000 x 2600 pixels, and carry 600-DPI metadata. Publication
overall, pocket, and 2D-interaction PNGs also carry 600-DPI metadata.

## Mandatory final check in PyMOL

Open the generated `.pse` in the real PyMOL GUI. Inspect both saved scenes and rotate the model.
Confirm:

1. the ligand is the intended pose and lies in the expected pocket;
2. residue numbering and chain context match the supplied map;
3. candidate contact dashes connect the reported atoms;
4. no relevant metal, cofactor, water, alternate location, or neighboring chain was hidden;
5. transparency and camera angle do not obscure steric clashes or misleadingly imply contacts.

Record that GUI inspection separately. A saved session proves reproducibility, not that a human
has completed the inspection.
