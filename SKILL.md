---
name: pymol-vina-docking
description: Plan, run, audit, and visualize small-molecule protein docking with AutoDock Vina, PyMOL, and RDKit. Use for choosing and documenting search-box and Vina parameters; validating receptor, ligand, config, log, and pose files; checking docked poses in PyMOL; assessing redocking and multi-seed convergence; producing reproducible docking evidence packages; and creating reference-matched 600-DPI transparent-protein, publication-style overall, pocket, 2D-interaction, triptych, and editable PyMOL-session outputs.
---

# PyMOL + AutoDock Vina docking

Use AutoDock Vina as the docking engine and PyMOL as a coordinate-level inspection and
visualization tool. Never describe PyMOL inspection or a favorable Vina score as proof of
binding, efficacy, or experimental free energy.

## Required integrity rules

- Preserve all supplied structures and logs. Write new work into a timestamped or named run
  directory.
- Record file hashes, software versions, the complete command/configuration, scoring function,
  random seed, and preparation decisions.
- Treat protonation, tautomer, stereochemistry, metals, cofactors, conserved waters, missing
  residues, alternate locations, and chain selection as scientific decisions.
- Compare scores only within the same preparation and scoring protocol. Do not compare Vina,
  Vinardo, and AD4 scores as if they shared one energy scale.
- Call a dashed line a candidate polar contact unless donor/acceptor identity, distance,
  protonation, and geometry have been checked.
- Use multiple seeds and, when a cognate ligand exists, redocking. One best-scoring pose is not
  a validated result.

## Read-on-demand references

- Read [references/input-files.md](references/input-files.md) before preparing or accepting
  docking inputs.
- Read [references/parameter-selection.md](references/parameter-selection.md) whenever choosing
  the search box, scoring function, exhaustiveness, number of poses, energy range, seed, or
  receptor flexibility.
- Read [references/validation.md](references/validation.md) before interpreting or reporting a
  docking result.
- Read [references/publication-figures.md](references/publication-figures.md) when asked for a
  polished, transparent-protein, publication-style, multi-panel, or editable PyMOL figure.
- Inspect [assets/reference-style/publication_triptych_v1_reference.png](assets/reference-style/publication_triptych_v1_reference.png)
  with an image viewer when matching or auditing the publication style. Use it only as a visual
  reference for layout, color, hierarchy, and annotation density; never as structural evidence.

## Default docking visualization contract

- Treat requests to display, visualize, render, export, or generate an image for a docking result
  as requests for the `publication_triptych_v1` workflow unless the user explicitly requests a
  ligand-only image or a basic verification image.
- Use a verification report whose stored SHA-256 matches the exact combined complex. If only
  separate receptor and ligand files exist, first create and verify a combined complex.
- Produce overall, pocket, and RDKit-derived 2D interaction panels; compose them into a
  4000 x 2600 RGB triptych with 600-DPI metadata; and retain the editable PyMOL session and
  manifests.
- Return the triptych as the primary image. Do not substitute a bare SDF render for a docking
  visualization.

## Workflow

### 1. Inventory and identify the scientific question

State whether the task is redocking, known-site docking, virtual screening, or exploratory blind
docking. Record:

- receptor source/PDB ID, biological assembly, chain(s), and preparation pH;
- ligand identity, stereochemistry, protonation/tautomer state, and preparation source;
- binding-site evidence: cognate ligand, experimentally supported residues, validated pocket
  predictor, or unknown;
- retained metals, cofactors, and structural waters with rationale.

Refuse to infer a binding site from the final docked pose itself; that is circular.

### 2. Audit the input contract

For a normal Vina run require:

- docking-ready rigid receptor: `receptor.pdbqt`;
- docking-ready ligand: `ligand.pdbqt`;
- source/prepared coordinate files for human inspection (`.pdb`/`.cif`, `.sdf`/`.mol2`);
- an explicit box center and size in Angstrom;
- an output directory that does not overwrite inputs.

Run the environment check:

```powershell
python scripts/check_environment.py `
  --vina-path "E:\Pymol\tools\autodock-vina-1.2.7\vina.exe" `
  --json environment.json
```

On this Windows installation, use:

```text
PyMOL environment:  conda run -n pymol ...
Preparation env:    conda run -n rdkit-vis ...
Vina executable:    E:\Pymol\tools\autodock-vina-1.2.7\vina.exe
```

The preparation environment contains RDKit and Meeko. Verify these paths and versions rather
than assuming they remain unchanged.

### 3. Choose parameters from evidence

Use the decision process in `references/parameter-selection.md`. Prefer a cognate-ligand or
known-residue box. Compute box dimensions as:

```text
size_axis = selected_site_extent_axis + 2 * padding_per_side
```

Start with 4–6 Å padding per side for a small-molecule known-site box, then inspect it in PyMOL.
Do not silently use a generic 20 or 25 Å cube.

Generate a config and a machine-readable plan:

```powershell
python scripts/make_vina_config.py `
  --receptor receptor.pdbqt `
  --ligand ligand.pdbqt `
  --reference cognate_ligand.pdb `
  --padding 5 `
  --scoring vina `
  --exhaustiveness 16 `
  --num-modes 10 `
  --energy-range 3 `
  --seed 2026 `
  --output config.txt `
  --plan docking_plan.json
```

If no reference ligand is available, pass explicit, independently justified values with
`--center X Y Z --size SX SY SZ --box-rationale "..."`. The config generator refuses an
explicit box without this rationale; never use the resulting docked pose as its own justification.

### 4. Run and preserve the exact execution record

For the Vina CLI:

```powershell
$vina = "E:\Pymol\tools\autodock-vina-1.2.7\vina.exe"
& $vina --config config.txt 2>&1 | Tee-Object -FilePath vina.log
```

Vina 1.2.7 writes the run report to standard output and does not expose a `--log` option; capture
it with `Tee-Object` (PowerShell) or `tee` (POSIX). For platforms where the Python bindings are
installed, set the seed in the constructor:

```python
from vina import Vina

v = Vina(sf_name="vina", cpu=0, seed=2026)
v.set_receptor("receptor.pdbqt")
v.set_ligand_from_file("ligand.pdbqt")
v.compute_vina_maps(center=[x, y, z], box_size=[sx, sy, sz])
v.dock(exhaustiveness=16, n_poses=10)
v.write_poses("docked_poses.pdbqt", n_poses=10, energy_range=3.0, overwrite=False)
```

Do not pass `seed` to `v.dock`; the official API defines it on `Vina(...)`.

Parse and preserve the log:

```powershell
python scripts/parse_vina_log.py vina.log --config config.txt `
  --poses docked_poses.pdbqt --output vina_parsed.json
```

### 5. Verify the actual coordinates in PyMOL

Create a coordinate audit, session, and ray-traced image:

```powershell
conda run -n pymol python scripts/pymol_verify_pose.py `
  --complex best_complex.pdb `
  --ligand-selection "resn LIG" `
  --center X Y Z `
  --size SX SY SZ `
  --output-dir pymol_verification `
  --ray
```

Alternatively pass `--receptor prepared_receptor.pdb --ligand best_pose.sdf`. The script reports:

- loaded atom/state counts and input hashes;
- ligand atoms outside the requested box;
- nearest receptor–ligand heavy-atom distance and short-contact counts;
- pocket residues within 4 Å;
- candidate polar-contact pairs;
- PyMOL version and generated `.pse`, `.png`, and `.json` paths.

Open the produced `.pse` in the real GUI and rotate the structure. Confirm the box, pocket,
ligand orientation, close contacts, residue numbering, metal/cofactor context, and alternate
poses. The saved image is evidence that the coordinates loaded; it is not biological validation.

### 6. Create a verified publication-style figure

Only render contact annotations from a verification report whose stored SHA-256 matches the exact
combined complex. Use a docked SDF for the 2D ligand depiction and an explicit residue map when
native numbering differs from the complex. Run the reusable cross-environment workflow:

```powershell
python scripts/make_publication_figure.py `
  --complex best_complex.pdb `
  --ligand-sdf best_pose.sdf `
  --verification-report pymol_verification\verification_report.json `
  --residue-map residue_map.json `
  --score-report vina_parsed.json `
  --protein-name Protein `
  --ligand-name Ligand `
  --output-dir publication_run_001
```

The default overall and pocket protein-cartoon transparencies are 0.32 and 0.55. The ligand,
displayed pocket residues, and candidate-contact dashes remain opaque. Final PNG backgrounds are
opaque and publication-safe; “transparent” here means protein-cartoon transparency. Use a new,
empty output directory for every run, then open the generated `.pse` in the PyMOL GUI and inspect
both saved scenes. Follow `references/publication-figures.md` for inputs, parameters, residue-map
schema, output files, and scientific labeling rules.
The final triptych is 4000 x 2600 with 600-DPI metadata. Compare its visual style with
`assets/reference-style/publication_triptych_v1_reference.png`.

### 7. Validate before interpretation

Use `references/validation.md`. At minimum:

1. Repeat docking with at least three recorded seeds.
2. Compare top-pose clusters, not only top scores.
3. If a cognate ligand exists, redock it and report symmetry-aware heavy-atom RMSD and method.
4. Inspect clashes, strained ligand geometry, buried unsatisfied polar groups, and expected
   catalytic/hinge interactions.
5. Report limitations and distinguish Vina score from experimental affinity.

## Expected evidence package

Keep these together:

```text
run/
  inputs/                 immutable source and prepared structures
  config.txt
  docking_plan.json
  environment.json
  vina.log
  vina_parsed.json
  docked_poses.pdbqt
  pymol_verification/
    verification_report.json
    verification.pse
    overview.png
    pocket.png
  publication_run_001/
    docking_publication_overall.png
    docking_publication_pocket.png
    docking_publication_contacts_2d.png
    docking_publication_triptych.png
    docking_publication.pse
    docking_publication_run_manifest.json
```

## Upstream and authoritative sources

This skill was adapted from the CC-BY-4.0
`jaechang-hits/SciAgent-Skills/autodock-vina-docking` skill. Parameter and API behavior must be
checked against the official AutoDock Vina documentation and the installed `vina --help`.
