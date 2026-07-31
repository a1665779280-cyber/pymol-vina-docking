# Parameter-selection process

## 1. Choose the docking mode

| Mode | Site evidence | Preferred validation |
|---|---|---|
| Redocking | Cognate ligand in the same structure | Symmetry-aware heavy-atom RMSD and pose recovery |
| Known-site docking | Cognate ligand, mutagenesis, or well-supported residues | Multi-seed convergence and known-interaction consistency |
| Virtual screening | Validated fixed protocol | Enrichment/decoy performance plus repeatability |
| Exploratory blind docking | No defensible site | Treat as hypothesis generation; confirm pockets independently |

Do not derive the box from the pose that the same docking run produced.

## 2. Select the center and size

Evidence priority:

1. Cognate-ligand coordinates from a compatible receptor conformation.
2. Experimentally supported pocket residues.
3. A pocket predictor validated or cross-checked for the target.
4. Whole-protein/blind search only when explicitly justified.

For a reference ligand or residue selection:

```text
center_axis = (minimum_axis + maximum_axis) / 2
size_axis   = maximum_axis - minimum_axis + 2 × padding_per_side
```

Use 4–6 Å padding per side as a starting point for a typical small molecule. Increase only when
needed for ligand length, uncertain side-chain conformations, or movement around the known site.
Inspect all six box faces in PyMOL. A large box increases search difficulty; compensate through
validation and search effort rather than assuming that "bigger is safer."

With an explicit center and size, pass the independent site evidence to the helper as
`--box-rationale`. The helper rejects an explicit box when that rationale is omitted.

## 3. Choose the scoring function

- `vina`: default starting point for conventional noncovalent small-molecule docking.
- `vinardo`: alternative; use only after target-specific redocking or benchmark evidence.
- `ad4`: requires the appropriate AutoDock maps/workflow. Its score is not numerically
  comparable with Vina/Vinardo.

Do not average raw scores from different scoring functions as a physical consensus energy.
Consensus ranking requires a declared normalization/rank method and independent validation.

## 4. Choose search and output parameters

| Parameter | Official behavior | Selection rule |
|---|---|---|
| `exhaustiveness` | Search effort; CLI/API default 8 | Determine by convergence. Pilot 8/16/32 across ≥3 seeds; increase until pose clusters stabilize enough for the question |
| `seed` | `0` selects a random seed; explicit value is reproducible | Record every seed; use multiple seeds, not one favorite seed |
| `num_modes` / API `n_poses` | Number of poses retained/generated | Usually retain 10–20 for pose-family analysis; do not report only pose 1 |
| `energy_range` | Output filter from best score; default 3 kcal/mol | Increase when scientifically relevant alternative pose families fall outside 3; it does not improve the search |
| `min_rmsd` | API pose-separation threshold; default 1 Å | Keep default initially; change only with a reason tied to ligand size/symmetry |
| `cpu` | Parallel resources | Performance parameter; record it but do not interpret it scientifically |
| `max_evals` | Overrides heuristic search evaluations | Leave at default unless benchmarking and justification are available |

There is no universal "publication-quality exhaustiveness." `32` is an official tutorial example
for a challenging ligand, not a general proof of convergence.

## 5. Receptor flexibility

Rigid-receptor docking is the default approximation. Make side chains flexible only when:

- alternate rotamers or induced fit are supported by structures/experiments;
- the selected residues are limited and documented;
- increased dimensionality is matched by stronger convergence testing.

For major backbone rearrangements, ensemble docking against multiple prepared receptor
conformations is usually more interpretable than declaring many flexible side chains.

## 6. EGFR–erlotinib example from the supplied run

The supplied `vina.log` records:

```text
scoring        = vina
center         = 16.2921, 34.8708, 92.0353 Å
size           = 18.000, 19.032, 22.828 Å
exhaustiveness = 32
num_modes      = 10
seed           = 2026
```

These values are provenance for that run, not automatically correct defaults for another EGFR
structure. Reuse them only if the receptor coordinates are in the identical frame and the same
binding-site rationale still applies.

## Official behavior references

- AutoDock Vina basic docking tutorial:
  https://autodock-vina.readthedocs.io/en/stable/docking_basic.html
- AutoDock Vina Python API:
  https://autodock-vina.readthedocs.io/en/stable/vina.html
