# Validation and PyMOL verification

## Four evidence levels

1. **File/provenance integrity**: hashes, versions, exact command, seed, log, input identity.
2. **Geometric plausibility**: ligand inside the intended box/pocket; no severe overlaps; sensible
   ligand geometry and local contacts. PyMOL helps here.
3. **Protocol performance**: cognate-ligand redocking, multi-seed convergence, and
   target-appropriate controls.
4. **Biological validity**: orthogonal computation and, ultimately, experimental evidence.

Never promote evidence from one level into a stronger claim.

## PyMOL inspection checklist

- Load the exact prepared receptor and exact docked pose used in the report.
- Confirm object/state counts, ligand residue identity, chain/residue numbering, and coordinate
  frame.
- Display the search box and confirm the ligand and intended pocket fit inside it.
- Inspect at least the top 3–5 distinct pose families.
- Look for heavy-atom overlaps, ligand strain, buried unsatisfied polar groups, and unreasonable
  exposure of hydrophobic groups.
- Inspect expected hinge/catalytic residues, metals, cofactors, and retained waters.
- Measure distances, but label them as candidate contacts until chemistry and geometry are
  verified.
- Save both the `.pse` session and a ray-traced PNG.

## Redocking

When a compatible cognate ligand exists:

1. Remove it from the receptor without changing the coordinate frame.
2. Prepare it independently with the same protocol as test ligands.
3. Dock it using the planned box and multiple seeds.
4. Align on the receptor if necessary.
5. Calculate a symmetry-aware heavy-atom RMSD with explicit atom mapping.

An RMSD near or below 2 Å is a common pose-recovery benchmark, not a universal pass certificate.
Report the exact RMSD method, receptor conformation, ligand symmetry handling, and whether the
crystal ligand was minimized or altered.

## Multi-seed convergence

Use at least three recorded seeds for a focused known-site run. Cluster poses by symmetry-aware
ligand RMSD and summarize:

- fraction of runs recovering each pose family;
- best and median scores within each family;
- whether the top-ranked family is stable as exhaustiveness increases.

Do not select only the seed that gives the expected picture.

## Interpretation language

Acceptable:

- "Vina ranked this pose first under the stated protocol."
- "The pose is geometrically compatible with the selected EGFR pocket."
- "PyMOL shows candidate polar contacts at the reported heavy-atom distances."

Avoid:

- "The docking proves binding."
- "The Vina score is the experimental binding free energy."
- "A short distance alone proves a hydrogen bond."
- "A visually attractive pose validates the model."

## Suggested report fields

- scientific question and docking mode;
- receptor/ligand provenance and preparation;
- search-box rationale and dimensions;
- all search/output parameters and software versions;
- score table and pose-cluster summary across seeds;
- redocking RMSD or explicit reason it was unavailable;
- PyMOL session/image and geometric warnings;
- limitations and proposed orthogonal validation.
