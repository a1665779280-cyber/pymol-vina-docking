# Input files and preparation decisions

## Minimal docking inputs

| Role | Required format | Requirements |
|---|---|---|
| Rigid receptor | PDBQT | Correct chains, atom types, hydrogens/protonation, charges, and retained site components |
| Ligand | PDBQT | One intended chemical species with 3D coordinates, charges, atom types, and torsion tree |
| Search space | Config/CLI | `center_x/y/z` and `size_x/y/z` in Å with a recorded rationale |
| Output | PDBQT + log | Ranked poses plus complete execution parameters and random seed |

PDB or mmCIF is not a docking-ready substitute for receptor PDBQT. SDF, MOL2, or SMILES is not
a docking-ready substitute for ligand PDBQT. Preserve these richer source formats because PDBQT
does not retain every bond-order and stereochemical detail needed for chemical auditing.

## Recommended provenance inputs

- Receptor source PDB/mmCIF, PDB ID, biological assembly, chain selection, resolution/method,
  mutation state, and preparation tool/version.
- Prepared receptor PDB used for PyMOL inspection.
- Ligand source SDF/MOL2 and canonical/isomeric SMILES; compound identifier and salt handling.
- Cognate/reference ligand coordinates or an experimentally supported pocket-residue list.
- Preparation pH and chosen protonation/tautomer states.
- Records for missing atoms/residues, alternate conformers, retained waters, metals, and
  cofactors.
- File hashes for every final input.

## Receptor decision checklist

1. Choose the biologically relevant chain/assembly.
2. Inspect missing residues and side-chain atoms, especially within 8–10 Å of the site.
3. Resolve alternate locations explicitly; do not keep mutually exclusive conformers together.
4. Assign protonation for the stated pH. Histidine identity and buried acidic/basic residues
   deserve manual review.
5. Retain catalytic metals/cofactors and bridging waters when justified; do not delete every
   heteroatom mechanically.
6. Remove crystallization additives and unrelated ligands only after identifying them.
7. Confirm the prepared PDB and receptor PDBQT represent the same coordinates and chain set.

## Ligand decision checklist

1. Verify identity, bond orders, formal charge, stereochemistry, and salt/mixture handling.
2. Enumerate protonation or tautomer states only when scientifically plausible for the target
   pH; dock and report them as separate species.
3. Generate a reasonable 3D conformer and minimize it without changing intended stereochemistry.
4. Inspect aromaticity, amide/nonrotatable bonds, macrocycles, and metal-binding groups.
5. Confirm the PDBQT torsion tree and rotatable-bond count.
6. Preserve the source SDF/MOL2 and preparation command.

## Reject or pause when

- the receptor/ligand identity is ambiguous;
- coordinates are all zero, non-finite, or clearly 2D;
- the intended binding site is unsupported and the user did not request exploratory blind
  docking;
- a catalytic metal/cofactor was removed without rationale;
- the prepared and source structures cannot be reconciled;
- a covalent ligand is being treated as ordinary noncovalent Vina docking without an explicit
  compatible protocol.
