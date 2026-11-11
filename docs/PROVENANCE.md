# Provenance

Downloaded on 2026-11-11:

- `data/raw/2W9H.pdb`: https://files.rcsb.org/download/2W9H.pdb
- `data/raw/TOP.cif`: https://files.rcsb.org/ligands/download/TOP.cif
- Local executable (not distributed): official AutoDock Vina 1.2.7 Windows release,
  https://github.com/forlilab/AutoDock-Vina/releases/tag/v1.2.7

The preparation audit records raw structure and protocol hashes. The run manifest
records prepared inputs, protocol, source scripts, executable, logs and pose hashes.
Original bytes are retained; every file under `data/` and `experiments/` is
exempt from Git newline conversion, including JSON and logs.

The project was created from the user's example description, not copied from
an existing AChE implementation. No AChE outputs or numerical results were reused.

## Preparation history

An initial exploratory receptor preparation omitted the incomplete terminal Lys158.
The user requested its completion before the final experiment. The released receptor
retains the residue. An unrefined PDBFixer completion produced long new bonds; a pilot
docking was stopped and excluded on that geometry criterion before evaluating its
pose metrics. The final input uses fixed-observed-atom Amber14 refinement. No pilot
results are represented as final evidence, and no historical baseline is claimed.

The user subsequently requested exhaustiveness 32. The run started at 16 was
stopped, and the final ten-seed experiment restarted at 32 with unchanged repaired
receptor, independent ligand, box and seeds. The preparation-time protocol is
archived separately in `data/prepared/preparation_protocol.json`; its docking
exhaustiveness field is historical and does not control the final experiment.
The active and run-local `protocol.json` files control docking.
The final 32-exhaustiveness run started after local commit `31f77d4`.

The local Git environment had no author identity configured. Commits created by
the assistant use `Codex <codex@localhost>` without changing global Git identity.

## Interpretation

TOP A1160 contains alternate locations A/B, including a methoxy difference.
Only A and shared atoms define the primary crystal reference. Ligand bond orders
come from the CCD by atom name; coordinates come from this specific PDB instance.
Neutral TMP is an explicit computational model, not a claim about its dominant
bound or solution microstate. Waters, EDO and Met0 are absent from prepared inputs
by the documented policy. Modeled Lys158 is not an experimentally validated rotamer.
