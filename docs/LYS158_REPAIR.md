# Lys158 completion and checks

PDB 2W9H REMARK 470 lists missing CA, C, O, CB, CG, CD, CE and NZ for Lys A158.
The deposited N is retained. PDBFixer 1.12.0 adds these eight atoms and terminal
OXT with seed 2026. OpenMM 8.6.1 / Amber14 initializes hydrogens at pH 7 and
minimizes added heavy atoms and hydrogens with all observed heavy atoms assigned
zero mass (fixed). Reference platform, tolerance 1 kJ/mol/nm, maximum 5000
iterations; no solvent or whole-protein relaxation. Temporary hydrogens are
removed, and Meeko 0.7.1 creates the final receptor using its residue templates.

Observed heavy-atom displacement after PDB export: **0.000 A**.
New heavy-atom bonds range from **1.245 to 1.549 A** before PDB rounding.
The existing Lys157 C–Lys158 N bond is approximately 1.233 A in the original
structure and remains unchanged; the procedure is not a global geometry repair.

The CA signed tetrahedral volume has the same sign as all ten observed lysines
(L stereochemistry). Among added heavy atoms, excluding atom pairs connected
by up to three covalent bonds, the nearest heavy-atom pair is Lys158 NZ to
Asp130 OD1 at approximately **2.668 A**. No pair in that check is below 1.8 A.
This is a severe-overlap screen, not a full crystallographic validation score.

The minimum Lys158-to-native-TOP heavy-atom distance is approximately **23.287 A**.
All checks refer to the archived `data/prepared/receptor_repaired.pdb` and
`preparation_audit.json`. Tests recompute added-bond lengths, handedness,
residue completeness and preservation of observed coordinates.

Meeko exchanges the names O and OXT of the two equivalent terminal carboxylate
oxygens. Their coordinate set is unchanged. Validation checks these two atoms
as an unordered pair and every other heavy atom by its name; no physical
heavy-atom coordinate displacement occurs during Meeko preparation.

One observed backbone atom cannot uniquely determine a terminal residue's
conformation. This is one physically relaxed completion, not an experimentally
resolved rotamer or an ensemble. Use the archived prepared receptor when
reproducing the docking calculation.
