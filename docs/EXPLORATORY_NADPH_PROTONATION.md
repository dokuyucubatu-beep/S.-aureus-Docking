# Exploratory add-on: NADPH cofactor and N1-protonated trimetoprim

**Status: exploratory, post hoc.** The primary, pre-specified result of this repository
(`experiments/independent/`, top-ranked RMSD < 2 A in 0/10 runs) is unchanged and is
not replaced by anything below. These runs were added after seeing that result.

## Motivation

PDB 2W9H contains trimetoprim (TOP) and ethylene glycol but no NADPH, although DHFR
binds trimetoprim together with the cofactor. The primary protocol also docks the
neutral ligand, whereas the N1 of the diaminopyrimidine is commonly protonated in
DHFR complexes. Two single-factor changes and their combination were tested.

## Design (everything else identical to the primary protocol)

Same repaired receptor (`data/prepared/receptor.pdbqt`), box, Vina 1.2.7,
exhaustiveness 32, 20 modes, energy range 5, seeds 2026-2035, same RDKit
ETKDGv3/MMFF94s embedding seed, symmetry-aware `rdMolAlign.CalcRMS` against
`reference.sdf` (formal charges removed from both sides before the identity check).

| Variant | Receptor | Ligand SMILES |
|---|---|---|
| neutral (primary) | repaired 2W9H | `COc1cc(Cc2cnc(N)nc2N)cc(OC)c1OC` |
| protonated | repaired 2W9H | `COc1cc(Cc2c[nH+]c(N)nc2N)cc(OC)c1OC` |
| neutral_NADPH | repaired 2W9H + grafted NADPH | neutral |
| protonated_NADPH | repaired 2W9H + grafted NADPH | N1-protonated |

NADPH source: PDB 3FRE (wild-type *S. aureus* DHFR with NADPH and trimetoprim, 2.2 A,
[doi:10.2210/pdb3FRE/pdb](https://doi.org/10.2210/pdb3FRE/pdb)). 157 matched C-alpha atoms
were superposed onto the repaired 2W9H receptor (RMSD 0.463 A, `scripts/graft_nadph.py`).
The 3FRE trimetoprim centroid lies about 0.7 A from the 2W9H reference after this
superposition. NADPH was converted with Open Babel (`-xr -p 7.4 -h`, rigid, Gasteiger
charges) and its atom lines were appended to `receptor.pdbqt`
(`inputs/receptor_NADPH.pdbqt`). Vina does not use partial charges, so charge sums are
not scored.

## Results

| Variant | Seeds | Top-1 < 2 A | Top-1 RMSD mean +/- SD (A) | Best sampled < 2 A | Median rank of best pose |
|---|---|---|---|---|---|
| neutral (primary, one seed re-run only) | 1 | 0/1 | 3.609 | 1/1 | 14 |
| protonated | 10 | 0/10 | 3.558 +/- 0.030 | 10/10 | 5 |
| neutral_NADPH | 10 | 0/10 | 6.620 +/- 0.881 | 10/10 | 5 |
| protonated_NADPH | 10 | **10/10** | **1.192 +/- 0.146** | 10/10 | 2 |

Per-seed values: `experiments/exploratory_nadph_protonation/results/per_seed.csv`;
raw poses and Vina logs are archived next to it. The primary ten-seed numbers
(3.599 +/- 0.023 A) are in `experiments/independent/analysis/`.

Neither change alone placed the correct pose first. NADPH with the neutral ligand
moved the top-ranked pose far from the reference (about 6.9 A in 9/10 runs); this was
not investigated. Only the combination recovered the crystallographic pose at rank 1.

## Limits

- Post hoc and unregistered: the choices were made after the primary result.
- The cofactor pose comes from a trimetoprim-bound complex (3FRE), so it carries
  holo information about the pocket; a cofactor taken from a ligand-free or different
  complex (for example 3FRD) would be a stricter control and was not run.
- NADPH preparation was a quick Open Babel conversion without manual inspection of
  atom typing or protonation of the phosphates.
- Cognate self-docking only; Vina has no electrostatics, so the N1 protonation acts
  through hydrogen-bond donor typing, not charge. This is not evidence of potency or
  of docking accuracy for other ligands.
- The neutral, no-NADPH variant was re-run for seed 2026 only (3.609 A, matching the
  archived primary run); the other nine seeds are the archived primary results.
- Archived dependencies: Vina 1.2.7, RDKit 2025.9.6, Meeko 0.7.1.

## Reproduce

```bash
python experiments/exploratory_nadph_protonation/scripts/graft_nadph.py
VINA=/path/to/vina RECEPTOR=experiments/exploratory_nadph_protonation/inputs/receptor_NADPH.pdbqt TAG=_NADPH \
  python experiments/exploratory_nadph_protonation/scripts/run_variant.py protonated
python experiments/exploratory_nadph_protonation/scripts/summarize.py   # from saved poses
```

## Turkce ozet

Bu bolum **kesifsel ve sonradan eklenmistir**; asil, onceden belirlenmis sonuc (top-1 < 2 A:
0/10) degismez. 2W9H'de NADPH yoktur; NADPH, 3FRE'den (Calpha RMSD 0.46 A) aktarildi ve
N1-protonlanmis trimetoprim denendi. Yalniz biri tek basina yetmedi (0/10); ikisi birlikte
top-1 poz 10/10 tohumda 1.19 +/- 0.15 A ile cozuldu. Kofaktor TMP'li bir kompleksten geldigi
icin sonuc hafifce iyimser olabilir.
