"""Summarize saved variant poses (symmetry-aware RMSD vs. the 2W9H reference). No Vina needed."""
import csv, re
from pathlib import Path
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolAlign
from meeko import PDBQTMolecule, RDKitMolCreate

HERE = Path(__file__).resolve().parents[1]
REPO = HERE.parents[1]
ref = Chem.SDMolSupplier(str(REPO / "data/prepared/reference.sdf"), removeHs=True)[0]

def neutral(m):
    m = Chem.RemoveHs(m)
    for a in m.GetAtoms():
        a.SetFormalCharge(0); a.SetNoImplicit(False); a.SetNumExplicitHs(0)
    Chem.SanitizeMol(m); return m

R = neutral(ref); rows = []
for variant in ("neutral", "protonated", "neutral_NADPH", "protonated_NADPH"):
    for d in sorted((HERE / "results" / variant).glob("seed_*")):
        text = (d / "poses.pdbqt").read_text()
        scores = [float(s) for s in re.findall(r"^REMARK VINA RESULT:\s+(\S+)", text, re.M)]
        m = neutral(RDKitMolCreate.from_pdbqt_mol(PDBQTMolecule(text, skip_typing=True))[0])
        assert Chem.MolToSmiles(m) == Chem.MolToSmiles(R)
        rm = [float(rdMolAlign.CalcRMS(m, R, prbId=c.GetId())) for c in m.GetConformers()]
        b = int(np.argmin(rm))
        rows.append(dict(variant=variant, seed=int(d.name[5:]), top1_score=scores[0], top1_rmsd_A=round(rm[0], 3),
                         best_rmsd_A=round(min(rm), 3), best_rank=b + 1, best_score=scores[b]))
with (HERE / "results/per_seed.csv").open("w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
out = []
for v in ("neutral", "protonated", "neutral_NADPH", "protonated_NADPH"):
    r = [x for x in rows if x["variant"] == v]; t = np.array([x["top1_rmsd_A"] for x in r])
    out.append(dict(variant=v, seeds=len(r), top1_below_2A=int((t < 2).sum()), top1_rmsd_mean=round(float(t.mean()), 3),
                    top1_rmsd_sd=round(float(t.std(ddof=1)), 3) if len(t) > 1 else "",
                    best_below_2A=sum(x["best_rmsd_A"] < 2 for x in r),
                    best_rank_median=float(np.median([x["best_rank"] for x in r]))))
with (HERE / "results/summary.csv").open("w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(out[0])); w.writeheader(); w.writerows(out)
for o in out: print(o)
