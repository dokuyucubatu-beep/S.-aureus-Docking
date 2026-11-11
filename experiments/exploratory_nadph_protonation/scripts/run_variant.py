"""Usage: VINA=/path/vina [RECEPTOR=inputs/receptor_NADPH.pdbqt TAG=_NADPH] python run_variant.py neutral|protonated [seeds]
Variant docking: neutral (baseline repro) or N1-protonated TMP, same receptor/box/protocol as the repo."""
import sys, json, subprocess, re
from pathlib import Path
import numpy as np
from concurrent.futures import ThreadPoolExecutor
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign
from meeko import MoleculePreparation, PDBQTWriterLegacy, PDBQTMolecule, RDKitMolCreate

import os
HERE = Path(__file__).resolve().parents[1]
REPO = HERE.parents[1]
VINA = os.environ.get("VINA", "vina")  # Vina 1.2.7 executable
variant = sys.argv[1]  # neutral | protonated
seeds = [int(s) for s in sys.argv[2].split(",")] if len(sys.argv) > 2 else list(range(2026, 2036))
RECEPTOR = os.environ.get("RECEPTOR", str(REPO / "data/prepared/receptor.pdbqt"))
out = HERE / "results" / (variant + os.environ.get("TAG", ""))
out.mkdir(parents=True, exist_ok=True)
smiles = {"neutral": "COc1cc(Cc2cnc(N)nc2N)cc(OC)c1OC",
          "protonated": "COc1cc(Cc2c[nH+]c(N)nc2N)cc(OC)c1OC"}[variant]

mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
ps = AllChem.ETKDGv3(); ps.randomSeed = 2026
assert AllChem.EmbedMolecule(mol, ps) == 0
AllChem.MMFFOptimizeMolecule(mol, mmffVariant="MMFF94s", maxIters=2000)
setup = MoleculePreparation().prepare(mol)[0]
pdbqt, ok, err = PDBQTWriterLegacy.write_string(setup)
assert ok, err
(out / "ligand.pdbqt").write_text(pdbqt)

box = json.loads((REPO / "data/prepared/box.json").read_text())
cfg = json.loads((REPO / "configs/protocol.json").read_text())
ref = Chem.SDMolSupplier(str(REPO / "data/prepared/reference.sdf"), removeHs=True)[0]

def dock(seed):
    d = out / f"seed_{seed}"; d.mkdir(exist_ok=True)
    cmd = [VINA, "--receptor", RECEPTOR, "--ligand", str(out / "ligand.pdbqt"),
           "--out", str(d / "poses.pdbqt"), "--seed", str(seed), "--cpu", "2", "--exhaustiveness", "32",
           "--num_modes", "20", "--energy_range", "5"]
    for i, a in enumerate("xyz"):
        cmd += [f"--center_{a}", str(box["center"][i]), f"--size_{a}", str(box["size"][i])]
    with (d / "vina.log").open("w") as log:
        subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=True)
    return seed

def neutralize(m):
    m = Chem.RemoveHs(m)
    for a in m.GetAtoms():
        a.SetFormalCharge(0); a.SetNoImplicit(False); a.SetNumExplicitHs(0)
    Chem.SanitizeMol(m)
    return m

def analyze(seed):
    text = (out / f"seed_{seed}/poses.pdbqt").read_text()
    scores = [float(s) for s in re.findall(r"^REMARK VINA RESULT:\s+(\S+)", text, re.M)]
    mols = RDKitMolCreate.from_pdbqt_mol(PDBQTMolecule(text, skip_typing=True))
    m = neutralize(mols[0])
    r = neutralize(ref)
    assert Chem.MolToSmiles(m) == Chem.MolToSmiles(r)
    rm = [float(rdMolAlign.CalcRMS(m, r, prbId=c.GetId())) for c in m.GetConformers()]
    return scores, rm

with ThreadPoolExecutor(max_workers=2) as pool:
    list(pool.map(dock, seeds))
rows = []
for s in seeds:
    sc, rm = analyze(s)
    b = int(np.argmin(rm))
    rows.append((s, sc[0], rm[0], min(rm), b + 1))
    print(f"seed {s}: top1 score {sc[0]:.2f} RMSD {rm[0]:.3f} | best RMSD {min(rm):.3f} at rank {b+1} (score {sc[b]:.2f})", flush=True)
t = np.array([r[2] for r in rows])
print(f"{variant}: top1<2A {int((t<2).sum())}/{len(rows)}, top1 RMSD mean {t.mean():.3f} +/- {t.std(ddof=1) if len(t)>1 else 0:.3f}, best<2A {sum(r[3]<2 for r in rows)}/{len(rows)}")
