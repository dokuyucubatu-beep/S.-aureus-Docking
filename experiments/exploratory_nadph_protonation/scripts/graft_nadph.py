"""Graft NADPH from 3FRE onto the repaired 2W9H receptor (CA superposition via RDKit).
numpy matmul/linalg are avoided: they crashed silently in the authors Windows env."""
import numpy as np
from pathlib import Path
HERE = Path(__file__).resolve().parents[1]
REPO = HERE.parents[1]
INP = HERE / "inputs"
def ca(path, chain=None):
    d = {}
    for l in Path(path).read_text().splitlines():
        if l.startswith("ATOM") and l[12:16].strip() == "CA" and l[16] in " A":
            d[int(l[22:26])] = (l[17:20], np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])]))
    return d
a = ca(REPO / "data/prepared/receptor_repaired.pdb"); b = ca(INP / "3FRE.pdb")
common = [k for k in a if k in b and a[k][0] == b[k][0]]
print("matched CA:", len(common), "of", len(a), len(b))
P = np.array([b[k][1] for k in common]); Q = np.array([a[k][1] for k in common])
from rdkit import Chem
from rdkit.Geometry import Point3D
from rdkit.Chem import rdMolAlign
def _pts(X):
    m = Chem.RWMol()
    for _ in X: m.AddAtom(Chem.Atom(6))
    c = Chem.Conformer(len(X))
    for i, x in enumerate(X): c.SetAtomPosition(i, Point3D(*map(float, x)))
    m.AddConformer(c); return m
def rot(X,R,t): return np.einsum('ij,nj->ni',R,X)+t
def kabsch(P, Q):
    rms, T = rdMolAlign.GetAlignmentTransform(_pts(P), _pts(Q))
    T = np.array(T); return T[:3, :3], T[:3, 3]
R, t = kabsch(P, Q)
for it in range(3):  # iterative outlier trimming
    dev = np.linalg.norm((rot(P,R,t)) - Q, axis=1)
    keep = dev < max(2.0, np.percentile(dev, 90)); R, t = kabsch(P[keep], Q[keep])
dev = np.linalg.norm((rot(P,R,t)) - Q, axis=1)
print(f"CA RMSD all {np.sqrt((dev**2).mean()):.3f}  core(<2A, n={int((dev<2).sum())}) {np.sqrt((dev[dev<2]**2).mean()):.3f}")
out = []
for l in (INP / "3FRE.pdb").read_text().splitlines():
    if l.startswith("HETATM") and l[17:20] == "NDP":
        x = rot(np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])]]),R,t)[0]
        out.append(l[:30] + "%8.3f%8.3f%8.3f" % tuple(x) + l[54:])
(INP / "NDP_grafted.pdb").write_text("\n".join(out) + "\nEND\n")
# TMP of 3FRE in 2W9H frame, for a cross-check of the superposition
tp = []
for l in (INP / "3FRE.pdb").read_text().splitlines():
    if l.startswith("HETATM") and l[17:20] == "TOP":
        tp.append(rot(np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])]]),R,t)[0])
ref = [np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])]) for l in (REPO/"data/prepared/native_ligand.pdb").read_text().splitlines() if l.startswith(("ATOM","HETATM")) and l[12:16].strip()[0] != "H"]
print("TMP centroid in 3FRE-grafted frame vs 2W9H native:", np.round(np.mean(tp,0),2), np.round(np.mean(ref,0),2))
