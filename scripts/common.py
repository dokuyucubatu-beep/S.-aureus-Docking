"""Shared validation and provenance primitives."""
import hashlib
import json
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolAlign

ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def protocol():
    return json.loads((ROOT / "configs/protocol.json").read_text())


def coordinates(mol, conf_id=-1):
    xyz = np.asarray(mol.GetConformer(conf_id).GetPositions())
    if not np.isfinite(xyz).all():
        raise ValueError("Nonfinite molecular coordinates")
    return xyz


def rmsd(probe, reference, conf_id=-1):
    a, b = Chem.RemoveHs(probe), Chem.RemoveHs(reference)
    if Chem.MolToSmiles(a) != Chem.MolToSmiles(b):
        raise ValueError("Chemical identity or formal charge mismatch")
    coordinates(a, conf_id)
    coordinates(b)
    return float(rdMolAlign.CalcRMS(a, b, prbId=conf_id))


def pdb_atoms(path):
    atoms = {}
    for line in Path(path).read_text().splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            key = (line[21], line[22:27].strip(), line[17:20].strip(), line[12:16].strip())
            if key in atoms:
                raise ValueError(f"Duplicate atom {key}")
            atoms[key] = np.array([float(line[i:i+8]) for i in (30, 38, 46)])
    return atoms


def contacts(receptor_atoms, ligand_xyz, cutoff=4.0):
    distances = {}
    for key, xyz in receptor_atoms.items():
        if key[3].startswith("H"):
            continue
        residue = ":".join(key[:3])
        distance = float(np.linalg.norm(ligand_xyz - xyz, axis=1).min())
        distances[residue] = min(distance, distances.get(residue, float("inf")))
    return {k for k, v in distances.items() if v <= cutoff}, distances
