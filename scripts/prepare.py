"""Prepare one fixed, audited 2W9H receptor and an independent SMILES conformer."""
import argparse
import importlib.metadata
import json
import subprocess
import sys
from pathlib import Path

import gemmi
import numpy as np
from meeko import MoleculePreparation, PDBQTWriterLegacy
from openmm import Platform, unit
from openmm.app import PDBFile
from pdbfixer import PDBFixer
from rdkit import Chem
from rdkit.Chem import AllChem

from common import ROOT, coordinates, digest, pdb_atoms, protocol, write_json


def write_sdf(path, mol):
    with Chem.SDWriter(str(path)) as writer:
        writer.write(mol)


def crystal_reference(lines, cfg, out):
    selected = [l[:16] + " " + l[17:] for l in lines
                if l.startswith("HETATM") and l[17:20] == cfg["ligand_resname"]
                and l[21] == cfg["chain"] and l[22:26].strip() == cfg["ligand_resid"]
                and l[16] in (" ", cfg["altloc"])]
    (out / "native_ligand.pdb").write_text("\n".join(selected) + "\nEND\n")
    atoms = {l[12:16].strip(): l for l in selected}
    if len(atoms) != len(selected):
        raise ValueError("Duplicate ligand atoms")
    block = gemmi.cif.read_file(str(ROOT / "data/raw/TOP.cif")).sole_block()
    cat = block.get_mmcif_category("_chem_comp_atom.")
    rw, index = Chem.RWMol(), {}
    for name, element, charge in zip(cat["atom_id"], cat["type_symbol"], cat["charge"]):
        if element == "H":
            continue
        atom = Chem.Atom(element)
        atom.SetFormalCharge(int(charge))
        atom.SetProp("ccd_atom_name", name)
        index[name] = rw.AddAtom(atom)
    if set(atoms) != set(index):
        raise ValueError("Crystal heavy atoms differ from CCD")
    bonds = block.get_mmcif_category("_chem_comp_bond.")
    orders = {"SING": Chem.BondType.SINGLE, "DOUB": Chem.BondType.DOUBLE, "TRIP": Chem.BondType.TRIPLE}
    for a, b, order in zip(bonds["atom_id_1"], bonds["atom_id_2"], bonds["value_order"]):
        if a in index and b in index:
            rw.AddBond(index[a], index[b], orders[order])
    mol = rw.GetMol()
    Chem.SanitizeMol(mol)
    conf = Chem.Conformer(mol.GetNumAtoms())
    for name, i in index.items():
        conf.SetAtomPosition(i, [float(atoms[name][j:j+8]) for j in (30, 38, 46)])
    mol.AddConformer(conf)
    expected = Chem.MolFromSmiles(cfg["smiles"])
    if Chem.MolToSmiles(mol) != Chem.MolToSmiles(expected):
        raise ValueError("CCD and configured SMILES disagree")
    write_sdf(out / "reference.sdf", mol)
    return mol


def prepare(out):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    cfg = protocol()
    lines = (ROOT / "data/raw/2W9H.pdb").read_text().splitlines()
    ref = crystal_reference(lines, cfg, out)
    selected = [l[:16] + " " + l[17:] for l in lines if l.startswith("ATOM  ")
                and l[21] == cfg["chain"] and l[16] in (" ", cfg["altloc"])]
    source = out / "receptor_observed.pdb"
    source.write_text("\n".join(selected) + "\nEND\n")
    fixer = PDBFixer(filename=str(source), platform=Platform.getPlatformByName("Reference"))
    fixer.findMissingResidues()
    fixer.missingResidues = {}  # Do not invent the absent N-terminal Met0.
    fixer.findMissingAtoms()
    missing = {f"{r.chain.id}:{r.id}:{r.name}": [a.name for a in aa] for r, aa in fixer.missingAtoms.items()}
    if set(missing) != {"A:158:LYS"}:
        raise ValueError(f"Unexpected missing residues/atoms: {missing}")
    fixer.addMissingAtoms(seed=cfg["repair_seed"])
    repaired = out / "receptor_repaired.pdb"
    with repaired.open("w") as handle:
        PDBFile.writeFile(fixer.topology, fixer.positions, handle, keepIds=True)
    before, after = pdb_atoms(source), pdb_atoms(repaired)
    if not set(before).issubset(after):
        raise ValueError("Repair removed observed atoms")
    max_shift = max(float(np.linalg.norm(after[k] - v)) for k, v in before.items())
    if max_shift > 0.002:
        raise ValueError(f"Observed atoms moved by {max_shift} A")
    added = {":".join(k): v.tolist() for k, v in after.items() if k not in before}
    xyz = fixer.positions.value_in_unit(unit.angstrom)
    bond_lengths = {}
    for a, b in fixer.topology.bonds():
        if a.residue.id == "158" or b.residue.id == "158":
            length = float(np.linalg.norm(np.array(xyz[a.index]) - np.array(xyz[b.index])))
            bond_lengths[f"{a.residue.id}:{a.name}-{b.residue.id}:{b.name}"] = length
            if not 1.1 < length < 1.8:
                raise ValueError(f"Implausible repaired bond: {length}")
    cmd = [sys.executable, "-m", "meeko.cli.mk_prepare_receptor", "--read_pdb", str(repaired),
           "-o", str(out / "receptor"), "-p", "-j", "--write_pdb", str(out / "receptorH.pdb")]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    (out / "receptor_preparation.log").write_text(proc.stdout + proc.stderr)
    proc.check_returncode()
    mol = Chem.AddHs(Chem.MolFromSmiles(cfg["smiles"]))
    params = AllChem.ETKDGv3()
    params.randomSeed = cfg["embedding_seed"]
    if AllChem.EmbedMolecule(mol, params) < 0:
        raise ValueError("Embedding failed")
    if AllChem.MMFFOptimizeMolecule(mol, mmffVariant="MMFF94s", maxIters=2000) != 0:
        raise ValueError("MMFF optimization failed to converge")
    write_sdf(out / "ligand_independent.sdf", mol)
    setups = MoleculePreparation().prepare(mol)
    if len(setups) != 1:
        raise ValueError("Expected one molecular setup")
    pdbqt, ok, error = PDBQTWriterLegacy.write_string(setups[0])
    if not ok:
        raise ValueError(error)
    (out / "ligand_independent.pdbqt").write_text(pdbqt)
    center = coordinates(ref).mean(axis=0).tolist()
    write_json(out / "box.json", {"center": center, "size": cfg["box_size"]})
    write_json(out / "preparation_audit.json", {
        "missing_atoms": missing, "added_atoms": added, "max_observed_atom_shift_A": max_shift,
        "repaired_bond_lengths_A": bond_lengths,
        "lys158_min_native_ligand_distance_A": min(float(np.linalg.norm(coordinates(ref)-v, axis=1).min()) for k,v in after.items() if k[1]=='158'),
        "policy": cfg["receptor_policy"], "ligand_formal_charge": Chem.GetFormalCharge(mol),
        "reference_heavy_atoms": ref.GetNumAtoms(),
        "versions": {x: importlib.metadata.version(x) for x in ("rdkit", "meeko", "pdbfixer", "openmm", "gemmi", "numpy")},
        "sha256": {str(p.relative_to(ROOT)): digest(p) for p in (ROOT/"data/raw/2W9H.pdb", ROOT/"data/raw/TOP.cif", ROOT/"configs/protocol.json")}
    })
    print(json.dumps({"prepared": str(out), "added_atoms": list(added), "max_shift_A": max_shift}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "data/prepared")
    prepare(parser.parse_args().output)
