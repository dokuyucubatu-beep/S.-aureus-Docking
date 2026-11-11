import json
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np
from rdkit import Chem

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from common import contacts, coordinates, pdb_atoms, rmsd
from analyze import read_poses


class ChemistryTests(unittest.TestCase):
    def setUp(self):
        self.ref = Chem.SDMolSupplier(str(ROOT / "data/prepared/reference.sdf"))[0]

    def test_translation_is_not_aligned_away(self):
        shifted = Chem.Mol(self.ref)
        conf = shifted.GetConformer()
        for i, p in enumerate(coordinates(shifted)):
            conf.SetAtomPosition(i, p + [10, 0, 0])
        self.assertAlmostEqual(rmsd(shifted, self.ref), 10.0, places=5)

    def test_atom_order_and_hydrogens(self):
        reordered = Chem.RenumberAtoms(self.ref, list(reversed(range(self.ref.GetNumAtoms()))))
        self.assertAlmostEqual(rmsd(reordered, self.ref), 0.0)
        self.assertAlmostEqual(rmsd(Chem.AddHs(self.ref, addCoords=True), self.ref), 0.0)

    def test_charge_mismatch_rejected(self):
        changed = Chem.Mol(self.ref)
        for atom in changed.GetAtoms():
            if atom.GetAtomicNum() == 7:
                atom.SetFormalCharge(1)
                break
        with self.assertRaises(ValueError):
            rmsd(changed, self.ref)

    def test_reference_matches_observed_crystal(self):
        native = np.array(list(pdb_atoms(ROOT / "data/prepared/native_ligand.pdb").values()))
        xyz = coordinates(self.ref)
        self.assertEqual(len(native), 21)
        distances = np.linalg.norm(xyz[:, None, :] - native[None, :, :], axis=2)
        self.assertLess(distances.min(axis=1).max(), 0.0002)

    def test_independent_embedding_not_native(self):
        mol = Chem.SDMolSupplier(str(ROOT / "data/prepared/ligand_independent.sdf"))[0]
        self.assertGreater(rmsd(mol, self.ref), 5)

    def test_lys158_completed_observed_atoms_preserved(self):
        source = pdb_atoms(ROOT / "data/prepared/receptor_observed.pdb")
        repaired = pdb_atoms(ROOT / "data/prepared/receptor_repaired.pdb")
        names = {k[3] for k in repaired if k[1] == "158"}
        self.assertEqual(names, {"N", "CA", "C", "O", "OXT", "CB", "CG", "CD", "CE", "NZ"})
        for k, xyz in source.items():
            np.testing.assert_allclose(xyz, repaired[k], atol=0.002)
        prepared = pdb_atoms(ROOT / "data/prepared/receptorH.pdb")
        self.assertTrue(all(k in prepared for k in repaired))
        for k, xyz in repaired.items():
            if k[1] == "158" and k[3] in ("O", "OXT"):
                continue  # Meeko exchanges equivalent carboxylate oxygen names.
            np.testing.assert_allclose(xyz, prepared[k], atol=0.002)
        before_oxygen = sorted(tuple(repaired[("A", "158", "LYS", n)]) for n in ("O", "OXT"))
        after_oxygen = sorted(tuple(prepared[("A", "158", "LYS", n)]) for n in ("O", "OXT"))
        np.testing.assert_allclose(before_oxygen, after_oxygen, atol=0.002)

    def test_contact_boundary(self):
        receptor = {("A", "1", "ALA", "CA"): np.array([0., 0., 0.])}
        self.assertEqual(contacts(receptor, np.array([[4., 0., 0.]]))[0], {"A:1:ALA"})
        self.assertEqual(contacts(receptor, np.array([[4.01, 0., 0.]]))[0], set())

    def test_repaired_lysine_geometry_and_handedness(self):
        atoms = pdb_atoms(ROOT / "data/prepared/receptor_repaired.pdb")
        def xyz(name):
            return atoms[("A", "158", "LYS", name)]
        for a, b in (("N", "CA"), ("CA", "C"), ("C", "O"), ("C", "OXT"),
                     ("CA", "CB"), ("CB", "CG"), ("CG", "CD"), ("CD", "CE"), ("CE", "NZ")):
            self.assertGreater(np.linalg.norm(xyz(a)-xyz(b)), 1.15)
            self.assertLess(np.linalg.norm(xyz(a)-xyz(b)), 1.65)
        # N/C/CB signed volume has the same sign as the observed L-lysines.
        signed = np.linalg.det(np.array([xyz(n)-xyz("CA") for n in ("N", "C", "CB")]))
        for residue in {k[1] for k in atoms if k[2] == "LYS" and k[1] != "158"}:
            reference = np.linalg.det(np.array([atoms[("A", residue, "LYS", n)] -
                                                atoms[("A", residue, "LYS", "CA")]
                                                for n in ("N", "C", "CB")]))
            self.assertGreater(signed * reference, 0)

    def test_malformed_pose_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            p = Path(directory) / "bad.pdbqt"
            for text in ("", "MODEL 1\nREMARK VINA RESULT: -5 0 0\n", "MODEL 1\nREMARK VINA RESULT: nan 0 0\nENDMDL\n"):
                p.write_text(text)
                with self.assertRaises(ValueError):
                    read_poses(p)


if __name__ == "__main__":
    unittest.main()
