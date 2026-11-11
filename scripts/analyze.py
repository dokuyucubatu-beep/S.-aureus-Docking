"""Analyze complete Vina multi-model files in the fixed receptor coordinate frame."""
import argparse
import csv
import json
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from meeko import PDBQTMolecule, RDKitMolCreate
from rdkit import Chem

from common import ROOT, contacts, coordinates, digest, pdb_atoms, protocol, rmsd, write_json


def read_poses(path):
    text = Path(path).read_text()
    starts = len(re.findall(r"^MODEL\s", text, re.M))
    ends = len(re.findall(r"^ENDMDL\s*$", text, re.M))
    scores = [float(s) for s in re.findall(r"^REMARK VINA RESULT:\s+(\S+)", text, re.M)]
    if not starts or starts != ends or starts != len(scores) or not np.isfinite(scores).all():
        raise ValueError("Incomplete poses or missing/nonfinite Vina scores")
    if "REMARK SMILES IDX" not in text:
        raise ValueError("Meeko atom mapping metadata required")
    mols = RDKitMolCreate.from_pdbqt_mol(PDBQTMolecule(text, skip_typing=True))
    if len(mols) != 1 or mols[0] is None or mols[0].GetNumConformers() != starts:
        raise ValueError("Expected one ligand with one conformer per MODEL")
    mol = Chem.RemoveHs(mols[0])
    for conf in mol.GetConformers():
        coordinates(mol, conf.GetId())
    return mol, scores


def csv_write(path, rows):
    with Path(path).open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def analyze(run_dir, output=None):
    run_dir = Path(run_dir)
    output = Path(output) if output else run_dir / "analysis"
    output.mkdir(parents=True, exist_ok=True)
    cfg = json.loads((run_dir / "protocol.json").read_text())
    manifest = json.loads((run_dir / "manifest.json").read_text())
    for rel, sha in manifest["sha256"].items():
        if digest(run_dir / rel) != sha:
            raise ValueError(f"Evidence hash mismatch: {rel}")
    reference = Chem.SDMolSupplier(str(run_dir / "inputs/reference.sdf"), removeHs=True)[0]
    if reference is None:
        raise ValueError("Invalid reference SDF")
    receptor = pdb_atoms(run_dir / "inputs/receptorH.pdb")
    native, native_dist = contacts(receptor, coordinates(reference), cfg["contact_cutoff_angstrom"])
    rows, runs, contact_rows = [], [], []
    for seed in cfg["seeds"]:
        mol, scores = read_poses(run_dir / f"seed_{seed}/poses.pdbqt")
        rr = []
        for rank, (conf, score) in enumerate(zip(mol.GetConformers(), scores), 1):
            row = {"seed": seed, "rank": rank, "score_kcal_mol": score,
                   "rmsd_A": rmsd(mol, reference, conf.GetId())}
            rows.append(row)
            rr.append(row)
        best = min(rr, key=lambda r: r["rmsd_A"])
        predicted, distances = contacts(receptor, coordinates(mol, 0), cfg["contact_cutoff_angstrom"])
        intersect = native & predicted
        runs.append({"seed": seed, "top1_rmsd_A": rr[0]["rmsd_A"], "best_rmsd_A": best["rmsd_A"],
                     "best_rank": best["rank"], "top1_score_kcal_mol": scores[0], "pose_count": len(rr),
                     "contact_recall": len(intersect)/len(native) if native else 0,
                     "contact_precision": len(intersect)/len(predicted) if predicted else 0,
                     "contact_jaccard": len(intersect)/len(native | predicted) if native | predicted else 0})
        for residue in sorted(native | predicted):
            contact_rows.append({"seed": seed, "residue": residue, "native_min_A": native_dist[residue],
                                 "top1_min_A": distances[residue], "native_contact": residue in native,
                                 "top1_contact": residue in predicted})
    top = np.array([r["top1_rmsd_A"] for r in runs])
    summary = {"runs": len(runs), "total_poses": len(rows),
               "top1_below_2A": sum(r["top1_rmsd_A"] < cfg["rmsd_cutoff_angstrom"] for r in runs),
               "any_pose_below_2A": sum(r["best_rmsd_A"] < cfg["rmsd_cutoff_angstrom"] for r in runs),
               "top1_mean_A": float(top.mean()), "top1_sample_sd_A": float(top.std(ddof=1)),
               "top1_median_A": float(np.median(top)), "top1_range_A": [float(top.min()), float(top.max())],
               "mean_top1_native_contact_recall": float(np.mean([r["contact_recall"] for r in runs]))}
    csv_write(output / "all_poses.csv", rows)
    csv_write(output / "per_run.csv", runs)
    csv_write(output / "contacts.csv", contact_rows)
    write_json(output / "summary.json", summary)
    fig, ax = plt.subplots(figsize=(8, 4.5), layout="constrained")
    ax.plot(cfg["seeds"], top, "o-", label="Top-ranked")
    ax.plot(cfg["seeds"], [r["best_rmsd_A"] for r in runs], "s--", label="Best sampled")
    ax.axhline(cfg["rmsd_cutoff_angstrom"], color="gray", linestyle=":", label="2 A cutoff")
    ax.set(xlabel="Vina seed", ylabel="Crystal-reference RMSD (A)", title="S. aureus DHFR / trimethoprim (2W9H)")
    ax.ticklabel_format(useOffset=False, style="plain")
    ax.legend()
    fig.savefig(output / "rmsd_by_seed.png", dpi=180)
    plt.close(fig)
    print(json.dumps(summary, indent=2))
    return summary


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--run-dir", type=Path, default=ROOT / "experiments/independent")
    p.add_argument("--output", type=Path)
    args = p.parse_args()
    analyze(args.run_dir, args.output)
