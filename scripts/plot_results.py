"""Additional plots from validated saved poses; does not rerun docking."""
import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from rdkit import Chem

from analyze import analyze, read_poses
from common import ROOT, coordinates


def plot(run_dir):
    run_dir = Path(run_dir)
    analyze(run_dir)
    out = run_dir / "analysis"
    with (out / "contacts.csv").open() as f:
        rows = list(csv.DictReader(f))
    seeds = sorted({int(r["seed"]) for r in rows})
    residues = sorted({r["residue"] for r in rows}, key=lambda s: int(s.split(":")[1]))
    values = np.zeros((len(residues), len(seeds)+1))
    for r in rows:
        i = residues.index(r["residue"])
        values[i, 0] = r["native_contact"] == "True"
        values[i, seeds.index(int(r["seed"]))+1] = r["top1_contact"] == "True"
    fig, ax = plt.subplots(figsize=(9, max(4, len(residues)*0.28)), layout="constrained")
    ax.imshow(values, aspect="auto", cmap="Blues", vmin=0, vmax=1, interpolation="nearest")
    ax.set_xticks(range(len(seeds)+1), ["Crystal"] + [str(s) for s in seeds], rotation=45)
    ax.set_yticks(range(len(residues)), residues)
    ax.set_title("Residue proximity contacts (heavy-atom distance <=4 A)")
    fig.savefig(out / "contacts.png", dpi=180)
    plt.close(fig)
    ref = Chem.SDMolSupplier(str(run_dir/"inputs/reference.sdf"))[0]
    poses, _ = read_poses(run_dir / f"seed_{seeds[0]}/poses.pdbqt")
    fig = plt.figure(figsize=(7, 6), layout="constrained")
    ax = fig.add_subplot(projection="3d")
    all_xyz = []
    for mol, color, label in ((ref, "#157f81", "Crystal altloc A"), (poses, "#e77935", "Top-ranked")):
        xyz = coordinates(mol, 0)
        all_xyz.append(xyz)
        ax.scatter(*xyz.T, color=color, s=22, label=label)
        for bond in mol.GetBonds():
            points = xyz[[bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()]]
            ax.plot(*points.T, color=color, linewidth=2)
    xyz = np.concatenate(all_xyz)
    center = (xyz.min(0)+xyz.max(0))/2
    radius = max(np.ptp(xyz, axis=0))/2+1
    ax.set(xlim=(center[0]-radius, center[0]+radius), ylim=(center[1]-radius, center[1]+radius),
           zlim=(center[2]-radius, center[2]+radius), xlabel="X (A)", ylabel="Y (A)", zlabel="Z (A)",
           title=f"2W9H / trimethoprim: seed {seeds[0]} (no alignment)")
    ax.set_box_aspect((1, 1, 1))
    ax.legend()
    fig.savefig(out / "overlay_seed2026.png", dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--run-dir", type=Path, default=ROOT / "experiments/independent")
    plot(p.parse_args().run_dir)
