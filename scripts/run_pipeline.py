"""Run the fixed ten-seed experiment; existing output directories are never overwritten."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import importlib.metadata
import json
import shutil
import subprocess
from pathlib import Path

from common import ROOT, digest, protocol, write_json


def run(output, vina, jobs=1):
    cfg = protocol()
    executable = shutil.which(str(vina))
    if not executable and Path(vina).is_file():
        executable = str(Path(vina).resolve())
    if not executable:
        raise ValueError("Vina executable not found; supply --vina PATH")
    version = subprocess.check_output([executable, "--version"], text=True).strip()
    if cfg["vina_version"] not in version:
        raise ValueError(f"Expected Vina {cfg['vina_version']}, found {version}")
    if not 1 <= jobs <= 10:
        raise ValueError("--jobs must be between 1 and 10")
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=False)
    inputs = output / "inputs"
    shutil.copytree(ROOT / "data/prepared", inputs)
    shutil.copy2(ROOT / "configs/protocol.json", output / "protocol.json")
    box = json.loads((inputs / "box.json").read_text())
    manifest = {"started_utc": datetime.now(timezone.utc).isoformat(), "vina_version": version,
                "vina_sha256": digest(executable), "status": "running", "jobs": jobs,
                "versions": {x: importlib.metadata.version(x) for x in ("rdkit", "meeko", "numpy", "matplotlib")},
                "sha256": {str(p.relative_to(output)).replace('\\', '/'): digest(p) for p in sorted(output.rglob('*')) if p.is_file()},
                "source_sha256": {str(p.relative_to(ROOT)).replace('\\', '/'): digest(p) for p in sorted((ROOT/'scripts').glob('*.py'))}}
    write_json(output / "manifest.json", manifest)

    def dock(seed):
        directory = output / f"seed_{seed}"
        directory.mkdir()
        cmd = [executable, "--receptor", str(inputs/"receptor.pdbqt"), "--ligand", str(inputs/"ligand_independent.pdbqt"),
               "--out", str(directory/"poses.pdbqt"), "--seed", str(seed)]
        for key in ("cpu", "exhaustiveness", "num_modes", "energy_range"):
            cmd += ["--" + key, str(cfg[key])]
        for i, axis in enumerate("xyz"):
            cmd += ["--center_" + axis, str(box["center"][i]), "--size_" + axis, str(box["size"][i])]
        write_json(directory / "command.json", cmd)
        with (directory / "vina.log").open("w") as log:
            subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=True)
        print(f"Finished seed {seed}", flush=True)

    try:
        with ThreadPoolExecutor(max_workers=jobs) as pool:
            list(pool.map(dock, cfg["seeds"]))
        for path in output.glob("seed_*/*"):
            manifest["sha256"][str(path.relative_to(output)).replace('\\', '/')] = digest(path)
        manifest["status"] = "complete"
        manifest["finished_utc"] = datetime.now(timezone.utc).isoformat()
        write_json(output / "manifest.json", manifest)
    except Exception:
        manifest["status"] = "failed"
        write_json(output / "manifest.json", manifest)
        raise
    from analyze import analyze
    analyze(output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "experiments/rerun1")
    parser.add_argument("--vina", default="vina")
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--analyze-only", action="store_true")
    args = parser.parse_args()
    if args.analyze_only:
        from analyze import analyze
        analyze(args.output)
    else:
        run(args.output, args.vina, args.jobs)
