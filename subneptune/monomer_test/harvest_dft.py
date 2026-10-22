#!/usr/bin/env python3
"""Read the 13 CP2K Wannier outputs into dft_dipoles.json, in Debye (A64).

    python harvest_dft.py --builder ~/scw/build_dipole_xyz_v3.py

The dipole is computed by the *same* ``dip()`` that built the training labels, so
the DFT side of this comparison is definitionally the same quantity the model was
fitted to. That function is obtained by reading ``build_dipole_xyz_v3.py`` and
executing it **only up to the line ``out=[]``**.

Why not import it. Everything after ``out=[]`` in that file is its harvest loop,
which rewrites the training set. Importing the module would run it. Truncating at
that marker takes the definitions and leaves the side effects behind — so this
script can never damage the training data, whatever else it does.

Needs numpy only, so it runs on a login node.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
TRUNCATE_AT = "out=[]"


def load_dip(builder: Path):
    """``dip`` from the label builder, with its harvest loop left unexecuted."""
    text = builder.read_text()
    marker = None
    for i, line in enumerate(text.splitlines()):
        if line.replace(" ", "").startswith(TRUNCATE_AT):
            marker = i
            break
    if marker is None:
        raise SystemExit(
            f"{builder} has no line starting '{TRUNCATE_AT}'. That marker is what "
            f"separates the definitions from the harvest loop that rewrites the "
            f"training set; without it this script will not execute the file at "
            f"all, rather than risk running the loop.")
    prefix = "\n".join(text.splitlines()[:marker])
    namespace: dict = {"__name__": "_dipole_defs"}
    exec(compile(prefix, str(builder), "exec"), namespace)  # noqa: S102
    if "dip" not in namespace:
        raise SystemExit(f"{builder}: no dip() defined before line {marker + 1}")
    print(f"[harvest] dip() taken from {builder} (first {marker} lines; the "
          f"harvest loop after '{TRUNCATE_AT}' was not executed)")
    return namespace["dip"]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--builder", type=Path,
                    default=Path("~/scw/build_dipole_xyz_v3.py"))
    ap.add_argument("--geoms", type=Path, default=HERE / "geoms.json")
    ap.add_argument("--run-dir", type=Path, default=HERE)
    ap.add_argument("--wannier-name", default="wannier")
    ap.add_argument("--out", type=Path, default=HERE / "dft_dipoles.json")
    args = ap.parse_args(argv)

    builder = args.builder.expanduser()
    if not builder.exists():
        raise SystemExit(
            f"{builder} not found. It is the label builder on ARC; its dip() is the "
            f"definition of the DFT dipole used for training and must be the one "
            f"used here.")
    dip = load_dip(builder)

    meta = json.loads(args.geoms.read_text())
    results, missing = [], []
    for entry in meta["geometries"]:
        name = entry["name"]
        wannier = args.run_dir / name / args.wannier_name
        if not wannier.exists():
            missing.append(name)
            continue
        mu = dip(str(wannier))
        mu = [float(v) for v in (mu.tolist() if hasattr(mu, "tolist") else mu)][:3]
        norm = float(sum(v * v for v in mu) ** 0.5)
        results.append({"name": name, "mode": entry["mode"], "sign": entry["sign"],
                        "reference": entry["reference"],
                        "displacement": entry["displacement"],
                        "displacement_unit": entry["displacement_unit"],
                        "mu_debye": mu, "mu_norm_debye": norm})
        print(f"  {name}  {entry['mode']:<13s} "
              f"mu = ({mu[0]:+.5f}, {mu[1]:+.5f}, {mu[2]:+.5f}) D  |mu| = {norm:.5f} D")

    if missing:
        print(f"\n  MISSING {len(missing)} of {len(meta['geometries'])}: "
              f"{' '.join(missing)}")
        print(f"  (no {args.wannier_name} file — those array elements have not run "
              f"or failed)")

    args.out.write_text(json.dumps({
        "meta": ("revPBE0/MLWF dipoles on the 13 monomer geometries (A64), via the "
                 "dip() of the label builder so this is the same quantity the model "
                 "was trained on. Units: Debye."),
        "builder": str(builder), "n_expected": len(meta["geometries"]),
        "n_found": len(results), "missing": missing, "results": results,
    }, indent=1) + "\n")
    print(f"\nwrote {args.out}  ({len(results)}/{len(meta['geometries'])} geometries)")
    return 0 if not missing else 1


if __name__ == "__main__":
    raise SystemExit(main())
