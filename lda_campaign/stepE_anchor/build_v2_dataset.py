#!/usr/bin/env python3
"""Build train_v2.xyz = base training file + harvest + anchor labelled frames.

Detects the energy/forces key convention from the base file's first frame
and writes the new frames identically, so training sees one homogeneous
format. Skips unconverged labels; prints per-frame max|F| so the
close-approach frames are visible (expected large, that is the point).

    python3 build_v2_dataset.py \
        --base-train <path to v1 training xyz> \
        --workdir /data/phys-plasma-dft/shug7909/lda_pipeline \
        --out /data/phys-plasma-dft/shug7909/lda_pipeline/train_v2.xyz
"""
import argparse
import pathlib
import re
import numpy as np
from ase.io import read, write
from ase import Atoms

H_EV, BOHR = 27.211386245988, 0.529177210903


def parse_label(outfile):
    t = pathlib.Path(outfile).read_text(errors="ignore")
    if "SCF run converged" not in t:
        return None, None
    Em = re.findall(r"ENERGY\|\s+Total FORCE_EVAL.*?:\s+([-\d.]+)", t)
    Fm = re.search(r"ATOMIC FORCES.*?\n(.*?)SUM OF ATOMIC FORCES", t, re.S)
    if not Em or not Fm:
        return None, None
    F = np.array([[float(x) for x in l.split()[3:6]]
                  for l in Fm.group(1).splitlines()
                  if l.split() and l.split()[0].isdigit()]) * H_EV / BOHR
    return float(Em[-1]) * H_EV, F


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base-train", required=True)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--subdirs", nargs="+", default=["harvest", "anchor"])
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    work = pathlib.Path(args.workdir)

    base = read(args.base_train, index=":")
    first = base[0]
    ekey = next((k for k in ("REF_energy", "energy", "dft_energy")
                 if k in first.info), None)
    fkey = next((k for k in ("REF_forces", "forces", "dft_forces")
                 if k in first.arrays), None)
    if ekey is None or fkey is None:
        raise SystemExit(f"could not detect keys in base file: "
                         f"info={list(first.info)} arrays={list(first.arrays)}")
    print(f"base: {len(base)} frames; keys: {ekey} / {fkey}")

    new, skipped = [], 0
    for sub in args.subdirs:
        for d in sorted((work / sub).glob("T*_rho*")):
            for w in sorted(p for p in d.glob("anchor_*") if p.is_dir()):
                E, F = parse_label(w / "label.out")
                if E is None:
                    skipped += 1
                    print(f"  skip (unconverged/missing): {sub}/{d.name}/{w.name}")
                    continue
                lines = (w / "frame.xyz").read_text().splitlines()[2:]
                sym = [l.split()[0] for l in lines]
                pos = np.array([l.split()[1:4] for l in lines], float)
                L = float(re.search(r"ABC\s+([\d.Ee+-]+)",
                                    (w / "label.inp").read_text()).group(1))
                a = Atoms(sym, positions=pos, cell=[L, L, L], pbc=True)
                a.info[ekey] = E
                a.info["config_type"] = f"{sub}_{d.name}"
                a.arrays[fkey] = F
                new.append(a)
                print(f"  + {sub}/{d.name}/{w.name}: max|F| = "
                      f"{np.abs(F).max():7.2f} eV/A")
    print(f"\nnew frames: {len(new)} (skipped {skipped})")
    if not new:
        raise SystemExit("nothing collected - check subdir names / label jobs")
    write(args.out, base + new)
    print(f"wrote {args.out}: {len(base) + len(new)} frames total")


if __name__ == "__main__":
    main()
