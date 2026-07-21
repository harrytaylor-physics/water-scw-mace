#!/usr/bin/env python3
"""Learning-curve protocol: make nested subsets + a runner script, then
analyse held-out parity per subset model.

Stage 1 (default): writes subset_XX.xyz + holdout.xyz + run_all.sh
Stage 2 (--analyze DIR): evaluates each trained model on holdout.xyz.
"""
import argparse, pathlib, subprocess, sys
import numpy as np

def stage1(args):
    from ase.io import read, write
    frames = read(args.train, index=":")
    rng = np.random.default_rng(args.seed)
    idx = rng.permutation(len(frames))
    n_hold = int(len(frames) * args.holdout)
    hold, pool = idx[:n_hold], idx[n_hold:]
    root = pathlib.Path(args.out); root.mkdir(exist_ok=True)
    write(root / "holdout.xyz", [frames[i] for i in hold])
    lines = ["#!/bin/bash", "set -e"]
    for f in args.fractions:
        n = int(len(pool) * f)
        sub = root / f"subset_{int(f*100):03d}.xyz"
        write(sub, [frames[i] for i in pool[:n]])
        name = f"lc_{int(f*100):03d}"
        lines.append(
            f"mace_run_train --name {name} --foundation_model {args.foundation} "
            f"--train_file {sub} --valid_fraction 0.1 "
            f"--energy_key REF_energy --forces_key REF_forces --E0s average "
            f"--r_max 6.0 --batch_size 8 --max_num_epochs {args.epochs} "
            f"--swa --start_swa {int(args.epochs*0.75)} --seed {args.seed} "
            f"--device cuda --model_dir {root}/models")
        print(f"{sub.name}: {n} frames")
    (root / "run_all.sh").write_text("\n".join(lines) + "\n")
    print(f"holdout: {n_hold} frames\nrun: bash {root}/run_all.sh  (GPU), "
          f"then: learning_curve.py {args.train} --analyze {root}")

def stage2(args):
    from ase.io import read
    from mace.calculators import MACECalculator
    root = pathlib.Path(args.analyze)
    hold = read(root / "holdout.xyz", index=":")
    xs, ys = [], []
    for m in sorted((root / "models").glob("lc_*_swa.model")) or \
             sorted((root / "models").glob("lc_*.model")):
        frac = int(m.name.split("_")[1].split(".")[0].replace("swa","").strip("_") or 0)
        calc = MACECalculator(model_paths=str(m), device=args.device,
                              default_dtype="float64")
        errs = []
        for a in hold:
            b = a.copy(); b.calc = calc
            errs.append((b.get_forces() - a.arrays["REF_forces"]).ravel())
        rmse = float(np.sqrt(np.mean(np.concatenate(errs) ** 2)) * 1000)
        xs.append(frac); ys.append(rmse)
        print(f"{m.name}: {rmse:.1f} meV/A")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.plot(xs, ys, "ko-")
    ax.set_xlabel("training fraction (%)"); ax.set_ylabel("holdout force RMSE (meV/A)")
    fig.tight_layout(); fig.savefig(root / "learning_curve.png", dpi=170)
    print(f"figure -> {root/'learning_curve.png'}")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("train")
    ap.add_argument("--fractions", nargs="+", type=float, default=[0.25, 0.5, 0.75, 1.0])
    ap.add_argument("--holdout", type=float, default=0.1)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--epochs", type=int, default=200)
    ap.add_argument("--foundation", default="mace-mpa-0-medium.model")
    ap.add_argument("--out", default="learning_curve_runs")
    ap.add_argument("--analyze", default=None)
    ap.add_argument("--device", default="cuda")
    args = ap.parse_args()
    stage2(args) if args.analyze else stage1(args)

if __name__ == "__main__":
    main()
