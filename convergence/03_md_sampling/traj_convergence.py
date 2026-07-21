#!/usr/bin/env python3
"""Running-average + block-error analysis of a thermo file.
Usage: traj_convergence.py thermo.dat --col 1 --time-col 0 --dt-fs 0.5
(col indices 0-based; --col is the observable, e.g. T or pressure)."""
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def block_error(x, nb=20):
    n = len(x) // nb
    if n < 2: return np.std(x) / max(len(x), 1) ** 0.5
    means = [x[i*n:(i+1)*n].mean() for i in range(nb)]
    return np.std(means, ddof=1) / nb ** 0.5

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("thermo"); ap.add_argument("--col", type=int, default=1)
    ap.add_argument("--time-col", type=int, default=0)
    ap.add_argument("--dt-fs", type=float, default=0.5)
    args = ap.parse_args()
    d = np.loadtxt(args.thermo)
    t_ps = d[:, args.time_col] * args.dt_fs / 1000
    x = d[:, args.col]
    run = np.cumsum(x) / np.arange(1, len(x) + 1)
    final = x[len(x)//3:]
    mu, err = final.mean(), block_error(final)
    tail = run[len(run)//3:]
    stable = np.all(np.abs(tail - mu) <= 3*err)
    # crude equilibration estimate: first index where remaining mean within 1 block-sigma
    eq = next((i for i in range(0, len(x), max(1, len(x)//200))
               if abs(x[i:].mean() - mu) < err*3), 0)
    print(f"mean (final 2/3): {mu:.4f} +/- {err:.4f} (block)")
    print(f"running mean stable within 3 sigma: {'YES' if stable else 'NO'}")
    print(f"equilibration estimate: ~{t_ps[min(eq, len(t_ps)-1)]:.1f} ps")
    fig, ax = plt.subplots(figsize=(7, 3.5))
    ax.plot(t_ps, x, "k-", lw=0.3, alpha=0.5)
    ax.plot(t_ps, run, "k-", lw=1.2)
    ax.axhspan(mu-err, mu+err, color="0.85")
    ax.set_xlabel("t (ps)"); ax.set_ylabel(f"column {args.col}")
    fig.tight_layout(); fig.savefig("traj_convergence.png", dpi=170)
    print("figure -> traj_convergence.png")

if __name__ == "__main__":
    main()
