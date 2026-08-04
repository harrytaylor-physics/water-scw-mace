"""P2 (deck slide ii): FP vs line-list opacity crossover with density.

    python figures/plot_p2_crossover.py config/fiducial.yaml

FP/HITEMP Planck- and Rosseland-mean ratios from the independent crossover
analysis (crossover_hitemp.json), plotted against density and coloured by
temperature. Shows the line-list vs first-principles disagreement increasing
with density across the grid.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # run as script or -m

import matplotlib.pyplot as plt
import numpy as np

from config import load_config
from figures._common import save, set_style


def main(argv=None):
    cfg = load_config((argv or sys.argv[1:])[0])
    set_style()
    data = json.load(open(cfg.fp_opacity_path / "crossover_hitemp.json"))["results"]
    rho = np.array([r["rho"] for r in data])
    T = np.array([r["T"] for r in data])
    pr = np.array([r["M2_planck_ratio"] for r in data])
    rr = np.array([r["M2_rosseland_ratio"] for r in data])

    fig, (axp, axr) = plt.subplots(1, 2, figsize=(10.6, 4.6), constrained_layout=True, sharey=False)
    for ax, y, lab in ((axp, pr, "Planck"), (axr, rr, "Rosseland")):
        sc = ax.scatter(rho, y, c=T, cmap="viridis", s=42, edgecolor="k", linewidth=0.3)
        ax.axhline(1.0, color="k", lw=0.8, ls=":")
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlabel("Density ρ [kg/m³]")
        ax.set_ylabel(f"FP / HITEMP {lab}-mean ratio")
        ax.set_title(f"{lab} mean")
    cb = fig.colorbar(sc, ax=axr, label="Temperature [K]")

    fig.suptitle("FP / HITEMP opacity ratio versus density (crossover_hitemp.json)\n"
                 "RADIS HITEMP2010, self-broadened Voigt; ratio rises with density; "
                 "Planck 1.4–18.6, Rosseland 1.7–39.6 across the 27-point grid", fontsize=9)
    path = save(fig, cfg, "P2_crossover")
    print(f"wrote {path}")
    # report the ranges for the slide bullets
    print(f"Planck ratio range: {pr.min():.2f}-{pr.max():.2f}; Rosseland: {rr.min():.2f}-{rr.max():.2f}")
    return path


if __name__ == "__main__":
    main()
