"""M4 (deck): planets in opacity space — where databases fail, where the
correction lives, where planets travel.

    python figures/plot_m4_opacity_space.py config/fiducial.yaml

The (T, ρ) plane showing: the line-list validity region (HITEMP grid), the FP
grid extent, the composite blend window, and the fiducial planet's T(ρ)
structure trajectory at two epochs (young T_int=75 K, old T_int=20 K) with the
RCBs marked. Background is shaded by the FP/line-list Planck-mean ratio (where
FP is valid, ρ≥109), unifying the argument in one frame.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # run as script or -m

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle

from config import load_config
from eos import AquaEOS
from figures._common import band_caption, save, set_style
from opacity import get_source
from structure import build_structure
from structure.model import FP_GRID_FLOOR_SI


def main(argv=None):
    cfg = load_config((argv or sys.argv[1:])[0])
    set_style()
    eos = AquaEOS.from_config(cfg)
    band = tuple(cfg.opacity.band_cm)
    fp = get_source("fp", cfg)
    ht = get_source("hitemp", cfg)

    fig, ax = plt.subplots(figsize=(8.2, 6.2))

    # Background: FP/HITEMP Planck ratio where FP is valid (ρ ≥ FP floor).
    rr = np.geomspace(FP_GRID_FLOOR_SI, 2200.0, 60)
    tt = np.linspace(1000.0, 3000.0, 60)
    R, T = np.meshgrid(rr, tt)
    ratio = fp.kappa_P(T.ravel(), R.ravel()) / ht.kappa_P(T.ravel(), R.ravel())
    ratio = ratio.reshape(R.shape)
    pcm = ax.pcolormesh(R, T, ratio, shading="gouraud", cmap="magma", vmin=1.0, vmax=3.0)
    cb = fig.colorbar(pcm, ax=ax, label="FP / line-list Planck-mean ratio")

    # Validity / grid regions.
    ax.add_patch(Rectangle((0.1, 500), 1000 - 0.1, 3000 - 500, fill=False,
                           ec="#009E73", lw=2.0, ls="-", label="line-list (HITEMP) grid"))
    ax.add_patch(Rectangle((FP_GRID_FLOOR_SI, 1000), 2229 - FP_GRID_FLOOR_SI, 3000 - 1000,
                           fill=False, ec="#0072B2", lw=2.0, ls="--", label="FP grid"))
    lo, hi = cfg.opacity.composite_blend_rho
    ax.axvspan(lo, hi, color="#CC79A7", alpha=0.18, label="blend window")

    # Fiducial T(ρ) trajectories at two epochs, RCBs marked.
    for t_int, style, lab in [(75.0, "-", "young: T_int=75 K"), (20.0, ":", "old: T_int=20 K")]:
        res = build_structure(cfg, 500.0, t_int, "composite", eos=eos, n_sample=1200)
        p = res.profile
        ax.plot(p["rho_si"], p["temperature_k"], style, color="k", lw=2.0, label=lab)
        for pr, rho_rcb in zip(res.rcb_pressures_pa, res.rcb_densities_si):
            t_rcb = float(np.interp(np.log(pr), np.log(p["pressure_pa"]), p["temperature_k"]))
            ax.plot(rho_rcb, t_rcb, "o", ms=8, mfc="white", mec="k", mew=1.5, zorder=5)

    ax.axvline(FP_GRID_FLOOR_SI, color="#0072B2", lw=0.8, ls=":")
    ax.set_xscale("log")
    ax.set_xlim(1e-3, 2300); ax.set_ylim(500, 3200)
    ax.set_xlabel("Density ρ [kg/m³]"); ax.set_ylabel("Temperature [K]")
    ax.set_title("Planets in opacity space — databases (line-list) at low ρ, the "
                 "FP correction at high ρ,\nand where planets travel. "
                 "White ○ = RCB. Background = FP/line-list Planck ratio (ρ≥109)\n"
                 f"sources: FP + HITEMP; {band_caption(band)}", fontsize=9)
    ax.legend(fontsize=8.3, loc="upper left", framealpha=0.9)
    path = save(fig, cfg, "M4_opacity_space")
    print(f"wrote {path}")
    return path


if __name__ == "__main__":
    main()
