"""D1: reachability of the FP-radiative regime versus core mass (A28/A30).

    python figures/plot_d1_fp_domain.py config/fiducial.yaml

Two curves per core mass, for the fiducial envelope. The lower curve is the
internal temperature at which the line-list radiative zone first reaches
rho = 109 kg/m^3 (the FP grid floor); a planet must cool below it to put the
first-principles opacity into the radiative interior. The upper curve is the
radiogenic internal-temperature floor set by an Earth-like core; a planet does
not cool below it. The floor lies above the reach threshold at every mass, so
the FP-radiative regime is not reached for these steam sub-Neptunes.
"""

from __future__ import annotations

import sys
from dataclasses import replace
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # run as script or -m

import matplotlib.pyplot as plt
import numpy as np

from constants import R_EARTH
from eos import AquaEOS
from evolution import t_int_floor_k
from figures._common import band_caption, load_cfg, save, set_style
from structure import build_structure
from structure.model import FP_GRID_FLOOR_SI

M_C_GRID = np.array([2.0, 3.0, 5.0, 8.0, 12.0, 20.0])
T_INT_PROBE = np.array([15.0, 20.0, 25.0, 30.0, 35.0, 40.0])
FLOOR_AGE_GYR = 10.0


def main(argv=None):
    cfg = load_cfg(argv if argv is not None else sys.argv[1:])
    set_style()
    eos = AquaEOS.from_config(cfg)
    t_eq = float(cfg.irradiation.t_eq_k[0])
    band = tuple(cfg.opacity.band_cm)

    reach, floor = [], []
    for m_c in M_C_GRID:
        c = replace(cfg, core=replace(cfg.core, mass_mearth=float(m_c)))
        rho, rp0 = [], None
        for ti in T_INT_PROBE:
            r = build_structure(c, t_eq, float(ti), "hitemp", eos=eos, n_sample=500, band=band)
            rho.append(r.max_radiative_rho_si)
            rp0 = rp0 or r.radius_p_rearth
        rho = np.array(rho)
        # max radiative density falls as T_int rises; invert for T_int at 109.
        reach.append(float(np.interp(FP_GRID_FLOOR_SI, rho[::-1], T_INT_PROBE[::-1])))
        floor.append(t_int_floor_k(float(m_c), rp0 * R_EARTH, FLOOR_AGE_GYR))
    reach, floor = np.array(reach), np.array(floor)

    fig, ax = plt.subplots(figsize=(7.6, 5.4))
    ax.fill_between(M_C_GRID, reach, floor, color="grey", alpha=0.18)
    ax.plot(M_C_GRID, reach, "o-", color="black", lw=2.0, ms=6,
            label="reach: $T_{int}$ for radiative zone to reach ρ = 109 kg/m³")
    ax.plot(M_C_GRID, floor, "s--", color="#D55E00", lw=2.0, ms=6,
            label="radiogenic floor: minimum $T_{int}$ at 10 Gyr")
    ax.plot(cfg.core.mass_mearth, np.interp(cfg.core.mass_mearth, M_C_GRID, floor),
            "*", color="black", ms=15, mec="white", mew=0.8)

    ax.set_xscale("log")
    ax.set_xticks(M_C_GRID)
    ax.set_xticklabels([f"{int(m)}" for m in M_C_GRID])
    ax.set_xlabel("Core mass $M_c$ [$M_⊕$]")
    ax.set_ylabel("Internal temperature $T_{int}$ [K]")
    ax.set_title("Radiogenic floor sits above the FP-reach threshold at every mass")
    ax.legend(fontsize=9, loc="upper left")
    gap = float(np.mean(floor - reach))
    fig.text(0.5, -0.01,
             f"$T_{{eq}}$={t_eq:.0f} K, $f_{{env}}$={cfg.envelope.f_env:.0%}, {band_caption(band)}; "
             f"floor above reach by {gap:.1f} K on average (grey band); star = fiducial",
             ha="center", va="top", fontsize=9)
    path = save(fig, cfg, "D1_fp_domain")
    print(f"wrote {path}; mean(floor - reach) = {gap:+.2f} K")
    return path


if __name__ == "__main__":
    main()
