"""D1: the FP-relevant domain in planet space (A28).

    python figures/plot_d1_fp_domain.py config/fiducial.yaml

Maps where the (pure line-list) radiative zone reaches ρ ≥ FP's grid floor
(109 kg/m³) — the region where first-principles opacity would change the deep
structure. Scanned over (M_c, f_env) at fixed T_eq for a cold "mature" T_int
(the RCB deepens as T_int drops, so a cold snapshot is the most FP-favourable).
The fiducial planet is marked. Also reports, per cell, the number of RCBs for
line-list vs composite (FP closes deep radiative windows, A28).
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
from figures._common import band_caption, comparison_band, load_cfg, save, set_style
from structure import build_structure
from structure.model import FP_GRID_FLOOR_SI

M_C_GRID = np.array([2.0, 3.0, 5.0, 8.0, 12.0, 20.0])       # Earth masses
F_ENV_GRID = np.array([0.01, 0.02, 0.05, 0.10, 0.25, 0.50])  # envelope fraction
T_INT_COLD = 20.0  # K, a mature/old planet (RCB deepest -> most FP-relevant)


def main(argv=None):
    cfg = load_cfg(argv if argv is not None else sys.argv[1:])
    set_style()
    eos = AquaEOS.from_config(cfg)
    t_eq = float(cfg.irradiation.t_eq_k[0])
    band = comparison_band(cfg, ["hitemp", "fp"])

    max_rad = np.full((F_ENV_GRID.size, M_C_GRID.size), np.nan)
    d_rp_pct = np.full_like(max_rad, np.nan)  # 100·|R_p(comp) − R_p(LL)|/R_p
    d_nrcb = np.zeros_like(max_rad)  # (#RCB line-list) − (#RCB composite)
    radiogenic_truncated = np.zeros_like(max_rad, dtype=bool)  # T_int=20 unreachable (A30)
    for i, f_env in enumerate(F_ENV_GRID):
        rp_guess = None
        for m_c in sorted(M_C_GRID):  # ascending mass warm-starts
            c = replace(cfg, core=replace(cfg.core, mass_mearth=float(m_c)),
                        envelope=replace(cfg.envelope, f_env=float(f_env)))
            j = list(M_C_GRID).index(m_c)
            try:
                ht = build_structure(c, t_eq, T_INT_COLD, "hitemp", eos=eos,
                                     n_sample=600, band=band, rp_guess_rearth=rp_guess)
                co = build_structure(c, t_eq, T_INT_COLD, "composite", eos=eos,
                                     n_sample=600, band=band, rp_guess_rearth=ht.radius_p_rearth)
                max_rad[i, j] = ht.max_radiative_rho_si
                d_rp_pct[i, j] = 100.0 * abs(co.radius_p_rearth - ht.radius_p_rearth) / ht.radius_p_rearth
                d_nrcb[i, j] = len(ht.rcb_pressures_pa) - len(co.rcb_pressures_pa)
                # Radiogenic floor (A30): is this cold T_int reachable at ~10 Gyr?
                floor = t_int_floor_k(float(m_c), ht.radius_p_rearth * R_EARTH, age_gyr=10.0)
                radiogenic_truncated[i, j] = floor > T_INT_COLD
                rp_guess = ht.radius_p_rearth
            except Exception as exc:  # noqa: BLE001
                print(f"  (M_c={m_c}, f_env={f_env}) failed: {exc}")

    fig, ax = plt.subplots(figsize=(7.2, 5.6))
    pcm = ax.pcolormesh(M_C_GRID, F_ENV_GRID * 100, d_rp_pct, shading="nearest", cmap="magma")
    fig.colorbar(pcm, ax=ax, label="|ΔR_p| composite vs line-list [%]  (FP structural impact)")
    # FP-relevance boundary: where the line-list radiative zone reaches ρ = 109.
    cs = ax.contour(M_C_GRID, F_ENV_GRID * 100, max_rad, levels=[FP_GRID_FLOOR_SI],
                    colors="cyan", linewidths=2.5)
    ax.clabel(cs, fmt={FP_GRID_FLOOR_SI: "ρ_rad=109"}, fontsize=9)
    # Cells where the line list opens a deep radiative window the composite closes.
    fi, fj = np.where(d_nrcb > 0)
    if fi.size:
        ax.plot(M_C_GRID[fj], F_ENV_GRID[fi] * 100, "w*", ms=13, mec="k", mew=0.5,
                ls="none", label="LL opens deep window")
    # Cells where radiogenic heating floors T_int above 20 K: this cold snapshot
    # is unreachable, so the FP-relevant regime is radiogenically truncated (A30).
    ti, tj = np.where(radiogenic_truncated)
    if ti.size:
        ax.plot(M_C_GRID[tj], F_ENV_GRID[ti] * 100, "x", color="white", ms=9, mew=1.8,
                ls="none", label="radiogenic floor > 20 K (unreachable)")
    ax.plot(cfg.core.mass_mearth, cfg.envelope.f_env * 100, "o", ms=14, ls="none",
            mfc="none", mec="lime", mew=2.5, label="fiducial")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Core mass $M_c$ [$M_⊕$]")
    ax.set_ylabel("Envelope fraction $f_{env}$ [%]")
    ax.set_title(f"FP structural impact across planet space — $T_{{eq}}$={t_eq:.0f} K, "
                 f"$T_{{int}}$={T_INT_COLD:.0f} K (mature)\n"
                 f"cyan = line-list radiative zone reaches FP floor (ρ≥109); {band_caption(band)}",
                 fontsize=9.5)
    ax.legend(loc="lower left", fontsize=8.5)
    path = save(fig, cfg, "D1_fp_domain")
    print(f"wrote {path}")
    n_relevant = int(np.nansum(max_rad >= FP_GRID_FLOOR_SI))
    print(f"FP-relevant cells (ρ_rad≥109): {n_relevant}/{max_rad.size}; "
          f"deep-window cells: {int((d_nrcb>0).sum())}; "
          f"max |ΔR_p| = {np.nanmax(d_rp_pct):.2f}%")
    return path


if __name__ == "__main__":
    main()
