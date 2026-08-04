"""M3 (deck): planetary-model anatomy — annotated T(P) with the governing physics.

    python figures/plot_m3_anatomy.py config/fiducial.yaml

A real fiducial composite profile, annotated with the model's anatomy: the
irradiated Guillot radiative outer layer, the radiative–convective boundary
(Schwarzschild), the convective adiabat on the AQUA EOS, the rocky core, and
where each opacity source applies (line-list dilute → blend window → FP dense).
This is the model-physics slide figure.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # run as script or -m

import matplotlib.pyplot as plt
import numpy as np

from config import load_config
from eos import AquaEOS
from figures._common import band_caption, set_style, save
from structure import build_structure


def main(argv=None):
    cfg = load_config((argv or sys.argv[1:])[0])
    set_style()
    eos = AquaEOS.from_config(cfg)
    res = build_structure(cfg, 600.0, 100.0, "composite", eos=eos, n_sample=1500)
    p = res.profile
    P_bar, T, rho = p["pressure_pa"] / 1e5, p["temperature_k"], p["rho_si"]
    conv = p["is_convective"]

    fig, ax = plt.subplots(figsize=(7.6, 6.4))
    # Colour the curve by convective (adiabat) vs radiative (Guillot).
    ax.plot(T[~conv], P_bar[~conv], ".", ms=2.4, color="#0072B2", label="radiative (Guillot T(τ))")
    ax.plot(T[conv], P_bar[conv], ".", ms=2.4, color="#D55E00", label="convective adiabat (AQUA EOS)")
    ax.set_yscale("log"); ax.invert_yaxis()
    ax.set_xlabel("Temperature [K]"); ax.set_ylabel("Pressure [bar]")

    # RCB marker(s) — the Schwarzschild boundary.
    for pr in res.rcb_pressures_pa:
        tr = float(np.interp(np.log(pr), np.log(p["pressure_pa"]), T))
        ax.axhline(pr / 1e5, color="k", ls="--", lw=1.0)
        ax.annotate("RCB — Schwarzschild:\n∇ = min(∇_ad, ∇_rad)", xy=(tr, pr / 1e5),
                    xytext=(tr + 800, pr / 1e5 / 12), fontsize=8.5,
                    arrowprops=dict(arrowstyle="->", lw=1.0))

    # Opacity-applicability bands by density (right-hand annotations).
    lo, hi = cfg.opacity.composite_blend_rho
    def p_at_rho(target):
        i = int(np.argmin(np.abs(rho - target)))
        return P_bar[i]
    ax.axhspan(p_at_rho(lo), p_at_rho(hi), color="#CC79A7", alpha=0.15)
    ax.text(0.97, 0.5, "line-list (dilute)\nρ < %g" % lo, transform=ax.transAxes,
            ha="right", va="center", fontsize=8.5, color="#333")
    ax.text(0.97, 0.30, "blend window\nρ = %g–%g kg/m³" % (lo, hi), transform=ax.transAxes,
            ha="right", va="center", fontsize=8.5, color="#CC79A7")
    ax.text(0.97, 0.13, "FP (dense)\nρ > %g" % hi, transform=ax.transAxes,
            ha="right", va="center", fontsize=8.5, color="#0072B2")

    # Photosphere + core annotations.
    ax.annotate("irradiated photosphere\n(Guillot two-stream, T_eq)",
                xy=(T[0], P_bar[0]), xytext=(T[0] + 500, P_bar[0] * 3),
                fontsize=8.5, arrowprops=dict(arrowstyle="->", lw=1.0))
    ax.annotate("rocky core (Zeng M–R)\nr(M_c) = R_c", xy=(T[-1], P_bar[-1]),
                xytext=(T[-1] - 2600, P_bar[-1] / 1.5), fontsize=8.5,
                arrowprops=dict(arrowstyle="->", lw=1.0))

    ax.text(0.03, 0.03,
            "hydrostatic: dP/dm = −Gm/4πr⁴   ·   dr/dm = 1/4πr²ρ\n"
            "L_int = 4πR_p²σT_int⁴ constant;  EOS: AQUA (Haldemann+ 2020)",
            transform=ax.transAxes, fontsize=8, color="#555", va="bottom",
            bbox=dict(boxstyle="round", fc="white", ec="#aaa", alpha=0.9))

    ax.set_title(f"Planetary-model anatomy — fiducial composite, $T_{{eq}}$=600 K, "
                 f"$T_{{int}}$=100 K\nsource: composite (HITEMP dilute + FP dense), "
                 f"{band_caption(res.band_cm)}; radiative clamp {res.clamp_fraction_radiative:.0%} "
                 f"({'reliable' if res.clamp_reliable else 'UNRELIABLE'})", fontsize=9)
    ax.legend(fontsize=8.5, loc="lower left", bbox_to_anchor=(0.0, 0.12))
    path = save(fig, cfg, "M3_anatomy")
    print(f"wrote {path}")
    return path


if __name__ == "__main__":
    main()
