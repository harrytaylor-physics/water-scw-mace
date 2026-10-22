"""Boundary sensitivity scans: γ and matching-pressure cap (SPEC §5, §10).

    python figures/plot_sensitivity.py config/fiducial.yaml

Documents that R_p is robust to both knobs (see docs/boundary_sensitivity.md).
Uses the config source if it builds, else the constant-κ reference.
"""

from __future__ import annotations

import sys
from dataclasses import replace
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # run as script or -m

import matplotlib.pyplot as plt

from eos import AquaEOS
from figures._common import available_sources, load_cfg, register_constant_reference, save, set_style
from structure import build_structure

GAMMAS = (0.25, 0.5, 1.0, 2.0, 4.0)
CAPS_BAR = (1.0, 3.0, 10.0, 30.0, 100.0)


def main(argv=None):
    cfg = load_cfg(argv if argv is not None else sys.argv[1:])
    set_style()
    eos = AquaEOS.from_config(cfg)
    register_constant_reference(cfg)
    source = (available_sources(cfg, [cfg.opacity.source]) or ["constant"])[0]
    t_eq = float(cfg.irradiation.t_eq_k[len(cfg.irradiation.t_eq_k) // 2])
    t_int = float(cfg.interior.t_int_k[len(cfg.interior.t_int_k) // 2])

    fig, (axg, axc) = plt.subplots(1, 2, figsize=(11.0, 4.6), constrained_layout=True)

    rp_g, rcb_g = [], []
    for g in GAMMAS:
        r = build_structure(replace(cfg, irradiation=replace(cfg.irradiation, gamma=g)),
                            t_eq, t_int, source, eos=eos)
        rp_g.append(r.radius_p_rearth)
        rcb_g.append(max(r.rcb_pressures_pa) / 1e5)
    axg.plot(GAMMAS, rp_g, "o-", color="#0072B2", label="$R_p$")
    axg.set_xscale("log")
    axg.set_xlabel(r"$\gamma = \kappa_{vis}/\kappa_{th}$")
    axg.set_ylabel("$R_p$ [$R_⊕$]", color="#0072B2")
    axg.tick_params(axis="y", colors="#0072B2")
    axg.set_title(f"γ sensitivity: ΔR_p = {100*(max(rp_g)-min(rp_g))/min(rp_g):.1f}%")
    axg2 = axg.twinx()
    axg2.plot(GAMMAS, rcb_g, "s--", color="#D55E00", alpha=0.85)
    axg2.set_ylabel("RCB [bar]", color="#D55E00")
    axg2.tick_params(axis="y", colors="#D55E00")
    axg2.set_yscale("log")

    rp_c, matched = [], []
    for cap in CAPS_BAR:
        r = build_structure(replace(cfg, boundary=replace(cfg.boundary, match_pressure_bar=cap)),
                            t_eq, t_int, source, eos=eos)
        rp_c.append(r.radius_p_rearth)
        matched.append(r.match_pressure_actual_bar)
    d_rp = 100 * (max(rp_c) - min(rp_c)) / min(rp_c)
    axc.plot(CAPS_BAR, rp_c, "o-", color="#009E73")
    mean = sum(rp_c) / len(rp_c)
    axc.set_ylim(mean * 0.985, mean * 1.015)  # fixed ~1.5% band so flatness reads as flat
    axc.set_xscale("log")
    axc.set_xlabel("matching-pressure cap [bar]")
    axc.set_ylabel("$R_p$ [$R_⊕$]")
    axc.set_title(f"matching-cap sensitivity: ΔR_p = {d_rp:.2f}%")
    note = ("adaptive match hits the RCB\n(cap does not bind)" if d_rp < 0.05
            else "cap binds below the RCB")
    axc.text(0.5, 0.1, note, transform=axc.transAxes, ha="center", fontsize=9,
             color="#555", style="italic")

    fig.suptitle(f"Boundary sensitivity — source={source}, $T_{{eq}}$={t_eq:.0f} K, "
                 f"$T_{{int}}$={t_int:.0f} K")
    path = save(fig, cfg, "sensitivity_gamma_matching")
    print(f"wrote {path}  (source={source})")
    return path


if __name__ == "__main__":
    main()
