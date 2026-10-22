"""Backup: Kelvin–Helmholtz cooling-time ratio at fixed T_int (diagnostic).

    python figures/plot_bf3_fixed_tint.py config/fiducial.yaml

At fixed T_int the internal luminosity L_int = 4πR_p²σT_int⁴ is imposed by the
boundary condition, so a more opaque, more extended planet has a larger L and a
shorter τ_KH = |E|/L. This ratio is therefore not a measure of the opacity effect
on cooling; the fixed-interior-entropy comparison (F3) is used instead. Shown for
completeness with the caveat stated.
"""

from __future__ import annotations

import sys
from dataclasses import replace
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # run as script or -m

import matplotlib.pyplot as plt

from eos import AquaEOS
from figures._common import SOURCE_COLOR, SOURCE_LABEL, available_sources, band_caption, \
    comparison_band, load_cfg, save, set_style
from structure import build_structure

F_ENV_GRID = (0.01, 0.02, 0.03, 0.05, 0.10)


def main(argv=None):
    cfg = load_cfg(argv if argv is not None else sys.argv[1:])
    set_style()
    eos = AquaEOS.from_config(cfg)
    t_int = float(cfg.interior.t_int_k[len(cfg.interior.t_int_k) // 2])
    t_eq = float(cfg.irradiation.t_eq_k[len(cfg.irradiation.t_eq_k) // 2])
    others = available_sources(cfg, ["composite", "pokazatel", "fp"])
    band = comparison_band(cfg, ["hitemp", *others])

    fig, ax = plt.subplots(figsize=(7.0, 5.0))
    for name in others:
        ratio, reliable = [], True
        for f_env in F_ENV_GRID:
            cf = replace(cfg, envelope=replace(cfg.envelope, f_env=f_env))
            base = build_structure(cf, t_eq, t_int, "hitemp", eos=eos, n_sample=400, band=band)
            src = build_structure(cf, t_eq, t_int, name, eos=eos, n_sample=400, band=band)
            ratio.append(src.tau_kh_s / base.tau_kh_s)
            reliable &= src.clamp_reliable
        ax.plot(F_ENV_GRID, ratio, marker="o", ms=5, color=SOURCE_COLOR[name],
                ls="-" if reliable else "--",
                label=f"{SOURCE_LABEL[name]} / HITEMP" + ("" if reliable else "  [unreliable]"))
    ax.axhline(1.0, color="k", lw=0.8, ls=":")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("$f_{env}$")
    ax.set_ylabel(r"$\tau_{KH}$(source) / $\tau_{KH}$(HITEMP) at fixed $T_{int}$")
    ax.set_title("τ_KH ratio to HITEMP at fixed T_int (diagnostic)\n"
                 f"L_int imposed by T_int, so this ratio is not the opacity cooling "
                 f"effect; $T_{{eq}}$={t_eq:.0f} K, $T_{{int}}$={t_int:.0f} K, {band_caption(band)}",
                 fontsize=9)
    ax.legend(fontsize=8.5)
    path = save(fig, cfg, "BF3_fixed_tint")
    print(f"wrote {path}")
    return path


if __name__ == "__main__":
    main()
