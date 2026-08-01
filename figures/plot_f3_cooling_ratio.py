"""F3: Kelvin–Helmholtz cooling-time ratio vs (f_env, T_eq) (SPEC §8).

    python figures/plot_f3_cooling_ratio.py config/fiducial.yaml

The headline is τ_FP / τ_LL ("planets stay puffy longer"). Until the line-list
tables land, the baseline is the constant-κ reference (τ_FP / τ_ref); the
baseline auto-switches to HITEMP when available. τ_KH = |E_env| / L_int.
"""

from __future__ import annotations

import sys
from dataclasses import replace
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # run as script or -m

import matplotlib.pyplot as plt
import numpy as np

from eos import AquaEOS
from figures._common import available_sources, load_cfg, register_constant_reference, save, set_style
from structure import build_structure

F_ENV_GRID = (0.01, 0.02, 0.03, 0.05, 0.10)


def _tau_kh(cfg, source, t_eq, t_int, f_env, eos):
    cfg_f = replace(cfg, envelope=replace(cfg.envelope, f_env=f_env))
    return build_structure(cfg_f, t_eq, t_int, source, eos=eos, n_sample=500).tau_kh_s


def main(argv=None):
    cfg = load_cfg(argv if argv is not None else sys.argv[1:])
    set_style()
    eos = AquaEOS.from_config(cfg)
    t_int = float(cfg.interior.t_int_k[len(cfg.interior.t_int_k) // 2])

    register_constant_reference(cfg)
    avail = available_sources(cfg, ["fp", "hitemp", "pokazatel", "composite"])
    numerator = "fp"
    # Prefer a real line-list baseline; HITEMP > POKAZATEL > constant placeholder.
    if "hitemp" in avail:
        baseline = "hitemp"
    elif "pokazatel" in avail:
        baseline = "pokazatel"
    else:
        baseline = "constant"

    fig, ax = plt.subplots(figsize=(6.4, 5.0))
    cmap = plt.get_cmap("viridis")
    t_eqs = cfg.irradiation.t_eq_k
    for k, t_eq in enumerate(t_eqs):
        ratio = []
        for f_env in F_ENV_GRID:
            tau_num = _tau_kh(cfg, numerator, float(t_eq), t_int, f_env, eos)
            tau_den = _tau_kh(cfg, baseline, float(t_eq), t_int, f_env, eos)
            ratio.append(tau_num / tau_den)
        ax.plot(F_ENV_GRID, ratio, marker="o", ms=5,
                color=cmap(k / max(1, len(t_eqs) - 1)), label=f"$T_{{eq}}$={t_eq:.0f} K")

    ax.axhline(1.0, color="k", lw=0.8, ls=":")
    ax.set_xscale("log")
    ax.set_xlabel("$f_{env}$")
    ax.set_ylabel(rf"$\tau_{{FP}} / \tau_{{{baseline.upper()}}}$")
    base_lbl = {"hitemp": "HITEMP", "pokazatel": "POKAZATEL (grid overlap; clamped)",
                "constant": "constant-κ ref (LL placeholder)"}[baseline]
    ax.set_title(f"KH cooling-time ratio vs baseline: {base_lbl}\n"
                 f"$T_{{int}}$={t_int:.0f} K, $M_c$={cfg.core.mass_mearth:.0f} $M_⊕$")
    # Honesty note: at FIXED T_int, L=4πR²σT_int⁴ makes the puffier (FP) planet
    # more luminous, so τ=|E|/L < 1 here does NOT show "stay puffy longer".
    # That headline needs a fixed-interior-entropy comparison or the Phase-2
    # time-integrated track (see docs / AMENDMENTS).
    ax.text(0.5, 0.03,
            "diagnostic at fixed $T_{int}$: L∝R² imposed by $T_{int}$,\n"
            "not opacity-throttled — see caveat (not the headline)",
            transform=ax.transAxes, ha="center", va="bottom", fontsize=8,
            color="#777", style="italic")
    ax.legend(fontsize=9, title="irradiation")
    path = save(fig, cfg, "F3_cooling_ratio")
    print(f"wrote {path}  (numerator=fp, baseline={baseline})")
    return path


if __name__ == "__main__":
    main()
