"""F3: Kelvin–Helmholtz cooling-time systematic vs (f_env, T_eq) (SPEC §8, revised).

    python figures/plot_f3_cooling_ratio.py config/fiducial.yaml

Reframed after the FP-validity finding (A27–A31). The reliable cooling-time
comparison is between the *line lists*: at fixed interior entropy, τ_source/τ_HITEMP
for the sources that actually cover the radiative zone (composite ≈ HITEMP; the
FP contribution does not reach the opacity-relevant zone for these planets, so
pure FP is shown flagged as an extrapolation artifact, not a headline). The
dominant, reliable opacity systematic on cooling is the **line-list choice**
(POKAZATEL vs HITEMP), ~2×.
"""

from __future__ import annotations

import sys
from dataclasses import replace
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # run as script or -m

import matplotlib.pyplot as plt

from eos import AquaEOS
from figures._common import (
    SOURCE_COLOR,
    SOURCE_LABEL,
    available_sources,
    band_caption,
    comparison_band,
    load_cfg,
    save,
    set_style,
)
from structure import build_at_entropy, build_structure

F_ENV_GRID = (0.01, 0.02, 0.03, 0.05, 0.10)


def main(argv=None):
    cfg = load_cfg(argv if argv is not None else sys.argv[1:])
    set_style()
    eos = AquaEOS.from_config(cfg)
    t_int_ref = float(cfg.interior.t_int_k[len(cfg.interior.t_int_k) // 2])
    t_eq = float(cfg.irradiation.t_eq_k[len(cfg.irradiation.t_eq_k) // 2])

    baseline = "hitemp"
    others = [s for s in available_sources(cfg, ["composite", "pokazatel", "fp"])]
    plotted = [baseline] + others
    band = comparison_band(cfg, plotted)

    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    for name in others:
        reliable_all = True
        ratio = []
        for f_env in F_ENV_GRID:
            cfg_f = replace(cfg, envelope=replace(cfg.envelope, f_env=f_env))
            base = build_structure(cfg_f, t_eq, t_int_ref, baseline, eos=eos, n_sample=400, band=band)
            src = build_at_entropy(cfg_f, base.interior_entropy_j_kg_k, t_eq, name,
                                   eos=eos, band=band, n_sample=400)
            ratio.append(src.tau_kh_s / base.tau_kh_s)
            reliable_all &= src.clamp_reliable
        ax.plot(F_ENV_GRID, ratio, marker="o", ms=5, color=SOURCE_COLOR[name],
                ls="-" if reliable_all else "--",
                label=f"{SOURCE_LABEL[name]} / HITEMP" + ("" if reliable_all else "  [UNRELIABLE]"))

    ax.axhline(1.0, color="k", lw=0.8, ls=":")
    ax.set_xscale("log")
    ax.set_yscale("log")  # FP artifact (~25×) and the reliable ~1-2× both readable
    ax.set_xlabel("$f_{env}$")
    ax.set_ylabel(r"$\tau_{KH}$(source) / $\tau_{KH}$(HITEMP)  at fixed interior entropy")
    ax.set_title("Cooling-time systematic vs HITEMP — the reliable spread is the "
                 "line-list choice\n"
                 f"composite ≈ HITEMP (FP doesn't reach the radiative zone); "
                 f"$T_{{eq}}$={t_eq:.0f} K, $T_{{int,ref}}$={t_int_ref:.0f} K, {band_caption(band)}",
                 fontsize=9)
    ax.legend(fontsize=8.5)
    path = save(fig, cfg, "F3_cooling_ratio")
    print(f"wrote {path}  (baseline=HITEMP; reliable sources: composite/pokazatel)")
    return path


if __name__ == "__main__":
    main()
