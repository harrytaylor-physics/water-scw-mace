"""F2: RCB pressure and R_p vs T_int, per opacity source (SPEC §8).

    python figures/plot_f2_rcb_radius.py config/fiducial.yaml
"""

from __future__ import annotations

import sys
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
    register_constant_reference,
    run_t_int_grid,
    save,
    set_style,
)


def main(argv=None):
    cfg = load_cfg(argv if argv is not None else sys.argv[1:])
    set_style()
    eos = AquaEOS.from_config(cfg)
    t_eq = float(cfg.irradiation.t_eq_k[0])

    register_constant_reference(cfg)
    # Reliable sources only: composite (which equals HITEMP at these T_int) and a
    # constant-opacity reference. Pure FP and POKAZATEL are opacity-extrapolated
    # in the radiative zone and are omitted.
    band = comparison_band(cfg, ["composite", "hitemp"])
    plotted = [("composite", "composite (= HITEMP at these $T_{int}$)", "black", "-"),
               ("constant", "constant opacity (reference)", "grey", "--")]

    fig, (ax_r, ax_p) = plt.subplots(1, 2, figsize=(10.5, 4.6))
    for name, lbl, color, ls in plotted:
        models = run_t_int_grid(cfg, name, t_eq, eos, band=band)
        t_int = [m.t_int_k for m in models]
        r_p = [m.radius_p_rearth for m in models]
        rcb = [max(m.rcb_pressures_pa) / 1e5 if m.rcb_pressures_pa else float("nan") for m in models]
        ax_r.plot(t_int, r_p, label=lbl, color=color, ls=ls, marker="o", ms=5)
        ax_p.plot(t_int, rcb, label=lbl, color=color, ls=ls, marker="o", ms=5)

    ax_r.set_xlabel("$T_{int}$ [K]")
    ax_r.set_ylabel("$R_p$ [$R_⊕$]")
    ax_r.set_title("Planet radius")
    ax_r.legend(fontsize=8.5)

    ax_p.set_xlabel("$T_{int}$ [K]")
    ax_p.set_ylabel("deepest RCB pressure [bar]")
    ax_p.set_yscale("log")
    ax_p.set_title("Radiative–convective boundary")

    fig.suptitle(f"$T_{{eq}}$={t_eq:.0f} K, $M_c$={cfg.core.mass_mearth:.0f} $M_⊕$, "
                 f"$f_{{env}}$={cfg.envelope.f_env:.0%}, {band_caption(band)}", y=1.02)
    path = save(fig, cfg, "F2_rcb_radius")
    print(f"wrote {path}")
    return path


if __name__ == "__main__":
    main()
