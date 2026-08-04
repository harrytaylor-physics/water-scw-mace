"""F1: T(P) profile and the Schwarzschild ratio at T_int = 20 K (SPEC §8).

    python figures/plot_f1_profiles.py config/fiducial.yaml

Plotted at T_int = 20 K, the cold regime where the composite (FP-informed) and
pure line-list interiors differ. Left: T(P) for composite (black) and HITEMP
(dashed). Right: grad_rad/grad_ad versus pressure; values below 1 are radiative.
Pure FP and POKAZATEL are opacity-extrapolated in the radiative zone and are not
shown. The shaded band marks the deep radiative window HITEMP opens (ratio < 1)
that the composite closes (ratio >= 1).
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # run as script or -m

import matplotlib.pyplot as plt
import numpy as np

from eos import AquaEOS
from figures._common import band_caption, comparison_band, load_cfg, save, set_style
from structure import build_structure

T_INT = 20.0
HITEMP_C = "#D55E00"  # single accent colour; composite is black


def main(argv=None):
    cfg = load_cfg(argv if argv is not None else sys.argv[1:])
    set_style()
    eos = AquaEOS.from_config(cfg)
    t_eq = float(cfg.irradiation.t_eq_k[0])
    band = comparison_band(cfg, ["composite", "hitemp"])

    comp = build_structure(cfg, t_eq, T_INT, "composite", eos=eos, n_sample=1500, band=band)
    ht = build_structure(cfg, t_eq, T_INT, "hitemp", eos=eos, n_sample=1500, band=band)

    fig, (ax_t, ax_g) = plt.subplots(1, 2, figsize=(10.6, 5.6), sharey=True,
                                     constrained_layout=True)

    for res, color, ls, lw, lab in [(comp, "black", "-", 2.4, "composite"),
                                    (ht, HITEMP_C, "--", 1.6, "HITEMP")]:
        p = res.profile
        pbar = p["pressure_pa"] / 1e5
        ax_t.plot(p["temperature_k"], pbar, color=color, ls=ls, lw=lw, label=lab)
        ratio = p["grad_rad"] / p["grad_ad"]
        ax_g.plot(ratio, pbar, color=color, ls=ls, lw=lw, label=lab)

    # Shade the deep radiative window: HITEMP radiative (ratio<1) at FP-grid
    # density (rho>=109 kg/m^3), which the composite keeps convective.
    ph = ht.profile
    rh = ph["grad_rad"] / ph["grad_ad"]
    window = (rh < 1.0) & (ph["rho_si"] >= 109.0)
    if window.any():
        pbar_h = ph["pressure_pa"] / 1e5
        lo, hi = pbar_h[window].min(), pbar_h[window].max()
        for ax in (ax_t, ax_g):
            ax.axhspan(lo, hi, color="grey", alpha=0.25)

    ax_t.set_yscale("log")
    ax_t.invert_yaxis()
    ax_t.set_xlabel("Temperature [K]")
    ax_t.set_ylabel("Pressure [bar]")
    ax_t.set_title("Temperature profile")
    ax_t.legend(fontsize=9, loc="upper right")

    ax_g.axvline(1.0, color="grey", lw=1.0, ls=":")
    ax_g.set_xscale("log")
    ax_g.set_xlabel(r"$\nabla_{rad}\,/\,\nabla_{ad}$  (< 1: radiative)")
    ax_g.set_title("Schwarzschild ratio; grey band = deep radiative window")

    fig.suptitle(f"Interior at $T_{{int}}$={T_INT:.0f} K — $T_{{eq}}$={t_eq:.0f} K, "
                 f"$M_c$={cfg.core.mass_mearth:.0f} $M_⊕$, $f_{{env}}$={cfg.envelope.f_env:.0%}, "
                 f"{band_caption(band)}", fontsize=11)
    path = save(fig, cfg, "F1_tp_profiles")
    print(f"wrote {path}  (composite RCBs: {len(comp.rcb_pressures_pa)}, "
          f"HITEMP RCBs: {len(ht.rcb_pressures_pa)})")
    return path


if __name__ == "__main__":
    main()
