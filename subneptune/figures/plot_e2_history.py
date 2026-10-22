"""E2 (Phase 2): cooling history under composite vs line-list opacity (A29).

    python figures/plot_e2_history.py config/fiducial.yaml

Evolves matched initial conditions under (a) composite (FP-informed) and
(b) pure line-list opacity and shows how the histories diverge: R_p(age), the
deepest RCB density(age) vs the FP floor, and the radius/entropy differences.
The atmospheric-escape window (first 100 Myr) and the FP-relevant epoch
(RCB ρ ≥ 109) are marked — a radius difference in the escape window would be an
irreversible imprint even though v1 does not compute escape.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # run as script or -m

import matplotlib.pyplot as plt
import numpy as np

from eos import AquaEOS
from evolution import compare_histories
from figures._common import available_sources, band_caption, load_cfg, save, set_style
from structure.model import FP_GRID_FLOOR_SI

ESCAPE_WINDOW_GYR = 0.1  # first 100 Myr


def main(argv=None):
    cfg = load_cfg(argv if argv is not None else sys.argv[1:])
    set_style()
    eos = AquaEOS.from_config(cfg)
    t_eq = float(cfg.irradiation.t_eq_k[0])

    line_list = "hitemp" if "hitemp" in available_sources(cfg, ["hitemp"]) else "pokazatel"
    # Reach cold T_int so the late, deep-RCB (FP-relevant) epoch is captured.
    t_int_grid = np.geomspace(12.0, 450.0, 44)
    h = compare_histories(cfg, "composite", line_list, t_eq, eos=eos, t_int_grid_k=t_int_grid)

    fig, (ax_r, ax_rho, ax_d) = plt.subplots(3, 1, figsize=(7.2, 9.2), sharex=True,
                                             constrained_layout=True)
    ep = h.fp_relevant_epoch_gyr()
    floor_age = h.radiogenic_floor_age_gyr(cfg.core.mass_mearth)  # A30

    def annotate_epochs(ax, label=True):
        ax.axvspan(h.age_gyr.min(), ESCAPE_WINDOW_GYR, color="orange", alpha=0.12,
                   label="first 100 Myr (escape window)" if label else None)
        if ep:
            ax.axvspan(ep[0], ep[1], color="green", alpha=0.10,
                       label="RCB density ≥ 109 kg/m³" if label else None)
        if floor_age is not None:
            ax.axvspan(floor_age, h.age_gyr.max(), color="grey", alpha=0.22,
                       label="radiogenically supported" if label else None)

    ax_r.plot(h.age_gyr, h.radius_a, color="#CC79A7", label="composite (FP-informed)")
    ax_r.plot(h.age_gyr, h.radius_b, color="#D55E00", ls="--", label=line_list.upper())
    annotate_epochs(ax_r)
    ax_r.set_ylabel("$R_p$ [$R_⊕$]")
    ax_r.set_title(f"Cooling history: composite vs {line_list.upper()}, "
                   f"$M_c$={cfg.core.mass_mearth:.0f} $M_⊕$, $f_{{env}}$={cfg.envelope.f_env:.0%}, "
                   f"{band_caption(cfg.opacity.band_cm or (100, 4000))}", fontsize=9.5)
    ax_r.legend(fontsize=8.5)

    ax_rho.plot(h.age_gyr, h.rcb_rho_a, color="#CC79A7", label="composite")
    ax_rho.plot(h.age_gyr, h.rcb_rho_b, color="#D55E00", ls="--", label=line_list.upper())
    ax_rho.axhline(FP_GRID_FLOOR_SI, color="k", lw=1, ls=":", label="FP floor 109 kg/m³")
    annotate_epochs(ax_rho, label=False)
    ax_rho.set_yscale("log")
    ax_rho.set_ylabel("deepest RCB ρ [kg/m³]")
    ax_rho.legend(fontsize=8.5)

    ax_d.plot(h.age_gyr, h.radius_frac_diff * 100, color="#0072B2", label="ΔR_p/R_p")
    ax_d.plot(h.age_gyr, h.entropy_frac_diff * 100, color="#009E73", ls="--", label="Δs/s (interior)")
    ax_d.axhline(0, color="k", lw=0.6)
    annotate_epochs(ax_d, label=False)
    ax_d.set_ylabel("composite − line-list [%]")
    ax_d.set_xlabel("Age [Gyr]")
    ax_d.set_xscale("log")
    ax_d.legend(fontsize=8.5)

    fig.text(0.5, -0.005,
             "Radiogenically supported epoch: internal temperature held at the radiogenic "
             "floor; core luminosity not evolved in this model.",
             ha="center", va="top", fontsize=8.5)
    path = save(fig, cfg, "E2_history")
    print(f"wrote {path}")
    print(f"peak |ΔR_p/R_p| = {h.peak_radius_frac_diff*100:.2f}%; "
          f"ΔR_p at 100 Myr = {h.diff_at_age(ESCAPE_WINDOW_GYR)*100:+.3f}%; "
          f"FP-relevant epoch = {ep} Gyr")

    # -- escape-window check and the FP-floor crossing ages (A40) --
    esc = h.escape_window_extrema(ESCAPE_WINDOW_GYR)
    cross = h.first_age_rcb_above(FP_GRID_FLOOR_SI)
    print(f"\nE2 — escape window (first {ESCAPE_WINDOW_GYR*1000:.0f} Myr), "
          f"composite − {line_list.upper()}, {esc['n_grid_points']} age-grid points "
          f"from {esc['age_min_gyr']:.2e} Gyr:")
    # These are null-level numbers, so print them with significant digits rather
    # than fixed decimals (0.0001% hides whether it is 1e-4 or 1e-7).
    print(f"  max |ΔR_p/R_p| = {esc['max_abs_radius_frac']*100:.3g}% "
          f"(at {esc['age_at_max_radius_gyr']:.4g} Gyr)")
    print(f"  max |Δs/s|     = {esc['max_abs_entropy_frac']*100:.3g}% "
          f"(at {esc['age_at_max_entropy_gyr']:.4g} Gyr)")
    print(f"  deepest RCB first reaches ρ ≥ {FP_GRID_FLOOR_SI:g} kg/m³: "
          + ", ".join(
              f"{name} at {('%.3f Gyr' % age) if age is not None else 'never (within the track)'}"
              for name, age in ((h.source_a, cross[h.source_a]),
                                (h.source_b, cross[h.source_b]))))
    return path


if __name__ == "__main__":
    main()
