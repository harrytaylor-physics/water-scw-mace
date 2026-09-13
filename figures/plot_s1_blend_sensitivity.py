"""S1: blend-window sensitivity of the RCB structure (A35).

    python figures/plot_s1_blend_sensitivity.py config/fiducial.yaml

The composite source blends HITEMP (dilute) into FP (dense) across a density
window. At the fiducial planet (T_int = 20 K, T_eq = 500 K) the composite keeps a
single RCB whose density lies *inside* the default window ρ = 100-300 kg/m³, so
its location is set partly by the blend prescription rather than by opacity data
alone. This scans the prescription:

  [100, 200]  narrower window        (linear blend)
  [100, 300]  default                (linear blend)
  [100, 500]  wider window           (linear blend)
  sharp @ 108.796 kg/m³              (hard switch at the FP grid floor)

and compares each with pure HITEMP. Reported per case: the number of RCBs, the
density of each RCB, the densest radiative point, and the deep radiative zones
(ρ > 300 kg/m³) that HITEMP opens and the composite is claimed to close (A28).

Outputs ``outputs/S1_blend_sensitivity.csv`` and ``outputs/S1_blend_sensitivity.png``.
"""

from __future__ import annotations

import csv
import sys
from dataclasses import replace
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # run as script or -m

import matplotlib.pyplot as plt
import numpy as np

from eos import AquaEOS
from figures._common import band_caption, comparison_band, load_cfg, save, set_style
from structure import FP_GRID_FLOOR_SI, build_structure, deep_radiative_zones

T_INT = 20.0  # cold fiducial: the regime where the composite and HITEMP differ
N_SAMPLE = 1500  # matches F1; the deep HITEMP radiative sliver is thin (A35)
DEEP_RHO_SI = 300.0  # closure threshold: radiative layers denser than this
# The prescription only acts near the seam; the dilute envelope (rho < 1) is
# identical in every case, so the panels zoom on the deep envelope.
RHO_XLIM = (1.0, 3.0e3)
HITEMP_C = "#D55E00"

# (label, blend_mode, (lo, hi)) — the default case is marked in the label.
CASES = [
    ("[100, 200] linear", "linear", (100.0, 200.0)),
    ("[100, 300] linear (default)", "linear", (100.0, 300.0)),
    ("[100, 500] linear", "linear", (100.0, 500.0)),
    (f"sharp @ {FP_GRID_FLOOR_SI:g} (FP floor)", "sharp", (FP_GRID_FLOOR_SI, FP_GRID_FLOOR_SI)),
]


def _row(label, mode, window, res):
    """One table row of RCB diagnostics for a converged model."""
    zones = deep_radiative_zones(res, DEEP_RHO_SI)
    return {
        "case": label,
        "blend_mode": mode,
        "blend_lo_kg_m3": "" if window is None else f"{window[0]:g}",
        "blend_hi_kg_m3": "" if window is None else f"{window[1]:g}",
        "n_rcb": len(res.rcb_densities_si),
        "rcb_rho_kg_m3": ";".join(f"{r:.1f}" for r in res.rcb_densities_si),
        "rcb_pressure_bar": ";".join(f"{p / 1e5:.3g}" for p in res.rcb_pressures_pa),
        "max_radiative_rho_kg_m3": f"{res.max_radiative_rho_si:.1f}",
        "n_deep_radiative_zones": len(zones),
        "deep_radiative_zones_kg_m3": ";".join(f"{a:.0f}-{b:.0f}" for a, b in zones),
        "radius_p_rearth": f"{res.radius_p_rearth:.5f}",
        "clamp_fraction_radiative": f"{res.clamp_fraction_radiative:.4f}",
        "clamp_reliable": res.clamp_reliable,
        "n_sample": N_SAMPLE,
    }


def run_cases(cfg, eos=None, n_sample: int = N_SAMPLE):
    """Build HITEMP + every composite blend case. Returns (rows, results)."""
    eos = eos or AquaEOS.from_config(cfg)
    t_eq = float(cfg.irradiation.t_eq_k[0])
    band = comparison_band(cfg, ["composite", "hitemp"])

    rows, results = [], {}
    ht = build_structure(cfg, t_eq, T_INT, "hitemp", eos=eos, n_sample=n_sample, band=band)
    results["HITEMP"] = ht
    rows.append(_row("HITEMP (line list only)", "n/a", None, ht))

    for label, mode, window in CASES:
        c = replace(cfg, opacity=replace(cfg.opacity, blend_window_kg_m3=window,
                                         blend_mode=mode))
        res = build_structure(c, t_eq, T_INT, "composite", eos=eos, n_sample=n_sample,
                              band=band)
        results[label] = res
        rows.append(_row(label, mode, window, res))
    return rows, results


def main(argv=None):
    cfg = load_cfg(argv if argv is not None else sys.argv[1:])
    set_style()
    eos = AquaEOS.from_config(cfg)
    t_eq = float(cfg.irradiation.t_eq_k[0])
    band = comparison_band(cfg, ["composite", "hitemp"])

    rows, results = run_cases(cfg, eos=eos)
    ht = results["HITEMP"]

    # -- figure: one panel per composite case, HITEMP overlaid in each --
    fig, axes = plt.subplots(2, 2, figsize=(11.0, 8.4), sharex=True, sharey=True,
                             constrained_layout=True)
    ph = ht.profile
    for ax, (label, mode, window) in zip(axes.ravel(), CASES):
        res = results[label]
        p = res.profile
        ax.plot(ph["rho_si"], ph["grad_rad"] / ph["grad_ad"], color=HITEMP_C, ls="--",
                lw=1.5, label="HITEMP")
        ax.plot(p["rho_si"], p["grad_rad"] / p["grad_ad"], color="black", lw=2.2,
                label="composite")
        lo, hi = window
        if mode == "sharp":
            ax.axvline(lo, color="#0072B2", lw=1.2, ls="-.", label="sharp switch")
        else:
            ax.axvspan(lo, hi, color="#0072B2", alpha=0.12, label="blend window")
        ax.axvline(DEEP_RHO_SI, color="grey", lw=1.0, ls=":")
        ax.axhline(1.0, color="grey", lw=1.0, ls=":")
        for rho_rcb in res.rcb_densities_si:
            ax.plot([rho_rcb], [1.0], marker="o", ms=6, color="black", zorder=5)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(RHO_XLIM)
        ax.set_title(f"{label}  —  {len(res.rcb_densities_si)} RCB, "
                     f"max radiative ρ = {res.max_radiative_rho_si:.0f}", fontsize=9.5)
        ax.legend(fontsize=8, loc="lower right")

    for ax in axes[-1]:
        ax.set_xlabel(r"$\rho$ [kg m$^{-3}$]")
    for ax in axes[:, 0]:
        ax.set_ylabel(r"$\nabla_{rad}\,/\,\nabla_{ad}$  (< 1: radiative)")

    fig.suptitle(
        f"S1 — blend-window sensitivity of the RCB structure: $T_{{int}}$={T_INT:.0f} K, "
        f"$T_{{eq}}$={t_eq:.0f} K, $M_c$={cfg.core.mass_mearth:.0f} $M_⊕$, "
        f"$f_{{env}}$={cfg.envelope.f_env:.0%}, {band_caption(band)}\n"
        f"grey dotted vertical line: ρ = {DEEP_RHO_SI:.0f} kg m⁻³ closure threshold; "
        f"panels zoom on ρ ≥ {RHO_XLIM[0]:g} kg m⁻³ (the dilute envelope is identical "
        f"in every case)",
        fontsize=10)
    path = save(fig, cfg, "S1_blend_sensitivity")

    csv_path = cfg.output_path / "S1_blend_sensitivity.csv"
    with csv_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    print(f"wrote {path}")
    print(f"wrote {csv_path}")
    print(f"\nS1 — blend-window sensitivity (T_int={T_INT:.0f} K, T_eq={t_eq:.0f} K, "
          f"n_sample={N_SAMPLE}, {band_caption(band)}):")
    head = f"{'case':30s} {'#RCB':>4s} {'RCB rho [kg/m3]':>26s} {'max rad rho':>11s} " \
           f"{'deep zones (rho>300)':>22s} {'R_p [Re]':>9s}"
    print(head)
    print("-" * len(head))
    for r in rows:
        print(f"{r['case']:30s} {r['n_rcb']:>4d} {r['rcb_rho_kg_m3'] or '-':>26s} "
              f"{r['max_radiative_rho_kg_m3']:>11s} "
              f"{(r['deep_radiative_zones_kg_m3'] or 'none') :>22s} "
              f"{r['radius_p_rearth']:>9s}")
    n_deep_ht = rows[0]["n_deep_radiative_zones"]
    worst = max(int(r["n_deep_radiative_zones"]) for r in rows[1:])
    print(f"\nHITEMP deep radiative zones (rho > {DEEP_RHO_SI:.0f}): {n_deep_ht}; "
          f"max over composite cases: {worst} "
          f"({'closed in every composite case' if worst == 0 else 'NOT closed'})")
    return path


if __name__ == "__main__":
    main()
