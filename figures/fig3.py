#!/usr/bin/env python3
"""Figure 3 — what the opacity difference does to a radiative-convective boundary.

    python figures/fig3.py config/fiducial.yaml

**a** The Schwarzschild ratio ∇_rad/∇_ad through the envelope of a 20 M⊕ water
world (M_c = 10 M⊕, water mass fraction 0.5, T_eq = 500 K, T_int = 30 K), for the
line list and for the composite opacity. The envelope is radiative where the ratio
falls below one. Both boundaries are marked.

**b** How far each boundary is from flipping, against how much the two opacities
disagree there. ∇_rad ∝ κ_R, so at fixed structure the factor that would flip a
boundary is exactly ∇_ad/∇_rad at that point — the bar. The paired point is the
first-principles / line-list κ_R ratio measured at the same place. Every boundary
shown has first-principles data interpolated at it, not clamped from elsewhere.

The point of the panel is the gap between bar and point: the disagreement exceeds
the flip threshold at every boundary, by one to two orders of magnitude.

Writes ``figs/fig3.pdf`` and merges into ``figs/fig_numbers.json``.
"""

from __future__ import annotations

import csv
import sys
from dataclasses import replace
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
from matplotlib.lines import Line2D

from config import load_config
from eos import AquaEOS
from figures._manuscript import (
    COMPOSITE_C,
    HITEMP_C,
    OKABE_ITO,
    FigureNumbers,
    figure,
    panel_label,
    save_pdf,
    set_manuscript_style,
)
from figures._common import comparison_band
from structure import build_structure

CASE = {"mass_core_mearth": 10.0, "wmf": 0.50, "t_eq_k": 500.0, "t_int_k": 30.0}
# Matches the domain scan and the boundary analysis (both 1500), so panel a's
# boundary densities are the same numbers panel b and the scan report. The located
# boundary moves by ~1.5% in density between 1500 and 2000 samples, which is
# resolution in finding the crossing, not physics — but the figure must agree with
# the table it sits next to.
N_SAMPLE = 1500
# Order matters: the composite uses the line list below the blend window, so the
# two curves coincide there. The line list is drawn last and dashed so both read.
SOURCES = (("composite", COMPOSITE_C, "-", "composite"),
           ("hitemp", HITEMP_C, "--", "line list"))


def _profile(cfg, source, band, eos):
    run_cfg = replace(
        cfg,
        core=replace(cfg.core, mass_mearth=CASE["mass_core_mearth"]),
        envelope=replace(cfg.envelope, f_env=CASE["wmf"]),
    )
    return build_structure(run_cfg, CASE["t_eq_k"], CASE["t_int_k"], source,
                           eos=eos, n_sample=N_SAMPLE, band=band)


def panel_a(ax, cfg, nums: FigureNumbers) -> dict:
    eos = AquaEOS.from_config(cfg)
    band = comparison_band(cfg, ["composite", "hitemp"])
    marked = {}
    for source, colour, style, label in SOURCES:
        res = _profile(cfg, source, band, eos)
        rho = np.asarray(res.profile["rho_si"], float)
        ratio = np.asarray(res.profile["grad_rad"], float) / np.asarray(
            res.profile["grad_ad"], float)
        order = np.argsort(rho)
        ax.plot(rho[order], ratio[order], color=colour, lw=1.1, ls=style,
                label=label, zorder=3 if source == "composite" else 4)
        rcbs = sorted(float(r) for r in res.rcb_densities_si)
        marked[source] = {
            "radius_rearth": float(res.radius_p_rearth),
            "rcb_densities_kg_m3": rcbs,
            "deepest_rcb_kg_m3": rcbs[-1] if rcbs else None,
            "max_radiative_rho_kg_m3": float(res.max_radiative_rho_si),
        }
        for rho_rcb in rcbs:
            ax.axvline(rho_rcb, color=colour, lw=0.7, ls=":", alpha=0.9, zorder=2)
            # at the foot of the axis, clear of the dex arrow at the top
            ax.annotate(f"{rho_rcb:.0f}", (rho_rcb, ax.get_ylim()[0]),
                        textcoords="offset points", xytext=(2, 3), fontsize=7.0,
                        color=colour, zorder=6)

    ax.axhline(1.0, color="black", lw=0.7, ls=":", zorder=1)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Density ρ (kg m$^{-3}$)")
    ax.set_ylabel("∇$_{rad}$ / ∇$_{ad}$")
    ax.legend(loc="upper left", fontsize=7.0)

    ht = marked["hitemp"]["deepest_rcb_kg_m3"]
    comp = marked["composite"]["deepest_rcb_kg_m3"]
    if ht and comp:
        d_dex = np.log10(comp) - np.log10(ht)
        y0, y1 = ax.get_ylim()
        ax.annotate("", xy=(ht, y1 / 40.0), xytext=(comp, y1 / 40.0),
                    arrowprops=dict(arrowstyle="<->", lw=0.7, color="black"), zorder=6)
        ax.annotate(f"{abs(d_dex):.2f} dex", xy=(np.sqrt(ht * comp), y1 / 22.0),
                    ha="center", fontsize=7.0, zorder=6)
        marked["boundary_shift_dex"] = float(d_dex)
    ax.text(0.03, 0.05, "radiative where the ratio < 1", transform=ax.transAxes,
            fontsize=7.0, color="#444444")

    nums.axes("a", ax, "Density (kg/m^3)", "grad_rad / grad_ad", xlog=True, ylog=True)
    nums.panel("a", case={**CASE, "mass_planet_mearth":
                          CASE["mass_core_mearth"] / (1.0 - CASE["wmf"]),
                          "n_profile_samples": N_SAMPLE},
               sampling_note=(
                   "the located boundary density shifts by ~1.5% between 1500 and "
                   "2000 profile samples; 1500 matches the domain scan and panel b"),
               **marked)
    return marked


def panel_b(ax, cfg, nums: FigureNumbers) -> None:  # noqa: C901
    path = cfg.output_path / "S8_rcb_flip.csv"
    if not path.exists():
        raise SystemExit(
            f"{path} not found — run figures/plot_s8_rcb_flip.py first (it needs the "
            f"extended domain scan)."
        )
    with path.open() as fh:
        rows = [r for r in csv.DictReader(fh) if r["fp_clamped_at_rcb"] == "False"]
    rows.sort(key=lambda r: float(r["fp_over_hitemp"]))

    x = np.arange(len(rows))
    ratio = np.array([float(r["fp_over_hitemp"]) for r in rows])
    thresh = np.array([float(r["threshold_factor"]) for r in rows])

    ax.bar(x, thresh - 1.0, width=0.55, bottom=1.0, color=OKABE_ITO["sky"],
           edgecolor="black", linewidth=0.4, zorder=2,
           label="flip threshold ∇$_{ad}$/∇$_{rad}$")
    ax.plot(x, ratio, ls="none", marker="o", ms=4.2, mfc=OKABE_ITO["blue"],
            mec="black", mew=0.4, zorder=4,
            label="first-principles / line-list κ$_R$ at the boundary")
    ax.vlines(x, thresh, ratio, color=OKABE_ITO["grey"], lw=0.8, zorder=3)

    ax.axhline(1.0, color="black", lw=0.7, ls=":", zorder=1)
    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels(
        [f"{float(r['t_k']):.0f}\n{float(r['rho_si']):.0f}" for r in rows],
        fontsize=7.0)
    ax.set_xlabel("Boundary:  T (K)  /  ρ (kg m$^{-3}$)")
    ax.set_ylabel("Factor")
    ax.set_ylim(0.85, ratio.max() * 3.2)
    for xi, rr, th in zip(x, ratio, thresh):
        ax.annotate(f"×{rr / th:.0f}", (xi, rr), textcoords="offset points",
                    xytext=(0, 6), ha="center", fontsize=7.0, zorder=5)
    ax.legend(loc="upper center", fontsize=7.0, ncol=1,
              bbox_to_anchor=(0.5, 1.0))

    # The stored boundary analysis predates the T = 2250/2750 K line-list rows
    # (A60). Those are far from these boundaries (726-1008 K), but a Delaunay
    # triangulation is not local, so the sensitivity is measured rather than
    # assumed: re-evaluate the ratio at each boundary with the current tables.
    from opacity.build_tables import get_source
    from opacity.table import CM2_PER_G_TO_SI

    band = comparison_band(cfg, ["fp", "hitemp"])
    fp_src = get_source("fp", cfg, band=band)
    ll_src = get_source("hitemp", cfg, band=band)
    drift = []
    for r in rows:
        t_k, rho_si = float(r["t_k"]), float(r["rho_si"])
        now = (fp_src.kappa_R(t_k, rho_si) / ll_src.kappa_R(t_k, rho_si))
        drift.append(abs(now / float(r["fp_over_hitemp"]) - 1.0))
    max_drift = max(drift)

    nums.axes("b", ax, "Radiative-convective boundary (ordered by ratio)", "Factor",
              ylog=True)
    nums.panel(
        "b",
        n_boundaries=len(rows),
        kappa_ratio_range=[float(ratio.min()), float(ratio.max())],
        flip_threshold_range=[float(thresh.min()), float(thresh.max())],
        margin_range=[float((ratio / thresh).min()), float((ratio / thresh).max())],
        boundaries=[{"T_K": float(r["t_k"]), "rho_kg_m3": float(r["rho_si"]),
                     "kappa_ratio": float(r["fp_over_hitemp"]),
                     "flip_threshold": float(r["threshold_factor"])} for r in rows],
        selection="boundaries with first-principles data interpolated, not clamped",
        input_vintage=(
            "boundary locations and kappa ratios come from the stored boundary "
            "analysis, which predates the 2250/2750 K line-list rows"),
        max_ratio_drift_since_stored=float(max_drift),
    )
    print(f"  stored-vs-current kappa ratio at these boundaries: max drift "
          f"{100 * max_drift:.2f}%")


def main(argv=None) -> int:
    import argparse

    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("config", type=Path)
    args = ap.parse_args(argv if argv is not None else sys.argv[1:])
    cfg = load_config(args.config)

    set_manuscript_style()
    nums = FigureNumbers(
        figure="fig3",
        caption_hint=(
            "a, Schwarzschild ratio through a 20 M_earth water world for the line "
            "list and the composite, both boundaries marked. b, per boundary, the "
            "first-principles/line-list kappa_R ratio against the factor that would "
            "flip it."),
    )
    fig, (ax_a, ax_b) = figure(82.0, ncols=2)
    marked = panel_a(ax_a, cfg, nums)
    panel_b(ax_b, cfg, nums)
    panel_label(ax_a, "a", dx=-0.15, dy=1.08)
    panel_label(ax_b, "b", dx=-0.15, dy=1.08)

    path = save_pdf(fig, "fig3", nums)
    print(f"wrote {path}")
    for source, _c, _ls, label in SOURCES:
        m = marked[source]
        print(f"  {label:<10s} R_p = {m['radius_rearth']:.5f} R_earth, "
              f"RCB rho = {', '.join(f'{v:.4g}' for v in m['rcb_densities_kg_m3'])} kg/m^3")
    if "boundary_shift_dex" in marked:
        print(f"  boundary shift: {marked['boundary_shift_dex']:+.3f} dex")
    b = nums.panels["b"]
    print(f"  {b['n_boundaries']} boundaries: kappa ratio "
          f"{b['kappa_ratio_range'][0]:.2f}-{b['kappa_ratio_range'][1]:.2f}, "
          f"flip threshold {b['flip_threshold_range'][0]:.3f}-"
          f"{b['flip_threshold_range'][1]:.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
