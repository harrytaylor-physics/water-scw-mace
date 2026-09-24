#!/usr/bin/env python3
"""Figure 4 — where the opacity disagreement lands relative to real envelopes.

    python figures/fig4.py config/fiducial.yaml

The (T, ρ) plane shaded by the first-principles / HITEMP Rosseland-mean ratio over
the 48 first-principles points, interpolated on the model's own triangulation. Over
it, the line-list radiative–convective boundary of each of the 108 scanned planets,
coloured by equilibrium temperature. The outline is the extent of the
first-principles grid — its convex hull, which is exactly where the model
interpolates rather than clamps.

The figure exists to answer one question: do the planets' boundaries sit where the
two opacities actually disagree? A boundary outside the outline is one where the
first-principles opacity is an extrapolation, and the shading there says nothing.

Writes ``figs/fig4.pdf`` and merges into ``figs/fig_numbers.json``.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
from matplotlib.colors import LogNorm
from matplotlib.lines import Line2D

from config import load_config
from figures._manuscript import (
    FigureNumbers,
    figure,
    save_pdf,
    set_manuscript_style,
)
from opacity.crossover import crossover_v3_rows

BAND_CM = (100.0, 4000.0)
NGRID = 340
SCAN_CSV = "S4_domain_scan_extended.csv"


def _hull(points: np.ndarray) -> np.ndarray:
    from scipy.spatial import ConvexHull

    h = ConvexHull(points)
    idx = list(h.vertices) + [h.vertices[0]]
    return points[idx]


def main(argv=None) -> int:
    import argparse

    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("config", type=Path)
    args = ap.parse_args(argv if argv is not None else sys.argv[1:])
    cfg = load_config(args.config)

    set_manuscript_style()
    nums = FigureNumbers(
        figure="fig4",
        caption_hint=(
            "(T, rho) plane shaded by the first-principles/HITEMP Rosseland ratio "
            "over the 48 first-principles points, with the line-list radiative-"
            "convective boundary of 108 scanned planets overlaid, coloured by T_eq. "
            "The outline is the first-principles grid extent."),
    )

    rows = crossover_v3_rows(cfg, BAND_CM)
    t_pts = np.array([r["T"] for r in rows])
    rho_pts = np.array([r["rho"] for r in rows])
    ratio_pts = np.array([r["ratio_rosseland_bare"] for r in rows])

    # Same scheme the model uses for opacity: linear on a Delaunay triangulation of
    # (T, log10 rho), here over log10 of the ratio. Blank outside the hull.
    from scipy.interpolate import LinearNDInterpolator

    pts = np.column_stack([t_pts, np.log10(rho_pts)])
    lin = LinearNDInterpolator(pts, np.log10(ratio_pts), rescale=True)

    t_lin = np.linspace(t_pts.min(), t_pts.max(), NGRID)
    lr_lin = np.linspace(np.log10(rho_pts.min()), np.log10(rho_pts.max()), NGRID)
    TT, LR = np.meshgrid(t_lin, lr_lin)
    surface = 10.0 ** lin(np.column_stack([TT.ravel(), LR.ravel()])).reshape(TT.shape)

    fig, ax = figure(108.0)
    finite = np.isfinite(surface)
    mesh = ax.pcolormesh(
        10.0**LR, TT, np.ma.masked_invalid(surface), cmap="viridis",
        norm=LogNorm(vmin=np.nanmin(surface[finite]), vmax=np.nanmax(surface[finite])),
        shading="auto", zorder=1, rasterized=True)
    cb = fig.colorbar(mesh, ax=ax, pad=0.015, fraction=0.045)
    cb.set_label("First-principles / HITEMP κ$_R$", fontsize=8.0)
    cb.ax.tick_params(labelsize=7.0)

    hull = _hull(np.column_stack([np.log10(rho_pts), t_pts]))
    ax.plot(10.0 ** hull[:, 0], hull[:, 1], color="black", lw=0.9, ls="-",
            zorder=4, label="first-principles grid extent")
    ax.scatter(rho_pts, t_pts, s=7.0, facecolors="none", edgecolors="white",
               linewidths=0.5, zorder=5)

    # the scanned planets' line-list boundaries
    scan = cfg.output_path / SCAN_CSV
    if not scan.exists():
        raise SystemExit(
            f"{scan} not found — run figures/plot_s4_domain_scan.py --extended first."
        )
    with scan.open() as fh:
        planets = [r for r in csv.DictReader(fh) if r.get("skipped") != "True"]
    rcb_rho, rcb_t, t_eq = [], [], []
    for r in planets:
        rho = r.get("hitemp_max_radiative_rho")
        t_k = r.get("hitemp_max_radiative_t_k")
        if not rho or not t_k:
            continue
        rho, t_k = float(rho), float(t_k)
        if rho <= 0 or t_k <= 0:
            continue
        rcb_rho.append(rho)
        rcb_t.append(t_k)
        t_eq.append(float(r["t_eq_k"]))
    rcb_rho, rcb_t, t_eq = np.array(rcb_rho), np.array(rcb_t), np.array(t_eq)

    teq_levels = sorted(set(t_eq.tolist()))
    markers = ["o", "s", "^", "D", "v", "P", "X", "*"]
    from matplotlib import colormaps

    cmap_teq = colormaps["plasma"]
    for i, lvl in enumerate(teq_levels):
        m = t_eq == lvl
        ax.scatter(rcb_rho[m], rcb_t[m], marker=markers[i % len(markers)], s=17.0,
                   facecolor=cmap_teq(i / max(1, len(teq_levels) - 1)),
                   edgecolor="black", linewidth=0.35, zorder=6,
                   label=f"T$_{{eq}}$ = {lvl:.0f} K")

    ax.set_xscale("log")
    ax.set_xlabel("Density ρ (kg m$^{-3}$)")
    ax.set_ylabel("Temperature T (K)")
    handles, labels = ax.get_legend_handles_labels()
    # Above the axes: every corner of this panel carries data.
    ax.legend(handles, labels, loc="lower center", bbox_to_anchor=(0.5, 1.005),
              fontsize=7.0, ncol=4, borderaxespad=0.0, columnspacing=1.2)

    nums.axes("", ax, "Density (kg/m^3)", "Temperature (K)", xlog=True)
    inside = np.isfinite(lin(np.column_stack([rcb_t, np.log10(rcb_rho)])))
    nums.panel(
        "",
        n_first_principles_points=int(t_pts.size),
        n_planets_plotted=int(rcb_rho.size),
        n_planets_in_scan=len(planets),
        ratio_range_over_grid=[float(ratio_pts.min()), float(ratio_pts.max())],
        colourbar_range=[float(np.nanmin(surface[finite])),
                         float(np.nanmax(surface[finite]))],
        boundary_T_range_K=[float(rcb_t.min()), float(rcb_t.max())],
        boundary_rho_range_kg_m3=[float(rcb_rho.min()), float(rcb_rho.max())],
        t_eq_levels_K=[float(v) for v in teq_levels],
        n_boundaries_inside_grid_extent=int(inside.sum()),
        grid_extent="convex hull of the 48 points in (T, log10 rho)",
        band_cm=list(BAND_CM),
    )

    path = save_pdf(fig, "fig4", nums)
    print(f"wrote {path}")
    p = nums.panels[""]
    print(f"  {p['n_first_principles_points']} first-principles points; ratio "
          f"{p['ratio_range_over_grid'][0]:.2f}-{p['ratio_range_over_grid'][1]:.2f}")
    print(f"  {p['n_planets_plotted']} planetary boundaries plotted "
          f"({p['n_boundaries_inside_grid_extent']} inside the grid extent), "
          f"T {p['boundary_T_range_K'][0]:.0f}-{p['boundary_T_range_K'][1]:.0f} K, "
          f"rho {p['boundary_rho_range_kg_m3'][0]:.3g}-"
          f"{p['boundary_rho_range_kg_m3'][1]:.3g} kg/m^3")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
