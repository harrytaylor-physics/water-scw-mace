#!/usr/bin/env python3
"""Figure 2 — the first-principles opacity surface, and how it differs from line lists.

    python figures/fig2.py config/fiducial.yaml

**a** The Rosseland mean κ_R(T, ρ) over the 48 first-principles points, band-limited
to 100-4,000 cm⁻¹ and interpolated exactly as the structure model reads it: linear
on a Delaunay triangulation of (T, log₁₀ρ) over log₁₀κ. Outside the convex hull the
model clamps to the nearest point; that is an extrapolation, not a surface, so it
is left blank here rather than shaded. The 48 computed points are drawn on top.

**b** First-principles / line-list ratio against density, coloured by temperature.
Three references: HITEMP at all 48 points, ExoMol POKAZATEL at the 27 points it
covers, and HITEMP plus the MT_CKD 4.3 water self-continuum — what a planetary
model would actually use. Filled marks are Planck means, open marks Rosseland.

Colour: marker shapes and categorical elements use Okabe-Ito; the two continuous
scales (κ_R in **a**, temperature in **b**) use viridis, which is perceptually
uniform and colourblind-safe — Okabe-Ito is a categorical palette and cannot
carry a continuous variable.

Writes ``figs/fig2.pdf`` and merges into ``figs/fig_numbers.json``.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
from matplotlib.lines import Line2D

from config import load_config
from figures._manuscript import (
    FigureNumbers,
    figure,
    panel_label,
    save_pdf,
    set_manuscript_style,
)
from opacity.build_tables import compute_point_means, get_source
from opacity.crossover import crossover_v3_rows
from opacity.table import CM2_PER_G_TO_SI

BAND_CM = (100.0, 4000.0)
T_VMIN, T_VMAX = 650.0, 3000.0
NGRID = 320

# shape carries the reference, fill carries the mean
REFERENCES = (
    ("hitemp", "o", 22.0, "HITEMP"),
    ("pokazatel", "s", 19.0, "POKAZATEL"),
    ("hitemp_mtckd", "^", 26.0, "HITEMP + MT_CKD continuum"),
)


def panel_a(ax, fig, cfg, nums: FigureNumbers) -> None:
    t_pts, rho_pts, kr, _ = compute_point_means("fp", cfg, BAND_CM)
    src = get_source("fp", cfg, band=BAND_CM)

    t_lin = np.linspace(t_pts.min(), t_pts.max(), NGRID)
    lr_lin = np.linspace(np.log10(rho_pts.min()), np.log10(rho_pts.max()), NGRID)
    TT, LR = np.meshgrid(t_lin, lr_lin)
    # The model's own triangulation; NaN outside the hull is left blank, because
    # the model clamps there and a clamp is not a surface.
    vals = src.table._lin(np.column_stack([TT.ravel(), LR.ravel()]))[:, 0]
    surface = (10.0 ** vals.reshape(TT.shape)) / CM2_PER_G_TO_SI  # SI -> cm^2/g

    from matplotlib.colors import LogNorm

    finite = np.isfinite(surface)
    norm = LogNorm(vmin=np.nanmin(surface[finite]), vmax=np.nanmax(surface[finite]))
    mesh = ax.pcolormesh(10.0**LR, TT, np.ma.masked_invalid(surface),
                         cmap="viridis", norm=norm, shading="auto", rasterized=True)
    ax.scatter(rho_pts, t_pts, s=9.0, facecolors="none", edgecolors="white",
               linewidths=0.6, zorder=3)
    ax.scatter(rho_pts, t_pts, s=9.0, facecolors="none", edgecolors="black",
               linewidths=0.35, zorder=4)
    ax.set_xscale("log")
    ax.set_xlabel("Density ρ (kg m$^{-3}$)")
    ax.set_ylabel("Temperature T (K)")
    cb = fig.colorbar(mesh, ax=ax, pad=0.012, fraction=0.05, aspect=28)
    cb.set_label("κ$_R$ (cm$^{2}$ g$^{-1}$)", fontsize=7.6)
    cb.ax.tick_params(labelsize=7.0)

    nums.axes("a", ax, "Density (kg/m^3)", "Temperature (K)", xlog=True)
    nums.panel(
        "a",
        n_computed_points=int(t_pts.size),
        T_range_K=[float(t_pts.min()), float(t_pts.max())],
        rho_range_kg_m3=[float(rho_pts.min()), float(rho_pts.max())],
        kappa_R_range_cm2_g=[float(kr.min()), float(kr.max())],
        colour_scale="viridis, log",
        colourbar_range_cm2_g=[float(np.nanmin(surface[finite])),
                               float(np.nanmax(surface[finite]))],
        band_cm=list(BAND_CM),
        blank_region="outside the convex hull, where the model clamps",
    )


def _pokazatel_rows(cfg, rows) -> list[dict]:
    """FP/POKAZATEL at the points POKAZATEL actually covers, matched exactly.

    POKAZATEL sits on the same 27 (T, ρ) states as the denser part of the
    first-principles grid, so no interpolation is needed or wanted: the ratio is
    computed between two spectra at one state.
    """
    from opacity.registry import get_spec

    spec = get_spec("pokazatel")
    t_pk, rho_pk, kr_pk, kp_pk = compute_point_means("pokazatel", cfg, BAND_CM)
    by_state = {(round(t, 3), round(r, 3)): (kr, kp)
                for t, r, kr, kp in zip(t_pk, rho_pk, kr_pk, kp_pk)}
    out = []
    for r in rows:
        key = (round(r["T"], 3), round(r["rho"], 3))
        if key not in by_state:
            continue
        kr, kp = by_state[key]
        out.append({"T": r["T"], "rho": r["rho"],
                    "planck": r["kappa_P_fp_cm2_g"] / kp,
                    "rosseland": r["kappa_R_fp_cm2_g"] / kr})
    return out


def panel_b(ax, fig, cfg, rows, nums: FigureNumbers) -> None:
    series = {
        "hitemp": [{"T": r["T"], "rho": r["rho"],
                    "planck": r["ratio_planck_bare"],
                    "rosseland": r["ratio_rosseland_bare"]} for r in rows],
        "hitemp_mtckd": [{"T": r["T"], "rho": r["rho"],
                          "planck": r["ratio_planck_mtckd"],
                          "rosseland": r["ratio_rosseland_mtckd"]} for r in rows],
        "pokazatel": _pokazatel_rows(cfg, rows),
    }

    sc = None
    for name, marker, size, _label in REFERENCES:
        data = series[name]
        if not data:
            continue
        x = [d["rho"] for d in data]
        c = [d["T"] for d in data]
        sc = ax.scatter(x, [d["planck"] for d in data], c=c, cmap="viridis",
                        marker=marker, s=size, vmin=T_VMIN, vmax=T_VMAX,
                        edgecolor="black", linewidth=0.25, zorder=3)
        ax.scatter(x, [d["rosseland"] for d in data], facecolors="none",
                   edgecolors=np.array([sc.cmap(sc.norm(v)) for v in c]),
                   marker=marker, s=size * 1.6, linewidth=0.9, zorder=2)

    ax.axhline(1.0, color="black", lw=0.6, ls=":", zorder=1)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Density ρ (kg m$^{-3}$)")
    ax.set_ylabel("First-principles / line-list mean opacity")
    cb = fig.colorbar(sc, ax=ax, pad=0.012, fraction=0.05, aspect=28)
    cb.set_label("Temperature T (K)", fontsize=7.6)
    cb.ax.tick_params(labelsize=7.0)

    handles = [Line2D([], [], ls="none", marker=m, mfc="#7f7f7f", mec="black",
                      ms=np.sqrt(s), label=lab)
               for _n, m, s, lab in REFERENCES]
    handles += [
        Line2D([], [], ls="none", marker="o", mfc="#7f7f7f", mec="black", ms=3.6,
               label="filled: Planck mean"),
        Line2D([], [], ls="none", marker="o", mfc="none", mec="#7f7f7f", ms=4.0,
               mew=0.9, label="open: Rosseland mean"),
    ]
    # below the axes: the upper left carries the dilute points
    ax.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 1.005),
              fontsize=7.0, ncol=2, borderaxespad=0.0, columnspacing=1.0)

    nums.axes("b", ax, "Density (kg/m^3)", "First-principles / line-list ratio",
              xlog=True, ylog=True)
    stats = {}
    for name, _m, _s, label in REFERENCES:
        data = series[name]
        if not data:
            continue
        stats[label] = {
            "n_points": len(data),
            "planck_range": [min(d["planck"] for d in data),
                             max(d["planck"] for d in data)],
            "rosseland_range": [min(d["rosseland"] for d in data),
                                max(d["rosseland"] for d in data)],
        }
    nums.panel("b", references=stats, unity_line=True,
               temperature_colour_range_K=[T_VMIN, T_VMAX],
               colour_scale="viridis, linear in T")


def main(argv=None) -> int:
    import argparse

    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("config", type=Path)
    args = ap.parse_args(argv if argv is not None else sys.argv[1:])
    cfg = load_config(args.config)

    set_manuscript_style()
    nums = FigureNumbers(
        figure="fig2",
        caption_hint=(
            "a, Rosseland-mean opacity over the 48 first-principles points, "
            "interpolated as the structure model reads it. b, first-principles / "
            "line-list ratio against density, coloured by temperature; filled "
            "Planck, open Rosseland."),
    )
    rows = crossover_v3_rows(cfg, BAND_CM)
    fig, (ax_a, ax_b) = figure(92.0, ncols=2, layout="constrained")
    panel_a(ax_a, fig, cfg, nums)
    panel_b(ax_b, fig, cfg, rows, nums)
    panel_label(ax_a, "a", dx=-0.17, dy=1.06)
    panel_label(ax_b, "b", dx=-0.17, dy=1.30)

    path = save_pdf(fig, "fig2", nums)
    print(f"wrote {path}")
    for label, s in nums.panels["b"]["references"].items():
        print(f"  {label:<32s} n={s['n_points']:>2d}  Planck "
              f"{s['planck_range'][0]:.3f}-{s['planck_range'][1]:.3f}  Rosseland "
              f"{s['rosseland_range'][0]:.3f}-{s['rosseland_range'][1]:.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
