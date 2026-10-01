#!/usr/bin/env python3
"""Figure 2 — the first-principles opacity surface, and how it differs from line lists.

    python figures/fig2.py config/fiducial.yaml

**a** The Rosseland mean κ_R(T, ρ) over the 52 first-principles points, band-limited
to 100-4,000 cm⁻¹ and interpolated exactly as the structure model reads it: linear
on a Delaunay triangulation of (T, log₁₀ρ) over log₁₀κ. Outside the convex hull the
model clamps to the nearest point; that is an extrapolation, not a surface, so it
is left blank here rather than shaded. The 52 computed points are drawn on top.

**b** First-principles / line-list ratio against density, coloured by temperature.
Three references: HITEMP at all 52 points, ExoMol POKAZATEL at the 27 points it
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
NGRID = 480  # A91: 1.5x, to soften the stepped hull edge

# shape carries the reference, fill carries the mean
REFERENCES = (
    ("hitemp", "o", 22.0, "HITEMP"),
    ("pokazatel", "s", 19.0, "POKAZATEL"),
    ("hitemp_mtckd", "^", 26.0, "HITEMP + MT_CKD continuum"),
)
# A91: Fig. 2b shows kappa_R against the two HITEMP references only; POKAZATEL
# and both Planck series moved to ED Fig. 2.
PANEL_B_REFERENCES = tuple(r for r in REFERENCES if r[0] != "pokazatel")


def panel_a(ax, fig, cfg, nums: FigureNumbers) -> None:
    t_pts, rho_pts, kr, _ = compute_point_means("fp", cfg, BAND_CM)
    src = get_source("fp", cfg, band=BAND_CM)

    # A91: NGRID 320 -> 480, and contourf rather than pcolormesh. The stepped
    # lower edge was the pcolormesh cell size at the hull boundary, not anything
    # in the model.
    t_lin = np.linspace(t_pts.min(), t_pts.max(), NGRID)
    lr_lin = np.linspace(np.log10(rho_pts.min()), np.log10(rho_pts.max()), NGRID)
    TT, LR = np.meshgrid(t_lin, lr_lin)
    # The model's own triangulation; NaN outside the hull is left blank, because
    # the model clamps there and a clamp is not a surface.
    vals = src.table._lin(np.column_stack([TT.ravel(), LR.ravel()]))[:, 0]
    surface = (10.0 ** vals.reshape(TT.shape)) / CM2_PER_G_TO_SI  # SI -> cm^2/g

    # A91: the model interpolates linearly on a Delaunay triangulation, so the
    # surface it defines is piecewise planar and the facets are real -- they are
    # the model, not a rendering artefact. For DISPLAY ONLY the same 52 points are
    # re-interpolated with a C1-continuous Clough-Tocher scheme, which removes the
    # faceting without moving the data. The hull mask is taken from the model's
    # own interpolant, so the blank region is unchanged, and the largest
    # display-vs-model difference is recorded below rather than left unstated.
    from scipy.interpolate import CloughTocher2DInterpolator

    smooth_fn = CloughTocher2DInterpolator(
        np.column_stack([t_pts, np.log10(rho_pts)]), np.log10(kr))
    smooth = 10.0 ** smooth_fn(TT, LR)
    inside = np.isfinite(surface)
    smooth = np.where(inside, smooth, np.nan)
    with np.errstate(invalid="ignore"):
        disp_dex = np.nanmax(np.abs(np.log10(smooth[inside] / surface[inside])))

    from matplotlib.colors import LogNorm

    finite = np.isfinite(surface)
    norm = LogNorm(vmin=np.nanmin(surface[finite]), vmax=np.nanmax(surface[finite]))
    levels = np.logspace(np.log10(norm.vmin), np.log10(norm.vmax), 48)
    mesh = ax.contourf(10.0**LR, TT, np.ma.masked_invalid(smooth), levels=levels,
                       cmap="viridis", norm=norm, extend="both")
    mesh.set_edgecolor("face")   # hairline seams between the filled bands
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
    # A91: contourf puts a tick on every level, which printed the two extreme
    # values in full precision and nothing else. Round decades inside the range.
    from matplotlib.ticker import FixedFormatter, FixedLocator

    decades = [10.0 ** e for e in range(-6, 7)
               if norm.vmin <= 10.0 ** e <= norm.vmax]
    cb.ax.yaxis.set_major_locator(FixedLocator(decades))
    cb.ax.yaxis.set_major_formatter(FixedFormatter([f"{d:g}" for d in decades]))
    cb.ax.yaxis.set_minor_locator(FixedLocator([]))

    nums.axes("a", ax, "Density (kg/m^3)", "Temperature (K)", xlog=True)
    nums.panel(
        "a",
        n_computed_points=int(t_pts.size),
        T_range_K=[float(t_pts.min()), float(t_pts.max())],
        rho_range_kg_m3=[float(rho_pts.min()), float(rho_pts.max())],
        kappa_R_range_cm2_g=[float(kr.min()), float(kr.max())],
        colour_scale="viridis, log",
        display_interpolant=("Clough-Tocher C1 on the same 52 points, for display "
                             "only; the model itself interpolates linearly on a "
                             "Delaunay triangulation"),
        display_vs_model_max_dex=float(disp_dex),
        hull_mask_from="the model's own linear interpolant",
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

    # A91: Rosseland only, and only the two HITEMP references. The panel carried
    # three references x two means = six overlapping series on one pair of axes,
    # and kappa_R is the quantity the structure model actually integrates. The
    # Planck ratios are ED Fig. 2b and the POKAZATEL comparison is ED Fig. 2c.
    sc = None
    for name, marker, size, _label in PANEL_B_REFERENCES:
        data = series[name]
        if not data:
            continue
        sc = ax.scatter([d["rho"] for d in data], [d["rosseland"] for d in data],
                        c=[d["T"] for d in data], cmap="viridis", marker=marker,
                        s=size, vmin=T_VMIN, vmax=T_VMAX, edgecolor="black",
                        linewidth=0.25, zorder=3)

    ax.axhline(1.0, color="black", lw=0.6, ls=":", zorder=1)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Density ρ (kg m$^{-3}$)")
    ax.set_ylabel("First principles / line list, κ$_R$")
    cb = fig.colorbar(sc, ax=ax, pad=0.012, fraction=0.05, aspect=28)
    cb.set_label("Temperature T (K)", fontsize=7.6)
    cb.ax.tick_params(labelsize=7.0)

    handles = [Line2D([], [], ls="none", marker=m, mfc="#7f7f7f", mec="black",
                      ms=np.sqrt(s), label=lab)
               for _n, m, s, lab in PANEL_B_REFERENCES]
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
    nums.panel(
        "b", references=stats, unity_line=True,
        # A91: the ledger keeps all three references and both means; the PANEL
        # draws only what this says, so a caption cannot claim the rest.
        shown_in_panel={"mean": "Rosseland only",
                        "references": [r[3] for r in PANEL_B_REFERENCES],
                        "moved_to_extended_data": ["POKAZATEL (ED Fig. 2c)",
                                                   "Planck ratios (ED Fig. 2b, 2c)"]},
        temperature_colour_range_K=[T_VMIN, T_VMAX],
        colour_scale="viridis, linear in T",
        line_list_asymmetry=(
            "the two line lists are not matched in isotopologue coverage: the "
            "HITEMP spectra carry isotopologues 1, 2 and 3, the POKAZATEL spectra "
            "the main one only. Both use the same ideal-gas broadening pressure "
            "and the same band. The minor isotopologues sit mostly in the windows, "
            "so they bear on the Rosseland ratio more than the Planck one; the "
            "size of that is not measured here and a caption should not imply the "
            "two references differ only in their line data"),
        pokazatel_matching=(
            "matched state by state, not interpolated: POKAZATEL sits on 27 of the "
            "52 first-principles states exactly"),
    )


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
            "a, Rosseland-mean opacity over the 52 first-principles points, "
            "interpolated as the structure model reads it. b, first-principles / "
            "line-list ratio against density, coloured by temperature; filled "
            "Planck, open Rosseland."),
    )
    rows = crossover_v3_rows(cfg, BAND_CM)
    fig, (ax_a, ax_b) = figure(92.0, ncols=2)
    panel_a(ax_a, fig, cfg, nums)
    panel_b(ax_b, fig, cfg, rows, nums)
    panel_label(ax_a, "a", dx=-0.17, dy=1.06)
    panel_label(ax_b, "b", dx=-0.17, dy=1.06)   # A91: aligned with a

    path = save_pdf(fig, "fig2", nums)
    print(f"wrote {path}")
    for label, s in nums.panels["b"]["references"].items():
        print(f"  {label:<32s} n={s['n_points']:>2d}  Planck "
              f"{s['planck_range'][0]:.3f}-{s['planck_range'][1]:.3f}  Rosseland "
              f"{s['rosseland_range'][0]:.3f}-{s['rosseland_range'][1]:.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
