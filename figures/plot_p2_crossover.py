"""P2 / manuscript Fig. 2b: FP vs line-list opacity ratio, against BOTH references.

    python figures/plot_p2_crossover.py config/fiducial.yaml

History. The panel began as a replot of the owner's delivered
``crossover_hitemp.json`` (27 points, RADIS). A53 recomputed the M2
Planck/Rosseland ratios through our own band-limited table builder over all 48 FP
points, with HITEMP evaluated by the same triangulated interpolation the structure
model uses, and kept the delivered JSON overplotted as the independent check.
A58 adds the **second reference** a reviewer will ask for: the same line list plus
the MT_CKD 4.3 water self-continuum, which is what every planetary model actually
uses. Both ratios are now reported everywhere, and the figure carries both.

Reading the figure. Each FP grid point appears twice, joined by a thin grey stem:

  * **filled, temperature-coloured** — FP / bare HITEMP (line absorption only);
  * **open vermillion** — FP / (HITEMP + MT_CKD self-continuum).

The marker *shape* carries which campaign delivered the point, because the caveats
differ by campaign:

  * **circle**   — one of the 27 original grid points;
  * **triangle** — one of the 12 cold-dense extension points (T ≤ 1000 K,
                   ρ = 100-500), the regime the structural claim lives in;
  * **diamond**  — one of the 9 dilute points (ρ ≤ 60), where A52 found the
                   transparency windows that carry κ_R sit at or below the dipole
                   model's noise floor;
  * **grey ×**   — HITEMP clamped (ρ > 1000, outside its grid). The continuum is
                   quadratic in density, so at these points it is held at the
                   ceiling and *under*-counted.

The horizontal line at 1 is where the two opacities agree; points that fall below
it under the continuum-inclusive reference are ones where the extrapolated
empirical continuum exceeds the first-principles absorption, and they are listed
explicitly in the printed summary and in ``crossover_v3.json``.

Writes ``outputs/crossover_hitemp_v2.json`` (bare, unchanged),
``outputs/crossover_v3.json`` (both references + per-group summary) and
``outputs/P2_crossover.png``.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # run as script or -m

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

from figures._common import band_caption, comparison_band, load_cfg, save, set_style
from opacity.crossover import (
    GROUPS,
    agreement_summary,
    crossover_rows,
    crossover_v3_rows,
    summarize_v3,
    write_json,
    write_v3_json,
)

LINE_LIST = "hitemp"
CONTINUUM_SOURCE = "hitemp_mtckd"
T_VMIN, T_VMAX = 650.0, 3000.0
CONT_COLOR = "#D55E00"  # vermillion (Okabe-Ito)

# marker, size, label — keyed by the group the point belongs to
GROUP_STYLE = {
    "original": ("o", 46, "original 27-point grid"),
    "cold_dense": ("^", 58, "cold-dense extension (T ≤ 1000 K, ρ = 100-500)"),
    "dilute": ("D", 40, "dilute (ρ ≤ 60): floor-limited κ_R (A52)"),
}
GROUP_LABEL = {
    "original": "27 original points",
    "cold_dense": "12 cold-dense points",
    "dilute": "9 dilute points",
}


def _scatter_both(ax, rows, bare_key, cont_key, cmap="viridis"):
    """Both references at every point, joined by a stem so the shift is the message."""
    # stems first, so markers sit on top
    for r in rows:
        ax.plot([r["rho"], r["rho"]], [r[bare_key], r[cont_key]],
                color="#999999", lw=0.7, zorder=1)
    sc = None
    for group, (marker, size, _lab) in GROUP_STYLE.items():
        these = [r for r in rows if r["group"] == group]
        if not these:
            continue
        x = [r["rho"] for r in these]
        sc = ax.scatter(x, [r[bare_key] for r in these], c=[r["T"] for r in these],
                        cmap=cmap, marker=marker, s=size, edgecolor="k", linewidth=0.35,
                        vmin=T_VMIN, vmax=T_VMAX, zorder=3)
        ax.scatter(x, [r[cont_key] for r in these], marker=marker, s=size,
                   facecolors="none", edgecolor=CONT_COLOR, linewidth=1.2, zorder=4)
    clamped = [r for r in rows if r["ll_clamped"]]
    if clamped:
        ax.scatter([r["rho"] for r in clamped], [r[bare_key] for r in clamped],
                   marker="x", s=22, color="#444444", linewidth=0.9, zorder=5)
    return sc


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    cfg = load_cfg(argv)
    set_style()
    band = comparison_band(cfg, ["fp", LINE_LIST])

    # the bare-reference product (A53) is unchanged and still written
    bare_rows = crossover_rows(cfg, band, line_list=LINE_LIST)
    v2_path = write_json(cfg, bare_rows, band, line_list=LINE_LIST)

    rows = crossover_v3_rows(cfg, band, line_list=LINE_LIST,
                             continuum_source=CONTINUUM_SOURCE)
    summary = summarize_v3(rows)
    v3_path = write_v3_json(cfg, rows, band, line_list=LINE_LIST,
                            continuum_source=CONTINUUM_SOURCE)

    # ------------------------------ figure ------------------------------
    fig, (axp, axr) = plt.subplots(1, 2, figsize=(11.6, 5.1), constrained_layout=True)
    sc = None
    for ax, bare_key, cont_key, lab in (
        (axp, "ratio_planck_bare", "ratio_planck_mtckd", "Planck"),
        (axr, "ratio_rosseland_bare", "ratio_rosseland_mtckd", "Rosseland"),
    ):
        sc = _scatter_both(ax, rows, bare_key, cont_key)
        ax.axhline(1.0, color="k", lw=0.8, ls=":")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("Density ρ [kg/m³]")
        ax.set_ylabel(f"FP / {LINE_LIST.upper()} {lab}-mean ratio")
        ax.set_title(f"{lab} mean")
    fig.colorbar(sc, ax=axr, label="Temperature [K]")

    handles = [
        Line2D([], [], ls="none", marker="o", mfc="#7f7f7f", mec="k", ms=7,
               label=f"filled: FP / bare {LINE_LIST.upper()}"),
        Line2D([], [], ls="none", marker="o", mfc="none", mec=CONT_COLOR, ms=7, mew=1.4,
               label=f"open: FP / ({LINE_LIST.upper()} + MT_CKD self-continuum)"),
    ] + [
        Line2D([], [], ls="none", marker=GROUP_STYLE[g][0], mfc="#7f7f7f", mec="k",
               ms=7, label=GROUP_STYLE[g][2])
        for g in GROUPS
    ] + [
        Line2D([], [], ls="none", marker="x", color="#444444", ms=7,
               label="HITEMP clamped (ρ > 1000): continuum under-counted"),
    ]
    axp.legend(handles=handles, fontsize=7.2, loc="upper left")

    n_sub = len(summary["sub_unity_mtckd"])
    n_illposed = sum(not r["decomposition_well_posed"] for r in rows)
    fig.suptitle(
        f"FP / {LINE_LIST.upper()} band-limited opacity ratio versus density, "
        f"{len(rows)}-point FP grid, {band_caption(band)}\n"
        f"against two references: bare line list, and line list + MT_CKD 4.3 self-"
        f"continuum (extrapolated in T and ρ; A58).  "
        f"{n_sub} point{'s' if n_sub != 1 else ''} fall below unity with the continuum; "
        f"{n_illposed}/{len(rows)} have Lorentz half-widths past the 25 cm⁻¹ bound the "
        f"continuum is defined against",
        fontsize=9.5)
    path = save(fig, cfg, "P2_crossover")

    # --------------------------- printed summary ---------------------------
    print(f"wrote {path}")
    print(f"wrote {v2_path}")
    print(f"wrote {v3_path}")
    _print_table(rows, band)
    _print_groups(summary)
    _print_sub_unity(summary)
    _print_cold_dense_shares(rows)
    _print_decomposition_caveat(rows)
    _print_agreement(bare_rows)
    return path


def _print_table(rows, band):
    print(f"\n  TASK 20 — all {len(rows)} FP grid points, both references, "
          f"{band_caption(band)}, sorted by T then ρ.")
    print(f"  share = 1 - κ_R(line)/κ_R(line+cont): the fraction of the combined "
          f"Rosseland mean\n  the continuum carries (= cont/(line+cont) in the grey "
          f"limit).")
    head = (f"\n  {'T [K]':>6s} {'rho':>9s} | {'FP/HITEMP':>19s} | "
            f"{'FP/(HITEMP+MT_CKD)':>19s} | {'MT_CKD':>7s} | {'note':<22s}")
    sub = (f"  {'':>6s} {'':>9s} | {'Planck':>9s} {'Ross':>9s} | "
           f"{'Planck':>9s} {'Ross':>9s} | {'share':>7s} |")
    print(head)
    print(sub)
    print("  " + "-" * (len(head) - 3))
    for r in rows:
        notes = [GROUP_LABEL[r["group"]].split()[1]]
        if r["ll_clamped"]:
            notes.append(f"CLAMPED /{r['continuum_underscale']:.1f}")
        flag = " <" if (r["ratio_planck_mtckd"] < 1.0
                        or r["ratio_rosseland_mtckd"] < 1.0) else "  "
        print(f"  {r['T']:>6.0f} {r['rho']:>9.4g} | "
              f"{r['ratio_planck_bare']:>9.3f} {r['ratio_rosseland_bare']:>9.3f} | "
              f"{r['ratio_planck_mtckd']:>9.3f} {r['ratio_rosseland_mtckd']:>9.3f}{flag}| "
              f"{r['mtckd_share_kappa_R']:>7.3f} | {', '.join(notes):<22s}")
    print("  ('<' marks a point where a continuum-inclusive ratio is below 1.)")


def _print_groups(summary):
    print(f"\n  per-group min/max of each ratio:")
    head = (f"  {'group':<28s} {'n':>3s} | {'FP/HITEMP Planck':>19s} "
            f"{'FP/HITEMP Ross':>19s} | {'FP/(H+MT) Planck':>19s} "
            f"{'FP/(H+MT) Ross':>19s} | {'MT_CKD share κ_R':>19s}")
    print(head)
    print("  " + "-" * (len(head) - 2))
    for group in GROUPS + ("all",):
        s = summary["groups"][group] if group in summary["groups"] else summary["all"]
        cells = []
        for key in ("ratio_planck_bare", "ratio_rosseland_bare",
                    "ratio_planck_mtckd", "ratio_rosseland_mtckd",
                    "mtckd_share_kappa_R"):
            cells.append(f"{s[key]['min']:>8.3f}-{s[key]['max']:<10.3f}")
        label = GROUP_LABEL.get(group, "all 48 points")
        print(f"  {label:<28s} {s['n']:>3d} | {cells[0]} {cells[1]} | "
              f"{cells[2]} {cells[3]} | {cells[4]}")


def _print_sub_unity(summary):
    sub = summary["sub_unity_mtckd"]
    print(f"\n  points where FP/(HITEMP+MT_CKD) < 1 in EITHER mean — the extrapolated\n"
          f"  empirical continuum exceeds the first-principles absorption there, and\n"
          f"  the manuscript has to say so ({len(sub)} of 48):")
    if not sub:
        print("    (none)")
        return
    head = (f"  {'T [K]':>6s} {'rho':>9s} {'group':<12s} | {'Planck':>8s} {'Ross':>8s} | "
            f"{'share κ_R':>10s}  which")
    print(head)
    print("  " + "-" * (len(head) - 2))
    for r in sub:
        which = []
        if r["ratio_planck_mtckd"] < 1.0:
            which.append("Planck")
        if r["ratio_rosseland_mtckd"] < 1.0:
            which.append("Rosseland")
        print(f"  {r['T']:>6.0f} {r['rho']:>9.4g} {r['group']:<12s} | "
              f"{r['ratio_planck_mtckd']:>8.3f} {r['ratio_rosseland_mtckd']:>8.3f} | "
              f"{r['mtckd_share_kappa_R']:>10.3f}  {' + '.join(which)}")


def _print_cold_dense_shares(rows):
    cold = [r for r in rows if r["group"] == "cold_dense"]
    print(f"\n  MT_CKD share of the line-list means at the {len(cold)} cold-dense "
          f"points — the regime\n  the structural claim lives in, and where the "
          f"continuum's (296/T)^n is least extrapolated:")
    head = (f"  {'T [K]':>6s} {'rho':>6s} | {'κ_R line':>10s} {'κ_R +cont':>10s} "
            f"{'share':>7s} | {'κ_P line':>10s} {'κ_P +cont':>10s} {'share':>7s}")
    print(head)
    print("  " + "-" * (len(head) - 2))
    for r in cold:
        print(f"  {r['T']:>6.0f} {r['rho']:>6.0f} | "
              f"{r['kappa_R_line_cm2_g']:>10.4g} {r['kappa_R_line_mtckd_cm2_g']:>10.4g} "
              f"{r['mtckd_share_kappa_R']:>7.3f} | "
              f"{r['kappa_P_line_cm2_g']:>10.4g} {r['kappa_P_line_mtckd_cm2_g']:>10.4g} "
              f"{r['mtckd_share_kappa_P']:>7.3f}")
    shares = [r["mtckd_share_kappa_R"] for r in cold]
    print(f"    κ_R share spans {min(shares):.1%}-{max(shares):.1%} across these 12 "
          f"points. At the top of that range the\n    line list's own absorption is a "
          f"minority of the opacity a planetary model would use —\n    and the "
          f"continuum supplying that majority is MT_CKD evaluated ~{500 / 0.7418:.0f}x "
          f"above its\n    reference density (0.742 kg/m³) and ~350 K above its "
          f"reference temperature (296 K),\n    far outside the range it was fitted "
          f"over. Treat it as what a planetary model would\n    do, not as a "
          f"measurement.")


def _print_decomposition_caveat(rows):
    """Where "line list + continuum" is a decomposition of anything, and where not."""
    from opacity.mtckd import LBLRTM_LINE_BOUND_CM

    ok = [r for r in rows if r["decomposition_well_posed"]]
    print(f"\n  IS THE SPLIT EVEN WELL POSED? MT_CKD is the residual beyond a "
          f"line-by-line calculation\n  truncated at {LBLRTM_LINE_BOUND_CM:g} cm⁻¹ "
          f"with the plinth removed, so it only means something while the Lorentz\n"
          f"  half-width is small against that bound. Using HITEMP's own broadening "
          f"parameters:")
    head = (f"  {'group':<28s} {'n':>3s} | {'Lorentz HWHM [cm^-1]':>22s} | "
            f"{'inside the bound':>17s}")
    print(head)
    print("  " + "-" * (len(head) - 2))
    for group in GROUPS:
        these = [r for r in rows if r["group"] == group]
        hw = [r["lorentz_hwhm_cm"] for r in these]
        n_ok = sum(r["decomposition_well_posed"] for r in these)
        print(f"  {GROUP_LABEL[group]:<28s} {len(these):>3d} | "
              f"{min(hw):>10.3g} - {max(hw):<9.4g} | {n_ok:>10d} / {len(these):<4d}")
    print(f"    {len(rows) - len(ok)} of {len(rows)} points have lines BROADER than the "
          f"truncation the continuum is defined\n    against, by one to two orders of "
          f"magnitude. The delivered spectra truncate at 50 cm⁻¹ with no\n    plinth "
          f"subtraction, so line + MT_CKD double-counts 25-50 cm⁻¹ — biasing the "
          f"combined\n    opacity HIGH and the FP ratios LOW, i.e. the sub-unity "
          f"points above are if anything\n    overstated. Read the continuum-inclusive "
          f"column as 'what a planetary model would\n    compute following standard "
          f"practice here', not as the physically correct opacity.")


def _print_agreement(bare_rows):
    summary = agreement_summary(bare_rows)
    print(f"\n  cross-check of the BARE ratios against the delivered "
          f"crossover_hitemp.json\n  ({sum(s['n'] for s in summary.values())} shared "
          f"points), by how κ_LL was obtained:")
    label = {"on_grid": "both coords on HITEMP's grid",
             "interpolated": "T interpolated (2250, 2750 K)",
             "clamped": "rho clamped (> 1000 kg/m^3)"}
    for key in ("on_grid", "interpolated", "clamped"):
        s = summary[key]
        if not s["n"]:
            continue
        print(f"    {label[key]:<32s} n={s['n']:>2d}  max |dev| planck "
              f"{s['max_abs_dev_planck_pct']:>6.2f}%  rosseland "
              f"{s['max_abs_dev_rosseland_pct']:>6.2f}%")
    print(f"    The first line is the real cross-check and it reproduces A25's ~1% "
          f"agreement.\n    The last is expected: RADIS computed those at the true "
          f"density, we hold κ at ρ = 1000.")


if __name__ == "__main__":
    main()
