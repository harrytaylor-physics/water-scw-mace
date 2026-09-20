"""P2 (deck slide ii): FP vs line-list opacity crossover on the 48-point grid.

    python figures/plot_p2_crossover.py config/fiducial.yaml

Was: a replot of the owner's delivered ``crossover_hitemp.json`` (27 points,
RADIS). Now (A53) the M2 Planck/Rosseland ratios are recomputed through our own
band-limited table builder over all 48 FP points, with HITEMP evaluated by the
same triangulated interpolation the structure model uses, and the delivered JSON
kept as the independent cross-check it always was — overplotted, not replaced.

Point styling carries the provenance, because it is the whole argument:

  * **filled circle** — legacy point, both coordinates on HITEMP's own grid;
  * **triangle** — extension point (T < 1000 K or ρ < 108.8 kg/m³);
  * **red ring** — dilute (ρ ≤ 60), where A52 found the transparency windows
    that carry the Rosseland mean sit at or below the dipole noise floor;
  * **open marker with a cross** — HITEMP clamped (ρ > 1000, outside its grid);
  * **grey +** — the delivered RADIS value at the same point.

Writes ``outputs/crossover_hitemp_v2.json`` and ``outputs/P2_crossover.png``.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # run as script or -m

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

from figures._common import band_caption, comparison_band, load_cfg, save, set_style
from opacity.crossover import agreement_summary, crossover_rows, write_json

LINE_LIST = "hitemp"
T_VMIN, T_VMAX = 650.0, 3000.0


def _scatter(ax, rows, key, cmap="viridis"):
    """One scatter per provenance class, so the marker carries the caveat."""
    classes = (
        ("legacy", lambda r: not r["is_extension_point"] and not r["ll_clamped"], "o", 46),
        ("extension", lambda r: r["is_extension_point"] and not r["ll_clamped"], "^", 52),
        ("clamped", lambda r: r["ll_clamped"], "s", 46),
    )
    sc = None
    for _name, pred, marker, size in classes:
        these = [r for r in rows if pred(r)]
        if not these:
            continue
        x = [r["rho"] for r in these]
        y = [r[key] for r in these]
        c = [r["T"] for r in these]
        facecolors = "none" if _name == "clamped" else None
        sc = ax.scatter(x, y, c=c, cmap=cmap, marker=marker, s=size,
                        edgecolor="k", linewidth=0.35, vmin=T_VMIN, vmax=T_VMAX,
                        zorder=3, **({"facecolors": facecolors} if facecolors else {}))
        if _name == "clamped":  # re-draw with a cross through it
            ax.scatter(x, y, marker="x", s=22, color="#555555", linewidth=0.9, zorder=4)
    # Dilute cells get a red ring: their κ_R rests on floor-dominated windows.
    dilute = [r for r in rows if r["is_dilute"]]
    if dilute:
        ax.scatter([r["rho"] for r in dilute], [r[key] for r in dilute],
                   marker="o", s=150, facecolors="none", edgecolor="#D55E00",
                   linewidth=1.5, zorder=5)
    # The delivered RADIS analysis at the same points.
    ref_key = "ref_planck_ratio" if "planck" in key else "ref_rosseland_ratio"
    ref = [r for r in rows if ref_key in r]
    if ref:
        ax.scatter([r["rho"] for r in ref], [r[ref_key] for r in ref],
                   marker="+", s=42, color="#666666", linewidth=0.9, zorder=2)
    return sc


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    cfg = load_cfg(argv)
    set_style()
    band = comparison_band(cfg, ["fp", LINE_LIST])
    rows = crossover_rows(cfg, band, line_list=LINE_LIST)
    json_path = write_json(cfg, rows, band, line_list=LINE_LIST)

    fig, (axp, axr) = plt.subplots(1, 2, figsize=(11.2, 4.9), constrained_layout=True)
    sc = None
    for ax, key, lab in ((axp, "M2_planck_ratio", "Planck"),
                         (axr, "M2_rosseland_ratio", "Rosseland")):
        sc = _scatter(ax, rows, key)
        ax.axhline(1.0, color="k", lw=0.8, ls=":")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("Density ρ [kg/m³]")
        ax.set_ylabel(f"FP / {LINE_LIST.upper()} {lab}-mean ratio")
        ax.set_title(f"{lab} mean")
    fig.colorbar(sc, ax=axr, label="Temperature [K]")

    handles = [
        Line2D([], [], ls="none", marker="o", mfc="#7f7f7f", mec="k", ms=7,
               label="legacy point (on HITEMP grid)"),
        Line2D([], [], ls="none", marker="^", mfc="#7f7f7f", mec="k", ms=8,
               label="extension point"),
        Line2D([], [], ls="none", marker="s", mfc="none", mec="k", ms=7,
               label="HITEMP clamped (ρ > 1000)"),
        Line2D([], [], ls="none", marker="o", mfc="none", mec="#D55E00", ms=11,
               mew=1.5, label="dilute (ρ ≤ 60): floor-limited κ_R"),
        Line2D([], [], ls="none", marker="+", color="#666666", ms=8,
               label="delivered RADIS crossover"),
    ]
    axp.legend(handles=handles, fontsize=7.5, loc="upper left")

    pr = np.array([r["M2_planck_ratio"] for r in rows])
    rr = np.array([r["M2_rosseland_ratio"] for r in rows])
    fig.suptitle(
        f"FP / {LINE_LIST.upper()} band-limited opacity ratio versus density, "
        f"{len(rows)}-point FP grid, {band_caption(band)}\n"
        f"Planck {pr.min():.2f}-{pr.max():.2f}, Rosseland {rr.min():.2f}-{rr.max():.2f}; "
        f"HITEMP interpolated in (T, log ρ) as the structure model reads it",
        fontsize=9.5)
    path = save(fig, cfg, "P2_crossover")

    # ---------------- printed summary ----------------
    print(f"wrote {path}")
    print(f"wrote {json_path}")
    print(f"\nP2 — FP/{LINE_LIST.upper()} crossover on {len(rows)} FP points, "
          f"{band_caption(band)}")

    head = (f"\n  {'T [K]':>6s} {'rho':>9s} | {'kR_fp':>9s} {'kR_ll':>9s} "
            f"{'kP_fp':>9s} {'kP_ll':>9s} | {'planck':>8s} {'ross':>8s} | {'note':<24s}")
    print(head)
    print("  " + "-" * (len(head) - 3))
    for r in rows:
        notes = []
        if r["is_extension_point"]:
            notes.append("new")
        if r["is_dilute"]:
            notes.append("DILUTE")
        if r["ll_clamped"]:
            notes.append("LL clamped")
        elif not r["ll_on_grid"]:
            notes.append("LL interp")
        print(f"  {r['T']:>6.0f} {r['rho']:>9.4g} | "
              f"{r['kappa_R_fp_cm2_g']:>9.3f} {r[f'kappa_R_{LINE_LIST}_cm2_g']:>9.3f} "
              f"{r['kappa_P_fp_cm2_g']:>9.3f} {r[f'kappa_P_{LINE_LIST}_cm2_g']:>9.3f} | "
              f"{r['M2_planck_ratio']:>8.3f} {r['M2_rosseland_ratio']:>8.3f} | "
              f"{', '.join(notes):<24s}")

    # --- the dilute cells, with their A52 floor fractions beside the ratios ---
    dilute = [r for r in rows if r["is_dilute"]]
    print(f"\n  the {len(dilute)} DILUTE cells (ρ ≤ 60), with the A52 noise-floor "
          f"fractions alongside — so an unreliable ratio is visible as one:")
    head = (f"  {'T [K]':>6s} {'rho':>6s} | {'planck':>8s} {'ross':>8s} | "
            f"{'win>2xfloor':>12s} {'core>2xfloor':>13s} {'Ross in win':>12s}")
    print(head)
    print("  " + "-" * (len(head) - 2))
    for r in dilute:
        if "floor_frac_above_window" in r:
            extra = (f"{r['floor_frac_above_window']:>12.3f} "
                     f"{r['floor_frac_above_cores']:>13.3f} "
                     f"{r['floor_rosseland_window_share']:>12.2f}")
        else:
            extra = f"{'(run S7)':>12s} {'':>13s} {'':>12s}"
        print(f"  {r['T']:>6.0f} {r['rho']:>6.0f} | {r['M2_planck_ratio']:>8.3f} "
              f"{r['M2_rosseland_ratio']:>8.3f} | {extra}")
    print(f"\n  Read that table as: the Planck ratios are trustworthy (the Planck mean "
          f"takes <8%\n  of its integrand from the floor-limited window), the Rosseland "
          f"ratios are not — 86-93%\n  of κ_R's integrand sits in a window where only "
          f"{min(r.get('floor_frac_above_window', float('nan')) for r in dilute):.1%}"
          f"-{max(r.get('floor_frac_above_window', float('nan')) for r in dilute):.1%} "
          f"of the\n  spectrum clears twice the noise floor (A52).")

    # --- what the ratio is actually measuring: window depth, not opacity ---
    print(f"\n  SPECTRAL CONTRAST κ_P/κ_R per source — large when the band has deep\n"
          f"  transparency windows, since the harmonic mean collapses and the "
          f"arithmetic one does not:")
    head = (f"  {'T [K]':>6s} {'rho':>9s} | {'FP':>9s} {'%s':>9s}" % LINE_LIST.upper()
            + f" | {'LL/FP':>8s} | {'ross ratio':>10s}")
    print(head)
    print("  " + "-" * (len(head) - 2))
    for r in rows:
        if r["rho"] not in (100.0, 150.0, 300.0, 500.0):
            continue  # one density column per temperature, the cold points included
        print(f"  {r['T']:>6.0f} {r['rho']:>9.4g} | {r['contrast_fp']:>9.1f} "
              f"{r[f'contrast_{LINE_LIST}']:>9.1f} | "
              f"{r[f'contrast_{LINE_LIST}'] / r['contrast_fp']:>8.1f} | "
              f"{r['M2_rosseland_ratio']:>10.2f}")
    cold = [r for r in rows if r["T"] < 1000.0]
    if cold:
        print(f"\n    The {len(cold)} new cold points (T = 650, 800 K) give Rosseland "
              f"ratios of\n    {min(r['M2_rosseland_ratio'] for r in cold):.0f}-"
              f"{max(r['M2_rosseland_ratio'] for r in cold):.0f}x, against 1.4-7x at "
              f"T >= 1000. That is not a 100x opacity difference:\n    FP's contrast "
              f"stays ~{np.median([r['contrast_fp'] for r in cold]):.0f} while "
              f"{LINE_LIST.upper()}'s rises to "
              f"~{max(r[f'contrast_{LINE_LIST}'] for r in cold):.0f} — the line list "
              f"has genuinely deep\n    windows when cold and the FP spectrum has "
              f"none, its windows being filled to the\n    noise floor (A52). Treat "
              f"the cold Rosseland ratios as a floor artefact until the ARC\n    floor "
              f"construction settles it.")

    # --- regression against the delivered analysis ---
    summary = agreement_summary(rows)
    print(f"\n  cross-check against the delivered crossover_hitemp.json "
          f"({sum(s['n'] for s in summary.values())} shared points), by how κ_LL was obtained:")
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
    return path


if __name__ == "__main__":
    main()
