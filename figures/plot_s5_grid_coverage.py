"""S5: does a planet's envelope ever cross the FP data? (A48)

    python figures/plot_s5_grid_coverage.py config/fiducial.yaml

A41 found that the headline closure at the old fiducial is made of FP opacity
clamped 300 K below the FP grid's floor, and A42/A47 found that almost no planet
puts a line-list radiative layer inside the FP box at all. This figure shows why,
directly: the 27 FP grid points and their interpolation hull on the (T, log ρ)
plane, with six envelope trajectories drawn through it, each segment coloured by
what the composite opacity call there actually was —

  * **interpolated** — inside the used sub-source's hull: real data;
  * **clamped in T** — outside the grid's temperature range (the FP grid starts
    at 1000 K), nearest-neighbour extrapolation;
  * **clamped in ρ** — outside its density range;
  * **in-box, outside the hull** — inside the grid's (T, ρ) box but outside the
    scattered hull (A6): the gap the next campaign can close *without* extending
    the grid's range at all;
  * **dilute-flat** — line-list-dominated below its ρ floor, κ held flat in ρ
    (A26): a justified extrapolation, not a clamp.

The printed summary gives the envelope mass fraction in each category per
trajectory, then the campaign specification: the (T, ρ) nodes — on a 100 K by
0.25 dex lattice — that would have to be added for the FP grid's hull to enclose
the in-box part of every trajectory.

Note: ``S5_g5_audit`` is a different figure (A43); this one writes
``outputs/S5_grid_coverage.png`` and ``outputs/S5_grid_coverage.csv``.
"""

from __future__ import annotations

import csv
import sys
from dataclasses import replace
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # run as script or -m

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

from eos import AquaEOS
from figures._common import band_caption, comparison_band, load_cfg, save, set_style
from opacity import get_source
from structure import FP_GRID_FLOOR_SI, build_structure
from structure.diagnostics import (
    CATEGORIES,
    campaign_nodes,
    category_mass_fractions,
    clamped_t_split,
    classify_evaluations,
)

FP_T_FLOOR_K = 1000.0
N_SAMPLE = 1500
T_STEP_K = 100.0       # campaign lattice spacing in temperature
LOG_RHO_STEP = 0.25    # ... and in log10 density
CAT_COLOR = {
    "interpolated": "#009E73",        # green: real data
    "in_box_out_of_hull": "#0072B2",  # blue: the fillable gap
    "clamped_t": "#D55E00",           # vermillion: outside the grid in T
    "clamped_rho": "#CC79A7",         # purple: outside the grid in ρ
    "dilute_flat": "#999999",         # grey: justified dilute extrapolation
}
CAT_LABEL = {
    "interpolated": "interpolated (inside hull)",
    "in_box_out_of_hull": "in-box, outside hull",
    "clamped_t": "clamped in T",
    "clamped_rho": "clamped in ρ",
    "dilute_flat": "dilute-flat (ρ < 0.1)",
}


def case(label, t_eq_k, t_int_k, wmf, m_core_mearth, note=""):
    return {"label": label, "t_eq_k": float(t_eq_k), "t_int_k": float(t_int_k),
            "wmf": float(wmf), "mass_core_mearth": float(m_core_mearth), "note": note}


# (a), (b), (d) are fixed; (c) is filled from the A47 scan CSV at run time.
FIXED_CASES = [
    case("(a) old fiducial", 500, 20, 0.05, 5, "A28/A35 headline case"),
    case("(b) candidate", 700, 20, 0.05, 5, "A42 in-domain case"),
]
HOT_CASE = case("(d) hot", 1100, 30, 0.20, 5, "hot irradiation, water-rich")


def scan_cases(csv_path: Path, n=3):
    """The ``n`` A47 cells with the largest HITEMP in-box radiative fraction."""
    import figures.plot_s4_domain_scan as s4

    if not Path(csv_path).exists():
        print(f"  (no {csv_path}; skipping the (c) cases — run the A47 scan first)")
        return []
    rows = s4.load_rows(csv_path)
    ranked = s4._rank(rows, "hitemp_in_box_radiative_mass_fraction", n)
    key = "hitemp_in_box_radiative_mass_fraction"
    if all(s4._num(r, key, 0.0) <= 0.0 for r in ranked):
        # A47: the in-box fraction is zero in every cell, so that ranking is a
        # 108-way tie. Break it on the binding constraint — how dense the
        # radiative zone gets, since every cell fails the box on ρ or on T — and
        # take one cell per T_eq, so the three panels are three different planets
        # rather than three envelope fractions of the same one.
        print("  (c): no extended-scan cell has ANY in-box radiative mass (A47); "
              "breaking the tie on max radiative ρ, one cell per T_eq")
        key = "hitemp_max_radiative_rho"
        ranked, seen = [], set()
        for r in sorted(rows, key=lambda r: s4._num(r, key, 0.0), reverse=True):
            if r["skipped"] or r["t_eq_k"] in seen:
                continue
            seen.add(r["t_eq_k"])
            ranked.append(r)
            if len(ranked) == n:
                break
    out = []
    for i, r in enumerate(ranked, 1):
        out.append(case(f"(c{i}) scan top-{i}", r["t_eq_k"], r["t_int_k"], r["wmf"],
                        r["mass_core_mearth"],
                        f"A47 rank {i} by {key.replace('hitemp_', '')} = "
                        f"{s4._num(r, key, 0.0):.6f}"))
    return out


def build_case(cfg, eos, band, spec, n_sample=N_SAMPLE):
    """Composite + HITEMP structures for one case, with the composite's categories."""
    cfg_c = replace(
        cfg,
        core=replace(cfg.core, mass_mearth=spec["mass_core_mearth"]),
        envelope=replace(cfg.envelope, f_env=spec["wmf"]),
    )
    out = dict(spec)
    for name in ("composite", "hitemp"):
        out[name] = build_structure(cfg_c, spec["t_eq_k"], spec["t_int_k"], name,
                                    eos=eos, n_sample=n_sample, band=band)
    source = get_source("composite", cfg_c, band=band)
    prof = out["composite"].profile
    out["categories"] = classify_evaluations(source, prof["temperature_k"], prof["rho_si"])
    out["mass_fractions"] = category_mass_fractions(out["composite"], out["categories"])
    return out


def fp_table(cfg, band):
    """The FP source's table (27 scattered points on (T, log ρ))."""
    return get_source("fp", cfg, band=band).table


def _plot_trajectory(ax, res, categories, label):
    t = np.asarray(res.profile["temperature_k"], float)
    rho = np.asarray(res.profile["rho_si"], float)
    # One scatter per category keeps the segment colours honest where the
    # category flips sample to sample (it does, near the blend window).
    for cat in CATEGORIES:
        sel = categories == cat
        if sel.any():
            ax.plot(t[sel], rho[sel], ".", ms=2.4, color=CAT_COLOR[cat], zorder=3)
    ax.plot(t, rho, "-", lw=0.5, color="#444444", alpha=0.45, zorder=2)
    # RCBs: where the Schwarzschild flip happens along this (composite) profile.
    conv = np.asarray(res.profile["is_convective"], bool)
    flips = np.flatnonzero(np.diff(conv.astype(int)) != 0)
    if flips.size:
        ax.plot(np.sqrt(t[flips] * t[flips + 1]), np.sqrt(rho[flips] * rho[flips + 1]),
                "o", ms=6, mfc="none", mec="black", mew=1.3, zorder=4)
    ax.set_title(label, fontsize=9.5)


def main(argv=None):
    import argparse

    argv = sys.argv[1:] if argv is None else argv
    ap = argparse.ArgumentParser()
    ap.add_argument("config", type=Path)
    ap.add_argument("--n-sample", type=int, default=N_SAMPLE)
    args = ap.parse_args(argv)
    cfg = load_cfg([str(args.config)])
    set_style()
    eos = AquaEOS.from_config(cfg)
    band = comparison_band(cfg, ["composite", "hitemp"])
    out_dir = cfg.output_path
    out_dir.mkdir(parents=True, exist_ok=True)

    specs = list(FIXED_CASES)
    specs += scan_cases(out_dir / "S4_domain_scan_extended.csv")
    specs.append(HOT_CASE)
    print(f"building {len(specs)} trajectories x 2 sources, n_sample={args.n_sample} ...")
    cases = [build_case(cfg, eos, band, s, args.n_sample) for s in specs]

    table = fp_table(cfg, band)
    t_grid, rho_grid = table.temperature_k, table.rho_si

    # ---------------- figure ----------------
    ncol = 3
    nrow = int(np.ceil(len(cases) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(5.0 * ncol, 4.3 * nrow),
                             sharex=True, sharey=True, constrained_layout=True)
    axes = np.atleast_1d(axes).ravel()
    from scipy.spatial import ConvexHull

    hull = ConvexHull(np.column_stack([t_grid, np.log10(rho_grid)]))
    verts = np.append(hull.vertices, hull.vertices[0])
    seg_t, seg_lr = [], []
    for a, b in zip(verts[:-1], verts[1:]):  # densify: straight in (T, log rho)
        seg_t.append(np.linspace(t_grid[a], t_grid[b], 40))
        seg_lr.append(np.linspace(np.log10(rho_grid[a]), np.log10(rho_grid[b]), 40))
    hull_t = np.concatenate(seg_t)
    hull_r = 10.0 ** np.concatenate(seg_lr)

    for ax, c in zip(axes, cases):
        ax.fill(hull_t, hull_r, color="#0072B2", alpha=0.10, zorder=0,
                label="FP interpolation hull")
        ax.plot(hull_t, hull_r, "-", color="#0072B2", lw=1.2, zorder=1)
        ax.plot(t_grid, rho_grid, "s", ms=3.5, color="#0072B2", zorder=5,
                label="FP grid points (27)")
        _plot_trajectory(ax, c["composite"], c["categories"],
                         f"{c['label']}\n$T_{{eq}}$={c['t_eq_k']:.0f} K, "
                         f"$T_{{int}}$={c['t_int_k']:.0f} K, wmf={c['wmf']:.2f}, "
                         f"$M_c$={c['mass_core_mearth']:.0f} $M_⊕$")
        ax.axhline(FP_GRID_FLOOR_SI, color="#0072B2", ls=":", lw=1.0)
        ax.axvline(FP_T_FLOOR_K, color="#0072B2", ls=":", lw=1.0)
        ax.set_yscale("log")
        ax.set_xscale("log")
        ax.set_xticks([300, 500, 1000, 2000, 4000, 6000])
        ax.get_xaxis().set_major_formatter(plt.matplotlib.ticker.ScalarFormatter())
        ax.set_xlabel("T [K]")
        ax.set_ylabel(r"$\rho$ [kg m$^{-3}$]")
    for ax in axes[len(cases):]:
        ax.set_visible(False)

    handles = [Line2D([], [], ls="none", marker="o", ms=5, color=CAT_COLOR[c],
                      label=CAT_LABEL[c]) for c in CATEGORIES]
    handles += [
        Line2D([], [], ls="none", marker="s", ms=5, color="#0072B2", label="FP grid point"),
        Line2D([], [], ls="none", marker="o", ms=7, mfc="none", mec="black",
               label="RCB (composite)"),
    ]
    axes[0].legend(handles=handles, fontsize=7.5, loc="lower right", ncol=1)
    fig.suptitle("S5 — FP grid coverage of envelope trajectories: what the composite "
                 f"opacity call actually was, sample by sample ({band_caption(band)})",
                 fontsize=10.5)
    path = save(fig, cfg, "S5_grid_coverage")

    # ---------------- campaign specification ----------------
    campaign = campaign_nodes(
        [(c["composite"], c["categories"]) for c in cases], table,
        t_step_k=T_STEP_K, log_rho_step=LOG_RHO_STEP, t_floor_k=FP_T_FLOOR_K,
        rho_floor_si=FP_GRID_FLOOR_SI,
    )

    csv_path = out_dir / "S5_grid_coverage.csv"
    with csv_path.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["kind", "temperature_k", "rho_si"])
        for t, r in campaign["fill_nodes"]:
            w.writerow(["fill_in_box", f"{t:.0f}", f"{r:.4f}"])
        for t, r in campaign["extend_nodes"]:
            w.writerow(["extend_below_t_floor", f"{t:.0f}", f"{r:.4f}"])

    # ---------------- printed summary ----------------
    print(f"wrote {path}")
    print(f"wrote {csv_path}")
    print(f"\nS5 — FP grid coverage, n_sample={args.n_sample}, {band_caption(band)}")
    print(f"  FP grid: {table.n_points} points, T [{t_grid.min():g}, {t_grid.max():g}] K, "
          f"rho [{rho_grid.min():g}, {rho_grid.max():g}] kg/m^3")

    head = (f"  {'trajectory':<22s} {'T_eq':>5s} {'T_int':>5s} {'wmf':>5s} {'M_c':>4s} | "
            + " ".join(f"{CAT_LABEL[c][:13]:>14s}" for c in CATEGORIES))
    print("\n  envelope mass fraction by opacity-call category (composite):")
    print(head)
    print("  " + "-" * (len(head) - 2))
    for c in cases:
        frac = c["mass_fractions"]
        print(f"  {c['label']:<22s} {c['t_eq_k']:>5.0f} {c['t_int_k']:>5.0f} "
              f"{c['wmf']:>5.2f} {c['mass_core_mearth']:>4.0f} | "
              + " ".join(f"{frac[cat]:>14.6f}" for cat in CATEGORIES))
    print("\n  the 'clamped in T' column splits into two different asks "
          f"(FP grid covers T [{t_grid.min():g}, {t_grid.max():g}] K):")
    for c in cases:
        split = clamped_t_split(c["composite"], c["categories"], table)
        print(f"    {c['label']:<22s} too COLD for the grid: {split['too_cold']:.6f}   "
              f"too HOT (deep convective interior): {split['too_hot']:.6f}")
    for c in cases:
        if c["note"]:
            print(f"    {c['label']}: {c['note']}")

    print(f"\n  CAMPAIGN SPECIFICATION — nodes on a {T_STEP_K:.0f} K x "
          f"{LOG_RHO_STEP:g} dex lattice needed so the FP hull encloses the IN-BOX part "
          f"(T >= {FP_T_FLOOR_K:.0f} K, rho >= {FP_GRID_FLOOR_SI:g}) of every trajectory:")
    if not campaign["fill_nodes"]:
        print("    none — the in-box part of every trajectory is already inside the hull "
              "(or no trajectory enters the box at all).")
    else:
        for t, r in campaign["fill_nodes"]:
            print(f"    T = {t:>6.0f} K   rho = {r:>10.3f} kg/m^3   "
                  f"(log10 rho = {np.log10(r):+.2f})")
    print(f"    in-box samples needing FP data: {campaign['n_in_box_needed']} over "
          f"{len(cases)} trajectories; enclosed by the listed nodes: "
          f"{'YES' if campaign['fill_encloses'] else 'NO'}")

    print(f"\n  (addendum, NOT part of the in-box specification) the same lattice, "
          f"extending the grid BELOW its {FP_T_FLOOR_K:.0f} K floor, would be needed to "
          f"cover the dense-but-cold segments where the composite currently clamps in T:")
    if not campaign["extend_nodes"]:
        print("    none.")
    else:
        t_lo = min(t for t, _ in campaign["extend_nodes"])
        t_hi = max(t for t, _ in campaign["extend_nodes"])
        r_lo = min(r for _, r in campaign["extend_nodes"])
        r_hi = max(r for _, r in campaign["extend_nodes"])
        print(f"    {len(campaign['extend_nodes'])} nodes spanning T [{t_lo:.0f}, {t_hi:.0f}] K, "
              f"rho [{r_lo:.1f}, {r_hi:.1f}] kg/m^3 "
              f"({campaign['n_extend_needed']} samples)")
        for t, r in campaign["extend_nodes"][:20]:
            print(f"    T = {t:>6.0f} K   rho = {r:>10.3f} kg/m^3")
        if len(campaign["extend_nodes"]) > 20:
            print(f"    ... and {len(campaign['extend_nodes']) - 20} more (see {csv_path})")

    return path


if __name__ == "__main__":
    main()
