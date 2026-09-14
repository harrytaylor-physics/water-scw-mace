"""S5: like-for-like audit of the G5 anchor definitions (A43).

    python figures/plot_s5_g5_audit.py config/fiducial.yaml

Before attributing the G5 residual to physics, three *definitions* have to match
Aguichine et al. (2021):

  1. irradiation — their Equation (9) defines T_irr = T_eff·√(R⋆/2a), i.e. the
     full-redistribution, zero-albedo equilibrium temperature (verified against
     their own Table 3: GJ 9827 b/c/d, 1184/820/686 K). Our gate divided it by
     √2, treating it as Guillot's irradiation temperature — one √2 too cold.
  2. radius — their reported radius is the *transiting* radius, interior
     boundary plus the atmosphere thickness taken at 0.1 Pa; ours stops at
     `boundary.p_top_bar` = 1e-5 bar = 1 Pa, one decade deeper.
  3. core — their x_core = CMF/(1−WMF) is the core fraction of the refractory
     part only, and Figure 8 has no metallic core, so x_core = 0 maps onto our
     cmf = 0 with M_c = (1−WMF)·M_p. This one already matched.

Each correction is applied separately so its contribution to the median |rel| is
visible, then together. Defaults are unchanged: the corrections are `anchor.*`
config flags, off by default.

Outputs ``outputs/S5_g5_audit.png`` and ``outputs/g5_residuals_conventions.json``.
"""

from __future__ import annotations

import csv
import json
import sys
from dataclasses import replace
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # run as script or -m

import matplotlib.pyplot as plt
import numpy as np

from eos import AquaEOS
from figures._common import load_cfg, save, set_style
from structure import build_structure
from validation.test_g5_literature import _REF, anchor_t_eq, our_radius

T_INT_G5 = 30.0
SOURCE = "hitemp"
THEIR_RADIUS_PRESSURE_BAR = 1.0e-6  # 0.1 Pa, the level their atmosphere thickness uses

# (label, anchor overrides) — applied one at a time, then together.
CONVENTIONS = [
    ("baseline (A32)", {}),
    ("+ irradiation (T_eq = T_irr)", {"t_irr_to_t_eq": "identity"}),
    ("+ radius level (0.1 Pa)", {"radius_pressure_bar": THEIR_RADIUS_PRESSURE_BAR}),
    ("+ both", {"t_irr_to_t_eq": "identity",
                "radius_pressure_bar": THEIR_RADIUS_PRESSURE_BAR}),
]
COLORS = ["#0072B2", "#D55E00", "#009E73", "#CC79A7"]


def _points():
    with _REF.open() as fh:
        return [r for r in csv.DictReader(l for l in fh if not l.startswith("#"))]


def residuals(cfg, eos, **anchor_overrides):
    """The eight G5 residuals under one set of anchor conventions."""
    c = replace(cfg, anchor=replace(cfg.anchor, **anchor_overrides))
    rows = []
    for r in _points():
        wmf, t_irr = float(r["wmf"]), float(r["t_irr_k"])
        m, r_ref = float(r["m_earth"]), float(r["r_earth"])
        r_ours = our_radius(c, m, wmf, t_irr, T_INT_G5, SOURCE, eos)
        rows.append({
            "wmf": wmf, "t_irr_k": t_irr, "m_earth": m, "r_aguichine": r_ref,
            "t_eq_used_k": anchor_t_eq(c, t_irr),
            "p_top_bar": c.anchor.radius_pressure_bar or c.boundary.p_top_bar,
            "r_ours": r_ours, "rel": (r_ours - r_ref) / r_ref,
        })
    return rows


def matching_pressure_sensitivity(cfg, eos, t_eq_k, t_int_k, source=SOURCE):
    """R_p across the configured 1-100 bar matching-pressure scan, at the fiducial."""
    out = []
    for cap in cfg.boundary.match_pressure_scan_bar:
        c = replace(cfg, boundary=replace(cfg.boundary, match_pressure_bar=float(cap)))
        res = build_structure(c, t_eq_k, t_int_k, source, eos=eos, n_sample=400)
        out.append({
            "cap_bar": float(cap), "radius_rearth": res.radius_p_rearth,
            "matched_at_bar": res.match_pressure_actual_bar,
            "matched_at_rcb": bool(res.matched_at_rcb),
        })
    return out


def main(argv=None):
    cfg = load_cfg(argv if argv is not None else sys.argv[1:])
    set_style()
    eos = AquaEOS.from_config(cfg)

    runs = [(label, residuals(cfg, eos, **kw), kw) for label, kw in CONVENTIONS]
    medians = {label: float(np.median([abs(x["rel"]) for x in rows]))
               for label, rows, _ in runs}

    t_eq_fid = float(cfg.irradiation.t_eq_k[0])
    scan = matching_pressure_sensitivity(cfg, eos, t_eq_fid, T_INT_G5)
    radii = [s["radius_rearth"] for s in scan]
    match_spread = (max(radii) - min(radii)) / min(radii)

    # ---------------- figure ----------------
    fig, (ax_p, ax_m) = plt.subplots(1, 2, figsize=(13.0, 5.2),
                                     gridspec_kw={"width_ratios": [2.1, 1.0]},
                                     constrained_layout=True)
    x = np.arange(8)
    width = 0.2
    labels = [f"{r['m_earth']:.0f} $M_⊕$\nwmf {r['wmf']:.1f}\n{r['t_irr_k']:.0f} K"
              for r in runs[0][1]]
    for k, (label, rows, _) in enumerate(runs):
        ax_p.bar(x + (k - 1.5) * width, [r["rel"] * 100 for r in rows], width=width,
                 color=COLORS[k], label=f"{label} (median |rel| {medians[label]:.1%})")
    ax_p.axhline(0.0, color="k", lw=0.8)
    ax_p.axhspan(-10, 10, color="grey", alpha=0.15, label="SPEC §7 G5 aim (±10%)")
    ax_p.set_xticks(x, labels, fontsize=7.5)
    ax_p.set_ylabel(r"$(R_{ours} - R_{Aguichine})\,/\,R_{Aguichine}$ [%]")
    ax_p.set_title("(a) G5 residuals per anchor point, by definition correction", fontsize=9.5)
    ax_p.legend(fontsize=7.5)

    ax_m.bar(np.arange(len(runs)), [medians[label] * 100 for label, _, _ in runs],
             color=COLORS)
    short = {"baseline (A32)": "baseline", "+ irradiation (T_eq = T_irr)": "+ irradiation",
             "+ radius level (0.1 Pa)": "+ radius\nlevel", "+ both": "+ both"}
    ax_m.set_xticks(np.arange(len(runs)), [short.get(label, label) for label, _, _ in runs],
                    fontsize=8.5)
    ax_m.set_ylabel("median |rel| [%]")
    ax_m.set_title("(b) median |ΔR_p/R_p| by correction", fontsize=9.5)
    for i, (label, _, _) in enumerate(runs):
        ax_m.text(i, medians[label] * 100 + 0.3, f"{medians[label] * 100:.1f}", ha="center",
                  fontsize=8.5)

    fig.suptitle("S5 — like-for-like audit of the G5 anchor definitions "
                 "(corrections are opt-in; defaults unchanged)", fontsize=10.5)
    path = save(fig, cfg, "S5_g5_audit")

    out_dir = cfg.output_path
    out_dir.mkdir(parents=True, exist_ok=True)
    json_path = out_dir / "g5_residuals_conventions.json"
    json_path.write_text(json.dumps({
        "anchor": "Aguichine et al. 2021 ApJ 914 84, Fig 8 (digitized); fit Eq. 29",
        "source": SOURCE, "t_int_k": T_INT_G5,
        "note": "Definition audit (A43); repository defaults are unchanged.",
        "conventions": [
            {"label": label, "anchor_overrides": kw,
             "median_abs_rel": medians[label], "residuals": rows}
            for label, rows, kw in runs
        ],
        "matching_pressure_scan": scan,
        "matching_pressure_rp_spread_frac": match_spread,
    }, indent=2))

    # ---------------- printed summary ----------------
    print(f"wrote {path}")
    print(f"wrote {json_path}")
    print(f"\nS5 — G5 definition audit (source={SOURCE}, T_int={T_INT_G5:g} K, "
          f"cmf={cfg.anchor.cmf:g}, n_sample=400)")
    print("  (1) irradiation: their Eq. (9) T_irr = T_eff*sqrt(R_star/2a) IS the "
          "full-redistribution T_eq;")
    print("      factor applied by the repository default: T_eq = T_irr/sqrt(2) = "
          f"{1 / np.sqrt(2):.4f}*T_irr; corrected factor: 1.0")
    print(f"  (2) radius: their reported radius is the transiting radius, atmosphere "
          f"thickness at 0.1 Pa = {THEIR_RADIUS_PRESSURE_BAR:g} bar; ours is the "
          f"top-of-atmosphere radius at p_top = {cfg.boundary.p_top_bar:g} bar "
          f"(NOT the {cfg.boundary.match_pressure_bar:g}-bar matching radius)")
    print(f"  (3) core: their x_core = CMF/(1-WMF) over the refractory part only; "
          f"Fig. 8 has no metallic core -> x_core = 0 = our cmf = {cfg.anchor.cmf:g}, "
          f"with M_c = (1-WMF)*M_p. Already matched: no correction.")

    print(f"\n  median |rel| by correction:")
    base = medians[CONVENTIONS[0][0]]
    for label, _, _ in runs:
        delta = medians[label] - base
        print(f"    {label:32s} {medians[label]:.4f}"
              + ("" if label == CONVENTIONS[0][0] else f"   ({delta:+.4f} vs baseline)"))

    head = (f"\n  {'wmf':>4s} {'T_irr':>6s} {'M':>5s} {'R_ref':>6s} "
            + " ".join(f"{lab[:14]:>15s}" for lab, _, _ in runs))
    print(head)
    for i in range(8):
        r0 = runs[0][1][i]
        cells = " ".join(f"{run[1][i]['r_ours']:>7.4f}/{run[1][i]['rel'] * 100:>6.2f}%"
                         for run in runs)
        print(f"  {r0['wmf']:>4.1f} {r0['t_irr_k']:>6.0f} {r0['m_earth']:>5.1f} "
              f"{r0['r_aguichine']:>6.2f} {cells}")
    print("     (each cell: R_p [R_earth] / relative residual)")

    print(f"\n  matching-pressure sensitivity at the fiducial "
          f"(M_c={cfg.core.mass_mearth:g} M_earth, f_env={cfg.envelope.f_env:.0%}, "
          f"T_eq={t_eq_fid:g} K, T_int={T_INT_G5:g} K, source={SOURCE}):")
    for s in scan:
        print(f"    cap {s['cap_bar']:>6.1f} bar -> R_p = {s['radius_rearth']:.5f} R_earth "
              f"(matched at {s['matched_at_bar']:.3g} bar"
              f"{', at the RCB' if s['matched_at_rcb'] else ''})")
    print(f"    total R_p spread over the {min(s['cap_bar'] for s in scan):g}-"
          f"{max(s['cap_bar'] for s in scan):g} bar scan: {match_spread * 100:.3f}%")
    return path


if __name__ == "__main__":
    main()
