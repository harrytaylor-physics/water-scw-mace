#!/usr/bin/env python3
"""Figure 1 — the first-principles pipeline, and what it produces (A62).

    python figures/fig1.py config/fiducial.yaml

**a** The route from configuration sampling to a monochromatic opacity, drawn as
vector art rather than replayed from a raster.

**b** κ(ν) at 2,000 K and 1,200 kg/m³, first-principles against HITEMP **at
identical conditions** — same temperature, same density, and therefore the same
ideal-gas broadening pressure of 11,077 bar. The 1,200 kg/m³ point is above
HITEMP's tabulated density ceiling, so the line-list curve here is the direct
line-by-line computation at the true density (A59), not the table's value held at
the ceiling. That distinction is the difference between comparing two opacities
and comparing one opacity with an extrapolation of the other.

A note for the caption. A draft caption marked the 2,550 cm⁻¹ window with 4 and
~200 cm² g⁻¹. At identical conditions the line-list value there is **9.5**, not 4;
**3.0** is what a 1 bar reference broadening would give, and A59 established that
the delivered tables do not use one. The panel is annotated with what it plots.

Writes ``figs/fig1.pdf`` and merges into ``figs/fig_numbers.json``.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

from config import load_config
from figures._manuscript import (
    FP_C,
    HITEMP_C,
    OKABE_ITO,
    FigureNumbers,
    figure,
    panel_label,
    save_pdf,
    set_manuscript_style,
)
from opacity.bands import BANDS

# A91: the 300-1000 cm^-1 window is "rotation" in the manuscript and "libration"
# in the code (opacity/bands.py). Fig. 1b is where a reader meets it first, so it
# carries both names; ED Fig. 4 uses BAND_DISPLAY's single "rotation".
BAND_DISPLAY = {"libration": "rotation/libration"}

POINT = (2000.0, 1200.0)
MARK_WN = 2550.0
BAND_CM = (100.0, 4000.0)

# A88: four stages, each carrying the number that validates it. The six-box
# version drew the spectral steps as separate stages, which gave the dipole
# autocorrelation the same visual weight as the whole DFT campaign; they are one
# stage of the route to kappa(nu). Captions come from tools/pipeline_facts.py,
# which is also where the ledger's pipeline.* entries come from.
STAGES = [
    ("Configuration\nsampling", "(T, ρ) states"),
    ("Hybrid DFT\nlabels", "energies, forces,\ndipoles"),
    ("MACE potential\n+ dipole model", "trained to the\nDFT labels"),
    ("NVE dynamics → κ(ν)", "dipole autocorrelation,\nmass opacity"),
]


def panel_a(ax, cfg, nums: FigureNumbers) -> None:
    from tools.pipeline_facts import pipeline_facts, stage_captions

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    facts = pipeline_facts(cfg)
    if not facts:
        raise SystemExit(
            f"{cfg.fp_opacity_path} holds no spectra, so the stage captions have "
            f"no source — Fig. 1a would be a schematic with blank validation")
    captions = stage_captions(facts)

    n = len(STAGES)
    gap = 0.030  # room for the arrow between boxes
    w = (1.0 - gap * (n - 1)) / n
    y0, h = 0.40, 0.56
    for i, ((title, sub), cap) in enumerate(zip(STAGES, captions)):
        x = i * (w + gap)
        face = OKABE_ITO["sky"] if i < n - 1 else OKABE_ITO["blue"]
        ax.add_patch(FancyBboxPatch(
            (x, y0), w, h, boxstyle="round,pad=0.0,rounding_size=0.02",
            linewidth=0.7, edgecolor="black",
            facecolor=face, alpha=0.28 if i < n - 1 else 0.55, zorder=2))
        ax.text(x + w / 2, y0 + h * 0.70, title, ha="center", va="center",
                fontsize=7.0, fontweight="bold", zorder=3, linespacing=1.25)
        ax.text(x + w / 2, y0 + h * 0.28, sub, ha="center", va="center",
                fontsize=7.0, color="#333333", zorder=3, linespacing=1.25)
        if i < n - 1:
            ax.add_patch(FancyArrowPatch(
                (x + w + 0.003, y0 + h / 2), (x + w + gap - 0.003, y0 + h / 2),
                arrowstyle="-|>", mutation_scale=6.0, linewidth=0.7,
                color="black", zorder=4))
        # What validates the stage, hung under it on a short leader so the reader
        # reads the check against the step rather than hunting a caption.
        ax.plot([x + w / 2] * 2, [y0 - 0.015, y0 - 0.075], color="#888888",
                lw=0.5, zorder=2)
        ax.text(x + w / 2, y0 - 0.095, cap, ha="center", va="top", fontsize=7.0,
                color="#333333", style="italic", zorder=3, linespacing=1.3)
    nums.panel("a", stages=[t.replace("\n", " ") for t, _ in STAGES],
               n_stages=n, kind="schematic; the captions are measured values",
               stage_captions=[c.replace("\n", " ") for c in captions],
               pipeline_facts={k: v["value"] for k, v in facts.items()},
               caption_correction=(
                   "A89 reported a draft caption's 44 meV/A and 7 meV/atom as "
                   "absent from the repository; A90 resolved it. They are the "
                   "HELD-OUT split error (44.4 and 7.0 on valid_Default) and come "
                   "from the training log train_v4_7924966.out, which is not in "
                   "the pull -- recorded in ed_data/mlip_heldout_errors.json. The "
                   "41/58 meV/A parity is a DIFFERENT quantity, training-set "
                   "parity by regime, and both are now shown. The NVE number "
                   "remains a drift rate, 0.017 meV/atom/ps against a 2.0 gate"))
    nums.panel("a", source_files=sorted({v["source_script"] for v in facts.values()}),
               source_per_fact={k: v["source_script"] for k, v in facts.items()},
               fact_notes={k: v["note"] for k, v in facts.items() if v["note"]})


def panel_b(ax, cfg, nums: FigureNumbers) -> None:
    t_k, rho = POINT
    fp_path = cfg.fp_opacity_path / f"T{t_k:g}_rho{rho:g}_opacity.npz"
    ll_path = cfg.output_path / "broadening" / f"T{t_k:g}_rho{rho:g}_ideal.npz"
    if not ll_path.exists():
        raise SystemExit(
            f"{ll_path} not found — it is the direct line-by-line spectrum at the "
            f"true density, produced by tools/broadening_sensitivity.py (A59). Run "
            f"that first."
        )
    with np.load(fp_path) as d:
        wn_fp, k_fp = np.asarray(d["wn"], float), np.asarray(d["kappa"], float)
    with np.load(ll_path) as d:
        wn_ll, k_ll = np.asarray(d["wn"], float), np.asarray(d["kappa"], float)
        p_bar = float(d["p_bar"])

    lo, hi = BAND_CM
    mf = (wn_fp >= lo) & (wn_fp <= hi)
    ml = (wn_ll >= lo) & (wn_ll <= hi)

    ax.plot(wn_ll[ml], k_ll[ml], color=HITEMP_C, lw=0.45, alpha=0.85,
            label="HITEMP line list", zorder=2)
    ax.plot(wn_fp[mf], k_fp[mf], color=FP_C, lw=1.0,
            label="first principles", zorder=3)

    v_fp = float(np.interp(MARK_WN, wn_fp, k_fp))
    v_ll = float(np.interp(MARK_WN, wn_ll, k_ll))
    ax.plot([MARK_WN, MARK_WN], [v_ll, v_fp], color="black", lw=0.6, ls=":", zorder=4)
    # Offsets in points, not data units: the lower marker sits near the axis floor
    # and a multiplicative offset would put its label outside the frame.
    for val, colour, off in ((v_fp, FP_C, (6, 5)), (v_ll, HITEMP_C, (6, -2))):
        ax.plot([MARK_WN], [val], marker="o", ms=3.6, color=colour,
                mec="black", mew=0.4, zorder=5)
        ax.annotate(f"{val:.0f}" if val >= 100 else f"{val:.1f}",
                    (MARK_WN, val), textcoords="offset points", xytext=off,
                    fontsize=7.0, color=colour, fontweight="bold", zorder=6)
    ax.annotate(f"{MARK_WN:.0f} cm$^{{-1}}$: ×{v_fp / v_ll:.0f}",
                (MARK_WN, np.sqrt(v_fp * v_ll)), xytext=(MARK_WN - 90.0, np.sqrt(v_fp * v_ll)),
                fontsize=7.0, ha="right", va="center", zorder=6)

    ax.set_yscale("log")  # x stays linear: the band structure is the point
    ax.set_xlim(lo, hi)
    kmin = min(k_fp[mf].min(), k_ll[ml].min())
    kmax = max(k_fp[mf].max(), k_ll[ml].max())
    ax.set_ylim(kmin / 2.2, kmax * 12.0)  # headroom for the band envelopes
    ax.set_xlabel("Wavenumber ν (cm$^{-1}$)")
    ax.set_ylabel("Opacity κ (cm$^{2}$ g$^{-1}$)")
    ax.legend(loc="lower left", fontsize=7.0, ncol=2)

    # band envelopes along the top, clear of the data
    ytop = ax.get_ylim()[1]
    for name, (blo, bhi) in BANDS.items():
        ax.plot([blo, bhi], [ytop / 2.6] * 2, color=OKABE_ITO["grey"], lw=2.2,
                solid_capstyle="butt", alpha=0.8, zorder=1)
        ax.text(0.5 * (blo + bhi), ytop / 2.1, BAND_DISPLAY.get(name, name),
                ha="center", va="bottom", fontsize=7.0, color="#444444", zorder=1)

    nums.axes("b", ax, "Wavenumber (cm^-1)", "Opacity (cm^2/g)", xlog=False, ylog=True)
    nums.panel(
        "b",
        state={"T_K": t_k, "rho_kg_m3": rho,
               "line_list_broadening_pressure_bar": p_bar},
        marked_wavenumber_cm=MARK_WN,
        marked_kappa_first_principles_cm2_g=round(v_fp, 4),
        marked_kappa_line_list_cm2_g=round(v_ll, 4),
        marked_ratio=round(v_fp / v_ll, 3),
        band_envelopes_cm={k: list(v) for k, v in BANDS.items()},
        band_display_names={k: BAND_DISPLAY.get(k, k) for k in BANDS},
        line_list_note=(
            "direct line-by-line computation at the true density; 1200 kg/m^3 is "
            "above the tabulated ceiling of 1000, so the table would clamp"),
        # A91. The marked channel sits on a RECTANGULAR PEDESTAL, not on line
        # structure. RADIS truncates each Lorentzian at `truncation` from line
        # centre and renormalises it to unit area, so where the true HWHM greatly
        # exceeds the cap -- here it is O(1e3) cm^-1 against a 50 cm^-1 cap --
        # every strong line becomes a ~2*cap-wide rectangle. The pedestal edges
        # pair at exactly 100.01 cm^-1 = 2 x 50 about lines at 2473.76 and
        # 2543.07, and 2550 lies on the second one.
        #
        # Raising the cap does NOT converge the comparison, so no cut-off here is
        # "the right one"; this is A71's systematic, measured at this state:
        line_list_truncation_cm=50.0,
        truncation_sensitivity={
            "note": ("the delivered convention is 50 cm^-1 and every line-list "
                     "table in the manuscript uses it, so this panel does too. "
                     "The ratio is NOT cut-off-independent and the caption must "
                     "name the cap."),
            "kappa_2550_cm2_g": {"50": 9.5213, "200": 12.1568, "500": 28.9124},
            "fp_over_line_list_at_2550": {"50": 35.0, "200": 27.4, "500": 11.5},
            "fp_over_line_list_band_median_2400_2600": {
                "50": 28.1, "200": 27.8, "500": 12.2},
            "pedestal_steps_above_2pct_in_2400_2600": {"50": 8, "200": 1, "500": 0},
            "producer": ("tools/broadening_sensitivity.py --only 2000,1200 "
                         "--conventions ideal --truncation <cap> --out-dir ..."),
        },
        caption_correction=(
            "a draft caption marked 4 and ~200 cm^2/g at 2550 cm^-1; at identical "
            "conditions the line-list value is 9.5, and 3.0 is what a 1 bar "
            "reference broadening gives"),
    )


def main(argv=None) -> int:
    import argparse

    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("config", type=Path)
    args = ap.parse_args(argv if argv is not None else sys.argv[1:])
    cfg = load_config(args.config)

    set_manuscript_style()
    nums = FigureNumbers(
        figure="fig1",
        caption_hint=(
            "a, pipeline from configuration sampling to kappa(nu). b, kappa(nu) at "
            "2000 K, 1200 kg/m^3, first principles vs HITEMP at identical "
            "conditions (same ideal-gas broadening pressure), with the line "
            "profiles truncated at 50 cm^-1 from line centre -- the delivered "
            "convention. The marked ratio depends on that cap (A71): it is 35 at "
            "50 cm^-1, 27 at 200 and 12 at 500."),
    )
    # the schematic grew a caption row under each box, so it needs more height
    fig, (ax_a, ax_b) = figure(132.0, nrows=2, gridspec_kw={"height_ratios": [1.0, 1.9]})
    panel_a(ax_a, cfg, nums)
    panel_b(ax_b, cfg, nums)
    panel_label(ax_a, "a", dx=-0.045, dy=1.14)
    panel_label(ax_b, "b", dx=-0.045, dy=1.06)

    path = save_pdf(fig, "fig1", nums)
    print(f"wrote {path}")
    p = nums.panels["b"]
    print(f"  at {MARK_WN:.0f} cm^-1: first principles {p['marked_kappa_first_principles_cm2_g']:.4g}, "
          f"line list {p['marked_kappa_line_list_cm2_g']:.4g} cm^2/g "
          f"(x{p['marked_ratio']:.3g})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
