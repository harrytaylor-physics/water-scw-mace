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
# A92: medians over a window, not a value at one wavenumber -- see panel_b.
WINDOW_CM = (2400.0, 2600.0)
BAND_CM = (100.0, 4000.0)
# The line list at two wing cut-offs: 50 cm^-1 is the delivered convention, 500
# the widest A71 scanned. Both are produced by tools/broadening_sensitivity.py.
LL_TRUNCATIONS = ((50.0, "broadening"), (500.0, "broadening_trunc500"))
LL_STYLE = {50.0: {"lw": 0.45, "ls": "-", "alpha": 0.85},
            500.0: {"lw": 0.6, "ls": "--", "alpha": 0.55}}

# A88: four stages, each carrying the number that validates it. The six-box
# version drew the spectral steps as separate stages, which gave the dipole
# autocorrelation the same visual weight as the whole DFT campaign; they are one
# stage of the route to kappa(nu). Captions come from tools/pipeline_facts.py,
# which is also where the ledger's pipeline.* entries come from.
STAGES = [
    ("Configuration\nsampling", "(T, ρ) states"),
    ("Hybrid DFT\nlabels", "energies, forces,\ndipoles"),
    ("MACE potential\n+ dipole model", "trained to the\nDFT labels"),
    ("NVE dynamics → κ(ν)", "dipole autocorrelation,\nκ(ν)"),
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
    """kappa(nu) at one state, with the line list at two wing cut-offs (A92).

    A91 showed the single channel this panel used to mark sits on a rectangular
    pedestal: RADIS truncates each Lorentzian at ``truncation`` from line centre
    and renormalises it to unit area, so where the true HWHM greatly exceeds the
    cap -- O(1e3) cm^-1 here against 50 -- every strong line becomes a
    ~2*cap-wide rectangle. Reading one wavenumber off that is reading the cap.

    So the panel now shows the line list at **both** 50 and 500 cm^-1 and quotes
    medians over a window rather than a value at a point. The two medians differ
    by more than a factor of two, which is the honest statement: this comparison
    is not cut-off-independent, and the figure says so instead of hiding it
    behind one number.
    """
    t_k, rho = POINT
    fp_path = cfg.fp_opacity_path / f"T{t_k:g}_rho{rho:g}_opacity.npz"
    stem = f"T{t_k:g}_rho{rho:g}_ideal.npz"
    ll_paths = {cap: cfg.output_path / sub / stem for cap, sub in LL_TRUNCATIONS}
    for cap, path in ll_paths.items():
        if not path.exists():
            raise SystemExit(
                f"{path} not found — the line-by-line spectrum at the true density "
                f"with a {cap:g} cm^-1 wing cut-off. Produce it with\n"
                f"  python tools/broadening_sensitivity.py config/fiducial.yaml "
                f"--only {t_k:g},{rho:g} --conventions ideal --truncation {cap:g} "
                f"--out-dir {path.parent}")

    with np.load(fp_path) as d:
        wn_fp, k_fp = np.asarray(d["wn"], float), np.asarray(d["kappa"], float)
    ll = {}
    for cap, path in ll_paths.items():
        with np.load(path) as d:
            ll[cap] = (np.asarray(d["wn"], float), np.asarray(d["kappa"], float),
                       float(d["p_bar"]), float(d["truncation_cm"]))

    lo, hi = BAND_CM
    wlo, whi = WINDOW_CM

    def _median(wn, k):
        return float(np.median(k[(wn >= wlo) & (wn <= whi)]))

    med_fp = _median(wn_fp, k_fp)
    med_ll = {cap: _median(w, k) for cap, (w, k, _p, _t) in ll.items()}

    # the window the medians are taken over, behind everything
    ax.axvspan(wlo, whi, color="#000000", alpha=0.07, lw=0, zorder=0)

    for cap, (w, k, _p, _t) in ll.items():
        style = LL_STYLE[cap]
        m = (w >= lo) & (w <= hi)
        ax.plot(w[m], k[m], color=HITEMP_C, lw=style["lw"], ls=style["ls"],
                alpha=style["alpha"], zorder=2,
                label=f"HITEMP, {cap:g} cm$^{{-1}}$ cut-off")
    mf = (wn_fp >= lo) & (wn_fp <= hi)
    ax.plot(wn_fp[mf], k_fp[mf], color=FP_C, lw=1.0, zorder=3,
            label="first principles")

    # median bars across the window, one per curve
    ax.plot([wlo, whi], [med_fp] * 2, color=FP_C, lw=1.6, zorder=6,
            solid_capstyle="butt")
    for cap, (_w, _k, _p, _t) in ll.items():
        ax.plot([wlo, whi], [med_ll[cap]] * 2, color=HITEMP_C, lw=1.6,
                ls=LL_STYLE[cap]["ls"], zorder=6, solid_capstyle="butt")

    ax.set_yscale("log")  # x stays linear: the band structure is the point
    ax.set_xlim(lo, hi)
    allk = [k_fp[mf]] + [k[(w >= lo) & (w <= hi)] for w, k, _p, _t in ll.values()]
    kmin = min(float(a.min()) for a in allk)
    kmax = max(float(a.max()) for a in allk)
    ax.set_ylim(kmin / 2.2, kmax * 12.0)  # headroom for the band envelopes
    ax.set_xlabel("Wavenumber ν (cm$^{-1}$)")
    ax.set_ylabel("Opacity κ (cm$^{2}$ g$^{-1}$)")
    ax.legend(loc="lower left", fontsize=7.0, ncol=1)

    # The medians and the two ratios, as a block in the clear lower-right corner.
    caps = sorted(ll)
    lines = [f"median over {wlo:.0f}–{whi:.0f} cm$^{{-1}}$:",
             f"  first principles  {med_fp:.0f}"]
    for cap in caps:
        lines.append(f"  HITEMP {cap:g} cm$^{{-1}}$  {med_ll[cap]:.0f}"
                     f"   ×{med_fp / med_ll[cap]:.0f}")
    ax.annotate("\n".join(lines), (0.985, 0.025), xycoords="axes fraction",
                ha="right", va="bottom", fontsize=7.0, linespacing=1.35,
                zorder=7)

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
               "line_list_broadening_pressure_bar": ll[caps[0]][2]},
        median_window_cm=list(WINDOW_CM),
        median_first_principles_cm2_g=round(med_fp, 4),
        median_line_list_cm2_g={f"{c:g}": round(med_ll[c], 4) for c in caps},
        median_ratio_fp_over_line_list={f"{c:g}": round(med_fp / med_ll[c], 3)
                                        for c in caps},
        line_list_truncation_cm=[float(c) for c in caps],
        band_envelopes_cm={k: list(v) for k, v in BANDS.items()},
        band_display_names={k: BAND_DISPLAY.get(k, k) for k in BANDS},
        line_list_note=(
            "direct line-by-line computation at the true density; 1200 kg/m^3 is "
            "above the tabulated ceiling of 1000, so the table would clamp"),
        why_two_cut_offs=(
            "A91: the line profiles are truncated at the cut-off and renormalised "
            "to unit area, so where the true HWHM exceeds it a line becomes a "
            "~2*cut-off-wide rectangle. The comparison does not converge in the "
            "cut-off -- these two medians differ by more than 2x -- so the panel "
            "shows both rather than quoting one. 50 cm^-1 is the delivered "
            "convention every line-list table in the manuscript is built at; "
            "500 cm^-1 is the widest A71 studied"),
        single_channel_marker_removed=(
            "the 2550 cm^-1 point marker is gone: A91 showed it sat on the "
            "pedestal of the line at 2543.07, so it measured the cut-off"),
        producer_for_the_500_spectrum=(
            "tools/broadening_sensitivity.py config/fiducial.yaml --only "
            "2000,1200 --conventions ideal --truncation 500 --out-dir "
            "outputs/broadening_trunc500"),
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
    w = p["median_window_cm"]
    print(f"  medians over {w[0]:.0f}-{w[1]:.0f} cm^-1: first principles "
          f"{p['median_first_principles_cm2_g']:.4g} cm^2/g")
    for cap, med in p["median_line_list_cm2_g"].items():
        print(f"    HITEMP at {cap} cm^-1 cut-off: {med:.4g} cm^2/g   "
              f"ratio x{p['median_ratio_fp_over_line_list'][cap]:.4g}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
