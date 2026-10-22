"""Assemble the reader-mode slide deck from the rendered figures.

    python figures/build_deck.py

Each slide is rendered as a 16:9 image (white background, black Arial text, a
single embedded figure occupying the majority of the slide, 2-3 bullets beneath).
The same slide images are assembled into outputs/deck/subneptune_deck.pdf and
into an image-per-slide outputs/deck/subneptune_deck.pptx, so the PDF, the PPTX,
and the inspection images are identical. A manifest lists the figure per slide.

Style: plain declarative text; numbers carry units and provenance; citations in
author-year form; reliability stated on the slide that uses the result.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"
DECK = OUT / "deck"
DECK.mkdir(parents=True, exist_ok=True)

# 16:9 at 150 dpi -> 2000 x 1125 px.
FIG_W, FIG_H, DPI = 13.333, 7.5, 150
plt.rcParams.update({"font.family": "Arial", "font.size": 13})

TITLE_FS, BULLET_FS, NUM_FS = 24, 14.5, 12


# --- slide specifications -------------------------------------------------
SLIDES = [
    dict(kind="title",
         title="First-principles water opacity in steam sub-Neptunes",
         lines=[
             "Interior structure and thermal evolution of pure-steam envelopes on rocky cores",
             "",
             "Model: hydrostatic structure; Schwarzschild convection; Guillot (2010) irradiated",
             "boundary; AQUA equation of state (Haldemann et al. 2020); Zeng et al. (2016) core",
             "",
             "subneptune model, v1 — reader-mode deck (self-contained)",
         ]),

    dict(figure="P2_crossover.png",
         title="Line-list and first-principles opacity disagree",
         bullets=[
             "FP/HITEMP Planck-mean ratio 1.4–18.6 and Rosseland-mean ratio 1.7–39.6 across the "
             "27-point (T,ρ) grid; the ratio increases with density (crossover_hitemp.json; "
             "Rothman et al. 2010).",
             "HITEMP and POKAZATEL agree to ~1.5× in the Rosseland mean and ~5% in the Planck "
             "mean; both differ from FP by the same factor (Polyansky et al. 2018).",
             "Line lists are built and validated for dilute gas; the deep envelope reaches "
             "ρ ~ 10³ kg/m³, T ~ 10³–10⁴ K.",
         ]),

    dict(figure="M1_mechanism.png",
         title="Line lists leave transparent windows the continuum fills",
         bullets=[
             "At T=2000 K, ρ=1200 kg/m³ the line-list window floor is ~10 cm²/g; the FP continuum "
             "floor is ~240 cm²/g (M1_mechanism).",
             "The Rosseland mean is a harmonic mean weighted toward low-opacity windows, so the "
             "window difference sets the mean.",
             "Far-wing and continuum absorption beyond 4000 cm⁻¹ is truncated in the line list. "
             "At this density the collision time approaches the vibrational period, where the "
             "discrete-line description does not apply.",
         ]),

    dict(figure="M2_pipeline.png",
         title="First-principles opacity pipeline with validation gates",
         bullets=[
             "Stages: DFT labels; MLIP fine-tune; molecular dynamics; dipole model; band-mean "
             "opacity (M2_pipeline).",
             "Four validation gates: energy/force RMSE; MD stability; dipole sum rules; "
             "plausibility and crossover against line lists.",
             "Pipeline stage and gate labels are indicative and require confirmation against the "
             "production code.",
         ]),

    dict(figure="Fx_opacity_compare.png",
         title="FP is 2–3× more opaque at the line-list overlap",
         bullets=[
             "Band-limited FP/POKAZATEL Rosseland ratio 1.5–45× and Planck ratio 1.4–19× on band "
             "[100,4000] cm⁻¹, 27 shared grid points; FP/HITEMP at 3 shared points "
             "(Fx_opacity_compare.csv).",
             "FP/HITEMP ratios reproduce the crossover analysis to ~1% on the same band: 3.43 vs "
             "3.42 at (1500 K, 300 kg/m³); 2.11 vs 2.11; 1.52 vs 1.52.",
             "Band choice changes the Rosseland ratio ~2×: [50,4500] cm⁻¹ gives 6.47 vs 3.43 at "
             "(1500,300). The 1e-6 relative floor is inactive for this data.",
         ]),

    dict(figure="M3_anatomy.png",
         title="Steam envelope on a rocky core, radiative outer layer",
         bullets=[
             "dP/dm = −Gm/4πr⁴; ∇ = min(∇_ad, ∇_rad) (Schwarzschild); L_int = 4πR_p²σT_int⁴ "
             "constant (M3_anatomy).",
             "Guillot (2010) two-stream outer boundary; AQUA EOS (Haldemann et al. 2020); Zeng "
             "et al. (2016) rock mass-radius core.",
             "Opacity acts only through ∇_rad, above the RCB. Composite: line-list ρ<100, blend "
             "ρ=100–300, FP ρ>300 kg/m³; radiative-zone clamp 0% (reliable).",
         ]),

    dict(figure="M4_opacity_space.png",
         title="Where planets reach the first-principles grid",
         bullets=[
             "Line-list grid ρ=0.1–1000 kg/m³, T=500–3000 K; FP grid ρ=109–2229 kg/m³, "
             "T=1000–3000 K; blend window ρ=100–300 kg/m³ (M4_opacity_space).",
             "At T_int=75 K the RCB is at ρ≈2 kg/m³ (line-list regime); at T_int=20 K the RCB is "
             "at ρ≈200 kg/m³ (FP grid).",
             "Background is the FP/line-list Planck ratio (1–3) for ρ≥109 kg/m³, band "
             "[100,4000] cm⁻¹.",
         ]),

    dict(figure="F1_tp_profiles.png",
         title="Line-list radiative window at low internal temperature",
         bullets=[
             "At T_int=20 K the composite and HITEMP temperature profiles coincide above "
             "~150 bar and diverge below it (F1_tp_profiles).",
             "HITEMP has 3 radiative-convective boundaries and a radiative zone at ρ=150–500 "
             "kg/m³ (∇_rad/∇_ad < 1); the composite has 1 boundary and no deep radiative zone.",
             "FP is 2–3× more opaque there, so the composite deep interior stays convective and "
             "warmer. Pure FP and POKAZATEL are opacity-extrapolated here and are not shown.",
         ]),

    dict(figure="D1_fp_domain.png",
         title="Radiogenic floor exceeds the FP-reach temperature",
         bullets=[
             "The line-list radiative zone reaches ρ=109 kg/m³ only below T_int=24–32 K (rising "
             "with M_c); a planet must cool below this to put FP in the radiative interior "
             "(D1_fp_domain).",
             "The radiogenic floor (Earth-like core ≈ M_c × 20 TW, 10 Gyr) is 25–33 K, above the "
             "reach threshold by ~1 K at every mass from 2 to 20 M⊕.",
             "The planet cools to the radiogenic floor and stops, so the FP-radiative regime is "
             "not reached; where it is reached (denser envelopes) the radius change is ≤1.7%.",
         ]),

    dict(figure="E2_history.png",
         title="Escape window unchanged; deep entropy biased when old",
         bullets=[
             "Composite and HITEMP tracks are identical through the first 100 Myr escape window: "
             "ΔR_p=0.00%, RCB at ρ=2–5 kg/m³ (E2_history).",
             "Interior entropy diverges up to ~8% at 2–13 Gyr, inside the radiogenically-supported "
             "epoch (>~0.5 Gyr) the planet does not reach.",
             "Peak radius divergence ~0.3%, reached late. Energy conservation |ΔE−∫L dt|/E < 1% "
             "(Lopez & Fortney 2014; Ginzburg et al. 2016).",
         ]),

    dict(kind="text",
         title="Limitations and v1 omissions",
         bullets=[
             "Composition: pure-steam envelope; no H/He; no mixtures. Heating: no atmospheric "
             "escape; no tidal heating; no radiogenic heating inside the envelope model.",
             "Rosseland mean is band-truncated to [100,4000] cm⁻¹; missing >4000 cm⁻¹ transparency "
             "windows bias κ_R high; band systematic ~2× (docs/band_limitation.md).",
             "FP grid covers ρ=109–2229 kg/m³, T=1000–3000 K; below ρ=0.1 kg/m³ line-list opacity "
             "is held flat in density (dilute limit). POKAZATEL is on the 27-point grid only and "
             "is not a reliable structural source.",
             "G5 residuals vs Aguichine et al. (2021) are 6–23%, from the pure-steam envelope vs "
             "their full water-phase layers and core-composition differences (g5_residuals.json).",
         ]),

    dict(kind="text",
         title="Summary",
         bullets=[
             "FP water opacity is 2–3× higher than HITEMP and POKAZATEL at ρ≥109 kg/m³, T≥1000 K, "
             "where the discrete-line description does not apply and DFT-MD is required.",
             "For steam sub-Neptunes this regime is convective while warm and below the radiogenic "
             "T_int floor while old, so FP changes radius by ≤1.7% (unreached) and does not change "
             "observable cooling; the escape window is unaffected.",
             "The applicability boundary of the FP correction is denser envelopes with the RCB at "
             "ρ≥109 kg/m³ while warm: massive water worlds and magma-ocean-coupled layers.",
         ]),

    dict(kind="text",
         title="References",
         bullets=[
             "Aguichine A., Mousis O., Deleuil M., Marcq E. 2021, ApJ, 914, 84.",
             "Barber R. J., Tennyson J., Harris G. J., Tolchenov R. N. 2006, MNRAS, 368, 1087.",
             "Ginzburg S., Schlichting H. E., Sari R. 2016, ApJ, 825, 29.",
             "Guillot T. 2010, A&A, 520, A27.",
             "Haldemann J., Alibert Y., Mordasini C., Benz W. 2020, A&A, 643, A105.",
             "Kippenhahn R., Weigert A., Weiss A. 2012, Stellar Structure and Evolution, 2nd ed.",
             "Lopez E. D., Fortney J. J. 2014, ApJ, 792, 1.",
             "Mlawer E. J., Payne V. H., Moncet J.-L., et al. 2012, Phil. Trans. R. Soc. A, 370, 2520.",
             "Pierrehumbert R. T. 2010, Principles of Planetary Climate (Cambridge Univ. Press).",
             "Polyansky O. L., Kyuberis A. A., Zobov N. F., et al. 2018, MNRAS, 480, 2597.",
             "Rothman L. S., Gordon I. E., Barber R. J., et al. 2010, JQSRT, 111, 2139.",
             "Zeng L., Sasselov D. D., Jacobsen S. B. 2016, ApJ, 819, 127.",
         ]),

    dict(kind="text",
         title="Backup: validation gates (48 pass)",
         bullets=[
             "G1 isothermal limit: κ→0 recovers the Guillot skin temperature.  "
             "G2 EOS: AQUA table values reproduced to ≤0.5%.",
             "G3 hydrostatic residual ~1e-7 (<1e-6).  G4 opacity-swap null: identical tables give "
             "identical structure.  G7 convergence: ΔR_p <0.1% under halved tolerances.",
             "G5 mass-radius vs Aguichine et al. (2021): residuals 6–23%, documented.  "
             "G6 energy bookkeeping: |ΔE−∫L dt|/E <1% over the cooling track.",
             "Opacity: factor-1000 plausibility; band auto-detection; dilute-flat; radiative-zone "
             "clamp reliability. Full suite: 48 pass, 0 fail.",
         ]),

    dict(kind="table",
         title="Backup: FP vs line-list opacity, selected points",
         table_csv="Fx_opacity_compare.csv",
         note="POKAZATEL, band [100,4000] cm⁻¹ (Fx_opacity_compare.csv). "
              "κ in cm²/g; ratio = FP / line-list."),

    dict(figure="sensitivity_gamma_matching.png",
         title="Backup: boundary sensitivity scans",
         bullets=[
             "γ = κ_vis/κ_th over 0.25–4: R_p varies ~1% (sensitivity_gamma_matching.png).",
             "Matching pressure: adaptive at the RCB; R_p independent of the 1–100 bar cap once "
             "the cap exceeds the RCB pressure.",
             "Rosseland band [50,4500] vs [100,4000] cm⁻¹: ~2× on κ_R (A25). Composite blend window "
             "default ρ=100–300 kg/m³; 30–100 documented as a sensitivity case.",
         ]),

    dict(figure="BF3_fixed_tint.png",
         title="Backup: cooling-time ratio at fixed T_int",
         bullets=[
             "At fixed T_int, L_int = 4πR_p²σT_int⁴ is imposed, so a more opaque, more extended "
             "planet has larger L and shorter τ_KH (BF3_fixed_tint.png).",
             "This ratio is therefore not a measure of the opacity effect on cooling; the "
             "fixed-interior-entropy comparison is used instead.",
             "composite/HITEMP = 1.0; FP/HITEMP ≈ 0.88 and POKAZATEL/HITEMP ≈ 0.98, both "
             "opacity-extrapolated (unreliable). Values below 1 reflect the imposed L.",
         ]),
]


# --- rendering ------------------------------------------------------------
# Arial lacks these four glyphs; map them to Arial-safe forms for slide text.
# (Figures render their own glyphs via matplotlib's font fallback.)
_GLYPH_FIX = {"⁻¹": "^-1", "⊕": "_E", "≲": "<~", "∇": "grad",
              "⁴": "^4", "⁰": "^0", "⁵": "^5", "⁶": "^6", "⁷": "^7",
              "⁸": "^8", "⁹": "^9", "⁻": "^-"}


def _sanitize(text: str) -> str:
    for bad, good in _GLYPH_FIX.items():
        text = text.replace(bad, good)
    return text


def _wrap(text, width):
    import textwrap
    return textwrap.fill(_sanitize(text), width)


def render_slide(spec, number, total):
    fig = plt.figure(figsize=(FIG_W, FIG_H), dpi=DPI)
    fig.patch.set_facecolor("white")
    kind = spec.get("kind", "figure")

    # Slide number (bottom-right) + rule.
    fig.text(0.97, 0.03, f"{number} / {total}", ha="right", va="bottom",
             fontsize=NUM_FS, color="black")
    # Title (baseline near top; rule sits clearly below the descenders).
    fig.text(0.05, 0.955, _sanitize(spec["title"]), ha="left", va="top", fontsize=TITLE_FS,
             fontweight="bold", color="black")
    fig.add_artist(plt.Line2D([0.05, 0.95], [0.905, 0.905], color="black", lw=1.2,
                              transform=fig.transFigure))

    if kind == "title":
        y = 0.66
        for ln in spec["lines"]:
            fig.text(0.05, y, _sanitize(ln), ha="left", va="top", fontsize=17, color="black")
            y -= 0.075
    elif kind == "text":
        y = 0.80
        for b in spec["bullets"]:
            wrapped = _wrap(b, 105)
            fig.text(0.06, y, "•  " + wrapped.replace("\n", "\n    "), ha="left", va="top",
                     fontsize=17, color="black")
            y -= 0.055 + 0.052 * wrapped.count("\n")
    elif kind == "table":
        _render_table(fig, spec)
    else:  # figure slide
        img = plt.imread(OUT / spec["figure"])
        h, w = img.shape[0], img.shape[1]
        # figure area: x 0.05-0.95, y 0.30-0.88 ; fit preserving aspect.
        ax = fig.add_axes([0.05, 0.30, 0.90, 0.58])
        ax.imshow(img)
        ax.axis("off")
        ax.set_aspect("equal")
        # bullets beneath.
        y = 0.255
        for b in spec["bullets"]:
            wrapped = _wrap(b, 128)
            fig.text(0.05, y, "•  " + wrapped.replace("\n", "\n    "), ha="left", va="top",
                     fontsize=BULLET_FS, color="black")
            y -= 0.045 + 0.043 * wrapped.count("\n")

    out = DECK / f"slide_{number:02d}.png"
    fig.savefig(out, dpi=DPI, facecolor="white")
    plt.close(fig)
    return out


def _render_table(fig, spec):
    import csv
    rows = list(csv.DictReader(open(OUT / spec["table_csv"])))
    # Select a spread of points.
    pick = rows[::4][:9]
    ax = fig.add_axes([0.06, 0.20, 0.88, 0.62]); ax.axis("off")
    cols = ["T_K", "rho_kg_m3", "kappaR_FP_cm2_g", "kappaR_LL_cm2_g", "kappaR_ratio",
            "kappaP_FP_cm2_g", "kappaP_LL_cm2_g", "kappaP_ratio"]
    hdr = ["T [K]", "ρ [kg/m³]", "κR FP", "κR LL", "κR ratio", "κP FP", "κP LL", "κP ratio"]
    cell = [[f"{float(r[c]):.4g}" for c in cols] for r in pick]
    tbl = ax.table(cellText=cell, colLabels=hdr, loc="center", cellLoc="center")
    tbl.auto_set_font_size(False); tbl.set_fontsize(12); tbl.scale(1, 1.6)
    for (rr, cc), c in tbl.get_celld().items():
        c.set_edgecolor("black")
        if rr == 0:
            c.set_text_props(fontweight="bold")
    fig.text(0.06, 0.15, _sanitize(spec["note"]), ha="left", va="top", fontsize=13, color="black")


def build_pptx(slide_pngs):
    from pptx import Presentation
    from pptx.util import Inches
    prs = Presentation()
    prs.slide_width = Inches(FIG_W); prs.slide_height = Inches(FIG_H)
    blank = prs.slide_layouts[6]
    for png in slide_pngs:
        s = prs.slides.add_slide(blank)
        s.shapes.add_picture(str(png), 0, 0, width=prs.slide_width, height=prs.slide_height)
    out = DECK / "subneptune_deck.pptx"
    prs.save(out)
    return out


def main():
    total = len(SLIDES)
    pngs = [render_slide(s, i + 1, total) for i, s in enumerate(SLIDES)]

    pdf_path = DECK / "subneptune_deck.pdf"
    with PdfPages(pdf_path) as pdf:
        for png in pngs:
            fig = plt.figure(figsize=(FIG_W, FIG_H))
            ax = fig.add_axes([0, 0, 1, 1]); ax.imshow(plt.imread(png)); ax.axis("off")
            pdf.savefig(fig, dpi=DPI); plt.close(fig)

    pptx_path = build_pptx(pngs)

    manifest = DECK / "manifest.txt"
    with manifest.open("w") as fh:
        for i, s in enumerate(SLIDES):
            fig_name = s.get("figure") or s.get("table_csv") or f"[{s.get('kind','figure')}]"
            fh.write(f"slide {i+1:02d}: {s['title']}  —  {fig_name}\n")

    print(f"wrote {pdf_path}")
    print(f"wrote {pptx_path}")
    print(f"wrote {manifest} ({total} slides)")


if __name__ == "__main__":
    main()
