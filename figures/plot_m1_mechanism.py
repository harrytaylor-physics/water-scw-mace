"""M1 (deck): why line lists fail in the warm-dense regime — κ(ν) mechanism.

    python figures/plot_m1_mechanism.py config/fiducial.yaml

Monochromatic opacity κ(ν) at a dense, hot grid point (T=2000 K, ρ=1200 kg/m³)
for a line list (POKAZATEL) and first-principles (FP), overlaid on log-y. FP is
drawn as a filled continuum to show it fills the transparent windows the discrete
Voigt lines leave. Three failure modes are annotated directly. This is the deck's
"why" figure.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # run as script or -m

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch

from config import load_config
from figures._common import save, set_style
from opacity.loaders import load_spectrum
from opacity.registry import get_spec

T_POINT, RHO_POINT = 2000.0, 1200.0
FP_C, LL_C = "#0072B2", "#D55E00"  # Okabe–Ito blue / vermillion


def _load(cfg, band):
    fp = load_spectrum(cfg.fp_opacity_path / "T2000_rho1200_opacity.npz", get_spec("fp"), band)
    pk = load_spectrum(cfg.data_root_path / "pokazatel_opacity" / "T2000_rho1200.npz",
                       get_spec("pokazatel"), band)
    return fp, pk


def main(argv=None):
    cfg = load_config((argv or sys.argv[1:])[0])
    set_style()
    band = (50.5, 4499.5)
    fp, pk = _load(cfg, band)

    fig, ax = plt.subplots(figsize=(8.4, 5.4))
    ax.fill_between(pk.wn_cm, 1e-1, pk.kappa_cm2_g, color=LL_C, lw=0.3, alpha=0.85,
                    label="POKAZATEL (line list): discrete Voigt lines")
    ax.plot(fp.wn_cm, fp.kappa_cm2_g, color=FP_C, lw=1.6,
            label="FP (DFT-MD): collective continuum")
    ax.set_yscale("log")
    ax.set_xlabel("Wavenumber ν [cm⁻¹]")
    ax.set_ylabel("κ(ν) [cm²/g]")
    ax.set_ylim(1.0, 8e3)
    ax.set_xlim(band)

    # (a) discrete lines + transparent windows vs continuum (point to the deep
    # ~2500 cm⁻¹ window where LL collapses but FP holds a continuum floor).
    ax.annotate(f"(a) discrete lines leave transparent windows;\n"
                f"     the collective continuum fills them\n"
                f"     (window floor here: FP ≈ {np.percentile(fp.kappa_cm2_g,10):.0f} "
                f"vs LL ≈ {np.percentile(pk.kappa_cm2_g,10):.0f} cm²/g)",
                xy=(2520, 6), xytext=(1350, 1.5), fontsize=8.5, color="k",
                arrowprops=dict(arrowstyle="->", color="k", lw=1.0))
    # (b) far-wing / continuum truncation at the band edges.
    for edge in (band[0], band[1]):
        ax.axvline(edge, color="grey", ls=":", lw=1.0)
    ax.annotate("(b) far-wing / continuum absorption\n     truncated beyond the line-list band",
                xy=(band[1], 20), xytext=(3050, 3), fontsize=8.5, color="grey",
                arrowprops=dict(arrowstyle="->", color="grey", lw=1.0))
    # (c) the line concept itself dissolves in the warm-dense regime (white gap
    # between the ν≈1500 and ν≈3700 line-list humps).
    ax.text(2350, 1500,
            "(c) at T=2000 K, ρ=1200 kg/m³ the collision time\n"
            "     approaches the vibrational period — the discrete-line\n"
            "     picture itself dissolves; only a first-principles\n"
            "     (DFT-MD) treatment captures it",
            fontsize=8.5, color=FP_C, style="italic", ha="center", va="center",
            bbox=dict(boxstyle="round", fc="white", ec=FP_C, alpha=0.9))

    ax.set_title(f"Monochromatic opacity κ(ν) at T={T_POINT:.0f} K, "
                 f"ρ={RHO_POINT:.0f} kg/m³\n"
                 f"sources: POKAZATEL (line list), FP (first-principles); "
                 f"native resolution, no band mean applied", fontsize=9)
    ax.legend(fontsize=9, loc="upper right")
    path = save(fig, cfg, "M1_mechanism")
    print(f"wrote {path}")
    return path


if __name__ == "__main__":
    main()
