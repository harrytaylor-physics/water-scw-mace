"""M2 (deck): the DFT-MD → opacity pipeline schematic (methods slide).

    python figures/plot_m2_pipeline.py config/fiducial.yaml

A single horizontal band: DFT labels → MLIP fine-tune → gated MD → dipole model
→ κ(ν), with validation gates drawn as first-class elements between stages.
Schematic only (no data); plain, colorblind-safe.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # run as script or -m

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrow, FancyBboxPatch

from config import load_config
from figures._common import save, set_style

STAGE_C = "#0072B2"
GATE_C = "#009E73"

STAGES = [
    ("DFT labels", "PBE/SCAN energies\n& forces on H₂O\nconfigurations"),
    ("MLIP fine-tune", "machine-learned\ninteratomic potential\ntrained to DFT"),
    ("Gated MD", "NVT/NPT dynamics\nat (T, ρ) grid points\n— gates below"),
    ("Dipole model", "polarization /\ndipole autocorrelation\n→ absorption"),
    ("κ(ν)", "monochromatic\nopacity → band\nRosseland & Planck"),
]
# Validation gates sit under the pipeline as first-class checkpoints.
GATES = [
    (0.5, "energy/force RMSE\nvs held-out DFT"),
    (1.5, "MD stability +\nequilibration"),
    (2.5, "sum rules /\nf-sum on dipole"),
    (3.5, "plausibility +\ncrossover vs line list"),
]


def main(argv=None):
    cfg = load_config((argv or sys.argv[1:])[0])
    set_style()
    fig, ax = plt.subplots(figsize=(11.5, 4.6))
    ax.set_xlim(-0.3, len(STAGES)); ax.set_ylim(-1.7, 1.5)
    ax.axis("off")

    w, y = 0.82, 0.35
    for i, (title, body) in enumerate(STAGES):
        box = FancyBboxPatch((i + 0.09, y - 0.55), w, 1.1, boxstyle="round,pad=0.02,rounding_size=0.06",
                             fc="white", ec=STAGE_C, lw=2.0)
        ax.add_patch(box)
        ax.text(i + 0.5, y + 0.34, title, ha="center", va="center", fontsize=11,
                fontweight="bold", color=STAGE_C)
        ax.text(i + 0.5, y - 0.16, body, ha="center", va="center", fontsize=8.3, color="#333")
        if i < len(STAGES) - 1:
            ax.add_patch(FancyArrow(i + 0.93, y, 0.14, 0, width=0.02, head_width=0.12,
                                    head_length=0.06, length_includes_head=True, color="#555"))

    # Gates as diamonds under the arrows.
    for x, label in GATES:
        ax.plot(x + 0.5, -0.55, "D", ms=13, color=GATE_C, mec="k", mew=0.5)
        ax.plot([x + 0.5, x + 0.5], [y - 0.55, -0.42], color=GATE_C, lw=1.0, ls=":")
        ax.text(x + 0.5, -1.05, label, ha="center", va="center", fontsize=7.8, color=GATE_C)
    ax.text(-0.12, -0.55, "gates:", ha="right", va="center", fontsize=9, color=GATE_C,
            fontweight="bold")

    ax.text(len(STAGES) / 2.0, 1.28,
            "First-principles water opacity pipeline with validation gates",
            ha="center", va="center", fontsize=10.5, fontweight="bold")
    ax.text(len(STAGES) / 2.0, -1.55,
            "Output is the FP opacity source; combined with HITEMP or POKAZATEL in "
            "the composite (line-list dilute, FP dense).",
            ha="center", va="center", fontsize=8.5, color="#555")
    path = save(fig, cfg, "M2_pipeline")
    print(f"wrote {path}")
    return path


if __name__ == "__main__":
    main()
