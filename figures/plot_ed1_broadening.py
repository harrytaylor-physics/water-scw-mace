"""Extended Data: what the broadening-pressure convention is worth (Task 21, A59).

    python tools/broadening_sensitivity.py config/fiducial.yaml   # produces the data
    python figures/plot_ed1_broadening.py config/fiducial.yaml    # this panel

The panel is the κ_R ratio at six points spanning the grid, one convention against
another. Three conventions are carried (see ``tools/broadening_sensitivity.py``
for why, and for the correction to the brief's premise that the delivered spectra
use a fixed 1 bar reference — they do not, they use P = ρRT/M per cell):

  ``ref1bar``  fixed 1 bar — the convention the brief assumed;
  ``ideal``    P = ρRT/M — the convention the delivered spectra actually use;
  ``aqua``     the tabulated AQUA P(T, ρ) — the structure model's own EOS.

Two things are plotted, because they answer two different questions:

  * **filled bars, left axis** — κ_R(convention) / κ_R(1 bar). How much the fixed
    1 bar reference the brief asked about would have cost, had it been used.
  * **open markers** — κ_R(AQUA) / κ_R(ideal). The residual systematic that is
    actually present in the delivered tables, since they are built at ideal-gas P
    and the structure model integrates on AQUA.

Mass opacity is convention-independent by construction (κ = σ·N_A/M), so every
difference here is a Voigt-width effect and nothing else. Lorentzian profiles are
truncated at 50 cm⁻¹ from line centre, so the highest-pressure points under-count
absorption and their ratios are lower bounds; the truncation is printed in the
caption.

Writes ``outputs/ED1_broadening.png``.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib.pyplot as plt
import numpy as np

from figures._common import band_caption, load_cfg, save, set_style

IDEAL_COLOR = "#0072B2"  # blue
AQUA_COLOR = "#D55E00"  # vermillion
RESIDUAL_COLOR = "#009E73"  # green

DATA_NAME = "broadening_sensitivity.json"


def load_results(cfg) -> dict:
    path = cfg.output_path / DATA_NAME
    if not path.exists():
        raise SystemExit(
            f"{path} not found — run\n"
            f"    python tools/broadening_sensitivity.py <config>\n"
            f"first (or, on ARC, sbatch tools/broadening_sensitivity.sbatch and copy "
            f"the outputs back)."
        )
    return json.loads(path.read_text())


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    cfg = load_cfg(argv)
    set_style()
    data = load_results(cfg)
    rows = sorted(data["results"], key=lambda r: (r["T"], r["rho"]))
    band = tuple(data["band_cm"])

    labels = [f"{r['T']:.0f} K\n{r['rho']:g} kg/m³" for r in rows]
    x = np.arange(len(rows))
    width = 0.38

    have = lambda key: all(key in r for r in rows)  # noqa: E731
    fig, (ax, ax2) = plt.subplots(2, 1, figsize=(10.0, 7.4), sharex=True,
                                  gridspec_kw={"height_ratios": [2.0, 1.0]},
                                  constrained_layout=True)

    # --- top: kappa_R and kappa_P against the 1 bar reference ---------------
    drawn = False
    for offset, conv, color, label in (
        (-width / 2, "ideal", IDEAL_COLOR, "ideal-gas P = ρRT/M  (the delivered convention)"),
        (+width / 2, "aqua", AQUA_COLOR, "AQUA P(T, ρ)  (the structure model's EOS)"),
    ):
        key = f"ratio_R_{conv}_over_ref1bar"
        if not have(key):
            continue
        drawn = True
        ax.bar(x + offset, [r[key] for r in rows], width, color=color, alpha=0.85,
               edgecolor="k", linewidth=0.4, label=f"κ_R: {label}", zorder=3)
        pkey = f"ratio_P_{conv}_over_ref1bar"
        if have(pkey):
            ax.plot(x + offset, [r[pkey] for r in rows], ls="none", marker="_",
                    ms=13, mew=2.0, color="k", zorder=5,
                    label="κ_P (same convention)" if offset < 0 else None)
    if not drawn:
        raise SystemExit(
            "no 1 bar-referenced ratios in the results — rerun the driver without "
            "restricting --conventions."
        )
    ax.axhline(1.0, color="k", lw=0.9, ls=":", zorder=2)
    ax.set_yscale("log")
    ax.set_ylabel("κ / κ (1 bar reference broadening)")
    ax.legend(fontsize=8.5, loc="best")
    ax.set_title("Broadening at the cell's own pressure, against the 1 bar reference",
                 fontsize=10.5)

    # --- bottom: the systematic actually present in the delivered tables ----
    key = "ratio_R_aqua_over_ideal"
    if have(key):
        ax2.bar(x, [r[key] for r in rows], 0.5, color=RESIDUAL_COLOR, alpha=0.85,
                edgecolor="k", linewidth=0.4, zorder=3)
        for xi, r in zip(x, rows):
            ax2.annotate(f"P: {r['p_aqua_pa'] / r['p_ideal_pa']:.2f}×",
                         (xi, r[key]), textcoords="offset points",
                         xytext=(0, 4 if r[key] >= 1 else -12), ha="center", fontsize=7.5)
    ax2.axhline(1.0, color="k", lw=0.9, ls=":", zorder=2)
    ax2.set_yscale("log")
    ax2.set_ylabel("κ_R(AQUA) / κ_R(ideal)")
    ax2.set_xticks(x)
    ax2.set_xticklabels(labels, fontsize=8.5)
    ax2.set_title("Residual systematic in the delivered tables: AQUA P versus the "
                  "ideal-gas P they were built at", fontsize=10.5)

    fig.suptitle(
        f"Broadening-convention sensitivity of the HITEMP means, {band_caption(band)}\n"
        f"wstep {data['wstep_cm']:g} cm⁻¹, lines loaded over "
        f"[{data['wavenum_range_cm'][0]:g}, {data['wavenum_range_cm'][1]:g}] cm⁻¹, "
        f"Lorentzian truncated at {data['truncation_cm']:g} cm⁻¹ from line centre "
        f"(so the high-P points under-count: these are lower bounds)",
        fontsize=9.5)
    path = save(fig, cfg, "ED1_broadening")

    # ------------------------------ summary -------------------------------
    print(f"wrote {path}")
    print(f"\n  ED1 — broadening-convention sensitivity, {band_caption(band)}")
    head = (f"\n  {'T [K]':>6s} {'rho':>9s} | {'P 1bar':>8s} {'P ideal':>11s} "
            f"{'P AQUA':>11s} | {'kR 1bar':>9s} {'kR ideal':>9s} {'kR aqua':>9s} | "
            f"{'id/1bar':>8s} {'aq/1bar':>8s} {'aq/id':>7s}")
    print(head)
    print("  " + "-" * (len(head) - 3))
    for r in rows:
        cell = lambda k, f="{:>9.4g}": (f.format(r[k]) if k in r else " " * 9)  # noqa: E731
        ratio = lambda k: (f"{r[k]:>8.3f}" if k in r else " " * 8)  # noqa: E731
        print(f"  {r['T']:>6.0f} {r['rho']:>9.4g} | {1.0:>8.3g} "
              f"{r['p_ideal_pa'] / 1e5:>11.5g} {r['p_aqua_pa'] / 1e5:>11.5g} | "
              f"{cell('kappa_R_ref1bar_cm2_g')} {cell('kappa_R_ideal_cm2_g')} "
              f"{cell('kappa_R_aqua_cm2_g')} | "
              f"{ratio('ratio_R_ideal_over_ref1bar')} "
              f"{ratio('ratio_R_aqua_over_ref1bar')} "
              f"{(f'{r[key]:>7.3f}' if key in r else ' ' * 7)}")
    for label, k in (("ideal / 1 bar", "ratio_R_ideal_over_ref1bar"),
                     ("AQUA  / 1 bar", "ratio_R_aqua_over_ref1bar"),
                     ("AQUA  / ideal", "ratio_R_aqua_over_ideal")):
        vals = [r[k] for r in rows if k in r]
        if vals:
            print(f"    kappa_R {label}: {min(vals):.3f}-{max(vals):.3f} "
                  f"over the {len(vals)} points")
    print(f"\n  {data['meta']}")
    return path


if __name__ == "__main__":
    main()
