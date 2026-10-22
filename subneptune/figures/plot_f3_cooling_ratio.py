"""F3: Kelvin–Helmholtz cooling-time systematic vs (f_env, T_eq) (SPEC §8, revised).

    python figures/plot_f3_cooling_ratio.py config/fiducial.yaml

Reframed after the FP-validity finding (A27–A31). The reliable cooling-time
comparison is between the *line lists*: at fixed interior entropy, τ_source/τ_HITEMP
for the sources that actually cover the radiative zone (composite ≈ HITEMP; the
FP contribution does not reach the opacity-relevant zone for these planets, so
pure FP is shown flagged as an extrapolation artifact, not a headline). The
dominant, reliable opacity systematic on cooling is the **line-list choice**
(POKAZATEL vs HITEMP), ~2×.

The printed summary reports the ratio's min/max and the grid extents so the
manuscript can quote the range (A40). The plotted curve is, as before, the f_env
scan at the middle T_eq of the config grid; ``--t-eq-scan`` additionally walks the
whole config T_eq grid for the *printed* range (SPEC §8's (f_env, T_eq) grid) and
leaves the figure unchanged.
"""

from __future__ import annotations

import math
import sys
from dataclasses import replace
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # run as script or -m

import matplotlib.pyplot as plt

from eos import AquaEOS
from figures._common import (
    SOURCE_COLOR,
    SOURCE_LABEL,
    available_sources,
    band_caption,
    comparison_band,
    load_cfg,
    save,
    set_style,
)
from structure import build_at_entropy, build_structure

F_ENV_GRID = (0.01, 0.02, 0.03, 0.05, 0.10)


def _ratio_row(cfg, eos, name, baseline, t_eq, t_int_ref, band):
    """τ_KH(source)/τ_KH(baseline) across the f_env grid at one T_eq."""
    ratio, reliable = [], True
    for f_env in F_ENV_GRID:
        cfg_f = replace(cfg, envelope=replace(cfg.envelope, f_env=f_env))
        base = build_structure(cfg_f, t_eq, t_int_ref, baseline, eos=eos, n_sample=400, band=band)
        src = build_at_entropy(cfg_f, base.interior_entropy_j_kg_k, t_eq, name,
                               eos=eos, band=band, n_sample=400)
        ratio.append(src.tau_kh_s / base.tau_kh_s)
        reliable &= src.clamp_reliable
    return ratio, reliable


def ratio_extremes(grids: dict, name: str, t_eqs) -> tuple[tuple | None, tuple | None]:
    """(min, max) of one source's ratio over the scanned (T_eq, f_env) grid.

    ``grids`` maps T_eq -> {source: [ratio per f_env]}. Each returned entry is
    (ratio, T_eq, f_env), so the manuscript can quote where the extreme sits.
    Cells that did not converge are stored as NaN and skipped; if every cell for
    a source is NaN the result is (None, None).
    """
    vals = [(grids[te][name][i], te, F_ENV_GRID[i])
            for te in t_eqs for i in range(len(F_ENV_GRID))
            if not math.isnan(grids[te][name][i])]
    if not vals:
        return None, None
    return min(vals), max(vals)


def main(argv=None):
    import argparse

    argv = sys.argv[1:] if argv is None else argv
    ap = argparse.ArgumentParser()
    ap.add_argument("config", type=Path, help="run config YAML")
    ap.add_argument("--t-eq-scan", action="store_true",
                    help="also print the ratio range over the whole config T_eq grid "
                         "(SPEC §8's (f_env, T_eq) grid); the figure is unchanged")
    args = ap.parse_args(argv)
    cfg = load_cfg([str(args.config)])
    set_style()
    eos = AquaEOS.from_config(cfg)
    t_int_ref = float(cfg.interior.t_int_k[len(cfg.interior.t_int_k) // 2])
    t_eq = float(cfg.irradiation.t_eq_k[len(cfg.irradiation.t_eq_k) // 2])

    baseline = "hitemp"
    others = [s for s in available_sources(cfg, ["composite", "pokazatel", "fp"])]
    plotted = [baseline] + others
    band = comparison_band(cfg, plotted)

    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    plotted_ratios = {}
    for name in others:
        ratio, reliable_all = _ratio_row(cfg, eos, name, baseline, t_eq, t_int_ref, band)
        plotted_ratios[name] = ratio
        ax.plot(F_ENV_GRID, ratio, marker="o", ms=5, color=SOURCE_COLOR[name],
                ls="-" if reliable_all else "--",
                label=f"{SOURCE_LABEL[name]} / HITEMP" + ("" if reliable_all else "  [UNRELIABLE]"))

    ax.axhline(1.0, color="k", lw=0.8, ls=":")
    ax.set_xscale("log")
    ax.set_yscale("log")  # FP artifact (~25×) and the reliable ~1-2× both readable
    ax.set_xlabel("$f_{env}$")
    ax.set_ylabel(r"$\tau_{KH}$(source) / $\tau_{KH}$(HITEMP)  at fixed interior entropy")
    ax.set_title("τ_KH ratio to HITEMP at fixed interior entropy\n"
                 f"composite equals HITEMP; FP and POKAZATEL opacity-extrapolated "
                 f"(unreliable); $T_{{eq}}$={t_eq:.0f} K, $T_{{int,ref}}$={t_int_ref:.0f} K, "
                 f"{band_caption(band)}", fontsize=9)
    ax.legend(fontsize=8.5)
    path = save(fig, cfg, "F3_cooling_ratio")
    print(f"wrote {path}  (baseline=HITEMP; reliable sources: composite/pokazatel)")

    # -- printed range for the manuscript (A40) --
    t_eq_grid = [float(x) for x in cfg.irradiation.t_eq_k] if args.t_eq_scan else [t_eq]
    grids = {t_eq: plotted_ratios}
    for te in t_eq_grid:
        if te in grids:
            continue
        row = {}
        for n in others:
            try:
                row[n] = _ratio_row(cfg, eos, n, baseline, te, t_int_ref, band)[0]
            except ValueError as exc:
                # A fixed-entropy model can be unreachable over the T_int bracket
                # for an extrapolated source; report the cell, do not drop the run.
                print(f"  [T_eq={te:g} K, {n}] no fixed-entropy model: {exc}")
                row[n] = [float("nan")] * len(F_ENV_GRID)
        grids[te] = row
    scanned = sorted(grids) if args.t_eq_scan else [t_eq]

    print(f"\nF3 — tau_KH(source)/tau_KH({baseline}) at fixed interior entropy, "
          f"{band_caption(band)}")
    print(f"  grid: f_env {min(F_ENV_GRID):g}-{max(F_ENV_GRID):g} "
          f"({len(F_ENV_GRID)} points: {', '.join(f'{f:g}' for f in F_ENV_GRID)})")
    print(f"        T_eq {min(scanned):g}-{max(scanned):g} K "
          f"({len(scanned)} point{'s' if len(scanned) > 1 else ''}: "
          f"{', '.join(f'{t:g}' for t in scanned)})"
          + ("" if args.t_eq_scan else "  [plotted T_eq only; --t-eq-scan for the full grid]"))
    print(f"        T_int,ref {t_int_ref:g} K, M_c {cfg.core.mass_mearth:g} M_earth")
    for name in others:
        lo, hi = ratio_extremes(grids, name, scanned)
        if lo is None:
            print(f"  {name:10s}/{baseline}: no converged cell on the scanned grid")
            continue
        print(f"  {name:10s}/{baseline}: min {lo[0]:.4f} (T_eq={lo[1]:g} K, f_env={lo[2]:g}), "
              f"max {hi[0]:.4f} (T_eq={hi[1]:g} K, f_env={hi[2]:g})")
    return path


if __name__ == "__main__":
    main()
