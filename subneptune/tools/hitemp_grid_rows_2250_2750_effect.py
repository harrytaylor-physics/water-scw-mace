#!/usr/bin/env python3
"""Task 23 / A60: what the 2,250 and 2,750 K HITEMP rows changed.

    python tools/hitemp_grid_rows_2250_2750_effect.py config/fiducial.yaml

The two temperatures were absent from the line-list grid, so first-principles
points at 2,250 and 2,750 K were compared against a line list *interpolated* in
temperature — which A59 identified as the source of the only 20-30% deviations
from the independent crossover analysis. This measures what adding the rows did.

Method. The added spectra are moved aside, the crossover is computed, they are
put back, and it is computed again. That is the only honest way to get a "before"
once the files exist: a Delaunay triangulation is not local, so a *partial* grid
gives numbers that belong to no grid at all — with one of the fourteen rows
present, a ratio at 650 K had already moved by 0.4%.

Writes ``outputs/crossover_v3_before_a60.json`` and
``outputs/crossover_v3_after_a60.json``, and prints the four affected points
before and after plus the grid-wide range of every ratio.
"""

from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from config import Config, load_config
from opacity.build_tables import clear_cache
from opacity.crossover import crossover_v3_rows, summarize_v3

BAND_CM = (100.0, 4000.0)
ADDED_TEMPERATURES = (2250.0, 2750.0)
RATIO_KEYS = ("ratio_planck_bare", "ratio_rosseland_bare",
              "ratio_planck_mtckd", "ratio_rosseland_mtckd")


def added_files(cfg: Config) -> list[Path]:
    """The spectra at the added temperatures, by their stored T (not by name)."""
    import numpy as np

    from opacity.registry import get_spec

    spec = get_spec("hitemp")
    out = []
    for path in sorted(cfg.spectra_dir(spec.subdir).glob(spec.glob)):
        with np.load(path) as d:
            if any(abs(float(d["T"]) - t) < 0.5 for t in ADDED_TEMPERATURES):
                out.append(path)
    return out


def _rows(cfg: Config) -> list[dict]:
    clear_cache()  # the tables are cached by (name, band, ...), not by file mtime
    return crossover_v3_rows(cfg, BAND_CM)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("config", type=Path)
    ap.add_argument("--expect", type=int, default=14,
                    help="number of added spectra expected (default 14)")
    args = ap.parse_args(argv)
    cfg = load_config(args.config)

    files = added_files(cfg)
    if len(files) != args.expect:
        raise SystemExit(
            f"found {len(files)} spectra at T = "
            f"{', '.join(f'{t:g}' for t in ADDED_TEMPERATURES)} K, expected "
            f"{args.expect}. A partial grid is not a grid — wait for the rebuild to "
            f"finish, or pass --expect if the target really has changed."
        )

    out_dir = cfg.output_path
    out_dir.mkdir(parents=True, exist_ok=True)
    stash = out_dir / "_hitemp_grid_rows_stash"
    stash.mkdir(parents=True, exist_ok=True)

    try:
        for path in files:
            shutil.move(str(path), str(stash / path.name))
        before = _rows(cfg)
    finally:  # put them back whatever happens
        for path in files:
            src = stash / path.name
            if src.exists():
                shutil.move(str(src), str(path))
        stash.rmdir()
    after = _rows(cfg)

    for name, rows in (("before", before), ("after", after)):
        (out_dir / f"crossover_v3_{name}_a60.json").write_text(json.dumps({
            "meta": (
                f"FP vs line-list crossover {name.upper()} the T = 2250 and 2750 K "
                f"HITEMP rows were added (A60). Band "
                f"[{BAND_CM[0]:g}, {BAND_CM[1]:g}] cm^-1."),
            "state": name,
            "added_temperatures_k": list(ADDED_TEMPERATURES),
            "n_added_spectra": len(files),
            "summary": summarize_v3(rows),
            "results": rows,
        }, indent=1) + "\n")

    key = lambda r: (r["T"], r["rho"])  # noqa: E731
    b = {key(r): r for r in before}
    a = {key(r): r for r in after}
    affected = sorted(k for k in a
                      if any(abs(k[0] - t) < 0.5 for t in ADDED_TEMPERATURES))

    print(f"\n  TASK 23 — the {len(files)} added line-list spectra at T = "
          f"{', '.join(f'{t:g}' for t in ADDED_TEMPERATURES)} K")
    print(f"  Before, those first-principles points were compared against a line "
          f"list interpolated\n  in temperature across the 2000-2500 and 2500-3000 K "
          f"gaps. Now they are not.")
    head = (f"\n  {'T [K]':>6s} {'rho':>8s} | {'FP/HITEMP Planck':>21s} | "
            f"{'FP/HITEMP Rosseland':>21s} | {'FP/(H+MT) Planck':>21s} | "
            f"{'FP/(H+MT) Rosseland':>21s}")
    sub = (f"  {'':>6s} {'':>8s} | {'before':>9s} {'after':>9s} | "
           f"{'before':>9s} {'after':>9s} | {'before':>9s} {'after':>9s} | "
           f"{'before':>9s} {'after':>9s}")
    print(head)
    print(sub)
    print("  " + "-" * (len(sub) - 2))
    for k in affected:
        cells = []
        for rk in RATIO_KEYS:
            cells.append(f"{b[k][rk]:>9.3f} {a[k][rk]:>9.3f}")
        print(f"  {k[0]:>6.0f} {k[1]:>8.4g} | " + " | ".join(cells))

    print(f"\n  change at those {len(affected)} points, as a percentage:")
    for k in affected:
        pct = "  ".join(f"{rk.replace('ratio_', ''):<17s} "
                        f"{100 * (a[k][rk] / b[k][rk] - 1):>+7.2f}%"
                        for rk in RATIO_KEYS)
        print(f"    T={k[0]:>5.0f} rho={k[1]:>8.4g}: {pct}")

    print(f"\n  grid-wide min/max of each ratio over all {len(after)} points:")
    head = f"  {'ratio':<24s} | {'before':>21s} | {'after':>21s}"
    print(head)
    print("  " + "-" * (len(head) - 2))
    for rk in RATIO_KEYS:
        bv = [r[rk] for r in before]
        av = [r[rk] for r in after]
        print(f"  {rk.replace('ratio_', ''):<24s} | "
              f"{min(bv):>9.3f}-{max(bv):<11.3f} | {min(av):>9.3f}-{max(av):<11.3f}")

    moved = sorted(
        ((max(abs(a[k][rk] / b[k][rk] - 1) for rk in RATIO_KEYS), k) for k in a),
        reverse=True)
    print(f"\n  the five points that moved most, including ones at other "
          f"temperatures — a\n  triangulation is not local, so adding rows moves "
          f"values away from them too:")
    for frac, k in moved[:5]:
        tag = " (added row)" if any(abs(k[0] - t) < 0.5 for t in ADDED_TEMPERATURES) else ""
        print(f"    T={k[0]:>5.0f} rho={k[1]:>8.4g}: {100 * frac:>6.2f}%{tag}")
    unaffected = [f for f, k in moved
                  if not any(abs(k[0] - t) < 0.5 for t in ADDED_TEMPERATURES)]
    print(f"    largest move at a point NOT at an added temperature: "
          f"{100 * max(unaffected):.2f}%")
    print(f"\n  wrote {out_dir / 'crossover_v3_before_a60.json'}")
    print(f"  wrote {out_dir / 'crossover_v3_after_a60.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
