#!/usr/bin/env python3
"""Line-list spectra at states off the production grid, for diagnostics (A65).

    python tools/extra_linelist_points.py --points 1000,60 1500,60 2000,60

Some diagnostics need the line list at a density the production grid does not
carry. ρ = 60 kg/m³ is one: it is a first-principles grid point, but the line-list
grid runs 0.1, 1, 10, 30, 100, 300, 1000, so a band-by-band comparison there has
no line-list spectrum to compare against, and band integrals cannot be
interpolated the way band-limited means can.

**These are written to a separate directory and are deliberately not part of the
production grid.** Adding a density column would change the Delaunay
triangulation and therefore every ratio, figure and gate downstream — for a
diagnostic that needs the spectra and not the table. The registry never sees
this directory; only the tool that asked for the points reads it.

Conventions are inherited from ``generate_linelist_tables.py`` by import, not by
copy, so these spectra are the same calculation as the delivered ones: ideal-gas
broadening pressure at the cell's density, truncation 50 cm⁻¹, wstep 0.01 cm⁻¹,
isotopologues 1,2,3.
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np

import generate_linelist_tables as gen

DEFAULT_OUT = Path("~/subneptune_data/hitemp_extra")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--points", nargs="+", required=True,
                    help="states as T,rho (e.g. 1000,60)")
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--database", default="hitemp")
    ap.add_argument("--wavenum-min", type=float, default=50.0,
                    help="line-loading lower edge; 50 matches A60's rows and "
                         "leaves the [100, 4000] science band unaffected")
    args = ap.parse_args(argv)

    out_dir = args.out.expanduser()
    out_dir.mkdir(parents=True, exist_ok=True)
    states = []
    for spec in args.points:
        t_s, rho_s = spec.split(",")
        states.append((float(t_s), float(rho_s)))

    todo = [(t, r) for t, r in states
            if not gen.out_path(out_dir, t, r).exists()]
    print(f"[extra] {len(states)} states, {len(todo)} to compute, out={out_dir}")
    print(f"[extra] conventions inherited from generate_linelist_tables: "
          f"pressure_mode=ideal, truncation {gen.WSTEP_CM and 50} cm^-1, "
          f"wstep {gen.WSTEP_CM} cm^-1, isotopes {gen.ISOTOPES[args.database]}")
    if not todo:
        print("[extra] nothing to do")
        return 0

    sf = gen.build_factory(args.database, args.wavenum_min)
    for t_k, rho in todo:
        p_pa = rho * gen.R_GAS * t_k / gen.M_H2O_KG_MOL  # ideal, as delivered
        t0 = time.time()
        wn, kappa = gen.compute_point(sf, t_k, p_pa)
        path = gen.out_path(out_dir, t_k, rho)
        tmp = path.with_name(path.stem + ".tmp.npz")
        np.savez_compressed(
            tmp, wn=wn, kappa=kappa, alpha=kappa * (rho / 1000.0), T=t_k, rho=rho,
            rho_kg_m3=rho, p_bar=p_pa / 1e5, source=f"RADIS/{args.database}",
            isotopes=gen.ISOTOPES[args.database], pressure_mode="ideal",
            wstep_cm=gen.WSTEP_CM, wavenum_min_cm=args.wavenum_min,
            note="off-grid diagnostic point; NOT part of the production grid")
        tmp.replace(path)
        print(f"  T={t_k:g} rho={rho:g}  P={p_pa / 1e5:.4g} bar  "
              f"kappa_max={np.nanmax(kappa):.3g} cm^2/g  ({time.time() - t0:.1f}s)",
              flush=True)
    print(f"[extra] done -> {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
