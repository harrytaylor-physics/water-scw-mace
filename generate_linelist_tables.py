#!/usr/bin/env python3
"""Generate line-list water opacity tables (HITEMP / ExoMol-POKAZATEL) with RADIS.

Standalone cluster driver for `subneptune` (see docs/linelist_tables.md). Writes
one npz per (T, rho) grid point in the subneptune schema, resumable, offline.

  python generate_linelist_tables.py --database hitemp --out ~/subneptune_data/hitemp_opacity
  python generate_linelist_tables.py --database pokazatel --out ~/subneptune_data/pokazatel_opacity

Assumes a fat compute node with the line-list databank cached locally and NO
internet (no downloads). The (T, rho) grid and the AQUA broadening pressures are
baked in below.

Opacity convention (matches the FP files, correct units, NO x1000):
  kappa[cm^2/g] = sigma(nu)[cm^2/molecule] * N_A / M_H2O[g/mol]
where the cross-section sigma is obtained from the RADIS absorption coefficient
abscoeff[cm^-1] and the number density n = P_broad/(k_B T):  sigma = abscoeff/n.
This makes the mass opacity density-consistent for ANY broadening pressure, so
the broadening-pressure choice (below) only affects Voigt line widths.

Broadening pressure (PRESSURE_MODE):
  "ideal" (default) -> P = rho R T / M   (ideal-gas pressure at the real number
      density; the width gamma_self ~ collision rate ~ real density, which is the
      physically correct impact-broadening pressure). Equals the AQUA pressure in
      the dilute regime where the line list is valid.
  "aqua"            -> the tabulated AQUA P(T, rho) (baked in). Differs from
      "ideal" only in the non-ideal / two-phase regime, where the line list is
      superseded by FP anyway.
Flip PRESSURE_MODE if you prefer the thermodynamic (AQUA) pressure for widths.

The RADIS setup (databank load, isotopes, broadening) is marked CONFIG below;
reconcile it with your ~/scw/hitemp_opacity.py HitempCalculator for your node.
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np

# --- physical constants (CODATA) -------------------------------------------
K_B = 1.380649e-23  # J/K
N_A = 6.02214076e23  # /mol
R_GAS = 8.314462618  # J/mol/K
M_H2O_KG_MOL = 18.01528e-3  # kg/mol
M_H2O_G_MOL = 18.01528  # g/mol

# --- baked grid (docs/linelist_tables.md) ----------------------------------
T_GRID_K = [500, 600, 700, 800, 900, 1000, 1250, 1500, 2000, 2500, 3000]
RHO_GRID_KGM3 = [1e-5, 1e-4, 1e-3, 1e-2, 0.1, 1.0, 10.0, 30.0, 100.0, 300.0, 1000.0]

# AQUA P(T, rho) in Pa; rows = RHO_GRID_KGM3, cols = T_GRID_K.
P_AQUA_PA = [
    [2.3077e00, 2.7692e00, 3.2307e00, 3.6922e00, 4.1538e00, 4.6155e00, 5.7706e00, 6.9407e00, 1.0196e01, 2.6016e01, 4.1074e01],
    [2.3077e01, 2.7692e01, 3.2307e01, 3.6923e01, 4.1538e01, 4.6155e01, 5.7700e01, 6.9313e01, 9.5833e01, 1.7236e02, 3.7486e02],
    [2.3076e02, 2.7692e02, 3.2307e02, 3.6922e02, 4.1538e02, 4.6155e02, 5.7697e02, 6.9270e02, 9.3726e02, 1.3468e03, 2.6358e03],
    [2.3074e03, 2.7690e03, 3.2306e03, 3.6922e03, 4.1537e03, 4.6154e03, 5.7695e03, 6.9249e03, 9.2918e03, 1.2230e04, 1.8341e04],
    [2.3055e04, 2.7676e04, 3.2295e04, 3.6913e04, 4.1531e04, 4.6152e04, 5.7696e04, 6.9242e04, 9.2585e04, 1.1812e05, 1.5407e05],
    [2.2858e05, 2.7538e05, 3.2191e05, 3.6832e05, 4.1468e05, 4.6125e05, 5.7694e05, 6.9237e05, 9.2435e05, 1.1654e06, 1.4435e06],
    [2.0778e06, 2.6139e06, 3.1149e06, 3.6028e06, 4.0840e06, 4.5863e06, 5.7695e06, 6.9236e06, 9.2369e06, 1.1590e07, 1.4088e07],
    [2.6416e06, 6.9092e06, 8.6625e06, 1.0284e07, 1.1843e07, 1.3588e07, 1.7309e07, 2.0771e07, 2.7706e07, 3.4723e07, 4.2025e07],
    [2.6416e06, 1.2354e07, 2.1958e07, 2.8898e07, 3.5282e07, 4.3557e07, 5.7695e07, 6.9234e07, 9.2338e07, 1.1562e08, 1.3952e08],
    [2.6416e06, 1.2354e07, 3.5817e07, 6.2126e07, 8.7856e07, 1.1633e08, 1.7247e08, 2.2279e08, 3.1794e08, 4.0698e08, 4.9238e08],
    [3.5614e08, 5.7727e08, 7.9708e08, 1.0103e09, 1.2155e09, 1.4118e09, 1.8652e09, 2.2773e09, 3.0120e09, 3.6646e09, 4.2617e09],
]

PRESSURE_MODE = "ideal"  # "ideal" | "aqua"

# --- CONFIG: reconcile with your HitempCalculator (~/scw/hitemp_opacity.py) --
WAVENUM_MIN_CM = 10.0
WAVENUM_MAX_CM = 20000.0  # > 6000 so the >band Rosseland truncation is quantifiable
WSTEP_CM = 0.01  # fine; Rosseland is sensitive to unresolved windows
ISOTOPES = "1,2,3"  # H2(16O), H2(18O), H2(17O)
BROADENING = "voigt"  # self-broadened Voigt (per crossover_hitemp.json meta)
# Local, offline databank identifiers. Set these to the databank names you have
# registered in ~/.radis (or adapt build_factory to load_databank(path=...)).
DATABANK = {
    "hitemp": "HITEMP-H2O",  # HITEMP2010 water, local cache
    "pokazatel": "POKAZATEL-H2O",  # ExoMol POKAZATEL water, local cache
}


def broadening_pressure_pa(i_rho: int, i_t: int) -> float:
    rho = RHO_GRID_KGM3[i_rho]
    temperature = T_GRID_K[i_t]
    if PRESSURE_MODE == "aqua":
        return P_AQUA_PA[i_rho][i_t]
    return rho * R_GAS * temperature / M_H2O_KG_MOL  # ideal-gas P at real density


def build_factory(database: str):
    """Create a RADIS SpectrumFactory and load the databank ONCE (offline).

    CONFIG: adapt the load to match your working HitempCalculator. The
    SpectrumFactory pressure/Tgas are set per point in compute_point().
    """
    from radis import SpectrumFactory

    sf = SpectrumFactory(
        wmin=WAVENUM_MIN_CM,
        wmax=WAVENUM_MAX_CM,
        wstep=WSTEP_CM,
        molecule="H2O",
        isotope=ISOTOPES,
        mole_fraction=1.0,  # pure water; broadening is self-broadening
        pressure=1.0,  # bar, overridden per point
        cutoff=0.0,
        verbose=False,
    )
    if database == "hitemp":
        sf.load_databank(DATABANK["hitemp"])  # local HITEMP2010 cache, offline
    elif database == "pokazatel":
        sf.load_databank(DATABANK["pokazatel"])  # local ExoMol POKAZATEL cache
    else:
        raise ValueError(f"unknown database {database!r}")
    return sf


def compute_point(sf, temperature_k: float, pressure_pa: float):
    """Return (wn[cm^-1], kappa[cm^2/g]) for one (T, P_broad) point.

    kappa = sigma * N_A / M, with sigma = abscoeff / n and n = P/(k_B T). This is
    density-consistent for any broadening pressure (see module docstring).
    """
    sf.pressure = pressure_pa / 1e5  # RADIS wants bar
    s = sf.eq_spectrum(Tgas=temperature_k, name=f"T{temperature_k}")
    wn, abscoeff = s.get("abscoeff", wunit="cm-1", Iunit="cm-1")  # abscoeff [cm^-1]
    n_cm3 = pressure_pa / (K_B * temperature_k) * 1e-6  # molecules / cm^3
    sigma_cm2 = np.asarray(abscoeff, float) / n_cm3  # cm^2 / molecule
    kappa_cm2_g = sigma_cm2 * N_A / M_H2O_G_MOL  # cm^2 / g
    return np.asarray(wn, float), kappa_cm2_g


def out_path(out_dir: Path, temperature_k: float, rho_kgm3: float) -> Path:
    return out_dir / f"T{temperature_k:g}_rho{rho_kgm3:g}_opacity.npz"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--database", required=True, choices=["hitemp", "pokazatel"])
    ap.add_argument("--out", required=True, type=Path, help="output directory")
    ap.add_argument("--pressure-mode", choices=["ideal", "aqua"], default=None,
                    help="override the module PRESSURE_MODE")
    ap.add_argument("--dry-run", action="store_true",
                    help="list points and pressures without running RADIS")
    args = ap.parse_args(argv)

    global PRESSURE_MODE
    if args.pressure_mode:
        PRESSURE_MODE = args.pressure_mode

    out_dir: Path = args.out
    out_dir.mkdir(parents=True, exist_ok=True)
    points = [(it, irho) for irho in range(len(RHO_GRID_KGM3)) for it in range(len(T_GRID_K))]
    total = len(points)
    print(f"[generate] database={args.database} out={out_dir} pressure_mode={PRESSURE_MODE}", flush=True)
    print(f"[generate] {total} grid points, wn=[{WAVENUM_MIN_CM},{WAVENUM_MAX_CM}] wstep={WSTEP_CM}", flush=True)

    if args.dry_run:
        for k, (it, irho) in enumerate(points, 1):
            t_k, rho = T_GRID_K[it], RHO_GRID_KGM3[irho]
            p = broadening_pressure_pa(irho, it)
            exists = out_path(out_dir, t_k, rho).exists()
            print(f"  [{k:3d}/{total}] T={t_k:5g} rho={rho:9.3g}  P_broad={p/1e5:11.4g} bar"
                  f"  {'(exists, skip)' if exists else ''}")
        return 0

    sf = build_factory(args.database)  # databank loaded once

    done = skipped = failed = 0
    for k, (it, irho) in enumerate(points, 1):
        t_k, rho = float(T_GRID_K[it]), float(RHO_GRID_KGM3[irho])
        path = out_path(out_dir, t_k, rho)
        if path.exists():  # resumable
            skipped += 1
            print(f"  [{k:3d}/{total}] T={t_k:5g} rho={rho:9.3g}  skip (exists)", flush=True)
            continue
        p_pa = broadening_pressure_pa(irho, it)
        t0 = time.time()
        try:
            wn, kappa = compute_point(sf, t_k, p_pa)
        except Exception as exc:  # keep going; the run is resumable
            failed += 1
            print(f"  [{k:3d}/{total}] T={t_k:5g} rho={rho:9.3g}  FAILED: {exc}", flush=True)
            continue
        rho_g_cm3 = rho / 1000.0
        tmp = path.with_suffix(".npz.tmp")  # atomic write
        np.savez_compressed(
            tmp,
            wn=wn,
            kappa=kappa,
            alpha=kappa * rho_g_cm3,  # cm^-1
            T=t_k,
            rho_kg_m3=rho,
            database=args.database,
            pressure_mode=PRESSURE_MODE,
            P_broadening_pa=p_pa,
            wstep_cm=WSTEP_CM,
        )
        tmp.replace(path)
        done += 1
        print(f"  [{k:3d}/{total}] T={t_k:5g} rho={rho:9.3g}  P={p_pa/1e5:.4g} bar"
              f"  kP~{np.nanmax(kappa):.3g} cm2/g  {time.time()-t0:.1f}s", flush=True)

    print(f"[generate] done={done} skipped={skipped} failed={failed} / {total}", flush=True)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
