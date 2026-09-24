#!/usr/bin/env python3
"""Generate line-list water opacity tables (HITEMP / ExoMol-POKAZATEL) with RADIS.

Standalone cluster driver for `subneptune` (see docs/linelist_tables.md). Writes
one npz per (T, rho) grid point in the subneptune schema, resumable, offline.

  python generate_linelist_tables.py --database hitemp --out ~/subneptune_data/hitemp_opacity
  python generate_linelist_tables.py --database pokazatel --out ~/subneptune_data/pokazatel_opacity

Assumes a fat compute node with the line-list databank cached locally and NO
internet (no downloads). The (T, rho) grid and the AQUA broadening pressures are
baked in below.

Opacity convention (matches ~/scw/hitemp_opacity.py HitempCalculator + the
delivered POKAZATEL files; correct units, NO x1000):
  kappa[cm^2/g] = abscoeff[cm^-1] / rho_g_cm3,  rho_g_cm3 = n_cm3 * M / N_A,
  n_cm3 = P_broad/(k_B T).  This is identically sigma * N_A / M with
  sigma = abscoeff/n, so mass opacity is density-consistent for ANY broadening
  pressure -- the broadening-pressure choice below only affects Voigt widths.

Broadening pressure (PRESSURE_MODE):
  "ideal" (default) -> P = rho R T / M   (ideal-gas pressure at the real number
      density; self-broadening width ~ collision rate ~ real density -- the
      physical impact-broadening pressure). The delivered POKAZATEL p_bar equals
      this to <0.1% at every grid point (AMENDMENTS A17), so "ideal" reproduces
      the production pipeline. Equals AQUA P in the dilute regime.
  "aqua"            -> the tabulated AQUA P(T, rho) (baked in). Differs from
      "ideal" only in the non-ideal / two-phase regime, where the line list is
      superseded by FP anyway.

The RADIS block below is locked to the working HitempCalculator (SpectrumFactory
+ fetch_databank, DEFAULT_DOWNLOAD_PATH=$DATA/radisdb, isotopes 1,2,3 for HITEMP
/ 1 for POKAZATEL, self-broadened; POKAZATEL via fetch_databank("exomol",
database="POKAZATEL") with a multi-signature fallback across RADIS versions).
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
# 2250 and 2750 were added by A60: the FP grid has points at those temperatures
# and interpolating HITEMP across the 2000-2500 and 2500-3000 gaps was worth
# 20-30% there (A59).
T_GRID_K = [500, 600, 700, 800, 900, 1000, 1250, 1500, 2000, 2250, 2500, 2750, 3000]
RHO_GRID_KGM3 = [1e-5, 1e-4, 1e-3, 1e-2, 0.1, 1.0, 10.0, 30.0, 100.0, 300.0, 1000.0]

# AQUA P(T, rho) in Pa; rows = RHO_GRID_KGM3, cols = T_GRID_K.
P_AQUA_PA = [
    [2.3077e+00, 2.7692e+00, 3.2307e+00, 3.6922e+00, 4.1538e+00, 4.6155e+00, 5.7706e+00, 6.9407e+00, 1.0196e+01, 1.5124e+01, 2.6016e+01, 3.5978e+01, 4.1074e+01],
    [2.3077e+01, 2.7692e+01, 3.2307e+01, 3.6923e+01, 4.1538e+01, 4.6155e+01, 5.7700e+01, 6.9313e+01, 9.5833e+01, 1.1976e+02, 1.7236e+02, 2.7180e+02, 3.7486e+02],
    [2.3076e+02, 2.7692e+02, 3.2307e+02, 3.6922e+02, 4.1538e+02, 4.6155e+02, 5.7697e+02, 6.9270e+02, 9.3726e+02, 1.0953e+03, 1.3468e+03, 1.8184e+03, 2.6358e+03],
    [2.3074e+03, 2.7690e+03, 3.2306e+03, 3.6922e+03, 4.1537e+03, 4.6154e+03, 5.7695e+03, 6.9249e+03, 9.2918e+03, 1.0609e+04, 1.2230e+04, 1.4567e+04, 1.8341e+04],
    [2.3055e+04, 2.7676e+04, 3.2295e+04, 3.6913e+04, 4.1531e+04, 4.6152e+04, 5.7696e+04, 6.9242e+04, 9.2585e+04, 1.0480e+05, 1.1812e+05, 1.3378e+05, 1.5407e+05],
    [2.2858e+05, 2.7538e+05, 3.2191e+05, 3.6832e+05, 4.1468e+05, 4.6125e+05, 5.7694e+05, 6.9237e+05, 9.2435e+05, 1.0427e+06, 1.1654e+06, 1.2967e+06, 1.4435e+06],
    [2.0778e+06, 2.6139e+06, 3.1149e+06, 3.6028e+06, 4.0840e+06, 4.5863e+06, 5.7695e+06, 6.9236e+06, 9.2369e+06, 1.0404e+07, 1.1590e+07, 1.2811e+07, 1.4088e+07],
    [2.6416e+06, 6.9092e+06, 8.6625e+06, 1.0284e+07, 1.1843e+07, 1.3588e+07, 1.7309e+07, 2.0771e+07, 2.7706e+07, 3.1195e+07, 3.4723e+07, 3.8318e+07, 4.2025e+07],
    [2.6416e+06, 1.2354e+07, 2.1958e+07, 2.8898e+07, 3.5282e+07, 4.3557e+07, 5.7695e+07, 6.9234e+07, 9.2338e+07, 1.0394e+08, 1.1562e+08, 1.2745e+08, 1.3952e+08],
    [2.6416e+06, 1.2354e+07, 3.5817e+07, 6.2126e+07, 8.7856e+07, 1.1633e+08, 1.7247e+08, 2.2279e+08, 3.1794e+08, 3.6307e+08, 4.0698e+08, 4.4999e+08, 4.9238e+08],
    [3.5614e+08, 5.7727e+08, 7.9708e+08, 1.0103e+09, 1.2155e+09, 1.4118e+09, 1.8652e+09, 2.2773e+09, 3.0120e+09, 3.3466e+09, 3.6646e+09, 3.9689e+09, 4.2617e+09],
]

PRESSURE_MODE = "ideal"  # "ideal" | "aqua"

# --- RADIS setup (locked to ~/scw/hitemp_opacity.py HitempCalculator) --------
WAVENUM_MIN_CM = 10.0
WAVENUM_MAX_CM = 11000.0  # HitempCalculator range; > 6000 lets us quantify the
#                           FP band-truncation vs the line list where they overlap
WSTEP_CM = 0.01  # HitempCalculator default; Rosseland is sensitive to windows
# Isotopes per database: HITEMP uses 1,2,3; the POKAZATEL production used 1.
ISOTOPES = {"hitemp": "1,2,3", "pokazatel": "1"}
# RADIS caches the line databases under $DATA/radisdb on the cluster.
RADISDB_ENV = "DATA"  # radis.config["DEFAULT_DOWNLOAD_PATH"] = $DATA/radisdb


def broadening_pressure_pa(i_rho: int, i_t: int) -> float:
    rho = RHO_GRID_KGM3[i_rho]
    temperature = T_GRID_K[i_t]
    if PRESSURE_MODE == "aqua":
        return P_AQUA_PA[i_rho][i_t]
    return rho * R_GAS * temperature / M_H2O_KG_MOL  # ideal-gas P at real density


def _fetch_databank(sf, database: str) -> None:
    """Load the line database ONCE (offline, from the local cache).

    Locked to HitempCalculator: HITEMP via fetch_databank("hitemp"); POKAZATEL
    via fetch_databank("exomol", database="POKAZATEL") with a multi-signature
    fallback (RADIS has moved this kwarg across versions).
    """
    if database == "hitemp":
        sf.fetch_databank("hitemp")
        return
    # POKAZATEL: try the signatures the working script cycles through.
    attempts = (
        lambda: sf.fetch_databank("exomol", database="POKAZATEL"),
        lambda: sf.fetch_databank("exomol", databank="POKAZATEL"),
        lambda: sf.fetch_databank("exomol", name="POKAZATEL"),
    )
    errors = []
    for attempt in attempts:
        try:
            attempt()
            return
        except Exception as exc:  # noqa: BLE001 - try the next signature
            errors.append(repr(exc))
    raise RuntimeError("could not load POKAZATEL databank; tried:\n  " + "\n  ".join(errors))


def build_factory(database: str, wavenum_min_cm: float = WAVENUM_MIN_CM):
    """Create a RADIS SpectrumFactory and load the databank ONCE (offline).

    Locked to HitempCalculator's SpectrumFactory signature; pressure/Tgas are
    set per point in compute_point().
    """
    import os

    import radis
    from radis import SpectrumFactory

    data_dir = os.environ.get(RADISDB_ENV)
    if data_dir:
        radis.config["DEFAULT_DOWNLOAD_PATH"] = os.path.join(data_dir, "radisdb")
    else:
        print(f"[generate] WARNING: ${RADISDB_ENV} unset; using RADIS default cache path",
              flush=True)

    skip_radis_registration()
    sf = SpectrumFactory(
        wavenum_min=wavenum_min_cm,
        wavenum_max=WAVENUM_MAX_CM,
        molecule="H2O",
        isotope=ISOTOPES[database],
        wstep=WSTEP_CM,
        pressure=1.0,  # bar; overridden per point in compute_point
        mole_fraction=1.0,  # pure water -> self-broadening
        verbose=False,
    )
    _fetch_databank(sf, database)
    return sf


def skip_radis_registration() -> None:
    """Let a fully-cached, offline databank load succeed (A59).

    ``fetch_hitemp`` calls ``HITEMPDatabaseManager.register()`` when the databank
    is absent from ``~/radis.json``, and that method formats ``self.wmin``, which
    is only set when files are actually downloaded. Registration is bookkeeping;
    the line data load regardless. Skipped when nothing was downloaded, so the
    user's home config is left alone.
    """
    try:
        from radis.api.hitempapi import HITEMPDatabaseManager as _Mgr
    except ImportError:  # older RADIS: no shim needed
        return
    if getattr(_Mgr.register, "_subneptune_shim", False):
        return
    original = _Mgr.register

    def register(self, download):
        if not download and (self.wmin is None or self.wmax is None):
            return
        return original(self, download)

    register._subneptune_shim = True
    _Mgr.register = register


def compute_point(sf, temperature_k: float, pressure_pa: float):
    """Return (wn[cm^-1], kappa[cm^2/g]) for one (T, P_broad) point.

    Locked to HitempCalculator.kappa: kappa = abscoeff / rho_g_cm3 with
    rho_g_cm3 = n_cm3 * M / N_A, n_cm3 = P/(k_B T). Identically sigma*N_A/M.
    """
    # RADIS >= 0.16 made SpectrumFactory slotted: `sf.pressure = ...` raises, and
    # the pressure is a per-call argument. The kwarg form works on older RADIS too,
    # so no version branch (A59).
    s = sf.eq_spectrum(Tgas=float(temperature_k), pressure=pressure_pa / 1e5)
    wn, abscoeff = s.get("abscoeff", wunit="cm-1")  # abscoeff [cm^-1]
    n_cm3 = (pressure_pa / (K_B * temperature_k)) * 1e-6  # molecules / cm^3
    rho_g_cm3 = n_cm3 * M_H2O_G_MOL / N_A
    kappa_cm2_g = np.asarray(abscoeff, float) / rho_g_cm3  # cm^2 / g
    return np.asarray(wn, float), kappa_cm2_g


def out_path(out_dir: Path, temperature_k: float, rho_kgm3: float) -> Path:
    # Matches the delivered POKAZATEL convention (registry glob T*_rho*.npz).
    return out_dir / f"T{temperature_k:g}_rho{rho_kgm3:g}.npz"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--database", required=True, choices=["hitemp", "pokazatel"])
    ap.add_argument("--out", required=True, type=Path, help="output directory")
    ap.add_argument("--pressure-mode", choices=["ideal", "aqua"], default=None,
                    help="override the module PRESSURE_MODE")
    ap.add_argument("--temperatures", default=None,
                    help="comma-separated subset of T_GRID_K to compute (default: all)")
    ap.add_argument("--densities", default=None,
                    help="comma-separated subset of RHO_GRID_KGM3 to compute. The "
                         "delivered grid covers only 0.1-1000; the four thinner rows "
                         "were never produced (default: all)")
    ap.add_argument("--wavenum-min", type=float, default=WAVENUM_MIN_CM,
                    help=f"line-loading lower edge, cm^-1 (default {WAVENUM_MIN_CM:g}). "
                         f"The local HITEMP cache starts at 50; with truncation 50 a "
                         f"line below 50 cannot reach past 100, so the science band "
                         f"[100, 4000] is unaffected either way (A60)")
    ap.add_argument("--dry-run", action="store_true",
                    help="list points and pressures without running RADIS")
    args = ap.parse_args(argv)

    global PRESSURE_MODE
    if args.pressure_mode:
        PRESSURE_MODE = args.pressure_mode

    out_dir: Path = args.out
    out_dir.mkdir(parents=True, exist_ok=True)
    wanted = None
    if args.temperatures:
        wanted = {float(x) for x in args.temperatures.split(",")}
        unknown = wanted - set(map(float, T_GRID_K))
        if unknown:
            raise SystemExit(f"--temperatures {sorted(unknown)} not in T_GRID_K={T_GRID_K}")
    wanted_rho = None
    if args.densities:
        wanted_rho = {float(x) for x in args.densities.split(",")}
        unknown = wanted_rho - set(map(float, RHO_GRID_KGM3))
        if unknown:
            raise SystemExit(f"--densities {sorted(unknown)} not in RHO_GRID_KGM3")
    points = [(it, irho) for irho in range(len(RHO_GRID_KGM3)) for it in range(len(T_GRID_K))
              if (wanted is None or float(T_GRID_K[it]) in wanted)
              and (wanted_rho is None or float(RHO_GRID_KGM3[irho]) in wanted_rho)]
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

    sf = build_factory(args.database, args.wavenum_min)  # databank loaded once

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
        # Atomic write: tmp must end in .npz or np.savez appends another .npz.
        tmp = path.with_name(path.stem + ".tmp.npz")
        np.savez_compressed(  # schema-consistent with the delivered POKAZATEL files
            tmp,
            wn=wn,
            kappa=kappa,  # cm^2/g
            alpha=kappa * rho_g_cm3,  # cm^-1
            T=t_k,
            rho=rho,
            rho_kg_m3=rho,
            p_bar=p_pa / 1e5,
            source=f"RADIS/{args.database}",
            isotopes=ISOTOPES[args.database],
            pressure_mode=PRESSURE_MODE,
            wstep_cm=WSTEP_CM,
            wavenum_min_cm=args.wavenum_min,
        )
        tmp.replace(path)
        done += 1
        print(f"  [{k:3d}/{total}] T={t_k:5g} rho={rho:9.3g}  P={p_pa/1e5:.4g} bar"
              f"  kP~{np.nanmax(kappa):.3g} cm2/g  {time.time()-t0:.1f}s", flush=True)

    print(f"[generate] done={done} skipped={skipped} failed={failed} / {total}", flush=True)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
