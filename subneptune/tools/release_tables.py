#!/usr/bin/env python3
"""Drop-in opacity table and power-law fit for Data availability (Task 41).

    python tools/release_tables.py config/fiducial.yaml

Writes, to ``data/release/``:

  ``subneptune_kappa_R_mesa.txt``   first-principles Rosseland mean on a regular
  ``subneptune_kappa_P_mesa.txt``   (log R, log T) grid in the plain-text format
                                    MESA's kap module reads;
  ``subneptune_kappa_grid.csv``     the same numbers with true blanks outside the
                                    hull, for readers who are not running MESA;
  ``subneptune_kappa_fit.json``     the two-regime power-law fit.

**The format is taken from the reader, not from memory.** ``kap/private/
load_kap.f90`` in MESA (``Read_Kap_X_Table``, and the header parse at lines
266-285) reads:

  line 1, 2      skipped
  line 3         form version X Z num_logRs logR_min logR_max num_logTs
                 logT_min logT_max      (>= 10 whitespace-separated numbers)
  then per table block:
  line 1, 2      skipped
  line 3         the num_logRs values of log R
  line 4         skipped
  next num_logTs lines   log T followed by num_logRs values of log10 kappa

with ``form = kap_table_fixed_metal_form = 1`` (``kap_def.f90``) and
``version >= min_version = 37`` (``load_co_kap.f90``). MESA's R is

    R = rho / T6^3,   rho in **g cm^-3**,   T6 = T / 10^6 K,

and kappa is cm^2 g^-1, both as log10.

**This is the format, not a drop-in replacement for a MESA opacity table.** On
that definition a steam envelope at 10-2000 kg m^-3 and 650-3000 K sits at
log R = 7-9, where MESA's own tables run -8 to +1: no MESA table covers this
regime, which is the point of releasing one. The table is written so that a MESA
reader can parse it; using it requires extending the kap module's R range.

**Nothing is extrapolated.** Cells outside the convex hull of the 52 computed
points carry the sentinel -99.0 in the MESA file (the OPAL "no data" value; MESA
has no blank convention because its parser requires a number in every column)
and are empty in the CSV.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from datetime import date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np

from config import Config, load_config

BAND_CM = (100.0, 4000.0)
DLOG_R = 0.25
DLOG_T = 0.05
MESA_FORM = 1
MESA_VERSION = 37
NO_DATA = -99.0
DENSE_RHO_MIN = 100.0
RHO_0_KG_M3 = 100.0   # round reference values, stated rather than fitted
T_0_K = 1000.0

CAVEATS = (
    ("Interaction-induced absorption",
     "These are first-principles molecular-dynamics opacities, so "
     "collision-induced and many-body absorption is present by construction and "
     "is not separable into a line list plus a pair continuum."),
    ("Dipole model",
     "The spectra come from a machine-learned dipole surface whose "
     "symmetric-stretch derivative is 2.8x the revPBE0 value at 0.957 A and 4.0x "
     "at 1.00 A, so the 3200-4000 cm^-1 stretch band is too strong and kappa_P is "
     "affected more than kappa_R."),
)


def log_r(rho_kg_m3, t_k):
    """MESA's log R = log10(rho[g/cm^3] / (T/1e6 K)^3)."""
    return np.log10((np.asarray(rho_kg_m3) / 1000.0)
                    / (np.asarray(t_k) / 1.0e6) ** 3)


def _grid(cfg: Config):
    from opacity.build_tables import compute_point_means

    # compute_point_means returns (T, rho, kappa_R, kappa_P) — Rosseland first.
    t_pt, rho_pt, kr, kp = compute_point_means("fp", cfg, BAND_CM)
    t_pt = np.asarray(t_pt, float)
    rho_pt = np.asarray(rho_pt, float)
    lr_pt = log_r(rho_pt, t_pt)
    lt_pt = np.log10(t_pt)

    lt_lo = np.floor(lt_pt.min() / DLOG_T) * DLOG_T
    lt_hi = np.ceil(lt_pt.max() / DLOG_T) * DLOG_T
    lr_lo = np.floor(lr_pt.min() / DLOG_R) * DLOG_R
    lr_hi = np.ceil(lr_pt.max() / DLOG_R) * DLOG_R
    log_ts = np.round(np.arange(lt_lo, lt_hi + 0.5 * DLOG_T, DLOG_T), 10)
    log_rs = np.round(np.arange(lr_lo, lr_hi + 0.5 * DLOG_R, DLOG_R), 10)
    return t_pt, rho_pt, np.asarray(kr, float), np.asarray(kp, float), \
        lr_pt, lt_pt, log_rs, log_ts


def _interpolate(t_pt, rho_pt, values, log_rs, log_ts):
    """log10 kappa on the (log R, log T) mesh; NaN outside the point hull.

    Interpolated in (T, log10 rho) over log10 kappa with ``rescale=True`` — the
    same interpolant, in the same space, that the structure model's opacity table
    uses — and then sampled at the mesh. Fitting a second interpolant in
    (log R, log T) would release numbers the paper's own model does not produce.
    ``LinearNDInterpolator`` returns NaN outside the convex hull of the input
    points, which is the "blank, never extrapolated" rule; the nearest-neighbour
    clamp the model applies outside the hull is deliberately NOT used here.
    """
    from scipy.interpolate import LinearNDInterpolator

    f = LinearNDInterpolator(np.column_stack([t_pt, np.log10(rho_pt)]),
                             np.log10(values), rescale=True)
    RR, TT = np.meshgrid(log_rs, log_ts, indexing="ij")
    t_k = 10.0 ** TT
    rho = 10.0 ** RR * (t_k / 1.0e6) ** 3 * 1000.0
    return f(np.column_stack([t_k.ravel(), np.log10(rho).ravel()])).reshape(RR.shape)


def write_mesa(path: Path, log_rs, log_ts, log_kappa, label: str, cfg_note: str):
    """log_kappa indexed [logR, logT]; NaN becomes the no-data sentinel."""
    k = np.where(np.isfinite(log_kappa), log_kappa, NO_DATA)
    # MESA's reader skips exactly two lines before the numeric header, so the
    # full statement of units, extent and caveats lives in README.md beside this
    # file; these two lines have to carry the essentials and no more.
    lines = [
        f"# {label}, pure H2O, log10 kappa [cm^2/g] vs log R = "
        f"log10(rho[g/cm^3]/(T/1e6 K)^3) and log10 T[K]; band "
        f"{BAND_CM[0]:g}-{BAND_CM[1]:g} cm^-1; {cfg_note}; see README.md",
        "#  form  version      X      Z  num_logR   logR_min   logR_max"
        "  num_logT   logT_min   logT_max",
        f"{MESA_FORM:6d}{MESA_VERSION:9d}{0.0:7.1f}{1.0:7.1f}"
        f"{log_rs.size:10d}{log_rs[0]:11.4f}{log_rs[-1]:11.4f}"
        f"{log_ts.size:10d}{log_ts[0]:11.4f}{log_ts[-1]:11.4f}",
        "",
        f"# pure H2O (X = 0, Z = 1); log10 kappa in cm^2/g; "
        f"{NO_DATA:g} = outside the computed hull, not extrapolated",
        "                 " + "".join(f"{v:9.4f}" for v in log_rs),
        "",
    ]
    for j, lt in enumerate(log_ts):
        lines.append(f"{lt:9.4f}   " + "".join(f"{v:9.4f}" for v in k[:, j]))
    path.write_text("\n".join(lines) + "\n")


def read_mesa(path: Path):
    """A strict reader of the same spec, for the round trip."""
    raw = [ln for ln in path.read_text().splitlines()]
    head = [float(x) for x in raw[2].split()]
    form, version = int(head[0]), int(head[1])
    n_r, n_t = int(head[4]), int(head[7])
    log_rs = np.array([float(x) for x in raw[5].split()])
    if log_rs.size != n_r:
        raise ValueError(f"{path.name}: header says {n_r} logR, row has {log_rs.size}")
    log_ts, k = [], []
    for ln in raw[7:7 + n_t]:
        v = [float(x) for x in ln.split()]
        if len(v) != 1 + n_r:
            raise ValueError(f"{path.name}: row has {len(v)} values, want {1 + n_r}")
        log_ts.append(v[0])
        k.append(v[1:])
    if len(log_ts) != n_t:
        raise ValueError(f"{path.name}: header says {n_t} logT, file has {len(log_ts)}")
    return {"form": form, "version": version, "log_rs": log_rs,
            "log_ts": np.array(log_ts), "log_kappa": np.array(k).T,
            "header": head}


def fit_power_law(t_pt, rho_pt, kr):
    """kappa_R = kappa_0 (rho/rho_0)^a (T/T_0)^b, per regime, by least squares."""
    out = {}
    for name, mask in (("dense", rho_pt >= DENSE_RHO_MIN),
                       ("dilute", rho_pt < DENSE_RHO_MIN)):
        n = int(mask.sum())
        if n < 4:
            out[name] = {"n_points": n, "blocked": "fewer than four points"}
            continue
        x = np.log(rho_pt[mask] / RHO_0_KG_M3)
        y = np.log(t_pt[mask] / T_0_K)
        z = np.log(kr[mask])
        A = np.column_stack([np.ones_like(x), x, y])
        coef, *_ = np.linalg.lstsq(A, z, rcond=None)
        resid = z - A @ coef
        out[name] = {
            "n_points": n,
            "kappa_0_cm2_g": float(np.exp(coef[0])),
            "a_rho_exponent": float(coef[1]),
            "b_T_exponent": float(coef[2]),
            "rho_0_kg_m3": RHO_0_KG_M3, "T_0_K": T_0_K,
            "rms_log10_residual": float(np.sqrt((resid ** 2).mean()) / np.log(10.0)),
            "max_abs_log10_residual": float(np.abs(resid).max() / np.log(10.0)),
            "validity_rho_kg_m3": [float(rho_pt[mask].min()), float(rho_pt[mask].max())],
            "validity_T_K": [float(t_pt[mask].min()), float(t_pt[mask].max())],
        }
    return out


def _readme(log_rs, log_ts, t_pt, rho_pt, k_r, note: str) -> str:
    return f"""# First-principles steam opacity — release tables

Rosseland and Planck mean opacity of pure H2O from first-principles molecular
dynamics, on a regular (log R, log T) grid.

| | |
|---|---|
| Units | kappa in cm^2 g^-1, tabulated as log10; T in K; rho in kg m^-3 (MESA's log R uses g cm^-3) |
| Band | {BAND_CM[0]:g}-{BAND_CM[1]:g} cm^-1, the widest band every source in the study shares |
| log R | {log_rs[0]:.2f} to {log_rs[-1]:.2f}, step {DLOG_R:g} dex, R = rho / (T/10^6 K)^3 with rho in g cm^-3 |
| log T | {log_ts[0]:.2f} to {log_ts[-1]:.2f}, step {DLOG_T:g} dex ({10 ** log_ts[0]:.0f}-{10 ** log_ts[-1]:.0f} K) |
| Computed points | {t_pt.size} states, T {t_pt.min():.0f}-{t_pt.max():.0f} K, rho {rho_pt.min():.0f}-{rho_pt.max():.0f} kg m^-3 |
| Coverage | {note} |
| Outside the hull | `{NO_DATA:g}` in the MESA files, blank in the CSV — never extrapolated |

## Files

- `subneptune_kappa_R_mesa.txt`, `subneptune_kappa_P_mesa.txt` — the plain-text
  format MESA's `kap` module reads (`kap/private/load_kap.f90`,
  `form = 1`, `version = {MESA_VERSION}`). **Format, not a drop-in table**: a steam
  envelope sits at log R {log_rs[0]:.1f}-{log_rs[-1]:.1f}, where MESA's own tables run -8 to +1.
- `subneptune_kappa_grid.csv` — the same numbers with true blanks, plus T and rho.
- `subneptune_kappa_fit.json` — the two-regime power-law fit and this metadata.

Values are the study's own opacity table sampled onto the regular grid: linear
on the Delaunay triangulation in (T, log10 rho) over log10 kappa.

## Two caveats, one line each

- **{CAVEATS[0][0]}.** {CAVEATS[0][1]}
- **{CAVEATS[1][0]}.** {CAVEATS[1][1]}
"""


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("config", type=Path)
    args = ap.parse_args(argv)
    cfg = load_config(args.config)
    # A87: redirectable for the same reason as FIGS_DIR -- this writes outside
    # cfg.output_path, and the A86 v3 pass clobbered the tracked v5 products.
    from layout import RELEASE_DIR as root
    root.mkdir(parents=True, exist_ok=True)

    t_pt, rho_pt, kr, kp, lr_pt, lt_pt, log_rs, log_ts = _grid(cfg)
    k_r = _interpolate(t_pt, rho_pt, kr, log_rs, log_ts)
    k_p = _interpolate(t_pt, rho_pt, kp, log_rs, log_ts)

    note = (f"grid {DLOG_R:g} dex in log R, {DLOG_T:g} dex in log T; "
            f"{int(np.isfinite(k_r).sum())} of {k_r.size} cells inside the "
            f"{t_pt.size}-point hull")
    paths = {}
    for label, arr, stem in (("Rosseland mean", k_r, "subneptune_kappa_R_mesa"),
                             ("Planck mean", k_p, "subneptune_kappa_P_mesa")):
        p = root / f"{stem}.txt"
        write_mesa(p, log_rs, log_ts, arr, label, note)
        back = read_mesa(p)
        got = np.where(back["log_kappa"] == NO_DATA, np.nan, back["log_kappa"])
        ok = np.isfinite(arr)
        err = float(np.max(np.abs(got[ok] - arr[ok]))) if ok.any() else 0.0
        if not np.array_equal(np.isfinite(got), ok):
            raise SystemExit(f"{p.name}: round trip lost the no-data mask")
        # the file carries four decimal places, MESA's own precision, so a round
        # trip agrees to half a unit in the last place and no better
        if err > 1.0e-4:
            raise SystemExit(f"{p.name}: round trip differs by {err:.2e}")
        paths[stem] = {"path": str(p.relative_to(root.parents[1])),
                       "round_trip_max_abs_dex": err,
                       "n_finite": int(ok.sum()), "n_cells": int(arr.size),
                       "written_decimals": 4,
                       "interpolant": ("linear on the Delaunay triangulation in "
                                       "(T, log10 rho) over log10 kappa, the "
                                       "model's own")}

    csv_path = root / "subneptune_kappa_grid.csv"
    with csv_path.open("w") as fh:
        fh.write("# first-principles H2O opacity, band "
                 f"{BAND_CM[0]:g}-{BAND_CM[1]:g} cm^-1; log R = log10(rho[g/cm^3] "
                 "/ (T/1e6 K)^3); blank = outside the computed hull\n")
        fh.write("log_R,log_T,T_K,rho_kg_m3,log10_kappa_R_cm2_g,"
                 "log10_kappa_P_cm2_g\n")
        for i, lr in enumerate(log_rs):
            for j, lt in enumerate(log_ts):
                t_k = 10.0 ** lt
                rho = 10.0 ** lr * (t_k / 1.0e6) ** 3 * 1000.0
                a = "" if not np.isfinite(k_r[i, j]) else f"{k_r[i, j]:.6f}"
                b = "" if not np.isfinite(k_p[i, j]) else f"{k_p[i, j]:.6f}"
                fh.write(f"{lr:.4f},{lt:.4f},{t_k:.4f},{rho:.6g},{a},{b}\n")

    readme = root / "README.md"
    readme.write_text(_readme(log_rs, log_ts, t_pt, rho_pt, k_r, note))

    fit = fit_power_law(t_pt, rho_pt, kr)
    doc = {
        "meta": __doc__.split("\n\n")[0],
        "generated": date.today().isoformat(),
        "band_cm": list(BAND_CM),
        "units": {"kappa": "cm^2/g", "rho": "kg/m^3 (MESA log R uses g/cm^3)",
                  "T": "K"},
        "grid": {"dlog_R": DLOG_R, "dlog_T": DLOG_T,
                 "log_R_range": [float(log_rs[0]), float(log_rs[-1])],
                 "log_T_range": [float(log_ts[0]), float(log_ts[-1])],
                 "n_logR": int(log_rs.size), "n_logT": int(log_ts.size),
                 "T_range_K": [float(10 ** log_ts[0]), float(10 ** log_ts[-1])],
                 "n_computed_points": int(t_pt.size),
                 "no_data_value": NO_DATA,
                 "outside_hull": "blank in the CSV, sentinel in the MESA files"},
        "mesa_format": {"source": "MESA kap/private/load_kap.f90 (Read_Kap_X_Table)",
                        "form": MESA_FORM, "version": MESA_VERSION,
                        "log_R_definition": "log10(rho[g/cm^3] / (T/1e6 K)^3)",
                        "note": ("this steam regime sits at log R 7-9; MESA's own "
                                 "tables run -8 to +1, so the files are in MESA's "
                                 "format but outside any MESA table's range")},
        "files": paths,
        "csv": str(csv_path.relative_to(root.parents[1])),
        "fit": fit,
        "caveats": [f"{a}: {b}" for a, b in CAVEATS],
    }
    (root / "subneptune_kappa_fit.json").write_text(json.dumps(doc, indent=1) + "\n")

    print(f"\n  RELEASE TABLES -> {root}")
    print(f"    grid: {log_rs.size} log R x {log_ts.size} log T, "
          f"log R {log_rs[0]:.2f} to {log_rs[-1]:.2f}, "
          f"log T {log_ts[0]:.2f} to {log_ts[-1]:.2f}")
    for stem, p in paths.items():
        print(f"    {stem}.txt  {p['n_finite']}/{p['n_cells']} cells inside the hull, "
              f"round trip {p['round_trip_max_abs_dex']:.1e} dex")
    print(f"    {csv_path.name}")
    print("\n  FIT  kappa_R = kappa_0 (rho/100)^a (T/1000)^b")
    for name, f in fit.items():
        if "blocked" in f:
            print(f"    {name:<7s} {f['blocked']}")
            continue
        print(f"    {name:<7s} kappa_0 = {f['kappa_0_cm2_g']:9.4g} cm^2/g   "
              f"a = {f['a_rho_exponent']:+.3f}   b = {f['b_T_exponent']:+.3f}   "
              f"RMS {f['rms_log10_residual']:.3f} dex (n = {f['n_points']})")
        print(f"            valid rho {f['validity_rho_kg_m3'][0]:g}-"
              f"{f['validity_rho_kg_m3'][1]:g} kg/m^3, "
              f"T {f['validity_T_K'][0]:g}-{f['validity_T_K'][1]:g} K")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
