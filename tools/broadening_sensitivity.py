#!/usr/bin/env python3
"""Broadening-convention sensitivity for the HITEMP spectra (Task 21, A59).

    python tools/broadening_sensitivity.py config/fiducial.yaml            # all
    python tools/broadening_sensitivity.py config/fiducial.yaml --only 1000,30
    python tools/broadening_sensitivity.py config/fiducial.yaml --dry-run  # P table

What is being tested, and a correction to the brief. The task was written on the
assumption that "the bare HITEMP spectra use a fixed 1 bar reference broadening".
They do not. Every delivered file carries ``pressure_mode = 'ideal'`` and a
``p_bar`` equal to ρRT/M at its own density (0.4615 bar at T = 1000 K,
ρ = 0.1 kg/m³; 4615 bar at ρ = 1000), and κ(ν) at fixed T differs by up to 143×
between ρ = 0.1 and ρ = 1000 — which cannot happen under a ρ-independent
broadening, since κ = α/ρ would then be identical. AMENDMENTS A17 recorded the
same thing for the delivered POKAZATEL set: the ``HitempCalculator`` *snippet*
used a 1 bar reference, the production run did not.

So the real convention question for this dataset is **ideal-gas P versus AQUA
P(T, ρ)**, which is what the structure model's EOS actually says the pressure is.
All three conventions are computed so the brief's literal question is answered
too:

  ``ref1bar``  — fixed 1 bar, the convention the brief assumed;
  ``ideal``    — P = ρRT/M, the convention the delivered spectra were built with
                 (so this run also *validates* the local pipeline against them);
  ``aqua``     — the tabulated AQUA P(T, ρ), the structure model's own EOS.

Mass opacity is convention-independent by construction: RADIS returns an
absorption coefficient at its own ideal number density n = P_broad/(k_B T), and
κ = α/(n·M/N_A) ≡ σ·N_A/M. The broadening pressure therefore moves *only* the
Voigt widths — which is exactly the systematic being quantified.

**The Lorentzian half-width cap.** RADIS truncates every line profile at
``truncation`` cm⁻¹ from line centre (default 50, used here). At these pressures
the Lorentzian HWHM runs far past that: self-broadened water is ~0.46 cm⁻¹/atm at
296 K, so at ρ = 2229 kg/m³ the nominal HWHM is O(10³) cm⁻¹.

**The truncation does not discard absorption — it redistributes it.** RADIS
renormalises each truncated profile to unit area
(``_calc_lorentzian_lineshape``: ``lineshape /= area``), so line strength is
conserved exactly. Verified on the delivered spectra: at T = 1000 K, ∫κ(ν)dν is
constant to 1 part in 10⁴ from ρ = 0.1 to 1000 kg/m³, i.e. across a 10⁴× range in
broadening pressure (0.46 → 4615 bar). What the truncation does instead is pile
the area that physically belongs beyond ±50 cm⁻¹ back into the window: where
HWHM ≫ 50 cm⁻¹ a line becomes, in effect, a ~100 cm⁻¹-wide rectangle carrying its
full strength. **The sign of that error on κ_R is not determined a priori** — it
depends on where the redistributed area lands relative to the transparency
windows — so the ratios reported here are *not* bounds in either direction. They
carry an unquantified line-shape systematic wherever HWHM > 50 cm⁻¹, which is 45
of the 48 FP grid points (A58). The cap is recorded per run in the output.

Writes ``outputs/broadening_sensitivity.json`` (+ per-point npz spectra). Run
``figures/plot_ed1_broadening.py`` afterwards for the Extended Data panel. If
RADIS or the line database is unavailable, see ``tools/broadening_sensitivity.sbatch``
for the ARC hand-off.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np

from config import Config, load_config
from constants import K_B, MMW_H2O, N_A
from opacity.means import planck_mean, rosseland_mean

R_GAS = 8.314462618  # J mol⁻¹ K⁻¹
M_H2O_G_MOL = MMW_H2O * 1.0e3

# The six points the brief names: four FP grid points spanning the grid, plus the
# cold-dense corner (650 K, 300) and a dilute one (1000 K, 30).
# ...ordered by ascending worst-case broadening pressure, so the cheap points
# land first and a failure in the extreme corner still leaves a usable run.
POINTS = (
    (1000.0, 30.0),
    (1000.0, 129.7),
    (650.0, 300.0),
    (1500.0, 471.3),
    (2000.0, 1200.0),
    (3000.0, 2228.7),
)

CONVENTIONS = ("ref1bar", "ideal", "aqua")

# Locked to generate_linelist_tables.py / the delivered files, so a run here is
# directly comparable with what is on disk.
# The delivered spectra were computed over [10, 11000] cm^-1. The local HITEMP
# cache starts at 50 cm^-1, and with truncation = 50 a line below 50 cm^-1 cannot
# reach past 100 cm^-1 — the band's lower edge — so the science band is
# unaffected. ``--wavenum-min`` overrides it on a machine with the 0-50 file
# (e.g. ARC), where 10.0 reproduces the delivered convention exactly.
WAVENUM_MIN_CM = 50.0
WAVENUM_MAX_CM = 11000.0
WSTEP_CM = 0.01
ISOTOPES = "1,2,3"
TRUNCATION_CM = 50.0  # RADIS default; the Lorentzian half-width cap (see above)


def aqua_pressure_pa(cfg: Config, temperature_k: float, rho_si: float) -> tuple[float, int]:
    """AQUA P(T, ρ) in Pa and the phase id, by bilinear interpolation in (log ρ, log T).

    The AQUA table is (ρ, T)-indexed with P as an output column (A2), so P(T, ρ)
    is a direct lookup — no inversion needed.
    """
    from eos.loader import load_aqua_tables

    tab = load_aqua_tables(str(cfg.aqua_table_path))
    lr, lt = np.log10(float(rho_si)), np.log10(float(temperature_k))
    i = int(np.clip(np.searchsorted(tab.logrho_grid, lr) - 1, 0, tab.n_rho - 2))
    j = int(np.clip(np.searchsorted(tab.logT_grid, lt) - 1, 0, tab.n_T - 2))
    fr = (lr - tab.logrho_grid[i]) / (tab.logrho_grid[i + 1] - tab.logrho_grid[i])
    ft = (lt - tab.logT_grid[j]) / (tab.logT_grid[j + 1] - tab.logT_grid[j])
    log_p = np.log10(tab.press[i : i + 2, j : j + 2])
    p = ((1 - fr) * (1 - ft) * log_p[0, 0] + fr * (1 - ft) * log_p[1, 0]
         + (1 - fr) * ft * log_p[0, 1] + fr * ft * log_p[1, 1])
    phase = int(tab.phase[i + int(fr > 0.5), j + int(ft > 0.5)])
    return float(10.0**p), phase


def ideal_pressure_pa(temperature_k: float, rho_si: float) -> float:
    return float(rho_si) * R_GAS * float(temperature_k) / MMW_H2O


def broadening_pressures(cfg: Config) -> list[dict]:
    """Per point, the three broadening pressures + the AQUA phase id."""
    out = []
    for t_k, rho in POINTS:
        p_aqua, phase = aqua_pressure_pa(cfg, t_k, rho)
        out.append({
            "T": t_k,
            "rho": rho,
            "p_ref1bar_pa": 1.0e5,
            "p_ideal_pa": ideal_pressure_pa(t_k, rho),
            "p_aqua_pa": p_aqua,
            "aqua_phase_id": phase,
        })
    return out


def _skip_radis_registration() -> None:
    """Don't let RADIS's bookkeeping step abort a fully-cached, offline load.

    ``fetch_hitemp`` calls ``HITEMPDatabaseManager.register()`` whenever the databank is not
    listed in ``~/radis.json``, and that method formats ``self.wmin``/``self.wmax``
    — which are only set when files are actually downloaded. On a machine whose
    cache is already complete it raises ``TypeError: unsupported format string
    passed to NoneType``. Registration is pure bookkeeping (it writes an entry to
    the user's home config); the line data are loaded afterwards regardless. So we
    skip it when nothing was downloaded, and leave the user's config untouched.
    """
    from radis.api.hitempapi import HITEMPDatabaseManager as _Mgr

    if getattr(_Mgr.register, "_subneptune_shim", False):
        return
    original = _Mgr.register

    def register(self, download):
        if not download and (self.wmin is None or self.wmax is None):
            return
        return original(self, download)

    register._subneptune_shim = True
    _Mgr.register = register


def build_factory(wavenum_min_cm: float = WAVENUM_MIN_CM):
    """RADIS SpectrumFactory, locked to the delivered-spectra conventions."""
    import radis
    from radis import SpectrumFactory

    import os

    data_dir = os.environ.get("DATA")
    if data_dir:
        radis.config["DEFAULT_DOWNLOAD_PATH"] = os.path.join(data_dir, "radisdb")
    _skip_radis_registration()
    sf = SpectrumFactory(
        wavenum_min=wavenum_min_cm,
        wavenum_max=WAVENUM_MAX_CM,
        molecule="H2O",
        isotope=ISOTOPES,
        wstep=WSTEP_CM,
        pressure=1.0,
        mole_fraction=1.0,  # pure water -> self-broadening
        truncation=TRUNCATION_CM,
        verbose=False,
    )
    sf.fetch_databank("hitemp")
    return sf


def compute_kappa(sf, temperature_k: float, pressure_pa: float):
    """(wn [cm⁻¹], κ [cm²/g]) at one (T, P_broad); κ = σ·N_A/M, convention-free."""
    # RADIS >= 0.16 made SpectrumFactory slotted: the pressure is a per-call
    # argument now, not a settable attribute (generate_linelist_tables.py was
    # written against the older API, which ARC still runs).
    s = sf.eq_spectrum(Tgas=float(temperature_k), pressure=pressure_pa / 1e5)
    wn, abscoeff = s.get("abscoeff", wunit="cm-1")
    n_cm3 = (pressure_pa / (K_B * temperature_k)) * 1e-6
    rho_g_cm3 = n_cm3 * M_H2O_G_MOL / N_A
    return np.asarray(wn, float), np.asarray(abscoeff, float) / rho_g_cm3


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("config", type=Path)
    ap.add_argument("--out-dir", type=Path, default=None)
    ap.add_argument("--only", default=None,
                    help="restrict to one point, as 'T,rho' (e.g. 1000,30)")
    ap.add_argument("--conventions", default=",".join(CONVENTIONS))
    ap.add_argument("--wavenum-min", type=float, default=WAVENUM_MIN_CM,
                    help=f"line-loading lower edge, cm^-1 (default {WAVENUM_MIN_CM:g}; "
                         f"use 10 where the 0-50 HITEMP file is cached)")
    ap.add_argument("--dry-run", action="store_true",
                    help="print the broadening pressures and stop (no RADIS)")
    args = ap.parse_args(argv)

    cfg = load_config(args.config)
    out_dir = args.out_dir or (cfg.output_path / "broadening")
    out_dir.mkdir(parents=True, exist_ok=True)
    band = tuple(cfg.opacity.band_cm or (100.0, 4000.0))
    conventions = [c for c in args.conventions.split(",") if c]

    table = broadening_pressures(cfg)
    if args.only:
        t_sel, rho_sel = (float(x) for x in args.only.split(","))
        table = [r for r in table
                 if abs(r["T"] - t_sel) < 0.5 and abs(r["rho"] - rho_sel) < 1e-3 * max(1, rho_sel)]
        if not table:
            raise SystemExit(f"--only {args.only} matches none of {POINTS}")

    print(f"[broadening] band {band[0]:g}-{band[1]:g} cm^-1, wstep {WSTEP_CM} cm^-1, "
          f"truncation (Lorentz half-width cap) {TRUNCATION_CM:g} cm^-1")
    head = (f"  {'T [K]':>6s} {'rho':>9s} | {'1 bar':>11s} {'ideal [bar]':>13s} "
            f"{'AQUA [bar]':>13s} {'AQUA/ideal':>11s}  phase")
    print(head)
    print("  " + "-" * (len(head) - 2))
    for r in table:
        print(f"  {r['T']:>6.0f} {r['rho']:>9.4g} | {1.0:>11.4g} "
              f"{r['p_ideal_pa'] / 1e5:>13.5g} {r['p_aqua_pa'] / 1e5:>13.5g} "
              f"{r['p_aqua_pa'] / r['p_ideal_pa']:>11.4g}  {r['aqua_phase_id']:>5d}")
    if args.dry_run:
        return 0

    sf = build_factory(args.wavenum_min)
    results = []
    for r in table:
        t_k, rho = r["T"], r["rho"]
        entry = dict(r)
        for conv in conventions:
            p_pa = r[f"p_{conv}_pa"]
            t0 = time.time()
            try:
                wn, kappa = compute_kappa(sf, t_k, p_pa)
            except Exception as exc:  # noqa: BLE001 - one bad corner must not end the run
                # At (3000 K, 2229) the AQUA broadening pressure is ~4.9e5 bar and
                # the nominal Lorentzian half-width is O(10^3) cm^-1. If RADIS
                # refuses that, the refusal text IS the answer to "what half-width
                # does RADIS accept" — so it is recorded, not swallowed.
                entry[f"error_{conv}"] = f"{type(exc).__name__}: {exc}"
                print(f"  T={t_k:g} rho={rho:g} {conv:<8s} P={p_pa / 1e5:>11.5g} bar  "
                      f"FAILED after {time.time() - t0:.1f}s: {type(exc).__name__}: "
                      f"{exc}", flush=True)
                continue
            kr = rosseland_mean(wn, kappa, t_k, band, cfg.opacity.rel_floor)
            kp = planck_mean(wn, kappa, t_k, band, cfg.opacity.rel_floor)
            npz = out_dir / f"T{t_k:g}_rho{rho:g}_{conv}.npz"
            np.savez_compressed(npz, wn=wn, kappa=kappa, T=t_k, rho_kg_m3=rho,
                                p_bar=p_pa / 1e5, convention=conv,
                                truncation_cm=TRUNCATION_CM, wstep_cm=WSTEP_CM,
                                isotopes=ISOTOPES, source="RADIS/hitemp")
            entry[f"kappa_R_{conv}_cm2_g"] = kr
            entry[f"kappa_P_{conv}_cm2_g"] = kp
            entry[f"spectrum_{conv}"] = npz.name
            print(f"  T={t_k:g} rho={rho:g} {conv:<8s} P={p_pa / 1e5:>11.5g} bar  "
                  f"kR={kr:>10.4g}  kP={kp:>10.4g} cm^2/g  ({time.time() - t0:.1f}s)",
                  flush=True)
        for conv in conventions:
            if (conv != "ref1bar" and "kappa_R_ref1bar_cm2_g" in entry
                    and f"kappa_R_{conv}_cm2_g" in entry):
                entry[f"ratio_R_{conv}_over_ref1bar"] = (
                    entry[f"kappa_R_{conv}_cm2_g"] / entry["kappa_R_ref1bar_cm2_g"])
                entry[f"ratio_P_{conv}_over_ref1bar"] = (
                    entry[f"kappa_P_{conv}_cm2_g"] / entry["kappa_P_ref1bar_cm2_g"])
        if {"kappa_R_aqua_cm2_g", "kappa_R_ideal_cm2_g"} <= set(entry):
            entry["ratio_R_aqua_over_ideal"] = (
                entry["kappa_R_aqua_cm2_g"] / entry["kappa_R_ideal_cm2_g"])
            entry["ratio_P_aqua_over_ideal"] = (
                entry["kappa_P_aqua_cm2_g"] / entry["kappa_P_ideal_cm2_g"])
        results.append(entry)
        _write(cfg, out_dir, band, conventions, results, args)  # crash-safe

    path = _write(cfg, out_dir, band, conventions, results, args)
    print(f"wrote {path}")
    return 0


def _write(cfg, out_dir: Path, band, conventions, results, args) -> Path:
    """(Re)write the results JSON, merging a single-point rerun into the full set."""
    payload = {
        "meta": (
            "HITEMP broadening-convention sensitivity (A59). Three conventions at each "
            "point: ref1bar (fixed 1 bar), ideal (P = rhoRT/M, the convention the "
            "DELIVERED spectra were built with — the brief's premise that they use a "
            "fixed 1 bar reference is incorrect, see the module docstring), and aqua "
            "(the structure model's own EOS). Mass opacity is convention-independent "
            "by construction (kappa = sigma*N_A/M), so the broadening pressure moves "
            "only the Voigt widths. Lorentzian profiles are truncated at "
            f"{TRUNCATION_CM:g} cm^-1 from line centre AND RENORMALISED to unit area, "
            "so line strength is conserved (verified: at T = 1000 K the integral of "
            "kappa over the delivered range is constant to 1e-4 across a 1e4 range in "
            "broadening pressure). The truncation redistributes far-wing absorption "
            "into the window rather than discarding it, and the sign of its effect on "
            "kappa_R is NOT determined a priori, so these ratios are not bounds in "
            "either direction."
        ),
        "band_cm": list(band),
        "rel_floor": cfg.opacity.rel_floor,
        "wstep_cm": WSTEP_CM,
        "truncation_cm": TRUNCATION_CM,
        "wavenum_range_cm": [args.wavenum_min, WAVENUM_MAX_CM],
        "points": [list(pt) for pt in POINTS],
        "isotopes": ISOTOPES,
        "conventions": conventions,
        "results": results,
    }
    path = out_dir.parent / "broadening_sensitivity.json"
    existing = json.loads(path.read_text()) if path.exists() else None
    if existing:  # keep points this invocation is not recomputing
        keep = [r for r in existing["results"]
                if not any(abs(r["T"] - n["T"]) < 0.5 and abs(r["rho"] - n["rho"]) < 1e-6
                           for n in results)]
        payload["results"] = keep + results
    payload["results"] = sorted(payload["results"], key=lambda r: (r["T"], r["rho"]))
    path.write_text(json.dumps(payload, indent=1) + "\n")
    return path


if __name__ == "__main__":
    raise SystemExit(main())
