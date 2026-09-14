"""Opacity-clamp audit for a converged structure (A41).

The composite closes the deep radiative window that the line list opens (A28,
A35). Whether that closure is *evidence* depends on where the FP opacity doing
the closing came from: interpolation inside the FP (T, ρ) grid, or a
nearest-neighbour clamp outside it. This module answers that question for a
given case, per density range, and splits the clamps by reason (T below the FP
grid's 1000 K floor, ρ outside its range, or inside the box but outside the
scattered grid's convex hull, A6).

CLI::

    python structure/diagnostics.py config/fiducial.yaml --t-int 20 --t-eq 500 \
        --rho-range 400,550
"""

from __future__ import annotations

if __name__ == "__main__":  # allow `python structure/diagnostics.py ...` as a CLI
    import sys as _sys
    from pathlib import Path as _Path

    _sys.path.insert(0, str(_Path(__file__).resolve().parents[1]))

import numpy as np

from config import Config
from eos import AquaEOS
from opacity import CompositeOpacitySource, get_source
from structure.model import build_structure


def used_subsource(source, rho_si):
    """The sub-source whose table actually supplies κ at each ρ (A12/A35).

    For the composite this is the FP ("dense") source where its blend weight is
    ≥ ½ and the line list below; for a plain grid source it is the source itself.
    """
    if not isinstance(source, CompositeOpacitySource):
        return np.full(np.shape(rho_si), source.name, dtype=object), {source.name: source}
    dense = source._weight(np.asarray(rho_si, float)) >= 0.5
    names = np.where(dense, source.dense.name, source.dilute.name)
    return names, {source.dense.name: source.dense, source.dilute.name: source.dilute}


def _grid_bounds(sub) -> tuple[tuple[float, float], tuple[float, float]] | None:
    table = getattr(sub, "table", None)
    if table is None:
        return None
    return (
        (float(table.temperature_k.min()), float(table.temperature_k.max())),
        (float(table.rho_si.min()), float(table.rho_si.max())),
    )


def clamp_reasons(source, temperature_k, rho_si) -> dict:
    """Split opacity evaluations into interpolated vs clamped, by reason.

    ``clamped`` is the authoritative hull test the run's reliability uses
    (``clamped_mask``); the reason counts are a box test against the used
    sub-source's grid, so a point can be in-box yet out-of-hull (the FP grid is
    scattered, A6). Line-list ρ-floor queries are dilute-flat, not clamps (A26).
    """
    t = np.atleast_1d(np.asarray(temperature_k, float))
    r = np.atleast_1d(np.asarray(rho_si, float))
    t, r = np.broadcast_arrays(t, r)
    clamped = np.asarray(source.clamped_mask(t, r)).ravel()
    names, subs = used_subsource(source, r)
    names = np.asarray(names).ravel()
    t, r = t.ravel(), r.ravel()

    out = {
        "n_points": int(t.size),
        "n_clamped": int(clamped.sum()),
        "n_interpolated": int((~clamped).sum()),
        "by_source": {},
    }
    for name, sub in subs.items():
        sel = names == name
        if not sel.any():
            continue
        bounds = _grid_bounds(sub)
        entry = {"n_points": int(sel.sum()), "n_clamped": int(clamped[sel].sum())}
        if bounds is not None:
            (t_lo, t_hi), (r_lo, r_hi) = bounds
            below_t = (t < t_lo) & sel
            above_t = (t > t_hi) & sel
            below_r = (r < r_lo) & sel
            above_r = (r > r_hi) & sel
            in_box = sel & ~(below_t | above_t | below_r | above_r)
            entry.update({
                "grid_t_k": [t_lo, t_hi],
                "grid_rho_si": [r_lo, r_hi],
                "n_below_t_floor": int(below_t.sum()),
                "n_above_t_ceiling": int(above_t.sum()),
                "n_below_rho_floor": int(below_r.sum()),
                "n_above_rho_ceiling": int(above_r.sum()),
                "n_in_box_but_clamped": int((in_box & clamped).sum()),
                "dilute_flat": bool(getattr(getattr(sub, "table", None), "dilute_flat", False)),
            })
        out["by_source"][name] = entry
    return out


def clamp_audit(
    cfg: Config,
    t_eq_k: float,
    t_int_k: float,
    source_name: str,
    rho_range: tuple[float, float],
    eos: AquaEOS | None = None,
    n_sample: int = 1500,
    band: tuple[float, float] | None = None,
) -> dict:
    """Build one structure and audit its opacity evaluations (whole profile and
    within ``rho_range``)."""
    eos = eos or AquaEOS.from_config(cfg)
    res = build_structure(cfg, t_eq_k, t_int_k, source_name, eos=eos, n_sample=n_sample,
                          band=band)
    source = get_source(source_name, cfg, band=band)
    prof = res.profile
    t, rho = prof["temperature_k"], prof["rho_si"]
    sel = (rho >= rho_range[0]) & (rho <= rho_range[1])

    return {
        "source": source_name,
        "t_eq_k": t_eq_k,
        "t_int_k": t_int_k,
        "n_sample": n_sample,
        "result": res,
        "rho_range": tuple(rho_range),
        "in_range_temperature_k": t[sel].tolist(),
        "in_range_rho_si": rho[sel].tolist(),
        "in_range_convective": prof["is_convective"][sel].tolist(),
        "in_range": clamp_reasons(source, t[sel], rho[sel]),
        "whole_profile": clamp_reasons(source, t, rho),
        "clamp_fraction_all": res.clamp_fraction_all,
        "clamp_fraction_radiative": res.clamp_fraction_radiative,
    }


def print_audit(audit: dict) -> None:
    """Human-readable audit report."""
    res = audit["result"]
    lo, hi = audit["rho_range"]
    print(f"\n--- {audit['source']}: T_int={audit['t_int_k']:g} K, T_eq={audit['t_eq_k']:g} K, "
          f"n_sample={audit['n_sample']} ---")
    print(f"  R_p = {res.radius_p_rearth:.5f} R_earth; {len(res.rcb_densities_si)} RCB at rho = "
          f"{', '.join(f'{x:.1f}' for x in res.rcb_densities_si) or '-'} kg/m^3; "
          f"max radiative rho = {res.max_radiative_rho_si:.1f}")
    print(f"  clamp_fraction_all = {audit['clamp_fraction_all']:.4f}  "
          f"clamp_fraction_radiative = {audit['clamp_fraction_radiative']:.4f}  "
          f"(profile samples; StructureResult.clamp_fraction={res.clamp_fraction:.4f} counts "
          f"ODE steps)")

    temps = audit["in_range_temperature_k"]
    rhos = audit["in_range_rho_si"]
    conv = audit["in_range_convective"]
    print(f"  samples with rho in [{lo:g}, {hi:g}] kg/m^3: {len(temps)}")
    for tt, rr, cc in zip(temps, rhos, conv):
        print(f"      rho = {rr:8.2f}   T = {tt:8.2f} K   {'convective' if cc else 'RADIATIVE'}")

    for label, key in (("in that rho range", "in_range"), ("whole profile", "whole_profile")):
        d = audit[key]
        if d["n_points"] == 0:
            print(f"  clamp attribution, {label}: no samples")
            continue
        print(f"  clamp attribution, {label}: {d['n_clamped']}/{d['n_points']} clamped, "
              f"{d['n_interpolated']} interpolated")
        for name, e in d["by_source"].items():
            line = (f"      via {name}: {e['n_clamped']}/{e['n_points']} clamped")
            if "grid_t_k" in e:
                line += (f"; grid T [{e['grid_t_k'][0]:g}, {e['grid_t_k'][1]:g}] K, "
                         f"rho [{e['grid_rho_si'][0]:g}, {e['grid_rho_si'][1]:g}] kg/m^3"
                         f"\n          below T floor: {e['n_below_t_floor']}, "
                         f"above T ceiling: {e['n_above_t_ceiling']}, "
                         f"below rho floor: {e['n_below_rho_floor']}"
                         f"{' (dilute-flat, not a clamp)' if e['dilute_flat'] else ''}, "
                         f"above rho ceiling: {e['n_above_rho_ceiling']}, "
                         f"in-box but out-of-hull: {e['n_in_box_but_clamped']}")
            print(line)


def _cli(argv=None) -> None:
    import argparse

    from config import load_config

    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("config", help="run config YAML")
    ap.add_argument("--t-int", type=float, default=20.0)
    ap.add_argument("--t-eq", type=float, default=None, help="default: first config T_eq")
    ap.add_argument("--sources", default="composite,hitemp")
    ap.add_argument("--rho-range", default="400,550", help="lo,hi in kg/m^3")
    ap.add_argument("--n-sample", type=int, default=1500)
    args = ap.parse_args(argv)

    cfg = load_config(args.config)
    t_eq = args.t_eq if args.t_eq is not None else float(cfg.irradiation.t_eq_k[0])
    lo, hi = (float(x) for x in args.rho_range.split(","))
    eos = AquaEOS.from_config(cfg)
    names = [s.strip() for s in args.sources.split(",")]
    from figures._common import comparison_band

    band = comparison_band(cfg, names)

    print(f"Opacity-clamp audit (A41): {cfg.name}, M_c={cfg.core.mass_mearth:g} M_earth, "
          f"f_env={cfg.envelope.f_env:.0%}, band=[{band[0]:g}, {band[1]:g}] cm^-1")
    audits = {}
    for name in names:
        audits[name] = clamp_audit(cfg, t_eq, args.t_int, name, (lo, hi), eos=eos,
                                   n_sample=args.n_sample, band=band)
        print_audit(audits[name])

    comp = audits.get("composite")
    if comp is not None:
        d = comp["in_range"]
        fp = d["by_source"].get("fp", {})
        n_fp = fp.get("n_points", 0)
        n_clamped = fp.get("n_clamped", 0)
        print(f"\nVERDICT for rho in [{lo:g}, {hi:g}] kg/m^3 at T_int={args.t_int:g} K, "
              f"T_eq={t_eq:g} K:")
        if n_fp == 0:
            print("  the composite does not use FP at these densities at all.")
        elif n_clamped == n_fp:
            print(f"  the FP opacity there is CLAMPED at every one of the {n_fp} samples "
                  f"({fp.get('n_below_t_floor', 0)} below the FP grid's "
                  f"{fp.get('grid_t_k', [float('nan')])[0]:g} K floor) — the closure of the "
                  f"deep line-list radiative zone at this case is decided by "
                  f"NEAREST-NEIGHBOUR EXTRAPOLATED FP opacity, not by interpolated FP data.")
        elif n_clamped == 0:
            print(f"  the FP opacity there is INTERPOLATED at all {n_fp} samples — the closure "
                  f"is decided by FP data inside its (T, rho) grid.")
        else:
            print(f"  mixed: {n_clamped} of {n_fp} FP samples are clamped, "
                  f"{n_fp - n_clamped} interpolated.")


if __name__ == "__main__":
    _cli()
