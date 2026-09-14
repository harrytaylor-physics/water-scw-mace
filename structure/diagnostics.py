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


def aqua_coverage(eos: AquaEOS, pressure_pa, temperature_k) -> dict:
    """Does a profile stay inside the AQUA table? (A47)

    ``AquaEOS`` interpolates with ``np.interp`` and ``np.clip``, so a query
    outside the table is silently *clamped* to its edge rather than raising. A
    scan over planet space has to test for that explicitly, or a cell whose
    envelope runs off the table is reported as though it were a model.

    The temperature test is against the table's log T grid; the pressure test is
    against the two bracketing T-columns' own log P ranges (the table is
    (ρ, T)-indexed, so each column covers a different pressure range, A2).
    """
    t = np.atleast_1d(np.asarray(temperature_k, float))
    p = np.atleast_1d(np.asarray(pressure_pa, float))
    t_lo, t_hi = eos.T_range_k
    below_t, above_t = t < t_lo, t > t_hi

    tg = eos._t.logT_grid
    lp = eos._t.logP_strict  # (n_rho, n_T), increasing down each column
    log_t = np.log10(np.clip(t, t_lo, t_hi))
    jt = np.clip(np.searchsorted(tg, log_t), 1, tg.size - 1)
    j0 = jt - 1
    p_lo = 10.0 ** np.maximum(lp[0, j0], lp[0, jt])
    p_hi = 10.0 ** np.minimum(lp[-1, j0], lp[-1, jt])
    below_p, above_p = p < p_lo, p > p_hi

    parts = []
    if below_t.any():
        parts.append(f"T down to {t.min():.1f} K below the table floor {t_lo:.1f} K")
    if above_t.any():
        parts.append(f"T up to {t.max():.1f} K above the table ceiling {t_hi:.1f} K")
    if below_p.any():
        parts.append(f"P below the column floor at {int(below_p.sum())} point(s)")
    if above_p.any():
        parts.append(f"P above the column ceiling at {int(above_p.sum())} point(s)")
    return {
        "inside": not (below_t.any() or above_t.any() or below_p.any() or above_p.any()),
        "n_below_t": int(below_t.sum()), "n_above_t": int(above_t.sum()),
        "n_below_p": int(below_p.sum()), "n_above_p": int(above_p.sum()),
        "t_range_k": (float(t.min()), float(t.max())),
        "table_t_range_k": (t_lo, t_hi),
        "summary": "; ".join(parts) or "inside the table",
    }


def _sample_mass_weights(mass_kg: np.ndarray) -> np.ndarray:
    """Per-sample mass weights (central differences) summing to the envelope mass."""
    m = np.asarray(mass_kg, float)
    edges = np.empty(m.size + 1)
    edges[1:-1] = 0.5 * (m[:-1] + m[1:])
    edges[0], edges[-1] = m[0], m[-1]
    return np.abs(np.diff(edges))


def radiative_layers(result, t_min_k: float = 0.0, rho_min_si: float = 0.0) -> list[dict]:
    """Contiguous radiative layers of a profile, with their extents and mass share.

    ``t_min_k``/``rho_min_si`` restrict the layers to a (T, ρ) box — with the FP
    grid's corner (1000 K, 108.796 kg/m³) this is the "in-box" radiative structure
    the first-principles data can actually speak to (A42).
    """
    prof = result.profile
    rho = np.asarray(prof["rho_si"], float)
    t = np.asarray(prof["temperature_k"], float)
    w = _sample_mass_weights(prof["mass_kg"])
    total = w.sum()
    mask = (~np.asarray(prof["is_convective"], bool)) & (t >= t_min_k) & (rho >= rho_min_si)

    layers: list[dict] = []
    idx = np.flatnonzero(mask)
    if idx.size == 0:
        return layers
    for run in np.split(idx, np.flatnonzero(np.diff(idx) > 1) + 1):
        layers.append({
            "t_min_k": float(t[run].min()), "t_max_k": float(t[run].max()),
            "rho_min_si": float(rho[run].min()), "rho_max_si": float(rho[run].max()),
            "mass_fraction": float(w[run].sum() / total),
            "n_samples": int(run.size),
        })
    return layers


def in_box_radiative_mass_fraction(result, t_min_k: float = 1000.0,
                                   rho_min_si: float = 108.796) -> float:
    """Envelope mass fraction that is radiative *and* inside the FP grid's corner."""
    prof = result.profile
    w = _sample_mass_weights(prof["mass_kg"])
    mask = (
        (~np.asarray(prof["is_convective"], bool))
        & (np.asarray(prof["temperature_k"], float) >= t_min_k)
        & (np.asarray(prof["rho_si"], float) >= rho_min_si)
    )
    return float(w[mask].sum() / w.sum())


CATEGORIES = ("interpolated", "in_box_out_of_hull", "clamped_t", "clamped_rho",
              "dilute_flat")


def classify_evaluations(source, temperature_k, rho_si) -> np.ndarray:
    """Label every opacity evaluation by what the call actually was (A48).

    One of :data:`CATEGORIES`, per point, judged against the sub-source the
    composite actually uses there (A12):

    ``dilute_flat``
        line-list-dominated below its ρ floor: κ held flat in ρ (A26), a
        justified dilute-limit extrapolation rather than a clamp — tested first,
        because such a point is also formally "below the grid's ρ floor".
    ``interpolated``
        inside the used sub-source's hull: real data.
    ``clamped_t`` / ``clamped_rho``
        clamped, and outside the grid's T / ρ range: the grid would have to be
        *extended* to cover it.
    ``in_box_out_of_hull``
        clamped, but inside the grid's (T, ρ) box — the scattered FP grid's own
        gaps (A6), which a campaign can fill without extending its range.
    """
    t = np.atleast_1d(np.asarray(temperature_k, float))
    r = np.atleast_1d(np.asarray(rho_si, float))
    t, r = np.broadcast_arrays(t, r)
    t, r = t.ravel(), r.ravel()
    clamped = np.asarray(source.clamped_mask(t, r)).ravel()
    names, subs = used_subsource(source, r)
    names = np.asarray(names).ravel()

    out = np.full(t.size, "interpolated", dtype=object)
    for name, sub in subs.items():
        sel = names == name
        if not sel.any():
            continue
        bounds = _grid_bounds(sub)
        if bounds is None:
            continue
        (t_lo, t_hi), (r_lo, r_hi) = bounds
        flat = bool(getattr(getattr(sub, "table", None), "dilute_flat", False))
        out_t = sel & ((t < t_lo) | (t > t_hi))
        out_r = sel & ((r < r_lo) | (r > r_hi))
        out[sel & clamped & ~out_t & ~out_r] = "in_box_out_of_hull"
        out[sel & clamped & out_r] = "clamped_rho"
        out[sel & clamped & out_t] = "clamped_t"  # T range is the FP grid's binding limit
        if flat:
            out[sel & (r < r_lo)] = "dilute_flat"
    return out


def category_mass_fractions(result, categories) -> dict:
    """Envelope mass fraction in each of :data:`CATEGORIES` (A48)."""
    w = _sample_mass_weights(result.profile["mass_kg"])
    cats = np.asarray(categories, dtype=object)
    total = w.sum()
    return {c: float(w[cats == c].sum() / total) for c in CATEGORIES}


def clamped_t_split(result, categories, table) -> dict:
    """Split the ``clamped_t`` mass into "too cold" and "too hot" (A48).

    Both are clamps in temperature, but they are different problems: below the
    grid's floor the FP data would have to be *extended downwards* (the old
    fiducial), while above its ceiling the deep convective interior is simply
    hotter than the table goes — where opacity does not shape the structure
    anyway.
    """
    w = _sample_mass_weights(result.profile["mass_kg"])
    t = np.asarray(result.profile["temperature_k"], float)
    cats = np.asarray(categories, dtype=object)
    sel = cats == "clamped_t"
    t_lo = float(np.min(table.temperature_k))
    t_hi = float(np.max(table.temperature_k))
    total = w.sum()
    return {
        "too_cold": float(w[sel & (t < t_lo)].sum() / total),
        "too_hot": float(w[sel & (t > t_hi)].sum() / total),
    }


def campaign_nodes(cases, table, t_step_k: float = 100.0, log_rho_step: float = 0.25,
                   t_floor_k: float = 1000.0, rho_floor_si: float = 108.796) -> dict:
    """(T, ρ) lattice nodes that would close the FP grid's gaps (A48).

    ``cases`` are ``(StructureResult, categories)`` pairs. A trajectory sample
    needs new FP data if it lies inside the FP grid's own (T, ρ) box yet outside
    its interpolation hull. Every such sample sits inside one cell of the
    ``t_step_k`` × ``log_rho_step`` lattice, so adding that cell's **four
    corners** puts the sample inside the convex hull of the enlarged grid — which
    is then verified, not assumed, with the same interpolator the model uses.

    ``extend_nodes`` is the separate, larger ask: the same lattice over the
    dense-but-cold segments (ρ inside the grid's range, T below its floor) where
    the composite currently clamps in T. Those need the grid *extended*, not
    filled, so they are reported apart from the in-box specification.

    Because the nodes are the *bracketing* corners of each sample's lattice cell,
    a node may sit up to one lattice step outside the grid's present (T, ρ)
    range: that is what enclosing a sample near the edge requires.
    """
    from scipy.interpolate import LinearNDInterpolator

    t_grid = np.asarray(table.temperature_k, float)
    r_grid = np.asarray(table.rho_si, float)
    t_lo, t_hi = float(t_grid.min()), float(t_grid.max())
    r_lo, r_hi = float(r_grid.min()), float(r_grid.max())
    pts = np.column_stack([t_grid, np.log10(r_grid)])
    lin = LinearNDInterpolator(pts, np.zeros(pts.shape[0]), rescale=True)

    need_fill, need_extend = [], []
    for result, _categories in cases:
        t = np.asarray(result.profile["temperature_k"], float)
        r = np.asarray(result.profile["rho_si"], float)
        in_box = (t >= max(t_lo, t_floor_k)) & (t <= t_hi) & (r >= max(r_lo, rho_floor_si)) & (r <= r_hi)
        outside_hull = np.isnan(lin(np.column_stack([t, np.log10(r)])))
        sel = in_box & outside_hull
        if sel.any():
            need_fill.append(np.column_stack([t[sel], r[sel]]))
        cold = (t < t_lo) & (r >= max(r_lo, rho_floor_si)) & (r <= r_hi)
        if cold.any():
            need_extend.append(np.column_stack([t[cold], r[cold]]))

    def corners(points):
        """The four lattice corners bracketing each point, deduplicated."""
        if not len(points):
            return []
        p = np.vstack(points)
        t_f = np.floor(p[:, 0] / t_step_k) * t_step_k
        lr_f = np.floor(np.log10(p[:, 1]) / log_rho_step) * log_rho_step
        nodes = set()
        for dt in (0.0, t_step_k):
            for dr in (0.0, log_rho_step):
                for tt, rr in zip(t_f + dt, lr_f + dr):
                    nodes.add((round(float(tt), 6), round(float(rr), 6)))
        return sorted((t, float(10.0**lr)) for t, lr in nodes)

    fill_nodes = corners(need_fill)
    extend_nodes = corners(need_extend)

    encloses = True
    if fill_nodes and need_fill:
        new = np.vstack([pts, np.column_stack([[t for t, _ in fill_nodes],
                                               np.log10([r for _, r in fill_nodes])])])
        lin2 = LinearNDInterpolator(new, np.zeros(new.shape[0]), rescale=True)
        probe = np.vstack(need_fill)
        encloses = bool(np.all(~np.isnan(
            lin2(np.column_stack([probe[:, 0], np.log10(probe[:, 1])])))))
    return {
        "fill_nodes": fill_nodes,
        "extend_nodes": extend_nodes,
        "n_in_box_needed": int(sum(len(x) for x in need_fill)),
        "n_extend_needed": int(sum(len(x) for x in need_extend)),
        "fill_encloses": encloses,
        "lattice": (t_step_k, log_rho_step),
    }


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
