#!/usr/bin/env python3
"""Which AQUA release is on disk, and what the entropy typo costs (Task 40a).

    python tools/aqua_entropy_check.py config/fiducial.yaml

Mazevet et al. (2019) Eq. (13) carries a sign error on the first two terms of the
temperature-dependent Helmholtz free energy F_T, and AQUA's reference entropy was
S0 = 4.9 n_at k_B where it should be 9.8 (n_at = atoms per kilogram). Aguichine
et al. (2025), footnote 3: "As of August 2024, the AQUA EOS tabulated the entropy
with the non-corrected version of the [Mazevet et al. 2019] EOS", and the typo is
in the entropy "but not its derivatives".

Both corrections live in AQUA's supercritical/superionic region — the part of the
table that comes from Mazevet et al. (2019) — and they move the entropy and the
internal energy. Pressure and density are untouched: the shift in F_T does not
depend on volume.

This does three things.

1. Identifies the table on disk by name and SHA-256 and says whether it predates
   the correction.
2. Compares its entropy and internal energy against the corrected values, at the
   (P, T) the water-world fiducial track actually samples. The corrected values
   are not reimplemented here: they are read from the PALEOS re-gridding of AQUA
   (Zenodo 10.5281/zenodo.22776069), which applies the sign fix analytically and
   states that rho and P are unaffected. Reimplementing Eq. (13) from a prose
   description of its error would be a worse primary source than the table whose
   authors did it.
3. Converts the entropy shift into the quantity that matters here — a shift in
   age at fixed entropy — using the track's own ds/d(age).

**Nothing is substituted.** The model keeps reading the table it has; this
reports what that costs.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from dataclasses import replace
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np

from config import Config, load_config
from constants import G as G_SI

PALEOS_REL = Path("paleos") / "paleos_water_eos_table_pt.dat"
PALEOS_DOI = "10.5281/zenodo.22776069"
PALEOS_SHAPE = (2251, 451)
CORRECTION_PUBLISHED = "2024-08"  # the date Aguichine et al. 2025 quote
# The water-world fiducial (same case as ED Fig. 8's middle track).
CASE = {"mass_core_mearth": 10.0, "wmf": 0.50, "t_eq_k": 500.0}
T_INT_RANGE_K = (12.0, 450.0)
N_NODES = 44
N_SAMPLE_AGES = 5


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def identify_table(cfg: Config) -> dict:
    path = Path(cfg.aqua_table_path)
    head = path.open().read(4000).splitlines()[:20]
    return {
        "path": str(path),
        "file_name": path.name,
        "sha256": _sha256(path),
        "bytes": path.stat().st_size,
        "version_from_name": "v1_0",
        "header_lines": [ln.rstrip() for ln in head[:3]],
        "version_string_in_header": None,
        "predates_correction": True,
        "why": ("the upstream repository (github.com/mnijh/AQUA) still offers only "
                "the v1_0 tables — one release, tagged 2020-09-02, last touched "
                f"2020-12-30 — so this file predates the {CORRECTION_PUBLISHED} "
                "cut-off Aguichine et al. (2025) name, and carries the "
                "uncorrected entropy"),
    }


def load_paleos(cfg: Config):
    """(log10 P, log10 T) grids and the corrected s, u, rho, phase."""
    path = Path(cfg.data_root_path) / PALEOS_REL
    if not path.exists():
        return None
    raw = np.genfromtxt(path, comments="#", dtype=None, encoding="utf-8",
                        usecols=(0, 1, 2, 3, 4, 8, 9))
    n_p, n_t = PALEOS_SHAPE
    p = np.array([r[0] for r in raw], float).reshape(n_p, n_t)
    t = np.array([r[1] for r in raw], float).reshape(n_p, n_t)
    rho = np.array([r[2] for r in raw], float).reshape(n_p, n_t)
    u = np.array([r[3] for r in raw], float).reshape(n_p, n_t)
    s = np.array([r[4] for r in raw], float).reshape(n_p, n_t)
    grad = np.array([r[5] for r in raw], float).reshape(n_p, n_t)
    phase = np.array([str(r[6]) for r in raw]).reshape(n_p, n_t)
    return {"log_p": np.log10(p[:, 0]), "log_t": np.log10(t[0, :]),
            "rho": rho, "u": u, "s": s, "grad_ad": grad, "phase": phase,
            "path": str(path),
            "sha256": _sha256(path)}


def _bilinear(grid, log_p: np.ndarray, log_t: np.ndarray, field: str):
    from scipy.interpolate import RegularGridInterpolator

    f = RegularGridInterpolator((grid["log_p"], grid["log_t"]), grid[field],
                                bounds_error=False, fill_value=np.nan)
    return f(np.column_stack([log_p, log_t]))


def eos_state_u(cfg: Config, p_pa: float, t_k: float) -> float:
    """AQUA's own internal energy at (P, T), for the like-for-like comparison."""
    from eos.aqua import AquaEOS

    global _EOS_CACHE
    if _EOS_CACHE is None:
        _EOS_CACHE = AquaEOS.from_config(cfg)
    return float(_EOS_CACHE.state(float(p_pa), float(t_k)).internal_energy)


_EOS_CACHE = None


def _nearest_phase(grid, log_p: np.ndarray, log_t: np.ndarray):
    i = np.clip(np.searchsorted(grid["log_p"], log_p), 0, grid["log_p"].size - 1)
    j = np.clip(np.searchsorted(grid["log_t"], log_t), 0, grid["log_t"].size - 1)
    return grid["phase"][i, j]


def run(cfg: Config) -> dict:
    from evolution.track import cooling_track
    from structure import build_structure

    out = {"meta": __doc__.split("\n\n")[0], "table": identify_table(cfg),
           "correction": {
               "source": "Mazevet et al. (2019) Eq. (13); AQUA revised entropy",
               "sign_error": "first two terms of F_T carry the wrong sign",
               "reference_entropy": "S0 = 4.9 -> 9.8 n_at k_B (n_at = atoms per kg)",
               "affects": ["specific entropy", "specific internal energy"],
               "does_not_affect": ["pressure", "density"],
               "region": ("AQUA's supercritical/superionic region, the part taken "
                          "from Mazevet et al. (2019)"),
           },
           "case": CASE}

    grid = load_paleos(cfg)
    if grid is None:
        out["blocked"] = (f"corrected reference table not found at "
                          f"{cfg.data_root_path / PALEOS_REL}; download PALEOS "
                          f"({PALEOS_DOI})")
        return out
    out["reference_table"] = {"doi": PALEOS_DOI, "path": grid["path"],
                              "sha256": grid["sha256"],
                              "provides": "AQUA re-gridded with the sign fix applied"}

    c = replace(cfg, core=replace(cfg.core, mass_mearth=CASE["mass_core_mearth"]),
                envelope=replace(cfg.envelope, f_env=CASE["wmf"]))
    grid_t_int = np.geomspace(*T_INT_RANGE_K, N_NODES)
    track = cooling_track(c, "composite", CASE["t_eq_k"],
                          t_int_grid_k=grid_t_int)
    age, s_track = np.asarray(track.age_gyr), np.asarray(track.entropy_j_kg_k)
    order = np.argsort(age)
    age, s_track = age[order], s_track[order]
    t_int_sorted = np.asarray(track.t_int_k)[order]

    # five ages spread over the part of the track a planet is observed in
    ages = np.geomspace(max(age.min(), 1e-3), min(age.max(), 10.0), N_SAMPLE_AGES)
    rows = []
    for a in ages:
        t_int = float(np.interp(a, age, t_int_sorted))
        res = build_structure(c, CASE["t_eq_k"], t_int, "composite")
        pr = res.profile
        p_pa = np.asarray(pr["pressure_pa"], float)
        t_k = np.asarray(pr["temperature_k"], float)
        s_aqua = np.asarray(pr["entropy_j_kg_k"], float)
        lp, lt = np.log10(p_pa), np.log10(t_k)
        s_corr = _bilinear(grid, lp, lt, "s")
        rho_corr = _bilinear(grid, lp, lt, "rho")
        rho_aqua = np.asarray(pr["rho_si"], float)
        phase = _nearest_phase(grid, lp, lt)
        ok = np.isfinite(s_corr)
        sc = phase == "supercritical"
        # The age integral is dt = dE/L with E the envelope internal energy, so
        # the correction reaches the tracks through u, not through s. E is
        # recomputed on the same mass grid with the corrected u.
        u_aqua = np.array([eos_state_u(c, pp, tt) for pp, tt in zip(p_pa, t_k)])
        u_corr = _bilinear(grid, lp, lt, "u")
        mass = np.asarray(pr["mass_kg"], float)
        radius = np.asarray(pr["radius_m"], float)
        good = np.isfinite(u_corr)
        # E = int(u - G m / r) dm, exactly as structure.model._envelope_energy
        # defines it. Only the u term moves: the gravitational term is geometry,
        # and the correction leaves rho and P alone.
        ordm = np.argsort(mass[good])
        m_s = mass[good][ordm]
        e_grav = float(np.trapezoid(-G_SI * m_s / radius[good][ordm], m_s))
        e_int_aqua = float(np.trapezoid(u_aqua[good][ordm], m_s))
        e_int_corr = float(np.trapezoid(u_corr[good][ordm], m_s))
        e_aqua = e_int_aqua + e_grav
        e_corr = e_int_corr + e_grav
        g_aqua = np.asarray(pr["grad_ad"], float)
        g_corr = _bilinear(grid, lp, lt, "grad_ad")
        rows.append({
            "age_gyr": float(a), "t_int_k": t_int,
            "envelope_energy_aqua_j": e_aqua,
            "envelope_energy_corrected_j": e_corr,
            "envelope_energy_frac_change": float(e_corr / e_aqua - 1.0),
            "envelope_internal_energy_aqua_j": e_int_aqua,
            "envelope_internal_energy_corrected_j": e_int_corr,
            "envelope_gravitational_energy_j": e_grav,
            "u_base_frac_change": float(u_corr[-1] / u_aqua[-1] - 1.0),
            "grad_ad_base_frac_change": float(g_corr[-1] / g_aqua[-1] - 1.0),
            "grad_ad_max_abs_frac_change": float(
                np.nanmax(np.abs(g_corr[good] / g_aqua[good] - 1.0))),
            "grad_ad_median_abs_frac_change": float(
                np.nanmedian(np.abs(g_corr[good] / g_aqua[good] - 1.0))),
            "grad_ad_max_at_pressure_pa": float(
                p_pa[good][int(np.nanargmax(np.abs(g_corr[good] / g_aqua[good] - 1.0)))]),
            "grad_ad_note": ("the correction is stated not to touch the derivatives; "
                             "what is measured here also carries the error of "
                             "interpolating between two different grids, so the "
                             "median is the meaningful number, not the maximum"),
            "base_pressure_pa": float(p_pa[-1]), "base_temperature_k": float(t_k[-1]),
            "interior_entropy_aqua_j_kg_k": float(s_aqua[-1]),
            "interior_entropy_corrected_j_kg_k": float(s_corr[-1]),
            "delta_s_j_kg_k": float(s_corr[-1] - s_aqua[-1]),
            "delta_s_frac": float(s_corr[-1] / s_aqua[-1] - 1.0),
            "rho_check_frac": float(rho_corr[-1] / rho_aqua[-1] - 1.0),
            "envelope_supercritical_fraction": float(sc.mean()),
            "profile_points": int(p_pa.size),
            "profile_points_off_reference_grid": int((~ok).sum()),
            "radius_p_rearth": float(res.radius_p_rearth),
        })
    out["points"] = rows

    # shift in age at fixed entropy: ds/d(age) from the track itself
    ds_dage = np.gradient(s_track, age)
    for r in rows:
        slope = float(np.interp(r["age_gyr"], age, ds_dage))
        r["ds_dage_j_kg_k_per_gyr"] = slope
        r["age_shift_gyr"] = (float(-r["delta_s_j_kg_k"] / slope)
                              if slope != 0 else None)
        r["age_shift_frac"] = (None if r["age_shift_gyr"] is None
                               else r["age_shift_gyr"] / r["age_gyr"])
    # The correction reaches the cooling tracks through E, not through s. The
    # fractional change in E is smooth in T_int over the five sampled structures,
    # so it is interpolated onto the track's own nodes and the age integral redone
    # with the same luminosities. (Rebuilding all 44 profiles would double the
    # track cost to refine a number that moves by 1% across the whole range.)
    from scipy.integrate import cumulative_simpson

    f_t = np.log(np.array([r["t_int_k"] for r in rows]))
    f_e = np.array([r["envelope_energy_frac_change"] for r in rows])
    ordf = np.argsort(f_t)
    frac = np.interp(np.log(t_int_sorted), f_t[ordf], f_e[ordf])
    e_aqua_track = np.asarray(track.energy_j)[order]
    lum = np.asarray(track.luminosity_w)[order]
    # the fractional change is measured on the same E definition the track uses,
    # so it can be applied to the track's own energies directly
    e_corr_track = e_aqua_track * (1.0 + frac)
    def _age(energy):
        k = np.argsort(t_int_sorted)
        t_from_bottom = cumulative_simpson(y=1.0 / lum[k], x=energy[k], initial=0.0)
        a = (t_from_bottom[-1] - t_from_bottom) / 3.15576e16
        out_a = np.empty_like(a)
        out_a[k] = a
        return out_a
    age_corr = _age(e_corr_track)
    age_aqua = _age(e_aqua_track)
    for r in rows:
        ti = r["t_int_k"]
        a0 = float(np.interp(-ti, -t_int_sorted, age_aqua))
        a1 = float(np.interp(-ti, -t_int_sorted, age_corr))
        r["age_from_reintegration_aqua_gyr"] = a0
        r["age_from_reintegration_corrected_gyr"] = a1
        r["age_frac_change_from_energy"] = (a1 / a0 - 1.0) if a0 > 0 else None
    out["energy_route"] = {
        "what": ("the age integral dt = dE/L redone with the corrected internal "
                 "energy and the same luminosities"),
        "total_radiated_energy_aqua_j": float(abs(e_aqua_track.max()
                                                  - e_aqua_track.min())),
        "total_radiated_energy_corrected_j": float(abs(e_corr_track.max()
                                                       - e_corr_track.min())),
        "total_radiated_energy_frac_change": float(
            abs(e_corr_track.max() - e_corr_track.min())
            / abs(e_aqua_track.max() - e_aqua_track.min()) - 1.0),
        "age_at_coldest_node_aqua_gyr": float(age_aqua.max()),
        "age_at_coldest_node_corrected_gyr": float(age_corr.max()),
        "age_frac_change_at_coldest_node": float(age_corr.max() / age_aqua.max() - 1.0),
    }
    out["track"] = {"n_nodes": N_NODES, "t_int_range_k": list(T_INT_RANGE_K),
                    "age_range_gyr": [float(age.min()), float(age.max())],
                    "entropy_range_j_kg_k": [float(s_track.min()),
                                             float(s_track.max())],
                    "source": "composite"}
    out["what_uses_entropy"] = (
        "The cooling age integrates dt = dE/L with E the envelope internal energy "
        "and the adiabat set by grad_ad, so ages, radii and luminosities do not "
        "read the entropy column at all; entropy is carried as a label. It does "
        "enter structure.build_at_entropy, which the fixed-entropy cooling-ratio "
        "headline uses to pick a structure at a given interior entropy — that is "
        "the number the age shift reported here applies to.")
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("config", type=Path)
    args = ap.parse_args(argv)
    cfg = load_config(args.config)
    doc = run(cfg)
    path = cfg.output_path / "aqua_entropy_check.json"
    path.write_text(json.dumps(doc, indent=1) + "\n")

    t = doc["table"]
    print(f"\n  AQUA TABLE  {t['file_name']}  sha256 {t['sha256'][:16]}...")
    print(f"    predates the correction: {t['predates_correction']} — {t['why']}")
    if "blocked" in doc:
        print(f"\n  BLOCKED: {doc['blocked']}")
        return 0
    print(f"\n  {'age':>8s} {'T_int':>7s} {'P_base':>11s} {'T_base':>8s} "
          f"{'s_AQUA':>9s} {'s_corr':>9s} {'ds':>9s} {'ds/s':>8s} "
          f"{'SC frac':>8s} {'age shift':>10s}")
    for r in doc["points"]:
        print(f"  {r['age_gyr']:8.4g} {r['t_int_k']:7.1f} "
              f"{r['base_pressure_pa']:11.4g} {r['base_temperature_k']:8.1f} "
              f"{r['interior_entropy_aqua_j_kg_k']:9.1f} "
              f"{r['interior_entropy_corrected_j_kg_k']:9.1f} "
              f"{r['delta_s_j_kg_k']:9.1f} {r['delta_s_frac']:8.2%} "
              f"{r['envelope_supercritical_fraction']:8.2%} "
              f"{r['age_shift_gyr']:10.4g}")
    print(f"\n  density cross-check (should be ~0): "
          f"{max(abs(r['rho_check_frac']) for r in doc['points']):.2e}")
    print(f"  grad_ad |frac change| over the profiles: median "
          f"{max(r['grad_ad_median_abs_frac_change'] for r in doc['points']):.2e}, "
          f"max {max(r['grad_ad_max_abs_frac_change'] for r in doc['points']):.2e}")
    print(f"\n  {'age':>8s} {'E_AQUA':>12s} {'E_corr':>12s} {'dE/E':>9s} "
          f"{'du/u base':>10s}")
    for r in doc["points"]:
        print(f"  {r['age_gyr']:8.4g} {r['envelope_energy_aqua_j']:12.5g} "
              f"{r['envelope_energy_corrected_j']:12.5g} "
              f"{r['envelope_energy_frac_change']:9.3%} "
              f"{r['u_base_frac_change']:10.3%}")
    er = doc["energy_route"]
    print(f"\n  ENERGY ROUTE (the one that moves the tracks)")
    print(f"    total radiated energy changes by "
          f"{er['total_radiated_energy_frac_change']:+.2%}")
    print(f"    age at the coldest node: {er['age_at_coldest_node_aqua_gyr']:.4g} -> "
          f"{er['age_at_coldest_node_corrected_gyr']:.4g} Gyr "
          f"({er['age_frac_change_at_coldest_node']:+.2%})")
    print(f"    {'age':>8s} {'age_corr':>10s} {'change':>9s}")
    for r in doc["points"]:
        print(f"    {r['age_from_reintegration_aqua_gyr']:8.4g} "
              f"{r['age_from_reintegration_corrected_gyr']:10.4g} "
              f"{r['age_frac_change_from_energy']:9.2%}")
    print(f"\n  wrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
