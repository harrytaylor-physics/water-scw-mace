"""S2: opacity-perturbation robustness of the structure and cooling history (A37).

    python figures/plot_s2_opacity_scaling.py config/fiducial.yaml

The FP table carries its own uncertainty, so the structural and evolutionary
conclusions must survive a constant perturbation of it. ``opacity.fp_scale``
multiplies the FP table's κ before blending; this repeats

  (a) the cold fiducial structure (T_int = 20 K, T_eq = 500 K), and
  (b) the E2 cooling-history comparison (composite vs HITEMP, matched initial
      conditions, E2's T_int grid),

for fp_scale ∈ {0.8, 1.0, 1.25}, and reports the RCB count and densities, the
interior-entropy difference between the composite and HITEMP tracks at 7 Gyr,
and the radius difference at 7 Gyr. The HITEMP track carries no FP, so it is
built once and reused as the common reference.

Outputs ``outputs/S2_opacity_scaling.csv`` and ``outputs/S2_opacity_scaling.png``.
"""

from __future__ import annotations

import csv
import sys
from dataclasses import replace
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # run as script or -m

import matplotlib.pyplot as plt
import numpy as np

from eos import AquaEOS
from evolution import cooling_track
from figures._common import band_caption, comparison_band, load_cfg, save, set_style
from structure import build_structure, deep_radiative_zones

FP_SCALES = (0.8, 1.0, 1.25)
T_INT = 20.0  # cold fiducial (S1/F1)
N_SAMPLE_STRUCT = 1500
DEEP_RHO_SI = 300.0  # closure threshold (A35)
RHO_XLIM = (1.0, 3.0e3)
AGE_REF_GYR = 7.0  # the manuscript's mature-age reference
T_INT_GRID_K = np.geomspace(12.0, 450.0, 44)  # E2's grid (A29), so fp_scale=1 reproduces E2
N_SAMPLE_TRACK = 400
SCALE_COLOR = {0.8: "#0072B2", 1.0: "black", 1.25: "#009E73"}
HITEMP_C = "#D55E00"


def _at_age(track, values, age_gyr: float) -> float:
    """Interpolate a track quantity at an age (tracks are ordered by T_int).

    Refuses to extrapolate: ``np.interp`` would silently clamp to the endpoint,
    which would quietly turn "the difference at 7 Gyr" into "the difference at
    the end of the track".
    """
    idx = np.argsort(track.age_gyr)
    ages = track.age_gyr[idx]
    if not ages[0] <= age_gyr <= ages[-1]:
        raise ValueError(
            f"age {age_gyr} Gyr is outside the {track.source} track "
            f"[{ages[0]:.3g}, {ages[-1]:.3g}] Gyr; extend the T_int grid"
        )
    return float(np.interp(age_gyr, ages, np.asarray(values, float)[idx]))


def run_scan(
    cfg,
    eos=None,
    scales=FP_SCALES,
    n_sample: int = N_SAMPLE_STRUCT,
    t_int_grid_k=None,
    n_sample_track: int = N_SAMPLE_TRACK,
    age_ref_gyr: float = AGE_REF_GYR,
    with_tracks: bool = True,
) -> dict:
    """Structures (and optionally cooling tracks) for every fp_scale.

    Returns ``{"rows", "structures", "hitemp", "tracks", "hitemp_track"}``; the
    gate calls it with coarse grids, the figure with the E2 grid.
    """
    eos = eos or AquaEOS.from_config(cfg)
    t_eq = float(cfg.irradiation.t_eq_k[0])
    band = comparison_band(cfg, ["composite", "hitemp"])
    t_int_grid_k = T_INT_GRID_K if t_int_grid_k is None else np.asarray(t_int_grid_k, float)

    hitemp = build_structure(cfg, t_eq, T_INT, "hitemp", eos=eos, n_sample=n_sample, band=band)
    hitemp_track = None
    if with_tracks:  # FP-free, so it is the same reference for every scale
        hitemp_track = cooling_track(cfg, "hitemp", t_eq, eos=eos, t_int_grid_k=t_int_grid_k,
                                     band=band, n_sample=n_sample_track)

    rows, structures, tracks = [], {}, {}
    for scale in scales:
        c = replace(cfg, opacity=replace(cfg.opacity, fp_scale=scale))
        res = build_structure(c, t_eq, T_INT, "composite", eos=eos, n_sample=n_sample, band=band)
        structures[scale] = res
        zones = deep_radiative_zones(res, DEEP_RHO_SI)
        row = {
            "fp_scale": scale,
            "n_rcb": len(res.rcb_densities_si),
            "rcb_rho_kg_m3": ";".join(f"{r:.1f}" for r in res.rcb_densities_si),
            "max_radiative_rho_kg_m3": f"{res.max_radiative_rho_si:.1f}",
            "n_deep_radiative_zones": len(zones),
            "deep_radiative_zones_kg_m3": ";".join(f"{a:.0f}-{b:.0f}" for a, b in zones),
            "radius_p_rearth": f"{res.radius_p_rearth:.5f}",
            "clamp_fraction_radiative": f"{res.clamp_fraction_radiative:.4f}",
            "n_sample": n_sample,
        }
        if with_tracks:
            tr = cooling_track(c, "composite", t_eq, eos=eos, t_int_grid_k=t_int_grid_k,
                               band=band, n_sample=n_sample_track)
            tracks[scale] = tr
            s_c = _at_age(tr, tr.entropy_j_kg_k, age_ref_gyr)
            s_h = _at_age(hitemp_track, hitemp_track.entropy_j_kg_k, age_ref_gyr)
            r_c = _at_age(tr, tr.radius_p_rearth, age_ref_gyr)
            r_h = _at_age(hitemp_track, hitemp_track.radius_p_rearth, age_ref_gyr)
            row.update({
                "age_ref_gyr": age_ref_gyr,
                "entropy_composite_j_kg_k": f"{s_c:.2f}",
                "entropy_hitemp_j_kg_k": f"{s_h:.2f}",
                "delta_entropy_j_kg_k": f"{s_c - s_h:+.3f}",
                "delta_entropy_frac": f"{(s_c - s_h) / s_h:+.3e}",
                "radius_composite_rearth": f"{r_c:.5f}",
                "radius_hitemp_rearth": f"{r_h:.5f}",
                "delta_radius_rearth": f"{r_c - r_h:+.5f}",
                "delta_radius_frac": f"{(r_c - r_h) / r_h:+.3e}",
                "track_age_max_gyr": f"{tr.age_gyr.max():.3f}",
                "energy_residual": f"{tr.energy_residual:.2e}",
            })
        rows.append(row)
    return {"rows": rows, "structures": structures, "hitemp": hitemp,
            "tracks": tracks, "hitemp_track": hitemp_track}


def _diff_vs_age(track_a, track_b, n: int = 150):
    """(age, ΔR/R, Δs/s) for composite (a) minus HITEMP (b) on a shared grid."""
    lo = max(track_a.age_gyr[track_a.age_gyr > 0].min(), track_b.age_gyr[track_b.age_gyr > 0].min())
    hi = min(track_a.age_gyr.max(), track_b.age_gyr.max())
    age = np.geomspace(lo, hi, n)
    ra = np.array([_at_age(track_a, track_a.radius_p_rearth, x) for x in age])
    rb = np.array([_at_age(track_b, track_b.radius_p_rearth, x) for x in age])
    sa = np.array([_at_age(track_a, track_a.entropy_j_kg_k, x) for x in age])
    sb = np.array([_at_age(track_b, track_b.entropy_j_kg_k, x) for x in age])
    return age, (ra - rb) / rb, (sa - sb) / sb


def main(argv=None):
    cfg = load_cfg(argv if argv is not None else sys.argv[1:])
    set_style()
    eos = AquaEOS.from_config(cfg)
    t_eq = float(cfg.irradiation.t_eq_k[0])
    band = comparison_band(cfg, ["composite", "hitemp"])

    scan = run_scan(cfg, eos=eos)
    rows, structures, hitemp = scan["rows"], scan["structures"], scan["hitemp"]

    fig, (ax_s, ax_h) = plt.subplots(1, 2, figsize=(12.4, 5.4), constrained_layout=True)

    # -- (a) cold-fiducial structure --
    ph = hitemp.profile
    ax_s.plot(ph["rho_si"], ph["grad_rad"] / ph["grad_ad"], color=HITEMP_C, ls="--", lw=1.6,
              label="HITEMP")
    for scale in FP_SCALES:
        p = structures[scale].profile
        ax_s.plot(p["rho_si"], p["grad_rad"] / p["grad_ad"], color=SCALE_COLOR[scale], lw=2.0,
                  label=f"composite, fp_scale={scale:g}")
    lo, hi = cfg.opacity.blend_window_kg_m3
    ax_s.axvspan(lo, hi, color="#0072B2", alpha=0.10)
    ax_s.axvline(DEEP_RHO_SI, color="grey", lw=1.0, ls=":")
    ax_s.axhline(1.0, color="grey", lw=1.0, ls=":")
    ax_s.set_xscale("log")
    ax_s.set_yscale("log")
    ax_s.set_xlim(RHO_XLIM)
    ax_s.set_xlabel(r"$\rho$ [kg m$^{-3}$]")
    ax_s.set_ylabel(r"$\nabla_{rad}\,/\,\nabla_{ad}$  (< 1: radiative)")
    ax_s.set_title(f"(a) structure at $T_{{int}}$={T_INT:.0f} K "
                   f"(shaded: blend window; dotted: ρ={DEEP_RHO_SI:.0f})", fontsize=9.5)
    ax_s.legend(fontsize=8.5, loc="lower right")

    # -- (b) cooling history, composite − HITEMP. Δs/s and ΔR/R differ by ~25×,
    #    so ΔR/R gets its own axis rather than being flattened onto zero. --
    ax_r = ax_h.twinx()
    ax_r.grid(False)
    for scale in FP_SCALES:
        age, dr, ds = _diff_vs_age(scan["tracks"][scale], scan["hitemp_track"])
        ax_h.plot(age, ds * 100, color=SCALE_COLOR[scale], lw=2.0,
                  label=f"$\\Delta s/s$, fp_scale={scale:g}  (left)")
        ax_r.plot(age, dr * 100, color=SCALE_COLOR[scale], lw=1.4, ls=":",
                  label=f"$\\Delta R_p/R_p$, fp_scale={scale:g}  (right)")
    ax_h.axvline(AGE_REF_GYR, color="grey", lw=1.0, ls="--")
    ax_h.axhline(0.0, color="k", lw=0.6)
    ax_h.set_xscale("log")
    ax_h.set_xlabel("Age [Gyr]")
    ax_h.set_ylabel(r"interior $\Delta s/s$, composite − HITEMP [%]")
    ax_r.set_ylabel(r"$\Delta R_p/R_p$, composite − HITEMP [%]")
    ax_h.set_title(f"(b) cooling history vs HITEMP (dashed line: {AGE_REF_GYR:.0f} Gyr)",
                   fontsize=9.5)
    handles = ax_h.get_legend_handles_labels()[0] + ax_r.get_legend_handles_labels()[0]
    labels_h = ax_h.get_legend_handles_labels()[1] + ax_r.get_legend_handles_labels()[1]
    ax_h.legend(handles, labels_h, fontsize=7.5, ncol=2, loc="upper left")

    fig.suptitle(
        f"S2 — opacity-perturbation robustness (fp_scale × FP table κ): "
        f"$T_{{eq}}$={t_eq:.0f} K, $M_c$={cfg.core.mass_mearth:.0f} $M_⊕$, "
        f"$f_{{env}}$={cfg.envelope.f_env:.0%}, {band_caption(band)}", fontsize=10.5)
    path = save(fig, cfg, "S2_opacity_scaling")

    csv_path = cfg.output_path / "S2_opacity_scaling.csv"
    with csv_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    print(f"wrote {path}")
    print(f"wrote {csv_path}")
    print(f"\nS2 — opacity-perturbation robustness (T_int={T_INT:.0f} K structure, "
          f"E2 cooling grid, {band_caption(band)}):")
    head = (f"{'fp_scale':>8s} {'#RCB':>4s} {'RCB rho [kg/m3]':>18s} {'max rad rho':>11s} "
            f"{'deep zones':>10s} {'ds(7Gyr) [J/kg/K]':>17s} {'ds/s':>10s} "
            f"{'dR(7Gyr) [Re]':>13s} {'dR/R':>10s}")
    print(head)
    print("-" * len(head))
    for r in rows:
        print(f"{r['fp_scale']:>8g} {r['n_rcb']:>4d} {r['rcb_rho_kg_m3']:>18s} "
              f"{r['max_radiative_rho_kg_m3']:>11s} "
              f"{(r['deep_radiative_zones_kg_m3'] or 'none'):>10s} "
              f"{r['delta_entropy_j_kg_k']:>17s} "
              f"{float(r['delta_entropy_frac']) * 100:>9.4f}% "
              f"{r['delta_radius_rearth']:>13s} "
              f"{float(r['delta_radius_frac']) * 100:>9.4f}%")
    signs = {np.sign(float(r["delta_entropy_j_kg_k"])) for r in rows}
    print(f"\nHITEMP reference track: R_p({AGE_REF_GYR:.0f} Gyr) = "
          f"{rows[0]['radius_hitemp_rearth']} R⊕, s = {rows[0]['entropy_hitemp_j_kg_k']} J/kg/K")
    print(f"entropy bias sign at {AGE_REF_GYR:.0f} Gyr: "
          f"{'same for all scales' if len(signs) == 1 else 'DIFFERS ACROSS SCALES'} "
          f"({', '.join(r['delta_entropy_j_kg_k'] for r in rows)})")
    worst = max(int(r["n_deep_radiative_zones"]) for r in rows)
    print(f"deep radiative zones (rho > {DEEP_RHO_SI:.0f}) over all scales: {worst} "
          f"({'closure holds' if worst == 0 else 'CLOSURE FAILS'})")
    return path


if __name__ == "__main__":
    main()
