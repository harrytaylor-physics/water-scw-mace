"""S3: band truncation of the HITEMP means and its effect on the G5 anchor (A39).

    python figures/plot_s3_band_truncation.py config/fiducial.yaml

SPEC §3 requires the absolute band-truncation systematic to be quantified by
comparing full-band with common-band means for the line-list source, where the
full band is available. HITEMP's spectra span [10, 11000] cm⁻¹ (the full band is
read from the files, never extrapolated); the science band is [100, 4000] cm⁻¹
(A25). This produces

  (a) κ_R(full)/κ_R(common) over HITEMP's (T, ρ) grid, and
  (b) the G5 residuals against Aguichine et al. (2021) recomputed with full-band
      HITEMP tables at the same eight points, alongside the common-band ones,

and writes ``outputs/g5_residuals_fullband.json``. Diagnostic only: the default
band is unchanged, and the FP-vs-line-list differential stays common-band (A21).

Because the Rosseland mean is harmonic it is dominated by transparency windows,
so the relative opacity floor (``opacity.rel_floor``) binds on the wide band
where it does not on the science band; both are reported.
"""

from __future__ import annotations

import csv
import json
import sys
from dataclasses import replace
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # run as script or -m

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm

from eos import AquaEOS
from figures._common import load_cfg, save, set_style
from opacity import compute_point_means, detect_common_band, resolve_band
from validation.test_g5_literature import _REF, our_radius

T_INT_G5 = 30.0  # the G5 gate's mature interior temperature
SOURCE = "hitemp"  # the only line list with full radiative-zone coverage (A33)


def _load_anchor_points():
    """The eight digitized Aguichine (2021) Fig. 8 points used by gate G5 (A32)."""
    with _REF.open() as fh:
        return [r for r in csv.DictReader(l for l in fh if not l.startswith("#"))]


def _means(cfg, band, rel_floor):
    t, rho, kr, kp = compute_point_means(SOURCE, cfg, band, rel_floor=rel_floor)
    order = np.lexsort((rho, t))
    return t[order], rho[order], kr[order], kp[order]


def _grid(t, rho, values):
    """Scatter (T, ρ) point values onto the source's tensor grid; NaN if absent."""
    t_ax, rho_ax = np.unique(t), np.unique(rho)
    out = np.full((t_ax.size, rho_ax.size), np.nan)
    for ti, ri, v in zip(t, rho, values):
        out[np.searchsorted(t_ax, ti), np.searchsorted(rho_ax, ri)] = v
    return t_ax, rho_ax, out


def g5_residuals(cfg, band, eos):
    """G5 residuals vs Aguichine (2021) with HITEMP built on ``band``.

    Reuses the gate's own matching (WMF→f_env, T_irr→T_eq=T_irr/√2, cmf=0,
    T_int=30 K) by pinning the config band, so both columns come from one path.
    """
    cfg_band = replace(cfg, opacity=replace(cfg.opacity, band_cm=tuple(band)))
    rows = []
    for r in _load_anchor_points():
        wmf, t_irr = float(r["wmf"]), float(r["t_irr_k"])
        m, r_ref = float(r["m_earth"]), float(r["r_earth"])
        r_ours = our_radius(cfg_band, m, wmf, t_irr, T_INT_G5, SOURCE, eos)
        rows.append({"wmf": wmf, "t_irr_k": t_irr, "m_earth": m,
                     "r_aguichine": r_ref, "r_ours": r_ours, "rel": (r_ours - r_ref) / r_ref})
    return rows


def main(argv=None):
    cfg = load_cfg(argv if argv is not None else sys.argv[1:])
    set_style()
    eos = AquaEOS.from_config(cfg)

    band_common = resolve_band(cfg, [SOURCE])  # science band (config-pinned, A25)
    band_full = detect_common_band(cfg, [SOURCE])  # HITEMP's own computed range
    floor = cfg.opacity.rel_floor

    t, rho, kr_c, kp_c = _means(cfg, band_common, floor)
    _, _, kr_f, kp_f = _means(cfg, band_full, floor)
    # Floor-free diagnostic: the harmonic Rosseland mean picks up far-band
    # transparency windows, so rel_floor binds on the wide band (A39).
    _, _, kr_c0, _ = _means(cfg, band_common, 0.0)
    _, _, kr_f0, _ = _means(cfg, band_full, 0.0)

    ratio_r = kr_f / kr_c
    ratio_p = kp_f / kp_c
    ratio_r0 = kr_f0 / kr_c0

    res_common = g5_residuals(cfg, band_common, eos)
    res_full = g5_residuals(cfg, band_full, eos)
    med_common = float(np.median([abs(x["rel"]) for x in res_common]))
    med_full = float(np.median([abs(x["rel"]) for x in res_full]))

    # ---------------- figure ----------------
    fig, (ax_k, ax_g) = plt.subplots(1, 2, figsize=(13.0, 5.4), constrained_layout=True)

    t_ax, rho_ax, grid = _grid(t, rho, ratio_r)
    # The ratio spans ~1.7 decades and never exceeds 1 (truncation biases κ_R
    # high), so a log-scaled sequential map reads better than a diverging one.
    norm = LogNorm(vmin=np.nanmin(grid), vmax=max(1.0, np.nanmax(grid)))
    mesh = ax_k.pcolormesh(np.arange(rho_ax.size + 1), np.arange(t_ax.size + 1), grid,
                           cmap="viridis", norm=norm, edgecolors="white", linewidth=0.5)
    for i in range(t_ax.size):
        for j in range(rho_ax.size):
            v = grid[i, j]
            ax_k.text(j + 0.5, i + 0.5, "—" if np.isnan(v) else f"{v:.2f}",
                      ha="center", va="center", fontsize=8,
                      color="black" if np.isnan(v) or v > 0.3 else "white")
    ax_k.set_xticks(np.arange(rho_ax.size) + 0.5, [f"{x:g}" for x in rho_ax])
    ax_k.set_yticks(np.arange(t_ax.size) + 0.5, [f"{x:g}" for x in t_ax])
    ax_k.set_xlabel(r"$\rho$ [kg m$^{-3}$]")
    ax_k.set_ylabel("T [K]")
    ax_k.grid(False)
    ax_k.set_title(f"(a) $\\kappa_R$(full)/$\\kappa_R$(common), {SOURCE.upper()}\n"
                   f"full [{band_full[0]:.0f}, {band_full[1]:.0f}] vs common "
                   f"[{band_common[0]:.0f}, {band_common[1]:.0f}] cm⁻¹, "
                   f"rel_floor={floor:g}", fontsize=9.5)
    fig.colorbar(mesh, ax=ax_k, label=r"$\kappa_R$(full) / $\kappa_R$(common)  (1 = no truncation)")

    x = np.arange(len(res_common))
    labels = [f"{r['m_earth']:.0f} $M_⊕$\nwmf {r['wmf']:.1f}\n{r['t_irr_k']:.0f} K"
              for r in res_common]
    ax_g.bar(x - 0.2, [r["rel"] * 100 for r in res_common], width=0.4,
             color="#0072B2", label=f"common band (median |rel| = {med_common:.1%})")
    ax_g.bar(x + 0.2, [r["rel"] * 100 for r in res_full], width=0.4,
             color="#D55E00", label=f"full band (median |rel| = {med_full:.1%})")
    ax_g.axhline(0.0, color="k", lw=0.8)
    ax_g.axhspan(-10, 10, color="grey", alpha=0.15, label="SPEC §7 G5 aim (±10%)")
    ax_g.set_xticks(x, labels, fontsize=7.5)
    ax_g.set_ylabel(r"$(R_{ours} - R_{Aguichine})\,/\,R_{Aguichine}$ [%]")
    ax_g.set_title("(b) G5 residuals vs Aguichine et al. (2021), Fig. 8 points\n"
                   f"{SOURCE.upper()} tables on each band, $T_{{int}}$={T_INT_G5:.0f} K, "
                   "cmf=0", fontsize=9.5)
    ax_g.legend(fontsize=8)

    fig.suptitle("S3 — band truncation of the line-list means (diagnostic; the default "
                 "band is unchanged)", fontsize=10.5)
    path = save(fig, cfg, "S3_band_truncation")

    # ---------------- outputs ----------------
    out_dir = cfg.output_path
    out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = out_dir / "S3_band_truncation.csv"
    with csv_path.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["temperature_k", "rho_kg_m3", "kappa_R_common_cm2_g", "kappa_R_full_cm2_g",
                    "kappa_R_ratio", "kappa_P_common_cm2_g", "kappa_P_full_cm2_g",
                    "kappa_P_ratio", "kappa_R_ratio_no_floor"])
        for i in range(t.size):
            w.writerow([f"{t[i]:g}", f"{rho[i]:g}", f"{kr_c[i]:.6e}", f"{kr_f[i]:.6e}",
                        f"{ratio_r[i]:.6f}", f"{kp_c[i]:.6e}", f"{kp_f[i]:.6e}",
                        f"{ratio_p[i]:.6f}", f"{ratio_r0[i]:.6f}"])

    json_path = out_dir / "g5_residuals_fullband.json"
    json_path.write_text(json.dumps({
        "anchor": "Aguichine et al. 2021 ApJ 914 84, Fig 8 (digitized); fit Eq. 29",
        "source": SOURCE, "t_int_k": T_INT_G5, "cmf": 0.0,
        "note": "Diagnostic band study (A39); the default science band is unchanged.",
        "band_common_cm": list(band_common), "band_full_cm": list(band_full),
        "rel_floor": floor,
        "median_abs_rel_common": med_common, "median_abs_rel_full": med_full,
        "kappa_R_ratio_full_over_common": {
            "min": float(ratio_r.min()), "median": float(np.median(ratio_r)),
            "max": float(ratio_r.max()),
        },
        "residuals_common": res_common, "residuals_full": res_full,
    }, indent=2))

    print(f"wrote {path}")
    print(f"wrote {csv_path}")
    print(f"wrote {json_path}")
    print(f"\nS3 — band truncation, source={SOURCE}, {t.size} grid points")
    print(f"  common band = [{band_common[0]:.12g}, {band_common[1]:.12g}] cm^-1 (science, A25)")
    print(f"  full band   = [{band_full[0]:.12g}, {band_full[1]:.12g}] cm^-1 "
          f"(HITEMP's own computed range)")
    print(f"  kappa_R(full)/kappa_R(common): min {ratio_r.min():.3f}, "
          f"median {np.median(ratio_r):.3f}, max {ratio_r.max():.3f}  "
          f"[rel_floor={floor:g}]")
    print(f"  same ratio with rel_floor=0:   min {ratio_r0.min():.3f}, "
          f"median {np.median(ratio_r0):.3f}, max {ratio_r0.max():.3f}")
    print(f"  kappa_P(full)/kappa_P(common): min {ratio_p.min():.3f}, "
          f"median {np.median(ratio_p):.3f}, max {ratio_p.max():.3f}")
    print(f"\n  G5 median |rel|: common band {med_common:.4f}, full band {med_full:.4f}")
    head = f"  {'wmf':>4s} {'T_irr':>6s} {'M [Me]':>7s} {'R_ref':>6s} " \
           f"{'R common':>9s} {'rel':>8s} {'R full':>8s} {'rel':>8s}"
    print(head)
    for a, b in zip(res_common, res_full):
        print(f"  {a['wmf']:>4.1f} {a['t_irr_k']:>6.0f} {a['m_earth']:>7.1f} "
              f"{a['r_aguichine']:>6.2f} {a['r_ours']:>9.4f} {a['rel'] * 100:>7.2f}% "
              f"{b['r_ours']:>8.4f} {b['rel'] * 100:>7.2f}%")
    return path


if __name__ == "__main__":
    main()
