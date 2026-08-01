"""Fx: direct opacity comparison FP vs line-list on the shared grid (SPEC §3).

    python figures/plot_fx_opacity_compare.py config/fiducial.yaml

This is the *exact, unclamped* FP-vs-line-list comparison available today:
band-limited κ_R and κ_P evaluated at the 27 shared (T, ρ) grid points on the
identical common band, so the differential is exact (SPEC §3). Unlike F1/F2 it
does not integrate a structure, so it is not affected by the FP/POKAZATEL grid
missing the radiative zone. Compares whatever line lists are available
(POKAZATEL now; HITEMP when delivered) to FP.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # run as script or -m

import matplotlib.pyplot as plt
import numpy as np

from figures._common import (
    SOURCE_COLOR,
    SOURCE_LABEL,
    available_sources,
    band_caption,
    comparison_band,
    load_cfg,
    save,
    set_style,
)
from opacity import compute_point_means


def _dump_csv(cfg, band, rows) -> Path:
    """Write the ratio table as a clean, DictReader-parseable CSV (header first).

    The band is carried in dedicated columns on every row (not comment lines),
    so the file is self-describing and quotable directly in the talk.
    """
    out = cfg.output_path / "Fx_opacity_compare.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    fields = ["line_list", "T_K", "rho_kg_m3", "band_cm_lo", "band_cm_hi",
              "kappaR_FP_cm2_g", "kappaR_LL_cm2_g", "kappaR_ratio",
              "kappaP_FP_cm2_g", "kappaP_LL_cm2_g", "kappaP_ratio"]
    with out.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(fields)
        for r in rows:
            w.writerow([r["ll"], f"{r['T']:.6g}", f"{r['rho']:.6g}",
                        f"{band[0]:.6g}", f"{band[1]:.6g}",
                        f"{r['krF']:.6g}", f"{r['krL']:.6g}", f"{r['ratio_r']:.6g}",
                        f"{r['kpF']:.6g}", f"{r['kpL']:.6g}", f"{r['ratio_p']:.6g}"])
    return out


def main(argv=None):
    cfg = load_cfg(argv if argv is not None else sys.argv[1:])
    set_style()

    line_lists = available_sources(cfg, ["hitemp", "pokazatel"])
    if not line_lists:
        print("no line-list source available yet; nothing to compare")
        return None
    # One exact band shared by FP and every line list plotted (A21).
    band = comparison_band(cfg, ["fp", *line_lists])

    # FP means on the shared grid (keyed by rounded (T, rho) for matching).
    t_fp, rho_fp, kr_fp, kp_fp = compute_point_means("fp", cfg, band)
    fp_key = {(round(t, 3), round(r, 3)): i for i, (t, r) in enumerate(zip(t_fp, rho_fp))}

    fig, (ax_r, ax_p) = plt.subplots(1, 2, figsize=(11.0, 4.8), constrained_layout=True)
    csv_rows = []

    for name in line_lists:
        t_ll, rho_ll, kr_ll, kp_ll = compute_point_means(name, cfg, band)
        idx = [fp_key.get((round(t, 3), round(r, 3))) for t, r in zip(t_ll, rho_ll)]
        keep = [i for i, j in enumerate(idx) if j is not None]
        if not keep:
            continue
        jj = [idx[i] for i in keep]
        t_shared, rho_shared = t_ll[keep], rho_ll[keep]
        ratio_r = kr_fp[jj] / kr_ll[keep]  # FP / line-list, Rosseland
        ratio_p = kp_fp[jj] / kp_ll[keep]  # FP / line-list, Planck
        c = SOURCE_COLOR[name]
        ax_r.scatter(t_shared, ratio_r, color=c, label=f"FP / {SOURCE_LABEL[name]}", s=36)
        ax_p.scatter(t_shared, ratio_p, color=c, label=f"FP / {SOURCE_LABEL[name]}", s=36)
        for k, i in enumerate(keep):
            # compute_point_means returns cm^2/g already (SI conversion happens
            # later in build_source), so these are quoted directly in cm^2/g.
            csv_rows.append({
                "ll": name, "T": t_shared[k], "rho": rho_shared[k],
                "krF": kr_fp[jj[k]], "krL": kr_ll[i], "ratio_r": ratio_r[k],
                "kpF": kp_fp[jj[k]], "kpL": kp_ll[i], "ratio_p": ratio_p[k],
            })

    for ax, lab in ((ax_r, "Rosseland"), (ax_p, "Planck")):
        ax.axhline(1.0, color="k", lw=0.8, ls=":")
        ax.set_xlabel("Temperature [K]")
        ax.set_ylabel(f"$\\kappa_{{{lab[0]}}}$(FP) / $\\kappa_{{{lab[0]}}}$(line list)")
        ax.set_yscale("log")
        ax.set_title(f"{lab} mean")
        ax.legend(fontsize=9)

    fig.suptitle(f"FP vs line-list opacity on the shared grid — {band_caption(band)} "
                 f"(exact differential, unclamped)")
    csv_rows.sort(key=lambda r: (r["ll"], r["T"], r["rho"]))
    csv_path = _dump_csv(cfg, band, csv_rows)
    path = save(fig, cfg, "Fx_opacity_compare")
    print(f"wrote {path}")
    print(f"wrote {csv_path}  ({len(csv_rows)} rows)")
    return path


if __name__ == "__main__":
    main()
