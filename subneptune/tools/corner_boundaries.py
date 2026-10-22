#!/usr/bin/env python3
"""The six boundaries the corner cells brought inside the hull (Task 44).

    python tools/corner_boundaries.py config/fiducial.yaml

Six of the 108 line-list radiative-convective boundaries fell between computed
points before A74 and were clamped to a 1000 K value. Two of them are among the
ten strongest and so already carry a flip threshold from ``plot_s8_rcb_flip.py``;
the other four have never been profiled, and a threshold needs grad_ad and
grad_rad at the point, which the scan CSV does not store. This profiles all six
by the same code path the S8 analysis uses — ``analyse_case`` imported, not
copied — so the six are comparable with the ten.

Writes ``outputs/corner_boundaries.json``.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import atomicio

import numpy as np

from config import load_config

BAND_CM = (100.0, 4000.0)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("config", type=Path)
    args = ap.parse_args(argv)
    cfg = load_config(args.config)

    import figures.plot_s4_domain_scan as s4
    from eos import AquaEOS
    from figures._common import comparison_band, fp_floors
    from figures.plot_s8_rcb_flip import SOURCES, _label, analyse_case
    from tools.corner_cells import triangle

    pre = json.loads((cfg.output_path / "pre52" / "state.json").read_text())
    want = {(round(b["t_eq_k"], 3), round(b["mass_core_mearth"], 3),
             round(b["wmf"], 3), round(b["t_int_k"], 3)): b
            for b in pre["between_computed_points"]}

    rows = [r for r in s4.load_rows(cfg.output_path / "S4_domain_scan_extended.csv")
            if not r["skipped"]]
    picked = []
    for r in rows:
        key = (round(float(r["t_eq_k"]), 3), round(float(r["mass_core_mearth"]), 3),
               round(float(r["wmf"]), 3), round(float(r["t_int_k"]), 3))
        if key in want:
            picked.append((r, want[key]))
    eos = AquaEOS.from_config(cfg)
    band = comparison_band(cfg, list(SOURCES))
    floors = fp_floors(cfg)

    out = {"meta": __doc__.split("\n\n")[0], "band_cm": list(band), "cases": []}
    for row, before in picked:
        case = analyse_case(cfg, eos, band, row, floors)
        # the deepest boundary, which is the one the scan reports and the one the
        # corner cells were computed for
        bs = [b for b in case["boundaries"] if np.isfinite(b["threshold_factor"])]
        if not bs:
            out["cases"].append({"case": _label(row),
                                 "skipped": "no boundary with a finite threshold"})
            continue
        a = max(bs, key=lambda b: b["rho_si"])
        a["case"] = _label(row)
        t_b, rho_b = float(a["t_k"]), float(a["rho_si"])
        rec = {
            "case": a["case"],
            "T_K": t_b, "rho_kg_m3": rho_b,
            "before": {
                "T_K": before["T"], "rho_kg_m3": before["rho"],
                "kappa_R_fp_cm2_g": before["kappa_R_cm2_g"],
                "provenance": before["provenance"],
                "clamped_to": before.get("clamped_to"),
            },
            "after": {
                "kappa_R_fp_cm2_g": float(a["kappa_r_fp_cm2_g"]),
                "kappa_R_hitemp_cm2_g": float(a["kappa_r_hitemp_cm2_g"]),
                "fp_over_hitemp": float(a["fp_over_hitemp"]),
                "provenance": "clamped" if a["fp_clamped_at_rcb"] else "interpolated",
                "triangle": triangle(cfg, t_b, rho_b),
            },
            "threshold_factor": float(a["threshold_factor"]),
            "exceeds_threshold": bool(float(a["fp_over_hitemp"])
                                      > float(a["threshold_factor"])),
            "d_log10_rho": float(a["d_log10_rho"]),
            "d_enclosed_mass_fraction": float(a["d_enclosed_mass_fraction"]),
        }
        # the ratio the clamp would have given, at the same boundary
        rec["before"]["fp_over_hitemp"] = (before["kappa_R_cm2_g"]
                                           / rec["after"]["kappa_R_hitemp_cm2_g"])
        rec["ratio_change"] = (rec["after"]["fp_over_hitemp"]
                               / rec["before"]["fp_over_hitemp"])
        out["cases"].append(rec)

    atomicio.write_text(cfg.output_path / "corner_boundaries.json", 
        json.dumps(out, indent=1) + "\n")

    print(f"\n  The six boundaries the corner cells brought inside "
          f"({len(out['cases'])} profiled)\n")
    hdr = (f"  {'case':<36s} {'T':>6s} {'rho':>7s} {'FP/LL before':>13s} "
           f"{'after':>8s} {'thr':>6s} {'>thr':>5s}")
    print(hdr); print("  " + "-" * (len(hdr) - 2))
    for c in out["cases"]:
        print(f"  {c['case']:<36s} {c['T_K']:6.1f} {c['rho_kg_m3']:7.2f} "
              f"{c['before']['fp_over_hitemp']:13.2f} "
              f"{c['after']['fp_over_hitemp']:8.2f} {c['threshold_factor']:6.3f} "
              f"{'yes' if c['exceeds_threshold'] else 'NO':>5s}")
    n = sum(1 for c in out["cases"] if c["exceeds_threshold"])
    print(f"\n  all six now interpolated: "
          f"{all(c['after']['provenance'] == 'interpolated' for c in out['cases'])}")
    print(f"  still exceeding the flip threshold: {n} of {len(out['cases'])}")
    print(f"\n  wrote {cfg.output_path / 'corner_boundaries.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
