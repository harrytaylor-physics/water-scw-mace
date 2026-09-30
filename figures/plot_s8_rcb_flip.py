"""S8: how far is each RCB from flipping, now that FP covers the boundary? (A55)

    python figures/plot_s8_rcb_flip.py config/fiducial.yaml

A49 asked how large a change in κ_R a radiative-convective boundary would even
notice, and answered it in the model's own terms: ∇_rad ∝ κ_R, so at fixed
structure the factor that flips a sample is exactly ∇_ad/∇_rad there. Most RCBs
turned out to sit within 3-8% of flipping, against an FP/line-list dispute of
2-3× — but at that time the FP κ_R *at the boundary* was a nearest-neighbour
clamp from 1000 K, so the dispute was not measured where it mattered.

The extension campaign covers those boundaries (T ≥ 650 K, ρ ≥ 10), so this
redoes A49 with the local FP/HITEMP κ_R ratio **interpolated at the RCB**, and
says per case whether the composite actually removes or moves the boundary —
in density and in enclosed mass, not just in principle.

Cases are the ``n`` strongest from the extended domain scan, ranked by the
structural change the composite makes, |ΔR_p|/R_p. Run the scan first
(``plot_s4_domain_scan.py --extended``); this reads its CSV.

Writes ``outputs/S8_rcb_flip.csv`` and ``outputs/S8_rcb_flip.png``.
"""

from __future__ import annotations

import csv

import sys
from dataclasses import replace
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # run as script or -m

import atomicio

import matplotlib.pyplot as plt
import numpy as np

from eos import AquaEOS
from figures._common import band_caption, comparison_band, fp_floors, load_cfg, save, set_style
from opacity import get_source
from opacity.table import CM2_PER_G_TO_SI
from structure import build_structure
from structure.diagnostics import rcb_opacity_sensitivity

N_CASES = 10
N_SAMPLE = 1500
SOURCES = ("hitemp", "composite")


def _strength(row) -> float:
    """|ΔR_p|/R_p between the two sources — the structural change itself."""
    try:
        a = float(row["composite_radius_rearth"])
        b = float(row["hitemp_radius_rearth"])
    except (TypeError, ValueError):
        return 0.0
    return abs(a - b) / b if b else 0.0


def strongest_cases(csv_path: Path, n: int = N_CASES) -> list[dict]:
    """The ``n`` converged scan cells the composite changes most."""
    import figures.plot_s4_domain_scan as s4

    rows = s4.load_rows(csv_path)
    live = [r for r in rows if not r["skipped"]]
    for r in live:
        r["_strength"] = _strength(r)
    return sorted(live, key=lambda r: r["_strength"], reverse=True)[:n]


def _enclosed_mass_fraction(result, rho_si: float) -> float:
    """Envelope mass fraction enclosed *below* the first sample at ``rho_si``.

    The profile runs top -> core, so this is the fraction of envelope mass
    interior to the boundary: what "the boundary moved" costs in mass terms.
    """
    prof = result.profile
    rho = np.asarray(prof["rho_si"], float)
    m = np.asarray(prof["mass_kg"], float)
    i = int(np.argmin(np.abs(rho - rho_si)))
    total = abs(m[0] - m[-1])
    return float(abs(m[i] - m[-1]) / total) if total else float("nan")


def analyse_case(cfg, eos, band, row, floors, n_sample: int = N_SAMPLE) -> dict:
    """One scan cell: HITEMP RCBs, the FP ratio there, and what the composite does."""
    t_floor, rho_floor = floors
    cfg_c = replace(
        cfg,
        core=replace(cfg.core, mass_mearth=float(row["mass_core_mearth"])),
        envelope=replace(cfg.envelope, f_env=float(row["wmf"])),
    )
    res = {name: build_structure(cfg_c, float(row["t_eq_k"]), float(row["t_int_k"]),
                                 name, eos=eos, n_sample=n_sample, band=band)
           for name in SOURCES}
    ht, comp = res["hitemp"], res["composite"]
    ht_src = get_source("hitemp", cfg_c, band=band)
    fp_src = get_source("fp", cfg_c, band=band)
    sens = rcb_opacity_sensitivity(ht, ht_src)

    # One entry per radiative LAYER, not per flip: an interior layer has two
    # flips but only one deepest sample, and it is that sample's ∇_ad/∇_rad that
    # says how far the layer is from closing. The boundary quoted is the RCB
    # nearest the layer's dense edge.
    boundaries = []
    rcb_rho = np.array([r["rho_si"] for r in sens["rcbs"]], float)
    for layer in sens["layers"]:
        if rcb_rho.size:
            k = int(np.argmin(np.abs(np.log10(rcb_rho) - np.log10(layer["rho_max_si"]))))
            rcb = sens["rcbs"][k]
        else:  # fully radiative profile: quote the layer's own dense edge
            rcb = {"t_k": layer["base"]["t_k"], "rho_si": layer["rho_max_si"],
                   "p_pa": layer["base"]["p_pa"]}
        t_k, rho = rcb["t_k"], rcb["rho_si"]
        kr_ht = float(ht_src.kappa_R(t_k, rho)) / CM2_PER_G_TO_SI
        kr_fp = float(fp_src.kappa_R(t_k, rho)) / CM2_PER_G_TO_SI
        fp_clamped = bool(np.asarray(fp_src.clamped_mask(t_k, rho)).ravel()[0])
        # Threshold: how much κ_R must change for the deepest radiative sample of
        # the layer this boundary caps to flip (A49's factor_close).
        threshold = layer["base"]["factor_close"]
        # Did the composite REMOVE this layer, or merely move its boundary?
        # Tested directly — "nearest composite RCB" is meaningless when the layer
        # is gone, since the nearest one is then some unrelated boundary.
        comp_prof = comp.profile
        comp_rho_all = np.asarray(comp_prof["rho_si"], float)
        comp_radiative = ~np.asarray(comp_prof["is_convective"], bool)
        within = ((comp_rho_all >= layer["rho_min_si"])
                  & (comp_rho_all <= layer["rho_max_si"]))
        survives = bool((within & comp_radiative).any())
        comp_rho = np.array(comp.rcb_densities_si, float)
        if not survives:
            # The whole layer went convective: that IS the enclosed-mass change.
            removed = True
            moved_rho, d_log_rho = float("nan"), float("nan")
            d_mass = -layer["mass_fraction"]
        elif comp_rho.size:
            j = int(np.argmin(np.abs(np.log10(comp_rho) - np.log10(rho))))
            moved_rho = float(comp_rho[j])
            d_log_rho = float(np.log10(moved_rho) - np.log10(rho))
            d_mass = (_enclosed_mass_fraction(comp, moved_rho)
                      - _enclosed_mass_fraction(ht, rho))
            removed = False
        else:
            removed = False
            moved_rho, d_log_rho, d_mass = float("nan"), float("nan"), float("nan")
        boundaries.append({
            "t_k": t_k, "rho_si": rho, "p_pa": rcb["p_pa"],
            "layer_mass_fraction": layer["mass_fraction"],
            "kappa_r_hitemp_cm2_g": kr_ht, "kappa_r_fp_cm2_g": kr_fp,
            "fp_over_hitemp": kr_fp / kr_ht,
            "fp_clamped_at_rcb": fp_clamped,
            "in_fp_box": bool(t_k >= t_floor and rho >= rho_floor),
            "threshold_factor": threshold,
            "removed": removed,
            "layer_survives_in_composite": survives,
            "composite_rho_si": moved_rho,
            "d_log10_rho": d_log_rho,
            "d_enclosed_mass_fraction": d_mass,
        })
    return {
        "row": row, "hitemp": ht, "composite": comp, "boundaries": boundaries,
        "n_rcb_hitemp": len(ht.rcb_densities_si),
        "n_rcb_composite": len(comp.rcb_densities_si),
        "radius_frac_change": _strength(row),
        "clamp_fraction_all_composite": comp.clamp_fraction_all,
    }


def _label(row) -> str:
    return (f"M_c={row['mass_core_mearth']:.0f} wmf={row['wmf']:.2f} "
            f"T_eq={row['t_eq_k']:.0f} T_int={row['t_int_k']:.0f}")


def plot_cases(cfg, cases, band, tag="S8_rcb_flip"):
    """Flip threshold vs the FP/HITEMP dispute, per boundary."""
    fig, ax = plt.subplots(figsize=(8.2, 5.6), constrained_layout=True)
    for case in cases:
        for b in case["boundaries"]:
            if not np.isfinite(b["threshold_factor"]):
                continue
            color = "#0072B2" if b["in_fp_box"] and not b["fp_clamped_at_rcb"] else "#999999"
            marker = "o" if not b["removed"] else "*"
            ax.plot(b["threshold_factor"], b["fp_over_hitemp"], marker, ms=11 if
                    b["removed"] else 7, color=color, mec="k", mew=0.4, zorder=3)
    lim = ax.get_xlim()
    xs = np.geomspace(max(lim[0], 1.0001), max(lim[1], 2.0), 50)
    ax.plot(xs, xs, "k--", lw=1.0, zorder=2)
    ax.fill_between(xs, xs, 1e4, color="#009E73", alpha=0.10, zorder=1)
    ax.text(0.97, 0.05, "above the line: the FP/HITEMP dispute\nexceeds what the "
            "boundary needs to flip", transform=ax.transAxes, ha="right", va="bottom",
            fontsize=8.5)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("flip threshold — factor on $\\kappa_R$ that moves the boundary")
    ax.set_ylabel("$\\kappa_R$(FP) / $\\kappa_R$(HITEMP) at the boundary")
    ax.set_title(f"S8 — is the FP/line-list opacity dispute big enough to move an RCB?\n"
                 f"{len(cases)} strongest scan cells, {band_caption(band)}; blue = FP "
                 f"interpolated at the boundary, grey = clamped; ★ = composite removes it",
                 fontsize=9.5)
    return save(fig, cfg, tag)


def main(argv=None):
    import argparse

    argv = sys.argv[1:] if argv is None else argv
    ap = argparse.ArgumentParser()
    ap.add_argument("config", type=Path)
    ap.add_argument("--scan-csv", type=Path, default=None)
    ap.add_argument("--n-cases", type=int, default=N_CASES)
    ap.add_argument("--n-sample", type=int, default=N_SAMPLE)
    args = ap.parse_args(argv)
    cfg = load_cfg([str(args.config)])
    set_style()
    eos = AquaEOS.from_config(cfg)
    band = comparison_band(cfg, list(SOURCES))
    floors = fp_floors(cfg)
    csv_path = args.scan_csv or (cfg.output_path / "S4_domain_scan_extended.csv")
    if not csv_path.exists():
        print(f"no {csv_path}; run plot_s4_domain_scan.py --extended first")
        return None

    rows = strongest_cases(csv_path, args.n_cases)
    print(f"S8 — RCB flip analysis (A55) on the {len(rows)} strongest scan cells, "
          f"ranked by |ΔR_p|/R_p, {band_caption(band)}")
    print(f"  FP box from config: T ≥ {floors[0]:g} K, ρ ≥ {floors[1]:g} kg/m³")
    cases = [analyse_case(cfg, eos, band, r, floors, args.n_sample) for r in rows]

    head = (f"\n  {'case':<40s} {'ΔR_p/R_p':>9s} {'nRCB HT':>8s} {'nRCB C':>7s} "
            f"{'C clamp_all':>12s}")
    print(head)
    print("  " + "-" * (len(head) - 3))
    for c in cases:
        print(f"  {_label(c['row']):<40s} {c['radius_frac_change']:>9.2%} "
              f"{c['n_rcb_hitemp']:>8d} {c['n_rcb_composite']:>7d} "
              f"{c['clamp_fraction_all_composite']:>12.4f}")

    print(f"\n  per HITEMP RCB — the dispute measured AT the boundary, interpolated "
          f"where the FP grid now reaches:")
    head = (f"  {'case':<40s} {'T_rcb':>7s} {'rho':>8s} {'P [bar]':>9s} | "
            f"{'kR HT':>8s} {'kR FP':>9s} {'FP/HT':>7s} {'thresh':>7s} | "
            f"{'verdict':<26s} {'Δlog10ρ':>8s} {'Δm_enc':>8s}")
    print(head)
    print("  " + "-" * (len(head) - 3))
    rows_out = []
    for c in cases:
        if not c["boundaries"]:
            print(f"  {_label(c['row']):<40s} (no RCB)")
            continue
        for b in c["boundaries"]:
            if b["removed"]:
                verdict = "composite CLOSES the layer"
            elif np.isfinite(b["threshold_factor"]) and \
                    b["fp_over_hitemp"] > b["threshold_factor"]:
                verdict = "dispute > threshold"
            else:
                verdict = "dispute < threshold"
            if b["fp_clamped_at_rcb"]:
                verdict += " (FP CLAMPED)"
            print(f"  {_label(c['row']):<40s} {b['t_k']:>7.0f} {b['rho_si']:>8.2f} "
                  f"{b['p_pa'] / 1e5:>9.2f} | {b['kappa_r_hitemp_cm2_g']:>8.3f} "
                  f"{b['kappa_r_fp_cm2_g']:>9.3f} {b['fp_over_hitemp']:>7.2f} "
                  f"{b['threshold_factor']:>7.3f} | {verdict:<26s} "
                  f"{b['d_log10_rho']:>8.3f} {b['d_enclosed_mass_fraction']:>8.4f}")
            rows_out.append({"case": _label(c["row"]), **b})

    interp = [b for c in cases for b in c["boundaries"] if not b["fp_clamped_at_rcb"]]
    clamped = [b for c in cases for b in c["boundaries"] if b["fp_clamped_at_rcb"]]
    print(f"\n  {len(interp)} of {len(interp) + len(clamped)} boundaries now have FP "
          f"INTERPOLATED at the RCB (A49 had none — every one was a 1000 K clamp).")
    if interp:
        rr = [b["fp_over_hitemp"] for b in interp]
        th = [b["threshold_factor"] for b in interp if np.isfinite(b["threshold_factor"])]
        print(f"    FP/HITEMP κ_R at those boundaries: {min(rr):.2f}-{max(rr):.2f}×")
        if th:
            print(f"    flip thresholds there:              {min(th):.3f}-{max(th):.3f}×")
        beats = [b for b in interp if np.isfinite(b["threshold_factor"])
                 and b["fp_over_hitemp"] > b["threshold_factor"]]
        print(f"    boundaries where the dispute exceeds the threshold: "
              f"{len(beats)} of {len(th)}")

    out_dir = cfg.output_path
    out_dir.mkdir(parents=True, exist_ok=True)
    csv_out = out_dir / "S8_rcb_flip.csv"
    if rows_out:
        with atomicio.open_write(csv_out, newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows_out[0]))
            w.writeheader()
            w.writerows(rows_out)
    path = plot_cases(cfg, cases, band)
    print(f"\nwrote {path}")
    print(f"wrote {csv_out}")
    return path


if __name__ == "__main__":
    main()
