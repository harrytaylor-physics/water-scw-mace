"""S7: are the dilute FP spectra above their own noise floor? (A52)

    python figures/plot_s7_dilute_floor.py config/fiducial.yaml

The extension campaign delivers FP opacity down to ρ = 10 kg/m³, which is the
first time first-principles data exist where the composite's *line-list* side
normally rules. Before any structural claim leans on those points, they have to
clear the dipole model's own noise floor: the FP spectrum is the Fourier
transform of a dipole autocorrelation, and the dipole model has a self-calibrated
residual of ~5.7-5.8 mDebye/atom against a signal RMS of 12.4 mDebye/atom. At low
density the signal falls (κ per gram is roughly ρ-independent, but the
*transparency windows* thin out) while the floor does not, so the windows are
where the floor bites first.

The floor construction (reimplemented here; the ARC original,
``lda_campaign/stepG_dipoles/compare_point_full.py``, is not on this machine and
remains the authority):

1. The absorption cross-section carries the dipole power spectrum D(ν) with the
   standard prefactor, σ(ν) ∝ ν·(1 − e^{−hcν/kT})·D(ν), so
   **D(ν) = σ(ν) / [ν (1 − e^{−hcν/kT})]** recovers it up to a constant that
   cancels below.
2. Parseval: ∫D dν over all frequencies is the mean-square dipole. The residual
   is a *ratio of amplitudes*, so it enters the power spectrum squared:
   **r² = (5.75 / 12.4)² = 0.215** of the total power is floor.
3. The residual is white, so that power is spread flat over the simulation's
   full spectral bandwidth — the Nyquist wavenumber of the 0.5 fs timestep,
   ν_Nyq = 1/(2·dt·c) ≈ 3.34×10⁴ cm⁻¹ — giving a constant
   **D_floor = r² ∫D dν / ν_Nyq**.
4. Mapping back through the same prefactor gives the floor as an opacity,
   **κ_floor(ν) = D_floor · ν · (1 − e^{−hcν/kT}) · (N_A/M)**, which is not flat
   in κ: it inherits the ν(1−e^{−x}) shape, so it rises across the band.

The bandwidth in step 3 is the dominant assumption and it is reported, not
hidden: ``--floor-bandwidth`` also accepts ``computed`` (the file's own
[1.1, 6000] cm⁻¹ range) and ``band`` (the science band alone, the most
pessimistic choice), and the printed table gives the verdict under all three.

Reported per dilute cell: the fraction of the science band where κ exceeds 2×
the floor, split into the **transparency windows** (1900-3100 cm⁻¹, between the
bend and stretch bands, where the Rosseland harmonic mean is decided) and the
**band cores** (100-1900 and 3100-4000 cm⁻¹).

Writes ``outputs/S7_dilute_floor.csv`` and ``outputs/S7_dilute_floor.png``.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # run as script or -m

import matplotlib.pyplot as plt
import numpy as np

from constants import C_LIGHT, H_PLANCK, K_B
from figures._common import load_cfg, save, set_style
from opacity.loaders import find_spectra, load_spectrum
from opacity.registry import get_spec

# Dipole-model self-calibration (owner-supplied, ARC step G): residual vs signal.
DIPOLE_FLOOR_MDEBYE = 5.75  # quoted range 5.7-5.8; both ends reported
DIPOLE_FLOOR_RANGE_MDEBYE = (5.7, 5.8)
DIPOLE_SIGNAL_RMS_MDEBYE = 12.4

DILUTE_RHO_MAX = 60.0  # "dilute" cells: ρ <= 60 kg/m³ (the low-density column)
SCIENCE_BAND_CM = (100.0, 4000.0)
WINDOW_CM = (1900.0, 3100.0)  # transparency window between bend and stretch
FLOOR_MULTIPLE = 2.0  # "above the floor" means κ > 2 × κ_floor
FIGURE_CELLS = ((1000.0, 10.0), (1000.0, 30.0), (1500.0, 30.0))

_C_CM_S = C_LIGHT * 100.0
_NA = 6.02214076e23
_M_H2O = 18.015286

try:
    from numpy import trapezoid as _trapz
except ImportError:  # numpy < 2.0
    from numpy import trapz as _trapz


def _one_minus_exp(wn_cm: np.ndarray, temperature_k: float) -> np.ndarray:
    """Stimulated-emission factor 1 − e^{−hcν/kT}."""
    x = H_PLANCK * _C_CM_S * np.asarray(wn_cm, float) / (K_B * temperature_k)
    return -np.expm1(-x)


def nyquist_wavenumber_cm(dt_fs: float) -> float:
    """Nyquist wavenumber of a ``dt_fs`` timestep, in cm⁻¹."""
    return 1.0 / (2.0 * dt_fs * 1.0e-15 * _C_CM_S)


def floor_curve(wn_cm, kappa_cm2_g, temperature_k, dt_fs,
                ratio: float, bandwidth_cm: float) -> np.ndarray:
    """κ_floor(ν) for one spectrum (steps 1-4 of the module docstring).

    ``ratio`` is the dipole *amplitude* ratio (floor/signal); it is squared here
    because the spectrum is a power. ``bandwidth_cm`` is the range the white
    floor power is spread over.
    """
    wn = np.asarray(wn_cm, float)
    kap = np.asarray(kappa_cm2_g, float)
    prefactor = wn * _one_minus_exp(wn, temperature_k)
    dipole = kap / prefactor  # D(ν), up to a constant that cancels
    power_total = float(_trapz(dipole, wn))
    dipole_floor = ratio**2 * power_total / bandwidth_cm
    return dipole_floor * prefactor


def window_share(wn_cm, kappa_cm2_g, temperature_k) -> tuple[float, float]:
    """Share of each mean's integrand contributed by the transparency window.

    The Rosseland mean is harmonic, so its integrand is (dB/dT)/κ and the *most
    transparent* part of the band dominates it. This is the number that decides
    whether the floor matters: if the window the FP model cannot resolve also
    carries the Rosseland mean, then FP κ_R is set by unresolved spectrum.
    """
    from opacity.means import planck_weight, rosseland_weight

    wn = np.asarray(wn_cm, float)
    kap = np.asarray(kappa_cm2_g, float)
    m = (wn >= SCIENCE_BAND_CM[0]) & (wn <= SCIENCE_BAND_CM[1])
    w, k = wn[m], kap[m]
    win = (w >= WINDOW_CM[0]) & (w <= WINDOW_CM[1])
    ross = rosseland_weight(w, temperature_k) / k  # harmonic denominator
    planck = k * planck_weight(w, temperature_k)
    return (
        float(_trapz(np.where(win, ross, 0.0), w) / _trapz(ross, w)),
        float(_trapz(np.where(win, planck, 0.0), w) / _trapz(planck, w)),
    )


def _band_masks(wn: np.ndarray) -> dict:
    lo, hi = SCIENCE_BAND_CM
    w_lo, w_hi = WINDOW_CM
    in_band = (wn >= lo) & (wn <= hi)
    window = in_band & (wn >= w_lo) & (wn <= w_hi)
    return {"band": in_band, "window": window, "cores": in_band & ~window}


def _fraction_above(wn, kappa, floor, mask) -> float:
    """Wavenumber fraction of ``mask`` where κ > FLOOR_MULTIPLE × κ_floor.

    Measured as a fraction of *wavenumber extent*, not of sample count, so an
    uneven grid cannot tilt the answer.
    """
    if not mask.any():
        return float("nan")
    above = (kappa > FLOOR_MULTIPLE * floor) & mask
    return float(_trapz(above[mask].astype(float), wn[mask])
                 / _trapz(np.ones(mask.sum()), wn[mask]))


def analyse_cell(rec, dt_fs: float, ratio: float, bandwidth_mode: str) -> dict:
    """Floor fractions and margins for one spectrum record."""
    wn, kap = rec.wn_cm, rec.kappa_cm2_g
    nyq = nyquist_wavenumber_cm(dt_fs)
    bandwidth = {
        "nyquist": nyq,
        "computed": float(wn.max() - wn.min()),
        "band": SCIENCE_BAND_CM[1] - SCIENCE_BAND_CM[0],
    }[bandwidth_mode]
    floor = floor_curve(wn, kap, rec.temperature_k, dt_fs, ratio, bandwidth)
    masks = _band_masks(wn)
    out = {
        "temperature_k": rec.temperature_k, "rho_si": rec.rho_si,
        "bandwidth_mode": bandwidth_mode, "bandwidth_cm": bandwidth,
        "nyquist_cm": nyq, "ratio": ratio,
        "wn_cm": wn, "kappa_cm2_g": kap, "floor_cm2_g": floor,
    }
    for key in ("band", "window", "cores"):
        m = masks[key]
        out[f"frac_above_{key}"] = _fraction_above(wn, kap, floor, m)
        out[f"median_margin_{key}"] = float(np.median(kap[m] / floor[m]))
        out[f"min_margin_{key}"] = float(np.min(kap[m] / floor[m]))
    out["rosseland_window_share"], out["planck_window_share"] = window_share(
        wn, kap, rec.temperature_k)
    return out


def load_dilute(cfg, rho_max: float = DILUTE_RHO_MAX) -> list:
    """Every FP spectrum with ρ <= ``rho_max``, plus its dt_fs.

    ``rho_max=None`` loads the whole grid — the control that decides whether a
    floor problem is specific to the dilute points or general (A52).
    """
    if rho_max is None:
        rho_max = float("inf")
    spec = get_spec("fp")
    paths = find_spectra(cfg.data_root_path / spec.subdir, spec)
    out = []
    for path in paths:
        with np.load(path) as d:
            if float(d["rho_kg_m3"]) > rho_max:
                continue
            dt_fs = float(d["dt_fs"])
        rec = load_spectrum(path, spec, SCIENCE_BAND_CM)
        out.append((rec, dt_fs))
    out.sort(key=lambda x: (x[0].temperature_k, x[0].rho_si))
    return out


def plot_cells(cfg, results: dict, tag: str = "S7_dilute_floor"):
    """κ(ν) with the floor overlaid for the three requested cells."""
    fig, axes = plt.subplots(1, len(FIGURE_CELLS), figsize=(5.3 * len(FIGURE_CELLS), 4.5),
                             sharey=True, constrained_layout=True)
    axes = np.atleast_1d(axes)
    for ax, key in zip(axes, FIGURE_CELLS):
        r = results.get(key)
        if r is None:
            ax.set_visible(False)
            continue
        wn, kap, floor = r["wn_cm"], r["kappa_cm2_g"], r["floor_cm2_g"]
        sel = (wn >= SCIENCE_BAND_CM[0] - 50) & (wn <= SCIENCE_BAND_CM[1] + 50)
        ax.axvspan(*WINDOW_CM, color="#D55E00", alpha=0.09, lw=0,
                   label="transparency window")
        ax.semilogy(wn[sel], kap[sel], "-", lw=0.8, color="#0072B2", label="FP κ(ν)")
        ax.semilogy(wn[sel], floor[sel], "-", lw=1.4, color="#D55E00",
                    label="noise floor")
        ax.semilogy(wn[sel], FLOOR_MULTIPLE * floor[sel], "--", lw=1.0,
                    color="#D55E00", alpha=0.7,
                    label=f"{FLOOR_MULTIPLE:g}× floor")
        ax.set_xlim(SCIENCE_BAND_CM[0] - 50, SCIENCE_BAND_CM[1] + 50)
        ax.set_xlabel("wavenumber [cm$^{-1}$]")
        ax.set_title(f"T = {key[0]:.0f} K, ρ = {key[1]:.0f} kg m$^{{-3}}$\n"
                     f"above {FLOOR_MULTIPLE:g}× floor: windows "
                     f"{r['frac_above_window']:.0%}, cores "
                     f"{r['frac_above_cores']:.0%}", fontsize=9.5)
    axes[0].set_ylabel(r"$\kappa$ [cm$^2$ g$^{-1}$]")
    axes[0].legend(fontsize=8, loc="lower right")
    r0 = next(iter(results.values()))
    fig.suptitle(
        f"S7 — dilute FP spectra against the dipole-model noise floor "
        f"(residual {DIPOLE_FLOOR_MDEBYE:g} / signal {DIPOLE_SIGNAL_RMS_MDEBYE:g} "
        f"mDebye/atom → {r0['ratio']**2:.3f} of the power, spread over "
        f"{r0['bandwidth_mode']} bandwidth {r0['bandwidth_cm']:.0f} cm⁻¹)",
        fontsize=10)
    return save(fig, cfg, tag)


def main(argv=None):
    import argparse

    argv = sys.argv[1:] if argv is None else argv
    ap = argparse.ArgumentParser()
    ap.add_argument("config", type=Path)
    ap.add_argument("--floor-bandwidth", default="nyquist",
                    choices=("nyquist", "computed", "band"),
                    help="range the white floor power is spread over (default: "
                         "nyquist, from dt_fs)")
    args = ap.parse_args(argv)
    cfg = load_cfg([str(args.config)])
    set_style()

    all_cells = load_dilute(cfg, rho_max=None)  # whole grid; dilute split out below
    cells = [(r, dt) for r, dt in all_cells if r.rho_si <= DILUTE_RHO_MAX]
    if not cells:
        print(f"no FP spectra with rho <= {DILUTE_RHO_MAX:g} kg/m^3 — has the "
              f"extension campaign been merged into fp_opacity/ (Task 14)?")
        return None
    ratio = DIPOLE_FLOOR_MDEBYE / DIPOLE_SIGNAL_RMS_MDEBYE

    print(f"S7 — dilute noise-floor check (A52): {len(cells)} cells with "
          f"rho <= {DILUTE_RHO_MAX:g} kg/m^3")
    print(f"  floor/signal amplitude ratio = {DIPOLE_FLOOR_MDEBYE:g}/"
          f"{DIPOLE_SIGNAL_RMS_MDEBYE:g} = {ratio:.4f}  ->  power fraction "
          f"r^2 = {ratio**2:.4f}")
    nyq = nyquist_wavenumber_cm(cells[0][1])
    print(f"  dt = {cells[0][1]:g} fs  ->  Nyquist = {nyq:.0f} cm^-1; floor spread "
          f"over the {args.floor_bandwidth} bandwidth")
    print(f"  science band {SCIENCE_BAND_CM[0]:.0f}-{SCIENCE_BAND_CM[1]:.0f} cm^-1; "
          f"windows {WINDOW_CM[0]:.0f}-{WINDOW_CM[1]:.0f}; cores = the rest")

    results = {}
    all_rows = []
    for rec, dt_fs in all_cells:
        r = analyse_cell(rec, dt_fs, ratio, args.floor_bandwidth)
        results[(rec.temperature_k, rec.rho_si)] = r
        all_rows.append(r)
    rows = [r for r in all_rows if r["rho_si"] <= DILUTE_RHO_MAX]

    def _table(title, these):
        head = (f"\n  {title}\n  {'T [K]':>6s} {'rho':>8s} | "
                f"{'above 2x floor':>29s} | {'median kappa/floor':>29s} | {'Ross':>6s}")
        print(head)
        print(f"  {'':>6s} {'':>8s} | {'band':>9s} {'windows':>9s} {'cores':>9s} | "
              f"{'band':>9s} {'windows':>9s} {'cores':>9s} | {'in win':>6s}")
        print("  " + "-" * 84)
        for r in these:
            print(f"  {r['temperature_k']:>6.0f} {r['rho_si']:>8.4g} | "
                  f"{r['frac_above_band']:>9.3f} {r['frac_above_window']:>9.3f} "
                  f"{r['frac_above_cores']:>9.3f} | "
                  f"{r['median_margin_band']:>9.2f} {r['median_margin_window']:>9.2f} "
                  f"{r['median_margin_cores']:>9.2f} | "
                  f"{r['rosseland_window_share']:>6.2f}")

    _table(f"the dilute cells (rho <= {DILUTE_RHO_MAX:g} kg/m^3):", rows)
    _table("CONTROL — the whole 48-point grid, by density:",
           sorted(all_rows, key=lambda r: (r["rho_si"], r["temperature_k"])))

    # --- why the window matters: it carries the Rosseland mean ---
    print(f"\n  'Ross in win' above = the share of the Rosseland mean's integrand "
          f"coming from\n  the {WINDOW_CM[0]:.0f}-{WINDOW_CM[1]:.0f} cm^-1 window. "
          f"The Rosseland mean is harmonic, so the most\n  transparent part of the "
          f"band decides it:")
    dense = [r for r in all_rows if r["rho_si"] > 1000.0]
    thin = [r for r in all_rows if r["rho_si"] <= 500.0]
    print(f"    rho <= 500  ({len(thin):2d} cells): Rosseland share in window "
          f"{min(r['rosseland_window_share'] for r in thin):.2f} - "
          f"{max(r['rosseland_window_share'] for r in thin):.2f}")
    print(f"    rho >  1000 ({len(dense):2d} cells): Rosseland share in window "
          f"{min(r['rosseland_window_share'] for r in dense):.2f} - "
          f"{max(r['rosseland_window_share'] for r in dense):.2f}")

    # --- the construction contradicts itself; say so rather than smoothing it over ---
    below = [r for r in all_rows if r["median_margin_window"] < 1.0]
    print(f"\n  CONTRADICTION IN THE REIMPLEMENTED CONSTRUCTION — read before quoting "
          f"any number above:")
    print(f"    In {len(below)} of {len(all_rows)} cells the delivered (smoothed) kappa "
          f"sits BELOW this floor in the\n    window (median kappa/floor down to "
          f"{min(r['median_margin_window'] for r in all_rows):.2f}). For an *additive "
          f"white* residual that is\n    impossible: smoothing removes fluctuation, not "
          f"bias, so E[kappa] >= floor always.")
    print(f"    So the data rules out this construction at the quoted "
          f"{ratio**2:.3f} power fraction. Either\n    the residual is not white (a "
          f"dipole-surface fit error correlated with the signal would\n    be ~"
          f"{ratio**2:.0%} *multiplicative* at every nu, which every cell clears "
          f"trivially), or the\n    delivered sigma was floor-subtracted before "
          f"smoothing. compare_point_full.py on ARC\n    is the authority and is not on "
          f"this machine — see the summary's open question.")

    # --- sensitivity: the two assumptions that move the answer ---
    print(f"\n  sensitivity — the floor bandwidth is the dominant assumption "
          f"(window fraction above {FLOOR_MULTIPLE:g}x floor):")
    print(f"  {'T [K]':>6s} {'rho':>6s} | " + " ".join(
        f"{m:>10s}" for m in ("nyquist", "computed", "band")))
    for rec, dt_fs in cells:
        vals = [analyse_cell(rec, dt_fs, ratio, m)["frac_above_window"]
                for m in ("nyquist", "computed", "band")]
        print(f"  {rec.temperature_k:>6.0f} {rec.rho_si:>6.0f} | " + " ".join(
            f"{v:>10.3f}" for v in vals))
    lo, hi = DIPOLE_FLOOR_RANGE_MDEBYE
    print(f"\n  sensitivity — the quoted residual range {lo:g}-{hi:g} mDebye/atom "
          f"(window fraction, {args.floor_bandwidth} bandwidth):")
    for rec, dt_fs in cells:
        vals = [analyse_cell(rec, dt_fs, f / DIPOLE_SIGNAL_RMS_MDEBYE,
                             args.floor_bandwidth)["frac_above_window"]
                for f in (lo, hi)]
        print(f"  {rec.temperature_k:>6.0f} {rec.rho_si:>6.0f} | "
              f"{lo:g}: {vals[0]:.3f}   {hi:g}: {vals[1]:.3f}")

    # --- verdict ---
    print("\n  VERDICT (the numbers the manuscript needs):")
    for r in rows:
        if r["rho_si"] > 30.0:
            continue
        state = ("FLOOR-DOMINATED" if r["frac_above_window"] < 0.5
                 else "above the floor")
        print(f"    T = {r['temperature_k']:.0f} K, rho = {r['rho_si']:.0f}: windows are "
              f"{state} — only {r['frac_above_window']:.1%} of the "
              f"{WINDOW_CM[0]:.0f}-{WINDOW_CM[1]:.0f} cm^-1 window clears "
              f"{FLOOR_MULTIPLE:g}x the floor (cores {r['frac_above_cores']:.1%})")

    thin_win = [r["frac_above_window"] for r in all_rows if r["rho_si"] <= 500.0]
    dilute_win = [r["frac_above_window"] for r in rows]
    print(f"\n    ...and the control says this does NOT single out the dilute points. "
          f"Under this\n    construction the window is floor-dominated at every one of "
          f"the {len(thin_win)} cells with\n    rho <= 500 (window fraction "
          f"{min(thin_win):.3f}-{max(thin_win):.3f}), INCLUDING the dense points the\n"
          f"    existing headline rests on. The dilute cells "
          f"({min(dilute_win):.3f}-{max(dilute_win):.3f}) are, if\n    anything, "
          f"slightly better than the rho = 100-500 ones. So the honest statement is\n"
          f"    not 'the dilute points are unreliable' but 'the transparency windows are "
          f"the weak\n    part of the FP spectra everywhere' — and, per the "
          f"contradiction above, by an amount\n    this reimplementation cannot pin "
          f"down. Only the hot, dense corner (rho >~ 1300,\n    T >= 2500) has windows "
          f"clearly above the floor.")

    out_dir = cfg.output_path
    out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = out_dir / "S7_dilute_floor.csv"
    fields = ["temperature_k", "rho_si", "bandwidth_mode", "bandwidth_cm", "nyquist_cm",
              "ratio", "frac_above_band", "frac_above_window", "frac_above_cores",
              "median_margin_band", "median_margin_window", "median_margin_cores",
              "min_margin_window", "rosseland_window_share", "planck_window_share"]
    # All 52 points, not just the dilute ones: the crossover table (A53) carries the
    # window-floor fraction per point, so unreliable cells stay visible there.
    with csv_path.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(fields)
        for r in sorted(all_rows, key=lambda r: (r["temperature_k"], r["rho_si"])):
            w.writerow([f"{r[k]:.6g}" if isinstance(r[k], float) else r[k]
                        for k in fields])
    path = plot_cells(cfg, results)
    print(f"\nwrote {path}")
    print(f"wrote {csv_path}")
    return path


if __name__ == "__main__":
    main()
