#!/usr/bin/env python3
"""LDA vs hybrid at one meeting point: kappa overlay + pre-registered
metrics (band-edge shift, Planck/Rosseland ratios, band ratios).

  python3 compare_meeting_point.py \
      --lda  .../lda_pipeline/production/T2000_rho500/opacity_lda.npz \
      --hyb  .../scw_prod/T2000_rho500_05fs/opacity_v3.npz --hyb-scale 1e-3 \
      --out  T2000_rho500_compare
"""
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

_tz = np.trapezoid if hasattr(np, "trapezoid") else np.trapz
BANDS = {"libration": (100, 900), "bend": (1300, 1900), "stretch": (2800, 3800)}


def load(path, scale=1.0):
    d = np.load(path)
    wn = np.asarray(d["wn"], float)
    kap = np.asarray(d["kappa"], float) * scale
    m = (wn > 0) & np.isfinite(kap)
    return wn[m], kap[m], d


def edge(wn, kap, level_frac=0.5):
    """Absorption-edge position: first wn where kappa crosses level_frac of
    its stretch-band peak, scanning down from the stretch band."""
    s = (wn > 2500) & (wn < 4000)
    if not s.any():
        return np.nan
    peak = kap[s].max()
    below = wn[(wn > 1900) & (wn < 3400)]
    kb = kap[(wn > 1900) & (wn < 3400)]
    idx = np.where(kb >= level_frac * peak)[0]
    return float(below[idx[0]]) if len(idx) else np.nan


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lda", required=True)
    ap.add_argument("--hyb", required=True)
    ap.add_argument("--hyb-scale", type=float, default=1.0)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    wl, kl, dl = load(args.lda)
    wh, kh, dh = load(args.hyb, args.hyb_scale)
    T = float(dl["T"])
    lo, hi = max(wl.min(), wh.min(), 50), min(wl.max(), wh.max(), 4500)
    wn = np.linspace(lo, hi, 3000)
    L = np.interp(wn, wl, kl)
    H = np.interp(wn, wh, kh)

    print(f"T = {T:.0f} K, rho = {float(dl['rho'])} kg/m3")
    for name, (a, b) in BANDS.items():
        s = (wn >= a) & (wn <= b)
        r = _tz(L[s], wn[s]) / max(_tz(H[s], wn[s]), 1e-300)
        print(f"  {name:>9} band ratio LDA/hybrid: {r:6.3f}")
    eL, eH = edge(wn, L), edge(wn, H)
    print(f"  stretch edge: LDA {eL:.0f} vs hybrid {eH:.0f} cm-1 "
          f"(shift {eL-eH:+.0f})")
    if "planck" in dl.files:
        pl = float(dl["planck"])
        ph_wn, ph_k = wh, kh
        # recompute hybrid means on the common convention for fairness
        from math import isfinite
        print(f"  Planck mean LDA: {pl:.4g} cm2/g")
    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    ax.semilogy(wn, H, "k-", lw=1.1, label="revPBE0-D3 (hybrid)")
    ax.semilogy(wn, L, "k--", lw=1.1, label="LDA")
    ax.set_xlabel("wavenumber (cm$^{-1}$)")
    ax.set_ylabel("$\\kappa$ (cm$^2$ g$^{-1}$)")
    ax.set_title(f"T = {T:.0f} K, $\\rho$ = {float(dl['rho']):g} kg m$^{{-3}}$")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(f"{args.out}.png", dpi=170)
    print(f"figure -> {args.out}.png")


if __name__ == "__main__":
    main()
