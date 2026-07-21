#!/usr/bin/env python3
"""Planck/Rosseland mean stability vs Rosseland window and floor.
Input: npz with wn (cm^-1), kappa (cm^2/g), T (K)."""
import argparse
import numpy as np

KB, H, C_CM, C = 1.380649e-23, 6.62607015e-34, 2.99792458e10, 2.99792458e8
_trapz = getattr(np, "trapezoid", None) or np.trapz

def means(wn, kap, T, lo, hi, floor):
    nu = wn * C_CM
    B = (2*H*nu**3/C**2) / np.maximum(np.exp(np.clip(H*nu/(KB*T), None, 700)) - 1, 1e-300)
    x = H*nu/(KB*T)
    ex = np.exp(np.clip(x, None, 700))
    dB = (2*H*nu**3/C**2) * ex/np.maximum((ex-1)**2, 1e-300) * (x/T)
    kp = _trapz(kap*B, nu) / _trapz(B, nu)
    m = (wn >= lo) & (wn <= hi)
    ks = np.maximum(kap[m], floor*np.nanmax(kap[m]))
    kr = _trapz(dB[m], nu[m]) / _trapz(dB[m]/ks, nu[m])
    return kp, kr

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("npz")
    ap.add_argument("--windows", nargs="+", default=["100:5000", "50:5500", "200:4500"])
    ap.add_argument("--floors", nargs="+", type=float, default=[1e-6, 1e-5, 1e-7])
    args = ap.parse_args()
    d = np.load(args.npz)
    wn, kap, T = np.asarray(d["wn"], float), np.asarray(d["kappa"], float), float(d["T"])
    m = wn > 0; wn, kap = wn[m], kap[m]
    lo0, hi0 = map(float, args.windows[0].split(":"))
    kp0, kr0 = means(wn, kap, T, lo0, hi0, args.floors[0])
    print(f"reference (window {args.windows[0]}, floor {args.floors[0]:g}): "
          f"kappa_P {kp0:.4g}  kappa_R {kr0:.4g}")
    print("\nwindow variation (floor fixed):")
    for w in args.windows[1:]:
        lo, hi = map(float, w.split(":"))
        kp, kr = means(wn, kap, T, lo, hi, args.floors[0])
        f1, f2 = 100*(kp/kp0-1), 100*(kr/kr0-1)
        print(f"  {w:>11}: dP {f1:+5.2f}%  dR {f2:+5.2f}%"
              f"{'' if max(abs(f1),abs(f2))<2 else '   <-- exceeds 2%'}")
    print("floor variation (window fixed):")
    for fl in args.floors[1:]:
        kp, kr = means(wn, kap, T, lo0, hi0, fl)
        f2 = 100*(kr/kr0-1)
        print(f"  {fl:>11g}: dR {f2:+5.2f}%"
              f"{'' if abs(f2)<2 else '   <-- exceeds 2% (opacity minima matter)'}")

if __name__ == "__main__":
    main()
