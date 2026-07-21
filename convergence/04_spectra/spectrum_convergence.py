#!/usr/bin/env python3
"""Spectrum-convergence scans from a saved dipole time series (dipoles.npy,
shape (n_frames, 3), same conventions as the production spectra.py)."""
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HBAR, KB, C_CM = 1.054571817e-34, 1.380649e-23, 2.99792458e10
BANDS = [(50, 900, "libration"), (1200, 2000, "bend"), (2800, 3800, "stretch")]

def ac(x):
    n = len(x)
    f = np.fft.rfft(x, n=2*n)
    a = np.fft.irfft(f*np.conj(f), n=2*n)[:n].real
    return a / (n - np.arange(n))

def spectrum(M, dt_s, T, maxlag, nfft=2**17):
    M = M - M.mean(0)
    C = sum(ac(M[:, k]) for k in range(3))
    m = min(maxlag, len(C)//2)
    Cw = C[:m] * np.hanning(2*m)[m:]
    I = np.fft.rfft(Cw, n=nfft).real
    freq = np.fft.rfftfreq(nfft, dt_s)
    wn = freq / C_CM
    w = 2*np.pi*freq
    alpha = w * (1 - np.exp(-HBAR*w/(KB*T))) * I
    return wn, np.abs(I), alpha

def band_metrics(wn, I):
    out = {}
    for lo, hi, name in BANDS:
        m = (wn > lo) & (wn < hi)
        if m.any():
            out[name] = (float(wn[m][np.argmax(I[m])]), float(np.trapezoid(I[m], wn[m])
                         if hasattr(np, "trapezoid") else np.trapz(I[m], wn[m])))
    return out

def overlay(fig_name, curves, labels, xmax):
    fig, ax = plt.subplots(figsize=(7, 3.5))
    styles = ["k-", "k--", "k:", "k-."]
    for (wn, I), lab, st in zip(curves, labels, styles):
        m = (wn > 5) & (wn < xmax)
        ax.semilogy(wn[m], I[m], st, lw=0.9, label=lab)
    ax.set_xlabel("wavenumber (cm$^{-1}$)"); ax.set_ylabel("|I| (arb)")
    ax.legend(frameon=False)
    fig.tight_layout(); fig.savefig(fig_name, dpi=170); plt.close(fig)

def report(tag, ref, test):
    for name in ref:
        if name in test:
            dpk = test[name][0] - ref[name][0]
            dint = 100*(test[name][1]/ref[name][1] - 1) if ref[name][1] else 0
            flag = "" if abs(dpk) < 10 and abs(dint) < 5 else "  <-- exceeds criterion"
            print(f"  {tag} {name:>9}: peak {dpk:+6.1f} cm-1, intensity {dint:+5.1f}%{flag}")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dipoles"); ap.add_argument("--dt-fs", type=float, required=True)
    ap.add_argument("--temp", type=float, required=True)
    ap.add_argument("--maxlag", type=int, default=8000)
    args = ap.parse_args()
    M = np.load(args.dipoles)
    dt = args.dt_fs * 1e-15
    wn0, I0, _ = spectrum(M, dt, args.temp, args.maxlag)
    ref = band_metrics(wn0, I0)
    nyq = wn0.max(); xmax = min(4500, 0.9*nyq)
    print(f"reference: {len(M)} frames, maxlag {args.maxlag}")

    # maxlag scan
    curves, labels = [(wn0, I0)], [f"maxlag {args.maxlag}"]
    print("\nACF window (maxlag):")
    for ml in (args.maxlag//2, args.maxlag*2):
        wn, I, _ = spectrum(M, dt, args.temp, ml)
        report(f"maxlag {ml:5d}", ref, band_metrics(wn, I))
        curves.append((wn, I)); labels.append(f"maxlag {ml}")
    overlay("conv_maxlag.png", curves, labels, xmax)

    # trajectory-length scan
    print("\ntrajectory length:")
    curves, labels = [], []
    for f in (0.25, 0.5, 0.75, 1.0):
        n = int(len(M)*f)
        wn, I, _ = spectrum(M[:n], dt, args.temp, min(args.maxlag, n//2))
        report(f"{int(f*100):3d}% traj ", ref, band_metrics(wn, I))
        curves.append((wn, I)); labels.append(f"{int(f*100)}%")
    overlay("conv_trajlen.png", curves, labels, xmax)
    print("\nfigures -> conv_maxlag.png, conv_trajlen.png")

if __name__ == "__main__":
    main()
