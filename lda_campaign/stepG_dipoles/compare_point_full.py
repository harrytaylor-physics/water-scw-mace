#!/usr/bin/env python3
"""Meeting figure per state point, v2:
(a) kappa overlay cropped to the LDA-valid band (<= Nyquist), with a
    SELF-CALIBRATED LDA noise-floor curve (floor/signal ratio from the
    spectrum file applied to kappa - unit-free);
(b) reflectance 1-e on log axis (emissivity differences visible even
    when e ~ 1), with Planck-weighted e_P annotated for both.
Usage unchanged from v1.
"""
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

KB, H, C_CM = 1.380649e-23, 6.62607015e-34, 2.99792458e10
_tz = np.trapezoid if hasattr(np, "trapezoid") else np.trapz


def load_opac(path, scale=1.0):
    d = np.load(path)
    wn = np.asarray(d["wn"], float)
    kap = np.asarray(d["kappa"], float) * scale
    m = (wn > 0) & np.isfinite(kap)
    return wn[m], kap[m], d


def load_emis(path):
    d = np.load(path)
    keys = list(d.files)
    wk = next((k for k in ("wn", "wn_cm", "wavenumber", "nu") if k in keys), None)
    ek = next((k for k in ("e", "emissivity", "eps", "e_nu") if k in keys), None)
    if wk is None or ek is None:
        arrs = {k: d[k] for k in keys if getattr(d[k], "ndim", 0) == 1}
        items = sorted(arrs.items(), key=lambda kv: -len(kv[1]))
        wk = wk or items[0][0]
        ek = ek or next(k for k, v in items if k != wk and len(v) == len(arrs[wk]))
    return np.asarray(d[wk], float), np.asarray(d[ek], float)


def planck_weight(wn, T):
    nu = wn * C_CM
    x = np.clip(H * nu / (KB * T), None, 700)
    return nu**3 / np.maximum(np.exp(x) - 1.0, 1e-300)


def smooth(y, n=15):
    k = np.ones(n) / n
    return np.convolve(y, k, mode="same")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lda-opac", required=True)
    ap.add_argument("--lda-spec", required=True)
    ap.add_argument("--hyb-opac", required=True)
    ap.add_argument("--hyb-scale", type=float, default=1.0)
    ap.add_argument("--lda-emis", default=None)
    ap.add_argument("--hyb-emis", default=None)
    ap.add_argument("--wn-max", type=float, default=3300.0,
                    help="common band cap (LDA Nyquist at 5 fs ~ 3335)")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    wl, kl, dl = load_opac(args.lda_opac)
    wh, kh, _ = load_opac(args.hyb_opac, args.hyb_scale)
    T, rho = float(dl["T"]), float(dl["rho"])

    # self-calibrated floor: kappa_floor = kappa * (floor / |I|), capped at 1
    s = np.load(args.lda_spec)
    wn_s, I, fl = np.asarray(s["wn"], float), np.abs(np.asarray(s["I"])), float(s["floor"])
    band = (wl > 30) & (wl <= args.wn_max)
    I_i = np.interp(wl[band], wn_s, smooth(I, 25))
    ratio = np.clip(fl / np.maximum(I_i, 1e-300), 0, 1.0)
    k_floor = smooth(kl[band] * ratio, 15)

    two = args.lda_emis and args.hyb_emis
    fig, ax = plt.subplots(1, 2 if two else 1,
                           figsize=(9.5 if two else 6.2, 3.8))
    ax = np.atleast_1d(ax)
    mh = (wh > 30) & (wh <= args.wn_max)
    ax[0].semilogy(wh[mh], kh[mh], "k-", lw=1.1, label="revPBE0-D3")
    ax[0].semilogy(wl[band], kl[band], "k--", lw=1.1, label="LDA")
    ax[0].semilogy(wl[band], k_floor, "k:", lw=0.9, label="LDA noise floor")
    ax[0].set_xlim(0, args.wn_max)
    ax[0].set_xlabel("wavenumber (cm$^{-1}$)")
    ax[0].set_ylabel(r"$\kappa$ (cm$^2$ g$^{-1}$)")
    ax[0].set_title(f"T = {T:.0f} K, $\\rho$ = {rho:g} kg m$^{{-3}}$")
    ax[0].legend(frameon=False, fontsize=8)

    if two:
        we_l, e_l = load_emis(args.lda_emis)
        we_h, e_h = load_emis(args.hyb_emis)
        for we, e, st, lab in ((we_h, e_h, "k-", "revPBE0-D3"),
                               (we_l, e_l, "k--", "LDA")):
            m = (we > 30) & (we <= args.wn_max)
            R = np.clip(1.0 - e[m], 1e-7, 1)
            B = planck_weight(we[m], T)
            eP = float(_tz(e[m] * B, we[m]) / _tz(B, we[m]))
            ax[1].semilogy(we[m], R, st, lw=1.1,
                           label=f"{lab}  ($e_P$ = {eP:.3f})")
        ax[1].set_xlim(0, args.wn_max)
        ax[1].set_xlabel("wavenumber (cm$^{-1}$)")
        ax[1].set_ylabel("reflectance $1-e$")
        ax[1].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(f"{args.out}.png", dpi=200)
    print(f"figure -> {args.out}.png")


if __name__ == "__main__":
    main()
