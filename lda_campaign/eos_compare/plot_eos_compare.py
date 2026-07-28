#!/usr/bin/env python3
"""Functional sensitivity of the EOS: isotherm overlay + dP/P vs (T, rho).
  python3 plot_eos_compare.py eos_hybrid.csv eos_lda.csv --out eos_functional
"""
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("hyb")
    ap.add_argument("lda")
    ap.add_argument("--out", default="eos_functional")
    args = ap.parse_args()
    H = np.genfromtxt(args.hyb, delimiter=",", names=True)
    L = np.genfromtxt(args.lda, delimiter=",", names=True)
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    Ts = sorted(set(H["T_K"]))
    markers = ["o", "s", "^", "v", "d", "*", "x"]
    lss = ["-", "--", ":", "-.", "-", "--", ":"]
    diffs = []
    for T, mk, ls in zip(Ts, markers, lss):
        h = np.atleast_1d(H[H["T_K"] == T])
        o = np.argsort(h["rho_kg_m3"])
        ax[0].plot(h["rho_kg_m3"][o], h["P_GPa"][o], marker=mk, ls=ls,
                   color="k", ms=4, lw=0.9, label=f"{T:.0f} K")
        l = np.atleast_1d(L[L["T_K"] == T])
        if l.size:
            o2 = np.argsort(l["rho_kg_m3"])
            ax[0].plot(l["rho_kg_m3"][o2], l["P_GPa"][o2], marker=mk, ls=ls,
                       color="0.55", ms=4, lw=0.9)
        for row in l:
            hm = np.atleast_1d(h[np.isclose(h["rho_kg_m3"],
                                            row["rho_kg_m3"])])
            if hm.size:
                dp = (row["P_GPa"] - hm[0]["P_GPa"]) / hm[0]["P_GPa"] * 100
                diffs.append((row["rho_kg_m3"], T, dp))
                ax[1].scatter(row["rho_kg_m3"], T,
                              s=float(np.clip(abs(dp) * 12, 20, 400)),
                              facecolors="none", edgecolors="k")
                ax[1].annotate(f"{dp:+.0f}%", (row["rho_kg_m3"], T),
                               textcoords="offset points", xytext=(6, 4),
                               fontsize=7)
    ax[0].set_xscale("log")
    ax[0].set_yscale("log")
    ax[0].set_xlabel("rho (kg/m3)")
    ax[0].set_ylabel("P (GPa)")
    ax[0].legend(frameon=False, fontsize=8, title="black: hybrid | grey: LDA")
    ax[1].set_xscale("log")
    ax[1].set_xlabel("rho (kg/m3)")
    ax[1].set_ylabel("T (K)")
    ax[1].set_title(r"$(P_{\rm LDA}-P_{\rm hyb})/P_{\rm hyb}$")
    fig.tight_layout()
    fig.savefig(f"{args.out}.png", dpi=170)
    print(f"figure -> {args.out}.png")
    if diffs:
        d = np.array([x[2] for x in diffs])
        print(f"matched points: {len(d)} | dP/P range {d.min():+.1f}% to "
              f"{d.max():+.1f}% | mean {d.mean():+.1f}%")


if __name__ == "__main__":
    main()
