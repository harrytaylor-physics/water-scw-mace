#!/usr/bin/env python3
"""Publication EOS figure: (a) representative isotherms, both functionals;
(b) full functional deviation dP/P vs density, one curve per isotherm.

  python3 paper_eos_figure.py eos_hybrid.csv eos_lda.csv \
      --isotherms 1500 2250 3000 --out eos_paper
"""
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({"font.size": 9.5, "axes.linewidth": 0.8,
                     "xtick.direction": "in", "ytick.direction": "in",
                     "xtick.top": True, "ytick.right": True})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("hyb")
    ap.add_argument("lda")
    ap.add_argument("--isotherms", nargs="+", type=float,
                    default=[1500, 3000])
    ap.add_argument("--out", default="eos_paper")
    args = ap.parse_args()
    H = np.genfromtxt(args.hyb, delimiter=",", names=True)
    L = np.genfromtxt(args.lda, delimiter=",", names=True)

    fig, ax = plt.subplots(1, 2, figsize=(7.0, 3.0))
    markers = {T: m for T, m in zip(sorted(set(H["T_K"])),
                                    ["o", "s", "^", "v", "D", "p", "x"])}

    # (a) two isotherms over matched densities; grey band = functional gap
    for T in args.isotherms:
        h = np.atleast_1d(H[H["T_K"] == T])
        l = np.atleast_1d(L[L["T_K"] == T])
        rho_m, ph, pl_, eh, el = [], [], [], [], []
        for row in np.atleast_1d(l):
            hm = np.atleast_1d(h[np.isclose(h["rho_kg_m3"], row["rho_kg_m3"])])
            if hm.size:
                rho_m.append(row["rho_kg_m3"])
                ph.append(hm[0]["P_GPa"]); eh.append(hm[0]["P_sem_GPa"])
                pl_.append(row["P_GPa"]); el.append(row["P_sem_GPa"])
        o = np.argsort(rho_m)
        rho_m = np.asarray(rho_m)[o]
        ph, pl_ = np.asarray(ph)[o], np.asarray(pl_)[o]
        eh, el = np.asarray(eh)[o], np.asarray(el)[o]
        ax[0].fill_between(rho_m, pl_, ph, color="0.88", lw=0)
        ax[0].errorbar(rho_m, ph, yerr=eh, marker="o", ms=3.5, color="k",
                       lw=1.0, ls="-", capsize=0)
        ax[0].errorbar(rho_m, pl_, yerr=el, marker="o", ms=3.8, mfc="white",
                       mec="k", color="k", lw=1.0, ls="--", capsize=0)
        mid = max(0, len(rho_m) // 2 - 1)
    ax[0].set_xscale("log")
    ax[0].set_yscale("log")
    ax[0].set_xlabel(r"$\rho$ (kg m$^{-3}$)")
    ax[0].set_ylabel(r"$P$ (GPa)")
    hyb_h = plt.Line2D([], [], color="k", ls="-", label="revPBE0-D3")
    lda_h = plt.Line2D([], [], color="k", ls="--", label="LDA")
    t_h = [plt.Line2D([], [], color="k", ls="none", marker=markers[T],
                      ms=4.0, label=f"{T:.0f} K") for T in args.isotherms]
    leg1 = ax[0].legend(handles=[hyb_h, lda_h], frameon=False, fontsize=9,
                        loc="upper left")
    ax[0].add_artist(leg1)
    ax[0].legend(handles=t_h, frameon=False, fontsize=9, loc="lower right")
    ax[0].text(0.0, 1.02, "(a)", transform=ax[0].transAxes, va="bottom",
               ha="left", fontsize=10)

    # (b) functional deviation as a grayscale contour map over (rho, T)
    import matplotlib.tri as mtri
    R, Tv, D = [], [], []
    for T in sorted(set(H["T_K"])):
        h = np.atleast_1d(H[H["T_K"] == T])
        for row in np.atleast_1d(L[L["T_K"] == T]):
            hm = np.atleast_1d(h[np.isclose(h["rho_kg_m3"], row["rho_kg_m3"])])
            if hm.size:
                R.append(row["rho_kg_m3"]); Tv.append(T)
                D.append((row["P_GPa"] / hm[0]["P_GPa"] - 1) * 100)
    R, Tv, D = map(np.asarray, (R, Tv, D))
    tri = mtri.Triangulation(np.log10(R), Tv)
    lo = 10 * np.floor(D.min() / 10)
    levels = np.arange(lo, 0.01, 10)
    cf = ax[1].tricontourf(tri, D, levels=levels, cmap="Greys_r")
    ax[1].tricontour(tri, D, levels=levels, colors="k", linewidths=0.4,
                     linestyles="solid")
    ax[1].plot(np.log10(R), Tv, "o", ms=3, mfc="white", mec="k", mew=0.7)
    xt = np.array([100, 200, 400, 800, 1600])
    xt = xt[(xt >= R.min() * 0.9) & (xt <= R.max() * 1.15)]
    ax[1].set_xticks(np.log10(xt))
    ax[1].set_xticklabels([str(v) for v in xt])
    ax[1].set_xlabel(r"$\rho$ (kg m$^{-3}$)")
    ax[1].set_ylabel(r"$T$ (K)")
    cb = fig.colorbar(cf, ax=ax[1], pad=0.02)
    cb.set_label(r"$(P_{\mathrm{LDA}}-P_{\mathrm{hyb}})/P_{\mathrm{hyb}}$ (%)")
    ax[1].text(0.0, 1.02, "(b)", transform=ax[1].transAxes, va="bottom",
               ha="left", fontsize=10)
    fig.tight_layout()
    fig.savefig(f"{args.out}.png", dpi=300)
    fig.savefig(f"{args.out}.pdf")
    print(f"figure -> {args.out}.png / .pdf")


if __name__ == "__main__":
    main()
