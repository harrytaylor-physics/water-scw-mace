#!/usr/bin/env python3
"""Waterfall of spectral emissivity: log10(1-e(nu)) ridgelines stacked
along density for one isotherm.

  python3 make_emissivity_waterfall.py --emis-root .../scw_prod/emissivity \
      --T 3000 --out emissivity_waterfall
"""
import argparse
import glob
import pathlib
import re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import matplotlib.cm as cm
from matplotlib.colors import Normalize


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--emis-root", required=True)
    ap.add_argument("--T", type=float, default=3000)
    ap.add_argument("--wn-max", type=float, default=4500)
    ap.add_argument("--out", default="emissivity_waterfall")
    ap.add_argument("--spacing", choices=["index", "log"], default="index",
                    help="index: even ridge spacing, true rho as labels")
    args = ap.parse_args()
    rows = []
    for d in sorted(glob.glob(str(pathlib.Path(args.emis_root) / "T*_rho*"))):
        name = pathlib.Path(d).name
        if name.endswith(("_5fs", "_gf2fail", "_gf2fail2")):
            continue
        m = re.match(r"T(\d+)_rho([\d.]+)", name)
        f = pathlib.Path(d) / "emissivity_v2.npz"
        if not (m and f.exists()) or float(m.group(1)) != args.T:
            continue
        dd = np.load(f)
        wn = np.asarray(dd["wn"], float)
        e = np.clip(np.asarray(dd["e"], float), 0, 1 - 1e-7)
        sel = (wn > 30) & (wn <= args.wn_max)
        rows.append((float(m.group(2)), wn[sel], np.log10(1 - e[sel])))
    if len(rows) < 3:
        raise SystemExit(f"only {len(rows)} points at T={args.T:g} - "
                         f"try another isotherm")
    rows.sort()
    print(f"T = {args.T:g} K: {len(rows)} densities "
          f"({rows[0][0]:g} - {rows[-1][0]:g} kg/m3)")

    zfloor = min(r[2].min() for r in rows) - 0.4
    fig = plt.figure(figsize=(8.5, 6))
    ax = fig.add_subplot(111, projection="3d")
    norm = Normalize(np.log10(rows[0][0]), np.log10(rows[-1][0]))
    verts, colors, ys = [], [], []
    for i, (rho, wn, z) in enumerate(rows):
        y = float(i) if args.spacing == "index" else np.log10(rho)
        ys.append(y)
        v = [(wn[0], zfloor)] + list(zip(wn, z)) + [(wn[-1], zfloor)]
        verts.append(v)
        colors.append(cm.magma(norm(y)))
    poly = PolyCollection(verts, facecolors=[(*c[:3], 0.85) for c in colors],
                          edgecolors="k", linewidths=0.7)
    ax.add_collection3d(poly, zs=ys, zdir="y")
    ax.set_xlim(0, args.wn_max)
    ax.set_ylim(min(ys) - 0.05, max(ys) + 0.05)
    ax.set_zlim(zfloor, max(r[2].max() for r in rows) + 0.2)
    if args.spacing == "index":
        ax.set_yticks(ys)
        ax.set_yticklabels([f"{r[0]:.0f}" for r in rows], fontsize=8)
    else:
        yt = np.array([100, 200, 400, 800, 1600])
        yt = yt[(yt >= rows[0][0] * 0.9) & (yt <= rows[-1][0] * 1.2)]
        ax.set_yticks(np.log10(yt))
        ax.set_yticklabels([str(v) for v in yt])
    ax.set_xlabel(r"wavenumber / cm$^{-1}$", labelpad=8)
    ax.set_ylabel(r"density / kg m$^{-3}$", labelpad=10)
    ax.set_zlabel(r"log$_{10}(1-e)$", labelpad=6)
    ax.set_title(f"spectral emissivity, T = {args.T:.0f} K", pad=2)
    ax.view_init(elev=25, azim=-58)
    fig.tight_layout()
    fig.savefig(f"{args.out}.png", dpi=220)
    print(f"figure -> {args.out}.png")


if __name__ == "__main__":
    main()
