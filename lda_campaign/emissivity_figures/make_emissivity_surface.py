#!/usr/bin/env python3
"""Planck-weighted emissivity surface over (T, rho) from per-point
emissivity npz files, in the style of the hybrid opacity surface.
Plots log10(1 - e_P) (Planck-weighted reflectance) for topography.

  python3 make_emissivity_surface.py \
      --emis-root /data/.../scw_prod/emissivity \
      --lda-root  /data/.../lda_pipeline/production/emissivity_lda \
      --out emissivity_surface
--style 3d (default, trisurf like the opacity surface) or contour.
"""
import argparse
import glob
import pathlib
import re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

KB, H, C_CM = 1.380649e-23, 6.62607015e-34, 2.99792458e10
_tz = np.trapezoid if hasattr(np, "trapezoid") else np.trapz


def e_planck(npz_path, T):
    d = np.load(npz_path)
    wn = np.asarray(d["wn"], float)
    e = np.asarray(d["e"], float)
    m = (wn > 30) & np.isfinite(e)
    wn, e = wn[m], np.clip(e[m], 0, 1)
    nu = wn * C_CM
    x = np.clip(H * nu / (KB * T), None, 700)
    B = nu**3 / np.maximum(np.exp(x) - 1.0, 1e-300)
    return float(_tz(e * B, wn) / _tz(B, wn))


def harvest(root):
    pts = []
    for d in sorted(glob.glob(str(pathlib.Path(root) / "T*_rho*"))):
        name = pathlib.Path(d).name
        if name.endswith(("_5fs", "_gf2fail", "_gf2fail2")):
            continue
        m = re.match(r"T(\d+)_rho([\d.]+)", name)
        f = pathlib.Path(d) / "emissivity_v2.npz"
        if m and f.exists():
            T, rho = float(m.group(1)), float(m.group(2))
            pts.append((T, rho, e_planck(f, T)))
    return pts


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--emis-root", required=True)
    ap.add_argument("--lda-root", default=None)
    ap.add_argument("--style", choices=["3d", "contour"], default="3d")
    ap.add_argument("--out", default="emissivity_surface")
    args = ap.parse_args()
    hyb = harvest(args.emis_root)
    if len(hyb) < 4:
        raise SystemExit(f"only {len(hyb)} hybrid points found under "
                         f"{args.emis_root} - check the path")
    T, R, E = map(np.asarray, zip(*hyb))
    Z = np.log10(np.clip(1.0 - E, 1e-6, 1))
    print(f"hybrid: {len(hyb)} points, e_P {E.min():.3f}-{E.max():.3f}")
    lda = harvest(args.lda_root) if args.lda_root else []
    for t_, r_, e_ in lda:
        print(f"  LDA T{t_:.0f} rho{r_:g}: e_P = {e_:.3f}")

    if args.style == "3d":
        import matplotlib.tri as mtri
        fig = plt.figure(figsize=(8.5, 6))
        ax = fig.add_subplot(111, projection="3d")
        tri0 = mtri.Triangulation(np.log10(R), T)
        refiner = mtri.UniformTriRefiner(tri0)
        tri_r, Z_r = refiner.refine_field(Z, subdiv=3)
        from matplotlib.colors import Normalize, LightSource
        import matplotlib.cm as cm
        s = ax.plot_trisurf(tri_r, Z_r, cmap="magma", linewidth=0.0,
                            edgecolor="none", antialiased=True)
        # Lambert relief shading: per-face normals x light direction
        xs, ys = tri_r.x, tri_r.y
        fx = (xs - xs.min()) / np.ptp(xs)
        fy = (ys - ys.min()) / np.ptp(ys)
        fz = (Z_r - Z_r.min()) / np.ptp(Z_r)
        trisF = tri_r.triangles
        v0 = np.c_[fx[trisF[:, 0]], fy[trisF[:, 0]], fz[trisF[:, 0]]]
        v1 = np.c_[fx[trisF[:, 1]], fy[trisF[:, 1]], fz[trisF[:, 1]]]
        v2 = np.c_[fx[trisF[:, 2]], fy[trisF[:, 2]], fz[trisF[:, 2]]]
        nrm = np.cross(v1 - v0, v2 - v0)
        nrm /= np.linalg.norm(nrm, axis=1, keepdims=True) + 1e-12
        nrm[nrm[:, 2] < 0] *= -1
        light = np.array([-0.5, 0.35, 0.79])
        light /= np.linalg.norm(light)
        inten = np.clip(nrm @ light, 0, 1)
        zface = Z_r[trisF].mean(axis=1)
        base = cm.magma(Normalize(Z.min(), Z.max())(zface))
        base[:, :3] *= (0.45 + 0.55 * inten)[:, None]
        s.set_fc(base)
        # draped iso-lines at their true heights
        for lev in np.linspace(Z.min(), Z.max(), 12)[1:-1]:
            ax.tricontour(tri_r, Z_r, levels=[lev], colors="k",
                          linewidths=0.35, alpha=0.55, zdir="z",
                          offset=lev + 0.008)
        zmin = Z.min() - 0.35 * (Z.max() - Z.min())
        ax.tricontourf(tri_r, Z_r, levels=14, cmap="magma", alpha=0.35,
                       zdir="z", offset=zmin)
        ax.tricontour(tri_r, Z_r, levels=14, colors="k", linewidths=0.4,
                      alpha=0.5, zdir="z", offset=zmin)
        ax.set_zlim(zmin, Z.max() + 0.05)
        if lda:
            Tl, Rl, El = map(np.asarray, zip(*lda))
            ax.scatter(np.log10(Rl), Tl,
                       np.log10(np.clip(1 - El, 1e-6, 1)) + 0.02,
                       c="cyan", s=45, ec="k", lw=0.8, depthshade=False,
                       label="LDA (this work)")
            ax.legend(loc="upper left", frameon=False, fontsize=9)
        xt = np.array([100, 200, 400, 1000, 2000])
        ax.set_xticks(np.log10(xt))
        ax.set_xticklabels([str(v) for v in xt])
        ax.set_xlabel(r"density / kg m$^{-3}$", labelpad=8)
        ax.set_ylabel("temperature / K", labelpad=8)
        ax.set_zlabel(r"log$_{10}(1-e_P)$", labelpad=6)
        cb = fig.colorbar(s, shrink=0.55, pad=0.08)
        cb.set_label(r"log$_{10}(1-e_P)$")
        ax.view_init(elev=22, azim=-120)
    else:
        import matplotlib.tri as mtri
        fig, ax = plt.subplots(figsize=(6.4, 4.6))
        tri = mtri.Triangulation(np.log10(R), T)
        cf = ax.tricontourf(tri, Z, levels=12, cmap="magma")
        ax.tricontour(tri, Z, levels=12, colors="k", linewidths=0.3)
        ax.plot(np.log10(R), T, "o", ms=3, mfc="white", mec="k", mew=0.6)
        if lda:
            Tl, Rl, _ = map(np.asarray, zip(*lda))
            ax.plot(np.log10(Rl), Tl, "s", ms=7, mfc="none", mec="cyan",
                    mew=1.6, label="LDA (this work)")
            ax.legend(frameon=False, fontsize=9)
        xt = np.array([100, 200, 400, 1000, 2000])
        ax.set_xticks(np.log10(xt))
        ax.set_xticklabels([str(v) for v in xt])
        ax.set_xlabel(r"$\rho$ (kg m$^{-3}$)")
        ax.set_ylabel(r"$T$ (K)")
        cb = fig.colorbar(cf, pad=0.02)
        cb.set_label(r"log$_{10}(1-e_P)$")
    fig.tight_layout()
    fig.savefig(f"{args.out}.png", dpi=220)
    print(f"figure -> {args.out}.png")


if __name__ == "__main__":
    main()
