#!/usr/bin/env python3
"""Absolute LDA opacity per point: dipoles.npy (written by spectra.py with
the LDA dipole model) -> opacity_lda.npz via the validated
opacity_pipeline backend. Fresh files, correct units - no legacy x1000.

  python3 make_opacity_lda.py --trajdir .../production/T1000_rho129.744 \
      --dt-fs 5.0
(T and rho parsed from the directory name; 64 molecules assumed.)
"""
import argparse
import pathlib
import re
import sys
import numpy as np

sys.path.insert(0, str(pathlib.Path.home() / "scw"))
from opacity_pipeline import spectrum_from_dipoles  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--trajdir", required=True)
    ap.add_argument("--dt-fs", type=float, default=5.0)
    ap.add_argument("--n-molecules", type=int, default=64)
    args = ap.parse_args()
    d = pathlib.Path(args.trajdir)
    m = re.match(r"T(\d+)_rho([\d.]+)", d.name)
    T, rho = float(m.group(1)), float(m.group(2))
    M = np.load(d / "dipoles.npy")
    res = spectrum_from_dipoles(M, T_K=T, dt_fs=args.dt_fs,
                                n_molecules=args.n_molecules,
                                density_g_cm3=rho / 1000.0,
                                units="eAng")
    np.savez(d / "opacity_lda.npz",
             wn=res["wn_cm-1"], kappa=res["kappa_cm2_g"],
             alpha=res["alpha_cm-1"],
             planck=res["planck_mean_cm2_g"],
             rosseland=res["rosseland_mean_cm2_g"],
             T=T, rho=rho, n_frames=res["n_frames"],
             source="LDA v1 MD + water_dipole_lda + opacity_pipeline")
    print(f"{d.name}: kappa_P = {res['planck_mean_cm2_g']:.4g} cm2/g, "
          f"kappa_R = {res['rosseland_mean_cm2_g']:.4g} cm2/g "
          f"({res['n_frames']} frames) -> opacity_lda.npz")


if __name__ == "__main__":
    main()
