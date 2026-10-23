#!/usr/bin/env python
"""
pyscf_monomer.py — revPBE0 dipoles of an isolated water molecule at the 13 monomer_test geometries,
dipole derivatives by central differences, and comparison with the water_dipole_v3 model.

Functional: revPBE0 = 0.25 HF + 0.75 revPBE exchange (libxc GGA_X_PBE_R) + PBE correlation, as in the
training labels. D3 is an energy-only correction and does not change the density or the dipole.
For an isolated neutral molecule the Wannier-centre dipole equals the dipole of the electron density,
which PySCF computes directly. Basis: all-electron Gaussian sets (the training labels used
DZVP-MOLOPT-SR-GTH with GTH pseudopotentials); three bases show the basis sensitivity.

Model dipoles are copied verbatim from the eval_model.py output (water_dipole_v3, same geometries, Debye).

usage: python pyscf_monomer.py [--dir .] [--bases cc-pvdz aug-cc-pvdz aug-cc-pvtz]
"""
import argparse, glob, json, os, numpy as np
from pyscf import gto, dft
from pyscf.dft import libxc

ap = argparse.ArgumentParser()
ap.add_argument("--dir", default=".")
ap.add_argument("--bases", nargs="+", default=["cc-pvdz", "aug-cc-pvdz", "aug-cc-pvtz"])
a = ap.parse_args()
XC = "0.25*HF + 0.75*GGA_X_PBE_R, GGA_C_PBE"
MU_EXP = 1.855  # D, gas-phase water

MODEL = {"g00": (0, 0, 1.82217), "g01": (0, 0, 1.83603), "g02": (0, 0, 1.80874),
         "g03": (0.01903, 0, 1.82248), "g04": (-0.01903, 0, 1.82248), "g05": (0, 0, 1.80527),
         "g06": (0, 0, 1.83875), "g07": (0, 0, 1.90323), "g08": (0, 0, 1.86991),
         "g09": (0.02516, 0, 1.88670), "g10": (-0.02516, 0, 1.88670), "g11": (0, 0, 1.86836),
         "g12": (0, 0, 1.90347)}

def read_xyz(p):
    lines = open(p).read().split("\n"); n = int(lines[0])
    return [(l.split()[0], [float(x) for x in l.split()[1:4]]) for l in lines[2:2 + n]]

def internals(atoms):
    O = np.array([c for s, c in atoms if s == "O"][0]); H = [np.array(c) for s, c in atoms if s == "H"]
    v1, v2 = H[0] - O, H[1] - O; r1, r2 = np.linalg.norm(v1), np.linalg.norm(v2)
    return r1, r2, np.arccos(np.dot(v1, v2) / (r1 * r2))

geo = {os.path.basename(p)[:-4]: read_xyz(p) for p in sorted(glob.glob(os.path.join(a.dir, "g??.xyz")))}
assert set(geo) == set(MODEL), f"expected g00..g12, found {sorted(geo)}"
print(f"13 geometries read; functional '{XC}', HF fraction {libxc.hybrid_coeff(XC):.2f}")

dft_mu = {}
for b in a.bases:
    dft_mu[b] = {}
    for k, atoms in geo.items():
        mol = gto.M(atom=[(s, c) for s, c in atoms], basis=b, unit="Angstrom", verbose=0)
        mf = dft.RKS(mol); mf.xc = XC; mf.grids.level = 5; mf.conv_tol = 1e-11
        mf.kernel(); assert mf.converged, f"SCF not converged: {b} {k}"
        dft_mu[b][k] = [float(x) for x in mf.dip_moment(unit="Debye", verbose=0)]
    print(f"  {b:12s} |mu(eq)| = {np.linalg.norm(dft_mu[b]['g00']):.4f} D")

pairs = [("g01", "g02"), ("g03", "g04"), ("g05", "g06"), ("g07", "g08"), ("g09", "g10"), ("g11", "g12")]
def deriv(mu, p, q):
    (r1a, r2a, ta), (r1b, r2b, tb) = internals(geo[p]), internals(geo[q])
    cands = [("sym", (r1a + r2a - r1b - r2b) / 2, "D/A"), ("asym", (r1a - r2a - (r1b - r2b)) / 2, "D/A"),
             ("bend", ta - tb, "D/rad")]
    mode, dq, unit = max(cands, key=lambda x: abs(x[1]))
    return mode, (r1a + r2a + r1b + r2b) / 4, np.linalg.norm(np.array(mu[p]) - np.array(mu[q])) / abs(dq), unit

big = a.bases[-1]
print(f"\nmu(eq): model {np.linalg.norm(MODEL['g00']):.4f} D | revPBE0/{big} {np.linalg.norm(dft_mu[big]['g00']):.4f} D | experiment {MU_EXP} D")
hdr = "".join(f"{b:>14s}" for b in a.bases)
print(f"\n{'mode':5s} {'r_ref':>6s}  {'model':>8s}{hdr}  model/DFT({big})")
out = {"xc": XC, "model_mu": MODEL, "dft_mu": dft_mu, "derivatives": []}
proxy = {}
for p, q in pairs:
    mode, ref, dm, unit = deriv(MODEL, p, q)
    dd = {b: deriv(dft_mu[b], p, q)[2] for b in a.bases}
    print(f"{mode:5s} {ref:6.3f}  {dm:8.3f}" + "".join(f"{dd[b]:14.3f}" for b in a.bases) + f"  {dm/dd[big]:8.2f}   ({unit})")
    out["derivatives"].append(dict(mode=mode, r_ref=ref, unit=unit, model=dm, dft=dd))
    proxy.setdefault(round(ref, 2), {})[mode] = (dm, dd[big])
print("\nstretch/bend intensity proxy (sym^2 + asym^2)/bend^2, model relative to DFT:")
for ref, d in sorted(proxy.items()):
    m = (d["sym"][0]**2 + d["asym"][0]**2) / d["bend"][0]**2; f = (d["sym"][1]**2 + d["asym"][1]**2) / d["bend"][1]**2
    print(f"  r_ref {ref:.2f} A: model/DFT = {m/f:.2f}   (dilute-limit spectra show 3.3-3.6)")
    out[f"proxy_ratio_{ref:.2f}"] = m / f
json.dump(out, open(os.path.join(a.dir, "pyscf_monomer.json"), "w"), indent=1)
print("\nReading: model/DFT near 1 for the stretch derivatives -> the dipole model is not the cause of the stretch excess;"
      "\n         1.5 or more -> it is. Wrote pyscf_monomer.json")
