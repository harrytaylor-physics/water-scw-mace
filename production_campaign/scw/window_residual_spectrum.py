#!/usr/bin/env python
"""
window_residual_spectrum.py — does the dipole model reproduce the DFT absorption, band by band?

Labels: DFT Wannier dipoles in <dir>/fNNNN (frames --dt fs apart), extracted with dip() taken verbatim from
build_dipole_xyz_v3.py: its source is exec'd only up to the line 'out=[]', so the extraction is identical to
the training labels and the file's harvest loop (which rewrites dipole_train_snip.xyz) never runs.
Predictions: the dipole model on the same coordinates. All frames are outside the training set.

Reports per band: A_pred/A_label and A_resid/A_label, where A = integral of nu(1-exp(-hc nu/kT)) |M(nu)|^2,
plus the same ratios on each half of the trajectory. Handles a partly finished label array: analyses the
longest leading block with >= 98% of frames labelled (gaps linearly interpolated, labels only).

usage: python window_residual_spectrum.py [--dir diplabels_window/T1000_rho130] [--dt 2.0] [--T 1000]
       [--model water_dipole_v3.model] [--smooth 60] [--out window_residual]
"""
import argparse, glob, os, re, sys, numpy as np
ap = argparse.ArgumentParser()
ap.add_argument("--dir", default="diplabels_window/T1000_rho130")
ap.add_argument("--dt", type=float, default=2.0); ap.add_argument("--T", type=float, default=1000.0)
ap.add_argument("--model", default="water_dipole_v3.model")
ap.add_argument("--smooth", type=float, default=60.0, help="Gaussian FWHM, cm^-1")
ap.add_argument("--out", default="window_residual"); ap.add_argument("--device", default="auto")
a = ap.parse_args()

# ---------- dipole extraction identical to the training labels ----------
src = open("build_dipole_xyz_v3.py").read()
assert "\nout=[]" in src, "harvest-loop marker not found in build_dipole_xyz_v3.py"
ns = {}; exec(src[:src.index("\nout=[]")], ns); dip = ns["dip"]

# ---------- labels ----------
dirs = sorted(glob.glob(os.path.join(a.dir, "f[0-9][0-9][0-9][0-9]")))
if not dirs: sys.exit(f"no frame directories under {a.dir}")
lab, bad = {}, dict(no_output=0, not_ended=0, assignment_failed=0)
for d in dirs:
    i = int(os.path.basename(d)[1:]); fo = f"{d}/dip.out"
    if not os.path.exists(fo): bad["no_output"] += 1; continue
    if "PROGRAM ENDED" not in open(fo).read(): bad["not_ended"] += 1; continue
    r = dip(d)
    if r is None: bad["assignment_failed"] += 1; continue
    lab[i] = r[0]
print(f"{len(dirs)} frames generated; {len(lab)} labelled; not labelled: {bad}")
avail = np.zeros(len(dirs), bool); avail[list(lab)] = True
N = max((k + 1 for k in range(len(avail)) if avail[:k + 1].mean() >= 0.98), default=0)
if N < 100: sys.exit(f"only {N} usable leading frames ({N*a.dt:.0f} fs); need >= 100 — wait for more labels")
print(f"analysing frames 0..{N-1} = {N*a.dt/1000:.3f} ps, resolution {1/(N*a.dt*2.99792458e-5):.0f} cm^-1, "
      f"{int((~avail[:N]).sum())} label gaps interpolated")

# ---------- predictions ----------
import torch
from ase.io import read
from mace import data
from mace.tools import torch_geometric, utils
dev = ("cuda" if torch.cuda.is_available() else "cpu") if a.device == "auto" else a.device
model = torch.load(a.model, map_location=dev, weights_only=False).to(dev).eval()
zt = utils.AtomicNumberTable([int(z) for z in model.atomic_numbers]); rmax = float(model.r_max)
L = float(re.search(r"ABC\s+([\d.]+)", open(f"{dirs[0]}/dip.inp").read()).group(1))
P = np.empty((N, 3)); M = np.full((N, 3), np.nan)
for i in range(N):
    at = read(os.path.join(a.dir, f"f{i:04d}", "coord.xyz")); at.set_cell([L]*3); at.set_pbc(True)
    ds = [data.AtomicData.from_config(data.config_from_atoms(at), z_table=zt, cutoff=rmax)]
    b = next(iter(torch_geometric.dataloader.DataLoader(ds, batch_size=1))).to(dev)
    with torch.no_grad():
        P[i] = model(b.to_dict(), compute_force=False)["dipole"].cpu().numpy().reshape(-1, 3)[0]
    if i in lab: M[i] = lab[i]
t = np.arange(N) * a.dt
for c in range(3):
    ok = ~np.isnan(M[:, c]); M[:, c] = np.interp(t, t[ok], M[ok, c])
jl, jp = np.abs(np.diff(M, axis=0)).max(), np.abs(np.diff(P, axis=0)).max()
print(f"largest frame-to-frame change: label {jl:.3f} e*A, prediction {jp:.3f} e*A")
if jl > 1.0: print("WARNING: label jump > 1 e*A (charge transfer or wrap) — inspect before using the spectrum")
nat = len(at); f = 4803.2 / nat
R = P - M
print(f"RMS per atom (mDebye): label fluctuation {np.sqrt(((M-M.mean(0))**2).sum(1).mean())*f:.2f}, "
      f"residual {np.sqrt(((R-R.mean(0))**2).sum(1).mean())*f:.2f}, mean offset {np.linalg.norm(R.mean(0))*f:.2f}")

# ---------- band-resolved absorption ratios ----------
C = 2.99792458e-5; HCK = 1.438777
bands = [("libration", 300, 1000), ("bend", 1400, 1800), ("window", 1900, 3100),
         ("stretch", 3200, 4000), ("overtone", 4200, 8000)]
def spec(x):
    x = x - x.mean(0); F = np.fft.rfft(x * np.hanning(len(x))[:, None], axis=0)
    Pw = (np.abs(F)**2).sum(1); nu = np.fft.rfftfreq(len(x), d=a.dt) / C
    s = a.smooth / 2.3548 / (nu[1] - nu[0]); k = int(4*s) + 1
    g = np.exp(-0.5*(np.arange(-k, k+1)/s)**2); g /= g.sum()
    return nu, np.convolve(Pw, g, mode="same")
def ratios(Mx, Px):
    nu, Ps = spec(Mx); _, Pp = spec(Px); _, Pr = spec(Px - Mx); w = nu * (1 - np.exp(-HCK*nu/a.T))
    out = {n: ((w*Pp)[(nu>=lo)&(nu<=hi)].sum() / (w*Ps)[(nu>=lo)&(nu<=hi)].sum(),
               (w*Pr)[(nu>=lo)&(nu<=hi)].sum() / (w*Ps)[(nu>=lo)&(nu<=hi)].sum()) for n, lo, hi in bands}
    return out, nu, w*Ps, w*Pp, w*Pr
full, nu, As, Ap, Ar = ratios(M, P); h = N // 2
h1 = ratios(M[:h], P[:h])[0]; h2 = ratios(M[h:2*h], P[h:2*h])[0]
tot = As[(nu >= 300) & (nu <= 8000)].sum()
print(f"\n{'band':10s} {'cm^-1':>11s} {'A_pred/A_label':>15s} {'A_resid/A_label':>16s} {'halves (pred/label)':>21s} {'share of label absorption':>26s}")
for n, lo, hi in bands:
    sh = As[(nu >= lo) & (nu <= hi)].sum() / tot
    print(f"{n:10s} {lo:5d}-{hi:<5d} {full[n][0]:15.3f} {full[n][1]:16.4f} {h1[n][0]:10.3f} / {h2[n][0]:.3f} {sh:26.3f}")
np.savez(a.out + ".npz", nu=nu, A_label=As, A_pred=Ap, A_resid=Ar, M_label=M, M_pred=P, N=N, dt=a.dt)
try:
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.semilogy(nu, As, lw=1.2, label="DFT Wannier (label)"); ax.semilogy(nu, Ap, lw=1, label="dipole model")
    ax.semilogy(nu, Ar, lw=1, label="residual (model − label)"); ax.axvspan(1900, 3100, alpha=0.12, label="window")
    ax.set_xlim(0, 8000); ax.set_xlabel("wavenumber (cm$^{-1}$)"); ax.set_ylabel(r"$\nu(1-e^{-hc\nu/kT})|M(\nu)|^2$ (arb.)")
    ax.legend(fontsize=8); fig.tight_layout(); fig.savefig(a.out + ".png", dpi=150); print(f"figure -> {a.out}.png")
except Exception as e:
    print("figure skipped:", e)
print("\nReading: window A_pred/A_label near 1 with A_resid/A_label << 1 -> the model reproduces the DFT window absorption at this state point.")
