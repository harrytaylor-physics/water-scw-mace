#!/usr/bin/env python
"""
residual_spectrum.py — measure the dipole-model noise floor directly.

Inputs:  CSV of Wannier total dipoles (step,time_fs,Mx_D,My_D,Mz_D,...) from a CP2K AIMD run,
         the matching CP2K positions file (*-pos-1.xyz, comment lines contain 'i = <step>'),
         the box length L [A], and the dipole model (.model, AtomicDipolesMACE).
Output:  per-band ratio of residual power to signal power, flatness test, and a figure.

usage: python residual_spectrum.py --csv FILE --xyz FILE --L 21.21 --model water_dipole_v3.model [--smooth 20] [--out residual_spectrum]
"""
import argparse, re, sys, numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("--csv", required=True); ap.add_argument("--xyz", required=True)
ap.add_argument("--L", type=float, required=True, help="cubic box length in Angstrom")
ap.add_argument("--model", default="water_dipole_v3.model")
ap.add_argument("--smooth", type=float, default=20.0, help="Gaussian smoothing FWHM in cm^-1")
ap.add_argument("--device", default="cuda"); ap.add_argument("--out", default="residual_spectrum")
ap.add_argument("--max_frames", type=int, default=0, help="0 = all")
a = ap.parse_args()

EA_TO_DEBYE = 4.80320  # 1 e*Angstrom = 4.80320 Debye
C_CM_PER_FS = 2.99792458e-5  # speed of light, cm/fs

# ---------- labels ----------
lab = np.genfromtxt(a.csv, delimiter=",", names=True, dtype=None, encoding="utf-8")
lab_step = lab["step"].astype(int); lab_t = lab["time_fs"].astype(float)
M_lab = np.stack([lab["Mx_D"], lab["My_D"], lab["Mz_D"]], 1).astype(float)
dt = float(np.median(np.diff(lab_t)))
print(f"labels: {len(lab_step)} frames, dt={dt:.3f} fs, steps {lab_step[0]}..{lab_step[-1]}")

# ---------- positions (CP2K xyz; comment line has 'i = N') ----------
frames = {}
with open(a.xyz) as f:
    while True:
        hdr = f.readline()
        if not hdr: break
        nat = int(hdr.split()[0]); cm = f.readline()
        m = re.search(r"i\s*=\s*(\d+)", cm); step = int(m.group(1)) if m else len(frames)
        sym = []; pos = np.empty((nat, 3))
        for k in range(nat):
            p = f.readline().split(); sym.append(p[0]); pos[k] = [float(x) for x in p[1:4]]
        frames[step] = (sym, pos)
common = [s for s in lab_step if s in frames]
if a.max_frames: common = common[:a.max_frames]
print(f"positions: {len(frames)} frames; matched to labels: {len(common)}")
if len(common) < 1000: sys.exit("fewer than 1000 matched frames — check step alignment / file paths")

# ---------- dipole model inference (standard MACE API) ----------
import torch
from ase import Atoms
from mace import data
from mace.tools import torch_geometric, utils
model = torch.load(a.model, map_location=a.device, weights_only=False).to(a.device).eval()
z_table = utils.AtomicNumberTable([int(z) for z in model.atomic_numbers]); r_max = float(model.r_max)
M_pred = np.empty((len(common), 3))
for i, s in enumerate(common):
    sym, pos = frames[s]
    at = Atoms(symbols=sym, positions=pos, cell=[a.L]*3, pbc=True)
    cfg = data.config_from_atoms(at)
    ds = [data.AtomicData.from_config(cfg, z_table=z_table, cutoff=r_max)]
    batch = next(iter(torch_geometric.dataloader.DataLoader(ds, batch_size=1, shuffle=False))).to(a.device)
    with torch.no_grad():
        out = model(batch.to_dict(), compute_force=False)
    d = out["dipole"].detach().cpu().numpy().reshape(-1, 3)[0] * EA_TO_DEBYE
    M_pred[i] = d
    if i % 2000 == 0: print(f"  {i}/{len(common)}", flush=True)
idx = np.array([np.where(lab_step == s)[0][0] for s in common]); M_l = M_lab[idx]

# ---------- residual, with Berry-phase branch check ----------
r = M_pred - M_l
r -= r.mean(0)  # constant offset is irrelevant to the spectrum
jumps = np.abs(np.diff(r, axis=0)).max()
print(f"residual RMS per component: {r.std(0)} D ; label fluctuation RMS: {(M_l-M_l.mean(0)).std(0)} D ; max frame-to-frame jump {jumps:.3f} D")
Mf = M_l - M_l.mean(0)
n_atoms = len(frames[common[0]][0])
print(f"residual RMS / atom: {np.sqrt((r**2).sum(1).mean())/n_atoms*1e3:.2f} mDebye/atom ; signal RMS / atom: {np.sqrt((Mf**2).sum(1).mean())/n_atoms*1e3:.2f} mDebye/atom")

# ---------- power spectra (same construction for both: Hann-windowed periodogram, Gaussian smooth) ----------
def power(x):
    x = x - x.mean(0); w = np.hanning(len(x))[:, None]
    F = np.fft.rfft(x * w, axis=0); P = (np.abs(F)**2).sum(1)
    nu = np.fft.rfftfreq(len(x), d=dt) / C_CM_PER_FS  # cm^-1
    if a.smooth > 0:
        sig = a.smooth / 2.3548 / (nu[1]-nu[0]); k = int(6*sig)+1
        g = np.exp(-0.5*(np.arange(-k, k+1)/sig)**2); g /= g.sum()
        P = np.convolve(P, g, mode="same")
    return nu, P
nu, P_sig = power(Mf); _, P_res = power(r)
bands = {"libration 300-1000": (300, 1000), "bend 1400-1800": (1400, 1800), "window 1900-3100": (1900, 3100),
         "stretch 3200-4000": (3200, 4000), "far-window 4200-6000": (4200, 6000)}
print("\nband                     residual/signal power   signal>2*residual fraction")
for name, (lo, hi) in bands.items():
    m = (nu >= lo) & (nu <= hi)
    ratio = P_res[m].sum() / P_sig[m].sum(); frac = np.mean(P_sig[m] > 2*P_res[m])
    print(f"{name:24s} {ratio:10.3f}              {frac:6.2f}")
# flatness: white noise gives residual power independent of frequency; correlated error tracks the signal
m = (nu >= 300) & (nu <= 6000)
corr = np.corrcoef(np.log(P_res[m]), np.log(P_sig[m]))[0, 1]
flat = P_res[m].std() / P_res[m].mean()
print(f"\nlog-log correlation(residual power, signal power) = {corr:.3f}   (white noise ~0, multiplicative error ~1)")
print(f"residual power coefficient of variation over 300-6000 cm^-1 = {flat:.2f}   (white noise ~small)")
np.savez(a.out + ".npz", nu=nu, P_sig=P_sig, P_res=P_res, M_pred=M_pred, M_lab=M_l, dt=dt)
try:
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.semilogy(nu, P_sig, lw=1, label="signal |M(ν)|²"); ax.semilogy(nu, P_res, lw=1, label="residual |M_pred−M_label|²")
    ax.axvspan(1900, 3100, alpha=0.15, label="window"); ax.set_xlim(0, 6000); ax.set_xlabel("wavenumber (cm⁻¹)")
    ax.set_ylabel("power (D², arb.)"); ax.legend(); fig.tight_layout(); fig.savefig(a.out + ".png", dpi=150)
    print(f"figure -> {a.out}.png")
except Exception as e:
    print("figure skipped:", e)
