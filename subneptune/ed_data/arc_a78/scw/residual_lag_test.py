#!/usr/bin/env python
"""
residual_lag_test.py — is the dipole-model error white noise or a smooth function of configuration?

Uses the Wannier-labelled snippets (25 segments x 8 frames @ 1 fs per cell). Within each segment
computes the normalised autocorrelation C(tau)/C(0), tau = 0..7 fs, for
  signal   : M_label - <M_label>
  residual : M_pred  - M_label  (mean-removed)
White noise gives residual C(1 fs)/C(0) ~ 0. An error that is a smooth function of the atomic
positions gives C(1 fs)/C(0) close to the signal's value and follows its oscillation.

usage: python residual_lag_test.py --xyz dipole_train_snip.xyz --model water_dipole_v3.model
       [--valid water_dipole_v3_valid_indices_1.txt] [--out lagtest]
"""
import argparse, json, re, numpy as np
ap = argparse.ArgumentParser()
ap.add_argument("--xyz", default="dipole_train_snip.xyz")
ap.add_argument("--model", default="water_dipole_v3.model")
ap.add_argument("--valid", default="")
ap.add_argument("--device", default="cuda"); ap.add_argument("--out", default="lagtest")
a = ap.parse_args()

# ---------- read extxyz: positions, cell, REF_dipole, config_type ----------
frames = []
with open(a.xyz) as f:
    while True:
        h = f.readline()
        if not h.strip(): break
        n = int(h); c = f.readline()
        L = float(re.search(r'Lattice="([^ ]+)', c).group(1))
        dip = np.array(json.loads(re.search(r'REF_dipole="_JSON (\[[^\]]*\])"', c).group(1)))
        ct = re.search(r"config_type=(\S+)", c).group(1)
        sym, pos = [], np.empty((n, 3))
        for k in range(n):
            p = f.readline().split(); sym.append(p[0]); pos[k] = [float(x) for x in p[1:4]]
        frames.append((ct, L, sym, pos, dip))
print(f"{len(frames)} labelled frames; config types: {sorted(set(fr[0] for fr in frames))}")
valid = set(int(x) for x in open(a.valid).read().split()) if a.valid else set()

# ---------- segment detection: consecutive frames 1 fs apart move < 0.5 A (min image) ----------
segs, cur = [], [0]
for i in range(1, len(frames)):
    same = frames[i][0] == frames[i-1][0]
    if same:
        d = frames[i][3] - frames[i-1][3]; Lb = frames[i][1]; d -= Lb*np.round(d/Lb)
        same = np.abs(d).max() < 0.5
    if same: cur.append(i)
    else: segs.append(cur); cur = [i]
segs.append(cur)
lens = np.bincount([len(s) for s in segs])
print(f"segments detected: {len(segs)}; length histogram (index=length): {dict((i,int(v)) for i,v in enumerate(lens) if v)}")

# ---------- model predictions (same units as REF_dipole, e*A) ----------
import torch
from ase import Atoms
from mace import data
from mace.tools import torch_geometric, utils
model = torch.load(a.model, map_location=a.device, weights_only=False).to(a.device).eval()
zt = utils.AtomicNumberTable([int(z) for z in model.atomic_numbers]); rmax = float(model.r_max)
pred = np.empty((len(frames), 3))
for i, (ct, L, sym, pos, dip) in enumerate(frames):
    at = Atoms(symbols=sym, positions=pos, cell=[L]*3, pbc=True)
    ds = [data.AtomicData.from_config(data.config_from_atoms(at), z_table=zt, cutoff=rmax)]
    b = next(iter(torch_geometric.dataloader.DataLoader(ds, batch_size=1))).to(a.device)
    with torch.no_grad():
        pred[i] = model(b.to_dict(), compute_force=False)["dipole"].cpu().numpy().reshape(-1, 3)[0]
lab = np.array([fr[4] for fr in frames]); res = pred - lab
nat = len(frames[0][2])

# ---------- per config type: RMS and lag autocorrelation over full 8-frame segments ----------
out = {}
print(f"\n{'cell':22s} {'resid RMS':>10s} {'signal RMS':>10s}  (mDebye/atom)   lag(fs):  " + " ".join(f"{t:5d}" for t in range(8)))
for ct in sorted(set(fr[0] for fr in frames)):
    S = [s for s in segs if frames[s[0]][0] == ct and len(s) == 8]
    if not S: continue
    idx = np.concatenate(S)
    sig = lab[idx] - lab[idx].mean(0); r = res[idx] - res[idx].mean(0)
    rms = lambda x: np.sqrt((x**2).sum(1).mean()) / nat * 4803.2  # e*A -> mDebye, per atom
    def acf(x):  # x: (nseg*8, 3) -> averaged normalised ACF over segments
        x = x.reshape(len(S), 8, 3); C = np.array([np.mean((x[:, :8-t] * x[:, t:]).sum(2)) for t in range(8)])
        return C / C[0]
    Cs, Cr = acf(sig), acf(r)
    print(f"{ct:22s} {rms(r):10.2f} {rms(sig):10.2f}   signal ACF: " + " ".join(f"{v:5.2f}" for v in Cs))
    print(f"{'':22s} {'':10s} {'':10s}   resid  ACF: " + " ".join(f"{v:5.2f}" for v in Cr))
    vi = [i for i in idx if i in valid]
    if vi: print(f"{'':22s} validation frames in segments: {len(vi)}, resid RMS on them {rms(res[vi]-res[vi].mean(0)):.2f}")
    out[ct] = dict(resid_rms=rms(r), signal_rms=rms(sig), acf_signal=Cs.tolist(), acf_resid=Cr.tolist())
json.dump(out, open(a.out + ".json", "w"), indent=1)
print("\nInterpretation: resid ACF(1 fs) ~ 0  -> white noise (Task-15 floor applies).")
print("                resid ACF(1 fs) ~ signal ACF(1 fs), same oscillation -> configuration-correlated error (floor construction invalid).")
