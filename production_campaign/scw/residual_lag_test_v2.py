#!/usr/bin/env python
"""
residual_lag_test_v2.py — is the dipole-model error white noise or a smooth function of configuration?

Labelled frames (dipole_train_snip.xyz, REF_dipole in e*A) are matched by coordinates to the complete
snippet trajectories (snippets/<config_type>.xyz, 25 segments x 8 frames @ 1 fs, in time order), so each
labelled frame gets its exact (segment, slot). Missing labels are handled exactly: lag tau uses only
pairs (seg, slot) and (seg, slot+tau) that both exist; pair counts are printed.

usage: python residual_lag_test_v2.py --xyz dipole_train_snip.xyz --snipdir snippets --model water_dipole_v3.model --out lagtest_v2
"""
import argparse, json, os, re, numpy as np
ap = argparse.ArgumentParser()
ap.add_argument("--xyz", default="dipole_train_snip.xyz"); ap.add_argument("--snipdir", default="snippets")
ap.add_argument("--model", default="water_dipole_v3.model"); ap.add_argument("--device", default="cuda")
ap.add_argument("--seglen", type=int, default=8); ap.add_argument("--out", default="lagtest_v2")
a = ap.parse_args()

def read_extxyz(path, labelled):
    out = []
    with open(path) as f:
        while True:
            h = f.readline()
            if not h.strip(): break
            n = int(h); c = f.readline()
            L = float(re.search(r'Lattice="([^ ]+)', c).group(1))
            dip = np.array(json.loads(re.search(r'REF_dipole="_JSON (\[[^\]]*\])"', c).group(1))) if labelled else None
            ct = re.search(r"config_type=(\S+)", c).group(1) if labelled else None
            sym, pos = [], np.empty((n, 3))
            for k in range(n):
                p = f.readline().split(); sym.append(p[0]); pos[k] = [float(x) for x in p[1:4]]
            out.append(dict(ct=ct, L=L, sym=sym, pos=pos, dip=dip))
    return out

lab = read_extxyz(a.xyz, True)
cells = sorted(set(fr["ct"] for fr in lab))
ref = {ct: read_extxyz(os.path.join(a.snipdir, ct + ".xyz"), False) for ct in cells if os.path.exists(os.path.join(a.snipdir, ct + ".xyz"))}
print(f"{len(lab)} labelled frames; cells: {cells}; reference trajectories found for: {sorted(ref)}")

# ---------- exact time index by coordinate match ----------
tag = []
for fr in lab:
    hit = None
    for j, rf in enumerate(ref.get(fr["ct"], [])):
        d = fr["pos"] - rf["pos"]; d -= rf["L"] * np.round(d / rf["L"])
        if np.abs(d).max() < 1e-3: hit = (j // a.seglen, j % a.seglen); break
    tag.append(hit)
for ct in cells:
    m = sum(1 for fr, t in zip(lab, tag) if fr["ct"] == ct and t is not None)
    tot = sum(1 for fr in lab if fr["ct"] == ct)
    print(f"  {ct:22s} labelled {tot:4d}  matched {m:4d}  (of {len(ref.get(ct, []))} generated)")

# ---------- model predictions (same units as labels, e*A) ----------
import torch
from ase import Atoms
from mace import data
from mace.tools import torch_geometric, utils
model = torch.load(a.model, map_location=a.device, weights_only=False).to(a.device).eval()
zt = utils.AtomicNumberTable([int(z) for z in model.atomic_numbers]); rmax = float(model.r_max)
pred = np.empty((len(lab), 3))
for i, fr in enumerate(lab):
    at = Atoms(symbols=fr["sym"], positions=fr["pos"], cell=[fr["L"]]*3, pbc=True)
    ds = [data.AtomicData.from_config(data.config_from_atoms(at), z_table=zt, cutoff=rmax)]
    b = next(iter(torch_geometric.dataloader.DataLoader(ds, batch_size=1))).to(a.device)
    with torch.no_grad():
        pred[i] = model(b.to_dict(), compute_force=False)["dipole"].cpu().numpy().reshape(-1, 3)[0]
M = np.array([fr["dip"] for fr in lab]); R = pred - M
nat = len(lab[0]["sym"]); to_mD = 4803.2 / nat  # e*A -> mDebye per atom

# ---------- lag autocorrelation with exact pairs ----------
res = {}
hdr = " ".join(f"{t:6d}" for t in range(a.seglen))
print(f"\nlag (fs):{'':34s}{hdr}")
for ct in cells:
    ii = [i for i, (fr, t) in enumerate(zip(lab, tag)) if fr["ct"] == ct and t is not None]
    if len(ii) < 20: print(f"{ct}: too few matched frames ({len(ii)})"); continue
    sig = {tag[i]: M[i] for i in ii}; rsd = {tag[i]: R[i] for i in ii}
    ms = np.mean(list(sig.values()), 0); mr = np.mean(list(rsd.values()), 0)
    def acf(d, mu):
        C, npair = [], []
        for t in range(a.seglen):
            prs = [(np.dot(d[(s, k)] - mu, d[(s, k + t)] - mu)) for (s, k) in d if (s, k + t) in d]
            C.append(np.mean(prs) if prs else np.nan); npair.append(len(prs))
        C = np.array(C); return C / C[0], npair
    Cs, ns = acf(sig, ms); Cr, _ = acf(rsd, mr)
    rs = np.sqrt(np.mean([np.sum((v - mr)**2) for v in rsd.values()])) * to_mD
    ss = np.sqrt(np.mean([np.sum((v - ms)**2) for v in sig.values()])) * to_mD
    print(f"{ct:22s} resid {rs:5.2f} / signal {ss:5.2f} mD/atom")
    print(f"{'':22s} signal ACF:          " + " ".join(f"{v:6.2f}" for v in Cs))
    print(f"{'':22s} resid  ACF:          " + " ".join(f"{v:6.2f}" for v in Cr))
    print(f"{'':22s} pairs:               " + " ".join(f"{v:6d}" for v in ns))
    res[ct] = dict(resid_mD=rs, signal_mD=ss, acf_signal=Cs.tolist(), acf_resid=Cr.tolist(), pairs=ns)
json.dump(res, open(a.out + ".json", "w"), indent=1)
print("\nresid ACF(1 fs) ~ 0 -> white noise.  resid ACF(1 fs) ~ signal ACF(1 fs), same oscillation -> configuration-correlated error.")
