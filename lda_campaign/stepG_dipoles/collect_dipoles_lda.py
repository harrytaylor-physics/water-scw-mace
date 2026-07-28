#!/usr/bin/env python3
"""Collect MLWF dipoles into dipole_train_lda.xyz.
Wannier-centre arithmetic cloned verbatim from build_dipole_xyz_v3.py
(proven on the hybrid campaign); only plumbing differs:
  - reads dip_lda.list (paths to dip.inp; dirs derived)
  - config_type from the snippet-stem directory name
  python3 collect_dipoles_lda.py \
      --list /data/.../diplabels/dip_lda.list \
      --out  /data/.../lda_pipeline/dipole_train_lda.xyz
"""
import argparse
import os
import pathlib
import re
import numpy as np
from ase.io import read, write


def dip(d, wannier_name="wannier"):
    fn = f"{d}/{wannier_name}"
    if not os.path.exists(fn):
        return None
    L = float(re.search(r'ABC\s+([\d.]+)', open(f"{d}/dip.inp").read()).group(1))
    raw = [l for l in open(fn).read().split('\n') if l.strip()]
    try:
        n = int(raw[0].split()[0]); body = raw[2:2+n]
    except Exception:
        return None
    sym = []; pos = []
    for ln in body:
        p = ln.split()
        if len(p) >= 4 and re.match(r'[A-Za-z]', p[0]):
            sym.append(p[0]); pos.append([float(x) for x in p[1:4]])
    sym = np.array(sym); pos = np.array(pos)
    Oi = np.where(sym == 'O')[0]; Hi = np.where(sym == 'H')[0]
    Xi = np.where(~np.isin(sym, ['O', 'H']))[0]
    if len(Xi) != 4*len(Oi) or len(Hi) != 2*len(Oi):
        return None
    Op = pos[Oi]
    mic = lambda a, b: (lambda x: x - L*np.round(x/L))(a - b)
    aX = np.array([int(np.argmin((mic(Op, pos[x])**2).sum(1))) for x in Xi])
    Hd = []
    for hi, h in enumerate(Hi):
        r = np.sqrt((mic(Op, pos[h])**2).sum(1))
        Hd.append((float(r.min()), hi, np.argsort(r)))
    Hd.sort(); cap = np.zeros(len(Oi), int); aH = np.full(len(Hi), -1)
    for _, hi, order in Hd:
        for j in order:
            if cap[j] < 2:
                aH[hi] = j; cap[j] += 1; break
    if (cap != 2).any():
        return None
    mu = np.zeros(3)
    for o in range(len(Oi)):
        rO = Op[o]
        for hi in np.where(aH == o)[0]:
            mu += mic(pos[Hi[hi]], rO)
        for xi in np.where(aX == o)[0]:
            mu += -2.0*mic(pos[Xi[xi]], rO)
    return mu, L


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--list", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--wannier-name", default="wannier")
    args = ap.parse_args()
    out, bad = [], 0
    for line in open(args.list):
        if not line.strip():
            continue
        d = str(pathlib.Path(line.strip()).parent)
        if not os.path.exists(f"{d}/dip.out") or \
                "PROGRAM ENDED" not in open(f"{d}/dip.out").read():
            bad += 1; continue
        r = dip(d, args.wannier_name)
        if r is None:
            bad += 1; continue
        mu, L = r
        at = read(f"{d}/coord.xyz")
        at.set_cell([L, L, L]); at.set_pbc(True)
        at.info['REF_dipole'] = mu.tolist()
        at.info['config_type'] = pathlib.Path(d).parent.name
        out.append(at)
    write(args.out, out)
    print(f"{args.out}: {len(out)} configs with MLWF dipoles "
          f"({bad} skipped: unfinished/unparseable)")
    from collections import Counter
    print(dict(Counter(a.info['config_type'] for a in out)))


if __name__ == "__main__":
    main()
