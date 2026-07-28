#!/usr/bin/env python3
"""Analyse the dispersion scan: binding curves, D3-vs-distance, bulk
intra/inter force projection, H-bond frequency shift."""
import re, sys, pathlib
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

H_EV, BOHR = 27.211386245988, 0.529177210903

def parse(f):
    p = pathlib.Path(f)
    if not p.exists():
        return None, None, False
    t = p.read_text(errors="ignore")
    Em = re.findall(r"ENERGY\|\s+Total FORCE_EVAL.*?:\s+([-\d.]+)", t)
    Fm = re.search(r"ATOMIC FORCES.*?\n(.*?)SUM OF ATOMIC FORCES", t, re.S)
    F = None
    if Fm:
        rows = [[float(x) for x in l.split()[3:6]] for l in Fm.group(1).splitlines()
                if l.split() and l.split()[0].isdigit()]
        F = np.array(rows) * H_EV / BOHR
    return (float(Em[-1]) * H_EV if Em else None), F, "SCF run converged" in t


def read_coords(d):
    for name in ("frame.xyz", "geom.xyz"):
        p = pathlib.Path(d) / name
        if p.exists():
            return p.read_text().splitlines()[2:]
    raise FileNotFoundError(f"no frame.xyz/geom.xyz in {d}")

def first_out(d):
    outs = sorted(pathlib.Path(d).glob("*.out"))
    return outs[0] if outs else pathlib.Path(d) / "missing.out"

def intra_inter_split(xyz, dF):
    intra2 = inter2 = 0.0
    for m in range(len(xyz) // 3):
        O = xyz[3*m]
        for Hi, fH in ((xyz[3*m+1], dF[3*m+1]), (xyz[3*m+2], dF[3*m+2])):
            u = Hi - O; u /= np.linalg.norm(u)
            a = float(np.dot(fH, u))
            intra2 += a*a
            inter2 += float(np.dot(fH, fH)) - a*a
        inter2 += float(np.dot(dF[3*m], dF[3*m]))
    return np.sqrt(intra2), np.sqrt(inter2)

def dimer_curve(root):
    rs, e1, e0 = [], [], []
    for d in sorted(root.glob("dimer/r*_d3")):
        r = float(d.name[1:-3])
        E1, _, c1 = parse(first_out(d))
        E0, _, c0 = parse(first_out(str(d).replace("_d3", "_nod3")))
        if c1 and c0:
            rs.append(r); e1.append(E1); e0.append(E0)
    if len(rs) < 5:
        print("dimer: not enough converged pairs yet"); return
    rs = np.array(rs); e1 = np.array(e1); e0 = np.array(e0)
    e1 -= e1[-1]; e0 -= e0[-1]
    fig, ax = plt.subplots(1, 2, figsize=(9.5, 3.6))
    ax[0].plot(rs, e1*1000, "k-", label="revPBE0-D3")
    ax[0].plot(rs, e0*1000, "k--", label="revPBE0 (no D3)")
    ax[0].set_xlabel("O-O separation (A)")
    ax[0].set_ylabel("binding energy (meV)")
    ax[0].legend(frameon=False)
    ax[1].plot(rs, (e1-e0)*1000, "k-")
    ax[1].set_xlabel("O-O separation (A)")
    ax[1].set_ylabel("D3 contribution (meV)")
    fig.tight_layout(); fig.savefig(root / "dimer_curve.png", dpi=170)
    if np.argmin(e1) == len(rs) - 1:
        print("dimer: NO interior minimum - stacked-parallel geometry is "
              "repulsive (parallel dipoles); D3-contribution panel remains "
              "valid; rerun with --hbond geometry for binding/frequency")
    for tag, E in (("+D3", e1), ("noD3", e0)):
        i = int(np.argmin(E))
        msg = f"  minimum {rs[i]:.2f} A, {E[i]*1000:.0f} meV"
        if 0 < i < len(rs)-1 and abs(rs[i+1]-rs[i]-(rs[i]-rs[i-1])) < 1e-6:
            k = (E[i-1]-2*E[i]+E[i+1])/(rs[i+1]-rs[i])**2
            mu = 18.015/2*1.6605e-27
            w = np.sqrt(max(k,0)*1.602e-19*1e20/mu)/(2*np.pi*2.9979e10)
            msg += f", H-bond stretch ~{w:.0f} cm-1"
        print(f"dimer {tag}:{msg}")
    print(f"D3 deepens the well by {((e0[np.argmin(e0)])-(e1[np.argmin(e1)]))*1000:.0f} meV")

def bulk_fractions(root):
    rows = []
    print(f"\n{'config':>10} {'|F_D3|/|F| %':>13} {'intra %':>9} {'inter %':>9}")
    for d in sorted(root.glob("bulk/cfg*_d3")):
        _, F1, c1 = parse(first_out(d))
        _, F0, c0 = parse(first_out(str(d).replace("_d3", "_nod3")))
        if not (c1 and c0) or F1 is None or F0 is None: continue
        xyz = np.array([l.split()[1:4] for l in
                        read_coords(d)], float)
        dF = F1 - F0
        frac = np.linalg.norm(dF)/np.linalg.norm(F1)*100
        intra, inter = intra_inter_split(xyz, dF)
        tot = np.hypot(intra, inter)
        print(f"{d.name[:-3]:>10} {frac:13.2f} {100*intra/tot:9.1f} {100*inter/tot:9.1f}")
        rows.append((frac, 100*intra/tot, 100*inter/tot))
    if rows:
        r = np.array(rows)
        print(f"{'MEAN':>10} {r[:,0].mean():13.2f} {r[:,1].mean():9.1f} {r[:,2].mean():9.1f}")
        print("hypothesis check: small intra% + dominant inter% => "
              "dispersion is structural, not vibrational")

if __name__ == "__main__":
    root = pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else "scan")
    dimer_curve(root)
    bulk_fractions(root)
