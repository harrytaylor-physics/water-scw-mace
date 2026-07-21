#!/usr/bin/env python3
"""Parse a cutoff scan; table + figure + recommendation vs the highest cutoff."""
import re, sys, pathlib
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HARTREE_EV, BOHR = 27.211386245988, 0.529177210903
E_TOL_MEV_ATOM, F_TOL_MEV_A = 1.0, 10.0

def parse(f):
    txt = f.read_text(errors="ignore")
    m = re.findall(r"ENERGY\|\s+Total FORCE_EVAL.*?:\s+([-\d.]+)", txt)
    E = float(m[-1]) if m else None
    m = re.search(r"ATOMIC FORCES.*?\n(.*?)SUM OF ATOMIC FORCES", txt, re.S)
    F = None
    if m:
        rows = [[float(x) for x in l.split()[3:6]] for l in m.group(1).splitlines()
                if l.split() and l.split()[0].isdigit()]
        F = np.array(rows) * HARTREE_EV / BOHR
    return E, F, "SCF run converged" in txt

def main(root):
    root = pathlib.Path(root)
    runs = {}
    for d in sorted(root.glob("c*")):
        outs = sorted(d.glob("*.out"))
        if outs:
            runs[int(d.name[1:])] = parse(outs[-1])
    if len(runs) < 2:
        sys.exit("need >=2 completed runs")
    ref_c = max(runs)
    Eref, Fref, _ = runs[ref_c]
    nat = len(Fref)
    rows, rec = [], None
    for c in sorted(runs):
        E, F, conv = runs[c]
        dE = abs(E - Eref) * HARTREE_EV * 1000 / nat
        dF = np.abs(F - Fref).max() * 1000
        ok = dE <= E_TOL_MEV_ATOM and dF <= F_TOL_MEV_A and conv
        if ok and rec is None and c < ref_c:
            rec = c
        rows.append((c, dE, dF, ok))
        print(f"{c:5d} Ry  dE {dE:8.3f} meV/atom  max dF {dF:8.2f} meV/A  "
              f"{'PASS' if ok else '----'}")
    print(f"\nreference: {ref_c} Ry;  recommended: {rec} Ry" if rec
          else f"\nno cutoff below {ref_c} converged — extend the scan")
    cs = [r[0] for r in rows[:-1]]
    fig, ax = plt.subplots(1, 2, figsize=(9, 3.5))
    ax[0].semilogy(cs, [r[1] for r in rows[:-1]], "ko-")
    ax[0].axhline(E_TOL_MEV_ATOM, color="k", ls="--", lw=0.7)
    ax[0].set_xlabel("CUTOFF (Ry)"); ax[0].set_ylabel("|dE| (meV/atom)")
    ax[1].semilogy(cs, [r[2] for r in rows[:-1]], "ko-")
    ax[1].axhline(F_TOL_MEV_A, color="k", ls="--", lw=0.7)
    ax[1].set_xlabel("CUTOFF (Ry)"); ax[1].set_ylabel("max |dF| (meV/A)")
    fig.tight_layout(); fig.savefig(root / "cutoff_convergence.png", dpi=170)
    print(f"figure -> {root/'cutoff_convergence.png'}")

if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "cutoff_scan")
