#!/usr/bin/env python3
"""Step E part 3: on-trajectory force RMSE vs the GE.1 gate."""
import argparse, json, pathlib, re
import numpy as np

H_EV, BOHR = 27.211386245988, 0.529177210903
PARITY = {"bulk": 12.3, "rcb": 20.2}
FACTOR = 1.5

def regime(T):
    return "bulk" if T <= 2250 else "rcb"

def parse_forces(f):
    t = pathlib.Path(f).read_text(errors="ignore")
    if "SCF run converged" not in t: return None
    m = re.search(r"ATOMIC FORCES.*?\n(.*?)SUM OF ATOMIC FORCES", t, re.S)
    rows = [[float(x) for x in l.split()[3:6]] for l in m.group(1).splitlines()
            if l.split() and l.split()[0].isdigit()]
    return np.array(rows) * H_EV / BOHR

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--model", required=True)
    ap.add_argument("--device", default="cuda")
    ap.add_argument("--subdir", default="anchor")
    args = ap.parse_args()
    from ase.io import read
    from mace.calculators import MACECalculator
    calc = MACECalculator(model_paths=args.model, device=args.device,
                          default_dtype="float64")
    work = pathlib.Path(args.workdir)
    errs = {"bulk": [], "rcb": []}
    for d in sorted((work / args.subdir).glob("T*_rho*")):
        T = int(re.match(r"T(\d+)_", d.name).group(1))
        for w in sorted(d.glob("anchor_*")):
            if not w.is_dir(): continue
            Fref = parse_forces(w / "label.out")
            if Fref is None:
                print(f"{w}: unconverged/missing - skipped"); continue
            a = read(w / "frame.xyz")
            a.cell = [float((w / "label.inp").read_text().split("ABC")[1].split()[0])]*3
            a.pbc = True; a.calc = calc
            dF = a.get_forces() - Fref
            errs[regime(T)].append(dF.ravel())
            errs.setdefault(d.name, []).append(dF.ravel())
    report, worst = {}, 0.0
    print("on-trajectory force RMSE vs GE.1 gate (<= 1.5x parity):")
    for reg, lst in errs.items():
        if not lst or reg not in PARITY: continue
        rmse = float(np.sqrt(np.mean(np.concatenate(lst)**2))*1000)
        gate = PARITY[reg]*FACTOR
        ratio = rmse/PARITY[reg]
        worst = max(worst, ratio)
        report[reg] = {"rmse_meV_A": rmse, "parity": PARITY[reg],
                       "gate": gate, "pass": rmse <= gate}
        print(f"  {reg:>5}: {rmse:6.1f}  (parity {PARITY[reg]}, gate {gate:.1f}) "
              f"{'PASS' if rmse <= gate else 'FAIL'}")
    print("\nper-point RMSE (meV/A):")
    for k, lst in errs.items():
        if k in PARITY: continue
        rmse = float(np.sqrt(np.mean(np.concatenate(lst)**2))*1000)
        print(f"  {k:>22}: {rmse:7.1f}  {'PASS' if rmse <= PARITY['rcb']*FACTOR else 'FAIL'} (gate {PARITY['rcb']*FACTOR:.1f})")
    verdict = "SKIP step D - proceed to production (F)" if worst <= FACTOR \
              else "RUN step D - on-trajectory ensemble differs"
    print(f"\nGE.1 verdict: {verdict}")
    (work / args.subdir / "anchor_report.json").write_text(
        json.dumps({"report": report, "verdict": verdict}, indent=2))

if __name__ == "__main__":
    main()
