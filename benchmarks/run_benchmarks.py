#!/usr/bin/env python3
"""Rerunnable benchmarks for water_scw_v4. CPU-capable (slow) or GPU."""
import argparse, json, math, pathlib, re, sys, time
import numpy as np

def regime_of(config_type):
    m = re.search(r"T(\d+)", str(config_type))
    if not m: return None
    T = int(m.group(1))
    if T > 3000: return None          # outside applied envelope
    return "bulk" if T <= 2250 else "rcb"

def parity(calc, labels_path):
    from ase.io import read
    frames = read(labels_path, index=":")
    errs = {}
    for a in frames:
        a = a.copy()
        Fref = a.arrays.get("REF_forces")
        if Fref is None:
            continue
        a.calc = calc
        dF = a.get_forces() - Fref
        regime = regime_of(a.info.get("config_type", ""))
        if regime is None: continue
        errs.setdefault(regime, []).append(dF.ravel())
    out = {}
    for k, v in errs.items():
        out[k] = float(np.sqrt(np.mean(np.concatenate(v) ** 2)) * 1000)  # meV/A
    return out

def eos_point(calc, T=673.0, rho=1000.0, n_mol=64, steps=4000):
    from ase import units, Atoms
    from ase.md.langevin import Langevin
    from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
    M, NA = 18.01528, 6.02214076e23
    L = (n_mol * M / NA / 1000.0 / rho) ** (1 / 3) * 1e10
    rng = np.random.default_rng(0)
    pos, sym = [], []
    ns = math.ceil(n_mol ** (1 / 3)); sp = L / ns
    r, ang = 0.9572, math.radians(104.52)
    base = np.array([[0, 0, 0],
                     [r * math.sin(ang / 2), 0, r * math.cos(ang / 2)],
                     [-r * math.sin(ang / 2), 0, r * math.cos(ang / 2)]])
    n = 0
    for i in range(ns):
        for j in range(ns):
            for k in range(ns):
                if n >= n_mol: continue
                c = np.array([(i + .5) * sp, (j + .5) * sp, (k + .5) * sp])
                q = rng.normal(size=4); q /= np.linalg.norm(q)
                w, x, y, z = q
                R = np.array([[1-2*(y*y+z*z), 2*(x*y-w*z), 2*(x*z+w*y)],
                              [2*(x*y+w*z), 1-2*(x*x+z*z), 2*(y*z-w*x)],
                              [2*(x*z-w*y), 2*(y*z+w*x), 1-2*(x*x+y*y)]])
                for el, v in zip("OHH", base):
                    sym.append(el); pos.append(c + R @ v)
                n += 1
    atoms = Atoms(sym, positions=pos, cell=[L, L, L], pbc=True)
    atoms.calc = calc
    MaxwellBoltzmannDistribution(atoms, temperature_K=T)
    dyn = Langevin(atoms, timestep=0.5 * units.fs, temperature_K=T, friction=0.02)
    ps = []
    def grab():
        s = atoms.get_stress(voigt=True)
        V = atoms.get_volume()
        p_kin = len(atoms) * units.kB * T / V
        ps.append((-(s[0] + s[1] + s[2]) / 3 + p_kin) / units.GPa)
    dyn.attach(grab, interval=50)
    dyn.run(steps)
    burn = len(ps) // 3
    return float(np.mean(ps[burn:]))  # GPa

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", required=True)
    ap.add_argument("--labels", default="reference/train_v4.xyz.gz")
    ap.add_argument("--device", default="cuda")
    ap.add_argument("--skip-eos", action="store_true")
    args = ap.parse_args()
    from mace.calculators import MACECalculator
    calc = MACECalculator(model_paths=args.model, device=args.device,
                          default_dtype="float32")
    ref = json.load(open(pathlib.Path(__file__).parent / "reference" / "reference.json"))
    results, failed = {}, False
    if pathlib.Path(args.labels).exists():
        t0 = time.time()
        results["parity_meV_A"] = parity(calc, args.labels)
        results["parity_time_s"] = time.time() - t0
        for regime, refval in ref["parity_meV_A"].items():
            got = results["parity_meV_A"].get(regime)
            ok = got is not None and got <= refval * 1.15
            print(f"parity[{regime:>6}]: {got and f'{got:.0f}'} meV/A "
                  f"(ref {refval}, +15% tol) {'PASS' if ok else 'FAIL'}")
            failed |= not ok
    else:
        print(f"parity: SKIPPED ({args.labels} not present)")
    if not args.skip_eos:
        p = eos_point(calc)
        ref_p = ref["eos_673K_1000kg_GPa"]
        ok = abs(p - ref_p) / ref_p <= 0.10
        print(f"EOS 673K/1000: {p:.3f} GPa (IAPWS {ref_p}) {'PASS' if ok else 'FAIL'}")
        results["eos_GPa"] = p; failed |= not ok
    json.dump(results, open("benchmark_results.json", "w"), indent=2)
    sys.exit(1 if failed else 0)

if __name__ == "__main__":
    main()
