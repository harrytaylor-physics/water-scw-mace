#!/usr/bin/env python3
"""EOS from existing trajectories: model stress on subsampled frames +
analytic ideal-gas kinetic term. One pass per (model, trajectory-set).

  python3 eval_eos.py --model <model> \
      --production "/data/.../scw_prod/T*_05fs" --traj-name prod.xyz \
      --stride 60 --out eos_hybrid.csv
  python3 eval_eos.py --model <water_lda_v1_swa.model> \
      --production "/data/.../lda_pipeline/production/T*" --traj-name prod.xyz \
      --stride 60 --out eos_lda.csv
GPU session required (torch).
"""
import argparse
import glob
import pathlib
import re
import numpy as np

KB = 1.380649e-23


def stream_frames(traj, stride):
    """Yield (symbols, positions) every `stride` frames, memory-safe."""
    with open(traj) as f:
        first = f.readline()
        nat = int(first.split()[0])
        f.seek(0)
        fi = -1
        while True:
            head = f.readline()
            if not head:
                return
            f.readline()
            fi += 1
            take = (fi % stride == 0)
            sym, pos = [], []
            for _ in range(nat):
                w = f.readline().split()
                if take:
                    sym.append(w[0])
                    pos.append([float(x) for x in w[1:4]])
            if take:
                yield sym, np.array(pos)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", required=True)
    ap.add_argument("--production", required=True, help="glob of point dirs")
    ap.add_argument("--traj-name", default="prod.xyz")
    ap.add_argument("--stride", type=int, default=60)
    ap.add_argument("--n-molecules", type=int, default=64)
    ap.add_argument("--device", default="cuda")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    from ase import Atoms
    from mace.calculators import MACECalculator
    calc = MACECalculator(model_paths=args.model, device=args.device,
                          default_dtype="float64")
    rows = []
    for d in sorted(glob.glob(args.production)):
        d = pathlib.Path(d)
        m = re.search(r"T(\d+)_rho([\d.]+)", d.name)
        if not m:
            continue
        T, rho = float(m.group(1)), float(m.group(2))
        traj = d / args.traj_name
        if not traj.exists():
            print(f"{d.name}: no {args.traj_name} - skip")
            continue
        L = (args.n_molecules * 18.01528 * 1e30 /
             (6.02214076e23 * 1000 * rho)) ** (1 / 3)
        V_m3 = (L * 1e-10) ** 3
        Ps = []
        for sym, pos in stream_frames(str(traj), args.stride):
            a = Atoms(sym, positions=pos, cell=[L, L, L], pbc=True)
            a.calc = calc
            s = a.get_stress(voigt=True)          # eV/A^3, ASE sign convention
            P_virial = -(s[0] + s[1] + s[2]) / 3.0 * 1.602176634e-19 / 1e-30
            P_kin = len(a) * KB * T / V_m3        # ideal-gas kinetic term (Pa)
            Ps.append(P_virial + P_kin)
        if not Ps:
            continue
        Ps = np.array(Ps) / 1e9                    # GPa
        sem = Ps.std(ddof=1) / np.sqrt(len(Ps))
        rows.append((T, rho, Ps.mean(), sem, len(Ps)))
        print(f"{d.name}: P = {Ps.mean():8.4f} +- {sem:.4f} GPa  "
              f"({len(Ps)} frames)")
    with open(args.out, "w") as f:
        f.write("T_K,rho_kg_m3,P_GPa,P_sem_GPa,n_frames\n")
        for r in rows:
            f.write(f"{r[0]:.0f},{r[1]},{r[2]:.6f},{r[3]:.6f},{r[4]}\n")
    print(f"wrote {args.out} ({len(rows)} points)")


if __name__ == "__main__":
    main()
