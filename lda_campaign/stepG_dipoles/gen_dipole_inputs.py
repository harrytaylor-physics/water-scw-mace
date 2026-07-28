#!/usr/bin/env python3
"""Per-frame CP2K MLWF inputs from snippets (clone of make_dipole_labels_v3
logic: ABC regex-substitution into the template, one dir per frame).

  python3 gen_dipole_inputs.py --template mlwf_lda.inp --snippets snippets \
      --out /data/.../lda_pipeline/diplabels
"""
import argparse
import glob
import os
import re
from ase.io import read, write


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--template", required=True)
    ap.add_argument("--snippets", default="snippets")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    tmpl = open(args.template).read()
    entries = []
    for pf in sorted(glob.glob(f"{args.snippets}/snip_T*_rho*.xyz")):
        stem = os.path.basename(pf)[:-4]
        frames = read(pf, index=":")
        L = float(frames[0].cell.lengths()[0]) if frames[0].cell.rank == 3 \
            else None
        if not L or L < 1:
            import re as _re
            m = _re.search(r"rho([\d.]+)", stem)
            rho = float(m.group(1))
            L = (64 * 18.01528 * 1e30 /
                 (6.02214076e23 * 1000 * rho)) ** (1 / 3)
        inp = re.sub(r"ABC\s+[\d. ]+", f"ABC {L:.4f} {L:.4f} {L:.4f}", tmpl)
        for i, at in enumerate(frames):
            d = f"{args.out}/{stem}/f{i:04d}"
            os.makedirs(d, exist_ok=True)
            write(f"{d}/coord.xyz", at, format="xyz")
            open(f"{d}/dip.inp", "w").write(inp)
            entries.append(f"{d}/dip.inp")
    lst = f"{args.out}/dip_lda.list"
    open(lst, "w").write("\n".join(entries) + "\n")
    print(f"{len(entries)} MLWF single points -> {lst}")
    print(f"submit: sbatch --array=0-{len(entries)-1}%24 dip_lda.slurm")


if __name__ == "__main__":
    main()
