#!/usr/bin/env python3
"""Generate CP2K single-point inputs across a cutoff scan.
The template must contain {cutoff}, {relcut}, {L}, {xyz} placeholders
(see the campaign label_template.inp for an example)."""
import argparse, pathlib, shutil
from ase.io import read

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("geometry", help="xyz with cell info (or use --box)")
    ap.add_argument("--template", required=True)
    ap.add_argument("--cutoffs", nargs="+", type=int, required=True)
    ap.add_argument("--rel-cutoff", type=int, default=60)
    ap.add_argument("--box", type=float, default=None, help="cubic L in A if xyz has no cell")
    ap.add_argument("--out", default="cutoff_scan")
    args = ap.parse_args()
    atoms = read(args.geometry)
    L = args.box or atoms.cell.lengths()[0]
    if L == 0:
        raise SystemExit("no cell in xyz — pass --box")
    tpl = pathlib.Path(args.template).read_text()
    root = pathlib.Path(args.out)
    for c in args.cutoffs:
        d = root / f"c{c:04d}"
        d.mkdir(parents=True, exist_ok=True)
        shutil.copy(args.geometry, d / "geom.xyz")
        (d / "sp.inp").write_text(tpl.format(cutoff=c, relcut=args.rel_cutoff,
                                             L=f"{L:.6f}", xyz="geom.xyz"))
    print(f"{len(args.cutoffs)} inputs under {root}/ — run each "
          f"(cp2k.psmp -i sp.inp -o sp.out), then analyze_cutoff_scan.py {root}")

if __name__ == "__main__":
    main()
