#!/usr/bin/env python3
"""Standalone H-bonded dimer +-D3 scan generator (no dependence on
make_dispersion_scan.py). Clones the concrete hybrid reference input,
substitutes ABC per geometry, strips &VDW_POTENTIAL for the no-D3 twin.

  python3 gen_hb_dimers.py --reference hybrid_label.inp --out scan_hb
"""
import argparse
import pathlib
import re
import numpy as np

SEPS = np.concatenate([np.arange(2.4, 3.61, 0.1), np.arange(3.8, 6.01, 0.2)])


def dimer_lines(d):
    return ["6", f"water dimer (H-bonded) O-O = {d:.2f} A",
            "O 0.000000 0.000000 0.000000",
            "H 0.757200 0.586500 0.000000",
            "H -0.757200 0.586500 0.000000",
            f"O 0.000000 0.000000 {d:.6f}",
            f"H 0.000000 0.000000 {d - 0.957200:.6f}",
            f"H 0.926500 0.000000 {d + 0.239900:.6f}"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reference", required=True)
    ap.add_argument("--out", default="scan_hb")
    args = ap.parse_args()
    tpl = pathlib.Path(args.reference).read_text()
    if "&VDW_POTENTIAL" not in tpl:
        raise SystemExit("reference lacks &VDW_POTENTIAL (D3)")
    cm = re.search(r"^\s*COORD_FILE_NAME\s+(\S+)", tpl, re.M)
    if cm is None:
        raise SystemExit("reference has no COORD_FILE_NAME")
    coord_name = cm.group(1)
    nod3 = re.sub(r"[ \t]*&VDW_POTENTIAL.*?&END VDW_POTENTIAL\n?", "",
                  tpl, flags=re.S)
    root = pathlib.Path(args.out)
    n = 0
    for d in SEPS:
        for tag, text in (("d3", tpl), ("nod3", nod3)):
            w = root / "dimer" / f"r{d:.2f}_{tag}"
            w.mkdir(parents=True, exist_ok=True)
            (w / coord_name).write_text("\n".join(dimer_lines(d)) + "\n")
            inp = re.sub(r"^(\s*ABC\s+)[\d.Ee+-]+\s+[\d.Ee+-]+\s+[\d.Ee+-]+",
                         r"\g<1>15.000000 15.000000 15.000000", text,
                         count=1, flags=re.M)
            (w / "sp.inp").write_text(inp)
            n += 1
    print(f"{n} inputs under {root}/ (coord file: {coord_name})")


if __name__ == "__main__":
    main()
