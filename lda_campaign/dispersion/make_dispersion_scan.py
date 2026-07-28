#!/usr/bin/env python3
"""Generate CP2K inputs for the dispersion test from a CONCRETE reference
input (a converged single-point .inp - e.g. the output of
make_hybrid_template.py run against a step-B label.inp, or any proven
hybrid SP input). Only the ABC cell line and the coordinate file are
substituted per geometry; a no-D3 twin is made by stripping the
&VDW_POTENTIAL block.

  python3 make_dispersion_scan.py --reference hybrid_label.inp \
      --pool-configs bulk_sample.xyz
"""
import argparse
import pathlib
import re
import numpy as np

DIMER_SEPS = np.concatenate([np.arange(2.4, 3.61, 0.1),
                             np.arange(3.8, 6.01, 0.2)])


def dimer_lines(d):
    """Canonical H-bonded dimer: donor O-H points at acceptor O along z."""
    return ["6", f"water dimer (H-bonded) O-O = {d:.2f} A",
            "O 0.000000 0.000000 0.000000",
            "H 0.757200 0.586500 0.000000",
            "H -0.757200 0.586500 0.000000",
            f"O 0.000000 0.000000 {d:.6f}",
            f"H 0.000000 0.000000 {d - 0.957200:.6f}",
            f"H 0.926500 0.000000 {d + 0.239900:.6f}"]


def strip_d3(t):
    return re.sub(r"[ \t]*&VDW_POTENTIAL.*?&END VDW_POTENTIAL\n?", "",
                  t, flags=re.S)


def set_abc(t, L):
    return re.sub(r"^(\s*ABC\s+)[\d.Ee+-]+\s+[\d.Ee+-]+\s+[\d.Ee+-]+",
                  rf"\g<1>{L:.6f} {L:.6f} {L:.6f}", t, count=1, flags=re.M)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reference", required=True,
                    help="concrete hybrid SP input with D3 block")
    ap.add_argument("--pool-configs", default=None)
    ap.add_argument("--out", default="scan")
    ap.add_argument("--dimer-only", action="store_true")
    args = ap.parse_args()
    tpl = pathlib.Path(args.reference).read_text()
    if "&VDW_POTENTIAL" not in tpl:
        raise SystemExit("reference input lacks &VDW_POTENTIAL (D3)")
    cm = re.search(r"^\s*COORD_FILE_NAME\s+(\S+)", tpl, re.M)
    if cm is None:
        raise SystemExit("reference input has no COORD_FILE_NAME")
    coord_name = cm.group(1)
    root = pathlib.Path(args.out)
    n = 0
    for d in DIMER_SEPS:
        for tag, text in (("d3", tpl), ("nod3", strip_d3(tpl))):
            w = root / "dimer" / f"r{d:.2f}_{tag}"
            w.mkdir(parents=True, exist_ok=True)
            (w / coord_name).write_text("\n".join(dimer_lines(d)) + "\n")
            (w / "sp.inp").write_text(set_abc(text, 15.0))
            n += 2
    if args.pool_configs and not args.dimer_only:
        from ase.io import read
        for i, a in enumerate(read(args.pool_configs, index=":")):
            L = a.cell.lengths()[0]
            lines = [f"{len(a)}", a.info.get("config_type", f"frame{i}")]
            for s, p in zip(a.get_chemical_symbols(), a.positions):
                lines.append(f"{s} {p[0]:.6f} {p[1]:.6f} {p[2]:.6f}")
            for tag, text in (("d3", tpl), ("nod3", strip_d3(tpl))):
                w = root / "bulk" / f"cfg{i:03d}_{tag}"
                w.mkdir(parents=True, exist_ok=True)
                (w / coord_name).write_text("\n".join(lines) + "\n")
                (w / "sp.inp").write_text(set_abc(text, L))
                n += 2
    print(f"{n} inputs under {root}/ (coord file per input: {coord_name})")


if __name__ == "__main__":
    main()
