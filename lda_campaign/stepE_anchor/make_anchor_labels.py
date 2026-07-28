#!/usr/bin/env python3
"""Step E part 2: CP2K inputs for the anchor frames.

Instead of re-formatting a template (whose placeholder set we do not
control), CLONE a concrete, converged step-B label.inp - identical
numerics by construction - and substitute only:
  - the cell length (ABC line)
  - the coordinates (written under whatever filename the reference input
    already points at)

Usage:
  python3 make_anchor_labels.py --workdir /data/.../lda_pipeline
  # optionally --reference-inp <path>; default: first converged step-B input
"""
import argparse
import pathlib
import re
from ase.io import read


def find_reference(work):
    for d in sorted((work / "labels").glob("T*_rho*/sel_*")):
        out = d / "label.out"
        if out.exists() and "SCF run converged" in out.read_text(errors="ignore"):
            return d / "label.inp"
    raise SystemExit("no converged step-B label.inp found under labels/")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--reference-inp", default=None)
    ap.add_argument("--subdir", default="anchor")
    args = ap.parse_args()
    work = pathlib.Path(args.workdir)
    ref = pathlib.Path(args.reference_inp) if args.reference_inp \
        else find_reference(work)
    tpl = ref.read_text()
    print(f"reference input: {ref}")

    if not re.search(r"^\s*ABC\s+[\d.Ee+-]+", tpl, re.M):
        raise SystemExit("no ABC line found in reference input")
    cm = re.search(r"^\s*COORD_FILE_NAME\s+(\S+)", tpl, re.M)
    if cm is None:
        raise SystemExit("reference input has no COORD_FILE_NAME - inline "
                         "&COORD not supported; pass a step-B input")
    coord_name = cm.group(1)
    print(f"coordinate file expected by input: {coord_name}")

    n = 0
    for d in sorted((work / args.subdir).glob("T*_rho*")):
        for f in sorted(d.glob("anchor_*.xyz")):
            a = read(f)
            w = d / f.stem
            w.mkdir(exist_ok=True)
            L = a.cell.lengths()[0]
            lines = [f"{len(a)}", d.name]
            for s, p in zip(a.get_chemical_symbols(), a.positions):
                lines.append(f"{s} {p[0]:.6f} {p[1]:.6f} {p[2]:.6f}")
            (w / coord_name).write_text("\n".join(lines) + "\n")
            new = re.sub(r"^(\s*ABC\s+)[\d.Ee+-]+\s+[\d.Ee+-]+\s+[\d.Ee+-]+",
                         rf"\g<1>{L:.6f} {L:.6f} {L:.6f}", tpl, count=1,
                         flags=re.M)
            (w / "label.inp").write_text(new)
            n += 1
    print(f"{n} anchor label inputs ready (numerics cloned from step B)")


if __name__ == "__main__":
    main()
