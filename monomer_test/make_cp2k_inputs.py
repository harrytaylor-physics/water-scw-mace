#!/usr/bin/env python3
"""Write g*/dip.inp and g*/coord.xyz for the 13 monomer CP2K/Wannier runs (A64).

    python make_cp2k_inputs.py --template ~/scw/mlwf_test/mlwf.inp

The template is taken as given and edited in exactly one place — the cell — using
the same substitution ``make_dipole_labels_v3.py`` uses:

    re.sub(r"ABC\\s+[\\d. ]+", "ABC 20.0000 20.0000 20.0000", text)

Nothing else about the template is parsed, assumed or rewritten. That is the
point: the DFT side of this test has to be the *same* calculation that produced
the training labels, or a model/DFT disagreement could just be a settings
difference. The template is not reproduced here and must be the one on ARC.

Prints the geometry list for the sbatch array.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
BOX_A = 20.0
ABC_PATTERN = r"ABC\s+[\d. ]+"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--template", type=Path, required=True,
                    help="CP2K MLWF input to start from (ARC: ~/scw/mlwf_test/mlwf.inp)")
    ap.add_argument("--geoms", type=Path, default=HERE / "geoms.json")
    ap.add_argument("--out-dir", type=Path, default=HERE)
    args = ap.parse_args(argv)

    template_path = args.template.expanduser()
    if not template_path.exists():
        raise SystemExit(
            f"{template_path} not found. This is the CP2K MLWF template on ARC; its "
            f"contents are not reproduced here and must not be guessed — a settings "
            f"difference would show up as a model/DFT disagreement.")
    template = template_path.read_text()
    if not re.search(ABC_PATTERN, template):
        raise SystemExit(
            f"{template_path} has no 'ABC <numbers>' line to substitute; the cell "
            f"cannot be set the way make_dipole_labels_v3.py sets it.")

    abc = f"ABC {BOX_A:.4f} {BOX_A:.4f} {BOX_A:.4f}"
    text = re.sub(ABC_PATTERN, abc, template)

    meta = json.loads(args.geoms.read_text())
    names = []
    for entry in meta["geometries"]:
        name = entry["name"]
        d = args.out_dir / name
        d.mkdir(parents=True, exist_ok=True)
        (d / "dip.inp").write_text(text)
        (d / "coord.xyz").write_text((args.geoms.parent / entry["file"]).read_text())
        names.append(name)

    print(f"template : {template_path}")
    print(f"cell set : {abc}")
    print(f"wrote    : {len(names)} x (dip.inp, coord.xyz) under {args.out_dir}")
    print(f"\nsbatch array indices: 0-{len(names) - 1}")
    print("geometries: " + " ".join(names))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
