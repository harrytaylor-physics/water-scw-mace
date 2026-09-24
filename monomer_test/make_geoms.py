#!/usr/bin/env python3
"""13 water-monomer geometries for the isolated-molecule dipole test (A64).

    python make_geoms.py [--out-dir .]

Why these thirteen. The dilute-limit band diagnostic (A61) found the
first-principles absorption sitting at 0.35-0.50 of the line list's in the
rotational band, 0.7-1.1 in the bend and 2.4-3.9 in the stretch, with a
stretch/bend intensity balance 3.3-3.6x the line list's — and the same at 1000,
1500 and 2000 K. A temperature-independent, mode-dependent ratio is what a wrong
dipole *derivative* looks like: ``water_dipole_v3`` was trained almost entirely
on dense configurations and may return a liquid-like dmu/dr for a stretching bond
in the gas. The isolated molecule is the direct test, and the direct test needs
only the permanent dipole and three derivatives, at two bond lengths.

Geometry set (g00..g12), all C2v-oriented with O at the box centre:

  g00              equilibrium, r_OH = 0.9572 A, HOH = 104.52 deg
  g01/g02          symmetric stretch, both bonds +/- 0.01 A
  g03/g04          asymmetric stretch, one bond +0.01, the other -0.01 A
  g05/g06          bend, HOH +/- 1.0 deg
  g07..g12         the same three pairs about a stretched reference, r_OH = 1.00 A

The stretched reference is not decoration: at 1000 K the bonds sample ~1.00 A, so
a derivative taken only at equilibrium need not be the one the spectrum uses. If
the model's error is anharmonic — right at 0.9572 A and wrong at 1.00 A — only
the second block shows it.

Orientation is fixed and identical across geometries (C2 axis along +z, molecule
in the xz plane), so the dipole *vectors* are directly comparable, not only their
magnitudes.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

BOX_A = 20.0
R_EQ_A = 0.9572
THETA_EQ_DEG = 104.52
R_STRETCHED_A = 1.00
D_R_A = 0.01
D_THETA_DEG = 1.0


def geometry(r1_a: float, r2_a: float, theta_deg: float) -> list[tuple[str, float, float, float]]:
    """O at the box centre, C2 axis along +z, molecule in the xz plane."""
    half = math.radians(theta_deg) / 2.0
    c = BOX_A / 2.0
    return [
        ("O", c, c, c),
        ("H", c + r1_a * math.sin(half), c, c + r1_a * math.cos(half)),
        ("H", c - r2_a * math.sin(half), c, c + r2_a * math.cos(half)),
    ]


def build() -> list[dict]:
    out = []
    for ref_name, r0 in (("equilibrium", R_EQ_A), ("stretched", R_STRETCHED_A)):
        if ref_name == "equilibrium":  # the undisplaced point exists only once
            out.append({"mode": "reference", "sign": 0, "reference": ref_name,
                        "r0_a": r0, "r1_a": r0, "r2_a": r0,
                        "theta_deg": THETA_EQ_DEG, "displacement": 0.0,
                        "displacement_unit": "-"})
        for mode in ("sym_stretch", "asym_stretch", "bend"):
            for sign in (+1, -1):
                r1, r2, theta = r0, r0, THETA_EQ_DEG
                if mode == "sym_stretch":
                    r1 = r2 = r0 + sign * D_R_A
                    disp, unit = sign * D_R_A, "A"
                elif mode == "asym_stretch":
                    r1, r2 = r0 + sign * D_R_A, r0 - sign * D_R_A
                    disp, unit = sign * D_R_A, "A"
                else:
                    theta = THETA_EQ_DEG + sign * D_THETA_DEG
                    disp, unit = sign * D_THETA_DEG, "deg"
                out.append({"mode": mode, "sign": sign, "reference": ref_name,
                            "r0_a": r0, "r1_a": r1, "r2_a": r2, "theta_deg": theta,
                            "displacement": disp, "displacement_unit": unit})
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out-dir", type=Path, default=Path(__file__).resolve().parent)
    args = ap.parse_args(argv)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    entries = build()
    assert len(entries) == 13, f"expected 13 geometries, built {len(entries)}"
    index = []
    for i, e in enumerate(entries):
        name = f"g{i:02d}"
        atoms = geometry(e["r1_a"], e["r2_a"], e["theta_deg"])
        lines = [f"{len(atoms)}",
                 f"{name} {e['reference']} {e['mode']} "
                 f"{e['displacement']:+g}{e['displacement_unit']} "
                 f"r1={e['r1_a']:.6f} r2={e['r2_a']:.6f} theta={e['theta_deg']:.4f} "
                 f"box={BOX_A:.4f}"]
        lines += [f"{s} {x:.10f} {y:.10f} {z:.10f}" for s, x, y, z in atoms]
        (args.out_dir / f"{name}.xyz").write_text("\n".join(lines) + "\n")
        index.append({"name": name, "file": f"{name}.xyz", "box_a": BOX_A,
                      "n_atoms": len(atoms), **e})
        print(f"  {name}  {e['reference']:<11s} {e['mode']:<13s} "
              f"{e['displacement']:+6g} {e['displacement_unit']:<3s} "
              f"r1={e['r1_a']:.4f} r2={e['r2_a']:.4f} theta={e['theta_deg']:.2f}")

    meta = {
        "meta": (
            "13 water-monomer geometries for the isolated-molecule dipole test "
            "(A64). Central-difference pairs for symmetric stretch, asymmetric "
            "stretch and bend, about two references: equilibrium (0.9572 A) and a "
            "stretched bond (1.00 A), the length the bonds sample at 1000 K. "
            "Orientation is identical across geometries, so dipole vectors compare "
            "directly."),
        "box_a": BOX_A, "r_eq_a": R_EQ_A, "theta_eq_deg": THETA_EQ_DEG,
        "r_stretched_a": R_STRETCHED_A, "d_r_a": D_R_A, "d_theta_deg": D_THETA_DEG,
        "geometries": index,
    }
    (args.out_dir / "geoms.json").write_text(json.dumps(meta, indent=1) + "\n")
    print(f"\nwrote {len(index)} xyz files and {args.out_dir / 'geoms.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
