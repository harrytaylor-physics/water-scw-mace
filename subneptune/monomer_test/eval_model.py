#!/usr/bin/env python3
"""Evaluate water_dipole_v3 on the 13 monomer geometries (A64).

    python eval_model.py                    # models/water_dipole_v3.model (layout.py)
    python eval_model.py --model ~/scw/water_dipole_v3.model
    python eval_model.py --dry-run          # data pipeline only, no model needed

**Must run on a compute node.** torch does not import on the ARC login nodes.

The inference path is the standard MACE dipole one, and is deliberately the same
call the production labelling used, so a difference here is a difference in the
molecule and not in how the model was asked:

    torch.load(..., weights_only=False)      # the checkpoint is a pickled module
    mace.data.config_from_atoms(atoms)
    mace.data.AtomicData.from_config(config, z_table, cutoff=model.r_max)
    torch_geometric.dataloader.DataLoader(..., batch_size=1)
    model(batch.to_dict(), compute_force=False)["dipole"]     # e*A

The model returns e*A; everything downstream of this script is in Debye
(1 e*A = 4.80320 D), so the conversion happens once, here.

``--dry-run`` builds the AtomicData and the batch and stops before the forward
pass. It needs no model file, so the whole data pipeline can be checked on a
machine that has mace but not the checkpoint — which is where this was written.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

EA_TO_DEBYE = 4.80320

HERE = Path(__file__).resolve().parent

sys.path.insert(0, str(HERE.parent))
from layout import DIPOLE_MODEL  # noqa: E402


def load_geometries(geoms_json: Path):
    from ase import Atoms

    meta = json.loads(geoms_json.read_text())
    out = []
    for entry in meta["geometries"]:
        lines = (geoms_json.parent / entry["file"]).read_text().splitlines()
        symbols, positions = [], []
        for line in lines[2:]:
            parts = line.split()
            if len(parts) != 4:
                continue
            symbols.append(parts[0])
            positions.append([float(v) for v in parts[1:]])
        box = float(entry["box_a"])
        atoms = Atoms(symbols=symbols, positions=positions,
                      cell=[box, box, box], pbc=True)
        out.append((entry, atoms))
    return meta, out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--model", type=Path, default=DIPOLE_MODEL)
    ap.add_argument("--geoms", type=Path, default=HERE / "geoms.json")
    ap.add_argument("--out", type=Path, default=HERE / "model_dipoles.json")
    ap.add_argument("--device", default="cpu")
    ap.add_argument("--dry-run", action="store_true",
                    help="build the batches and stop; no model file needed")
    args = ap.parse_args(argv)

    meta, geoms = load_geometries(args.geoms)
    print(f"[eval_model] {len(geoms)} geometries from {args.geoms}")

    import torch
    from mace import data as mace_data
    from mace.tools import torch_geometric, utils as mace_utils

    if args.dry_run:
        z_table = mace_utils.AtomicNumberTable([1, 8])
        cutoff, model = 5.0, None
        print(f"[eval_model] DRY RUN: z_table={z_table.zs}, cutoff={cutoff} A, "
              f"no model loaded")
    else:
        model_path = args.model.expanduser()
        if not model_path.exists():
            raise SystemExit(
                f"{model_path} not found. It is in git-LFS (git lfs pull --include "
                f"'models/*'); or pass --model, or use --dry-run to check the "
                f"data pipeline only.")
        model = torch.load(model_path, map_location=args.device, weights_only=False)
        model = model.to(args.device).eval()
        z_table = mace_utils.AtomicNumberTable(
            [int(z) for z in model.atomic_numbers])
        cutoff = float(model.r_max)
        print(f"[eval_model] loaded {model_path.name}: z={z_table.zs}, "
              f"r_max={cutoff} A, device={args.device}")

    results = []
    for entry, atoms in geoms:
        config = mace_data.config_from_atoms(atoms)
        atomic = mace_data.AtomicData.from_config(config, z_table=z_table,
                                                  cutoff=cutoff)
        loader = torch_geometric.dataloader.DataLoader(
            dataset=[atomic], batch_size=1, shuffle=False, drop_last=False)
        batch = next(iter(loader))
        if args.dry_run:
            keys = sorted(batch.to_dict())
            print(f"  {entry['name']}: batch keys = {len(keys)} "
                  f"({', '.join(keys[:6])}, ...)")
            continue
        with torch.no_grad():
            out = model(batch.to(args.device).to_dict(), compute_force=False)
        mu_ea = out["dipole"].detach().cpu().numpy().reshape(-1)[:3]
        mu_d = [float(v) * EA_TO_DEBYE for v in mu_ea]
        norm = float(sum(v * v for v in mu_d) ** 0.5)
        results.append({"name": entry["name"], "mode": entry["mode"],
                        "sign": entry["sign"], "reference": entry["reference"],
                        "displacement": entry["displacement"],
                        "displacement_unit": entry["displacement_unit"],
                        "mu_debye": mu_d, "mu_norm_debye": norm})
        print(f"  {entry['name']}  {entry['mode']:<13s} "
              f"mu = ({mu_d[0]:+.5f}, {mu_d[1]:+.5f}, {mu_d[2]:+.5f}) D  "
              f"|mu| = {norm:.5f} D")

    if args.dry_run:
        print("[eval_model] dry run complete: geometries, AtomicData and batches "
              "all build. Only the forward pass is untested.")
        return 0

    args.out.write_text(json.dumps({
        "meta": ("water_dipole_v3 dipoles on the 13 monomer geometries (A64). "
                 "Model returns e*A; converted here with 1 e*A = 4.80320 D."),
        "model": str(args.model), "device": args.device,
        "ea_to_debye": EA_TO_DEBYE, "cutoff_a": cutoff,
        "results": results,
    }, indent=1) + "\n")
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
