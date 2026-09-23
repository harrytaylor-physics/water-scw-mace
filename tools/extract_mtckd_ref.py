#!/usr/bin/env python3
"""Extract the MT_CKD water-vapour self-continuum reference table to an npz.

    python tools/extract_mtckd_ref.py \
        --nc ~/software/MT_CKD/data/absco-ref_wv-mt-ckd.nc \
        --out ~/subneptune_data/mtckd/mtckd_v4.3_self_ref.npz

Why this exists. ``opacity/mtckd.py`` needs three arrays from AER's distributed
reference file — the self-continuum coefficients at the reference state, their
per-wavenumber temperature exponents, and the wavenumber grid — plus the two
reference-state scalars. That file is netCDF, and netCDF4 is *not* a dependency
of this project (SPEC §1: numpy/scipy/matplotlib, PyYAML for config only). So it
is converted once, here, into the npz the runtime reads, and this tool is the
committed record of how (AMENDMENTS A58).

Nothing is transformed: the arrays are copied verbatim, in AER's own units
(cm² molecule⁻¹ / (molecule cm⁻³ / molecule cm⁻³), i.e. the coefficient at the
reference density, before the radiation term). The provenance strings written
into the npz are what ``opacity/mtckd.py`` reports and what the gates pin.
"""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

import numpy as np

DEFAULT_NC = "~/software/MT_CKD/data/absco-ref_wv-mt-ckd.nc"
DEFAULT_OUT = "~/subneptune_data/mtckd/mtckd_v4.3_self_ref.npz"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--nc", default=DEFAULT_NC, help="AER absco-ref_wv-mt-ckd.nc")
    ap.add_argument("--out", default=DEFAULT_OUT, help="output npz")
    args = ap.parse_args(argv)

    from netCDF4 import Dataset  # import here: not a runtime dependency

    nc_path = Path(args.nc).expanduser()
    out_path = Path(args.out).expanduser()
    out_path.parent.mkdir(parents=True, exist_ok=True)

    digest = hashlib.sha256(nc_path.read_bytes()).hexdigest()
    with Dataset(nc_path) as d:
        wavenumber_cm = np.asarray(d["wavenumbers"][:], dtype=float)
        self_absco_ref = np.asarray(d["self_absco_ref"][:], dtype=float)
        self_texp = np.asarray(d["self_texp"][:], dtype=float)
        ref_press_mb = float(d["ref_press"][:])
        ref_temp_k = float(d["ref_temp"][:])
        title = str(d.Title).strip()
        version_description = str(d.Version_description).strip()

    np.savez_compressed(
        out_path,
        wavenumber_cm=wavenumber_cm,
        self_absco_ref=self_absco_ref,
        self_texp=self_texp,
        ref_press_mb=ref_press_mb,
        ref_temp_k=ref_temp_k,
        title=title,
        version_description=version_description,
        source_file=str(nc_path),
        source_sha256=digest,
    )
    print(f"read  {nc_path}\n      sha256 {digest}")
    print(f"      {title}")
    print(f"      {wavenumber_cm.size} nodes, "
          f"{wavenumber_cm.min():g}-{wavenumber_cm.max():g} cm^-1, "
          f"d(wn) = {wavenumber_cm[1] - wavenumber_cm[0]:g} cm^-1")
    print(f"      reference state {ref_press_mb:g} mb, {ref_temp_k:g} K; "
          f"self_texp {self_texp.min():g}-{self_texp.max():g}")
    print(f"wrote {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
