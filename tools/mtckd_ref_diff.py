#!/usr/bin/env python3
"""Task 24 / A63: diff the ARC MT_CKD reference against the one built in Task 19.

    python tools/mtckd_ref_diff.py config/fiducial.yaml

A58 could not reach ARC (``scp`` needs an interactive password), so the
self-continuum was built from AER's own release 4.3 instead — the upstream
primary source rather than a reimplementation of a reimplementation. That left
exactly one thing open: whether ARC's ``mtckd_ref.npz`` is the same coefficients.
This answers it.

Three levels of comparison, because agreement at one does not imply the others:

1. **keys and shapes** — do the two files describe the same object at all;
2. **the reference arrays** — wavenumber grid, ``self_absco_ref``, ``self_texp``,
   and the two reference-state scalars, compared element by element;
3. **the derived coefficients** C(ν,T)·R(ν,T) at 296, 650, 1000 and 2000 K over
   the science band — what the opacity actually depends on, so a difference here
   is a difference in κ and a difference in level 2 that cancels here does not
   matter.

Writes ``outputs/mtckd_ref_diff.json``.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import atomicio

import numpy as np

from config import load_config
from opacity.mtckd import MTCKDReference, load_reference, self_continuum_coefficient

BAND_CM = (100.0, 4000.0)
TEMPERATURES = (296.0, 650.0, 1000.0, 2000.0)
ARC_RELPATH = "mtckd/mtckd_ref.npz"

# ARC's key names differ from AER's netCDF; this is the mapping, not a conversion.
KEY_MAP = {
    "wavenumbers": "wavenumber_cm",
    "self_absco_ref": "self_absco_ref",
    "self_texp": "self_texp",
    "ref_press_mbar": "ref_press_mb",
    "ref_temp_K": "ref_temp_k",
}


def _as_reference(path: Path) -> tuple[MTCKDReference, list[str]]:
    with np.load(path, allow_pickle=False) as d:
        keys = list(d.files)
        return MTCKDReference(
            wavenumber_cm=np.asarray(d["wavenumbers"], float),
            self_absco_ref=np.asarray(d["self_absco_ref"], float),
            self_texp=np.asarray(d["self_texp"], float),
            ref_press_mb=float(d["ref_press_mbar"]),
            ref_temp_k=float(d["ref_temp_K"]),
            title="(ARC mtckd_ref.npz — carries no version string)",
            source_sha256="",
        ), keys


def _rel(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Relative difference where b is non-zero, absolute difference where it is."""
    out = np.zeros_like(b, dtype=float)
    nz = b != 0.0
    out[nz] = a[nz] / b[nz] - 1.0
    out[~nz] = a[~nz] - b[~nz]
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("config", type=Path)
    args = ap.parse_args(argv)
    cfg = load_config(args.config)

    arc_path = cfg.data_root_path / ARC_RELPATH
    if not arc_path.exists():
        print(f"\n  TASK 24 — SKIPPED. {arc_path} is not present.")
        print(f"  Copy it first:\n"
              f"    scp shug7909@htc-login.arc.ox.ac.uk:'~/scw/mtckd_ref.npz' "
              f"{arc_path.parent}/")
        return 0

    arc, arc_keys = _as_reference(arc_path)
    ours = load_reference(cfg.mtckd_ref_path)

    print(f"\n  TASK 24 — ARC {arc_path.name} vs {Path(cfg.mtckd_ref_path).name}")
    print(f"\n  1. KEYS AND SHAPES")
    print(f"     ARC   : {arc_keys}")
    print(f"     ours  : wavenumber_cm, self_absco_ref, self_texp, ref_press_mb, "
          f"ref_temp_k, title, version_description, source_file, source_sha256")
    print(f"     The names differ; the content maps one to one "
          f"({', '.join(f'{a}->{b}' for a, b in KEY_MAP.items())}).")
    print(f"     ARC carries no version or provenance string; ours carries AER's "
          f"title and the\n     sha256 of the netCDF it came from.")

    same_shape = arc.wavenumber_cm.shape == ours.wavenumber_cm.shape
    print(f"\n  2. REFERENCE ARRAYS AND SCALARS")
    print(f"     grid points     ARC {arc.wavenumber_cm.size}   ours "
          f"{ours.wavenumber_cm.size}   {'same' if same_shape else 'DIFFERENT'}")
    if not same_shape:
        raise SystemExit("grids differ in length; a point-by-point diff is not defined")
    arrays = {}
    for label, a, b in (
        ("wavenumber_cm", arc.wavenumber_cm, ours.wavenumber_cm),
        ("self_absco_ref", arc.self_absco_ref, ours.self_absco_ref),
        ("self_texp", arc.self_texp, ours.self_texp),
    ):
        dev = np.abs(_rel(a, b))
        arrays[label] = float(dev.max())
        print(f"     {label:<16s} max |rel diff| = {dev.max():.3e}"
              f"{'   (identical)' if dev.max() == 0.0 else ''}")
    print(f"     ref_press_mb     ARC {arc.ref_press_mb:g}   ours {ours.ref_press_mb:g}")
    print(f"     ref_temp_k       ARC {arc.ref_temp_k:g}   ours {ours.ref_temp_k:g}")

    print(f"\n  3. DERIVED COEFFICIENTS C(ν,T)·R(ν,T) ON "
          f"[{BAND_CM[0]:g}, {BAND_CM[1]:g}] cm⁻¹ — what κ depends on")
    wn = ours.wavenumber_cm
    in_band = (wn >= BAND_CM[0]) & (wn <= BAND_CM[1])
    per_t = {}
    print(f"     {'T [K]':>7s} | {'max |rel diff|':>15s} | {'median':>12s}")
    print("     " + "-" * 40)
    for t_k in TEMPERATURES:
        a = self_continuum_coefficient(wn[in_band], t_k, arc)
        b = self_continuum_coefficient(wn[in_band], t_k, ours)
        dev = np.abs(_rel(a, b))
        per_t[t_k] = float(dev.max())
        print(f"     {t_k:>7.0f} | {dev.max():>15.3e} | {np.median(dev):>12.3e}")

    worst = max(list(arrays.values()) + list(per_t.values()))
    identical = worst == 0.0
    if identical:
        verdict = (
            "IDENTICAL. The ARC file is the same MT_CKD 4.3 self-continuum "
            "reference, bit for bit on every array and scalar, so A58's model was "
            "built on exactly the coefficients ARC uses. Nothing in Tasks 19-23 "
            "changes.")
    elif worst < 1e-6:
        verdict = (
            f"THE SAME DATA to {worst:.1e} — round-trip precision, not a different "
            f"release. Nothing in Tasks 19-23 changes.")
    else:
        verdict = (
            f"DIFFERENT at the {worst:.1e} level. That is above round-trip noise, "
            f"so the two files are not the same release of the coefficients and "
            f"the difference has to be chased before the continuum numbers are "
            f"quoted.")
    print(f"\n  VERDICT: {verdict}")

    # Which release is it? self_texp is the fingerprint: it first appears in
    # MT_CKD 3.6 (the power-law temperature dependence) and the self continuum was
    # last revised in 4.0; 4.1-4.3 changed only the foreign term, which is absent
    # from both files. So the self arrays cannot distinguish 4.0 from 4.3.
    print(f"\n  WHICH VERSION: ours is AER release 4.3 ({ours.title}), sha256 "
          f"{ours.source_sha256[:16]}…")
    if identical or worst < 1e-6:
        print(f"     The ARC file matches it, so it is the same coefficients. Note "
              f"that the SELF\n     arrays alone cannot separate MT_CKD 4.0 from "
              f"4.1-4.3: the self continuum was\n     last revised in 4.0 and "
              f"4.1-4.3 changed only the FOREIGN term, which neither\n     file "
              f"carries. The identification as 4.3 rests on our file's provenance, "
              f"not on\n     the numbers — and for a pure-steam envelope the "
              f"distinction is empty anyway,\n     since the foreign term carries "
              f"(1 − x_H2O) = 0.")

    out = cfg.output_path
    out.mkdir(parents=True, exist_ok=True)
    path = out / "mtckd_ref_diff.json"
    atomicio.write_text(path, json.dumps({
        "meta": (
            "Diff of ARC's mtckd_ref.npz against the MT_CKD 4.3 self-continuum "
            "reference extracted from AER's netCDF in Task 19 (A58/A63)."),
        "arc_file": str(arc_path),
        "ours_file": str(cfg.mtckd_ref_path),
        "arc_keys": arc_keys,
        "key_map": KEY_MAP,
        "array_max_rel_diff": arrays,
        "scalars": {
            "ref_press_mb": {"arc": arc.ref_press_mb, "ours": ours.ref_press_mb},
            "ref_temp_k": {"arc": arc.ref_temp_k, "ours": ours.ref_temp_k},
        },
        "coefficient_max_rel_diff_by_temperature_k": per_t,
        "band_cm": list(BAND_CM),
        "worst_max_rel_diff": worst,
        "identical": identical,
        "verdict": verdict,
    }, indent=1) + "\n")
    print(f"\n  wrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
