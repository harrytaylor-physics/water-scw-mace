"""One-shot: fix the FP opacity units at source and merge the extension grid (A50).

Until now the 27 original ``fp_opacity/`` spectra carried a legacy factor-1000
error and the loader divided by 1000 on every load (SPEC §2, A4). The extension
campaign (21 new points) is delivered in correct units, so a single on-load
scale can no longer serve both. This script fixes the units *in the data*:

1. ``fp_opacity/`` is copied to ``fp_opacity_legacy_x1000/`` — an untouched
   archive of exactly what was delivered.
2. Each of the 27 originals is rewritten with ``alpha``, ``kappa``, ``planck``
   and ``rosseland`` divided by 1000.
3. The 21 extension spectra are copied in unchanged.
4. Every file in ``fp_opacity/`` gains a :data:`UNITS_MARKER` key, and
   ``grid.json`` is regenerated over all 48 points.

**Which arrays are scaled, and why not ``sigma`` (A50).** The factor-1000 entered
at the ``alpha`` step, not at the cross-sections: in all 27 originals
``alpha = 1000 · n · sigma`` with ``n = rho·N_A/M`` the number density, while in
all 21 extension files ``alpha = n · sigma`` exactly. ``kappa = alpha/rho`` holds
in both. So ``sigma``/``sigma_raw`` are already correct cm²/molecule and are left
alone — scaling them would break ``alpha = n·sigma`` in the rewritten files and
corrupt the cross-sections the dilute noise-floor check reads (Task 15).
:func:`verify_alpha_convention` asserts this before anything is written.

The stored ``planck``/``rosseland`` scalars are FULL-BAND means and are scaled
only so the file is self-consistent; no science path reads them (band-limited
means are recomputed from ``kappa(wn)`` by ``opacity.build_tables``, A4).

Run once::

    python tools/fix_fp_units.py            # --dry-run to preview
"""

from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path

import numpy as np

# Marker key written into every npz in fp_opacity/. The loader REFUSES any FP
# spectrum without it, so a legacy x1000 file can never be silently re-ingested.
UNITS_MARKER = "units_fixed"
MARKER_VALUE = "cm2_per_g:A50"

LEGACY_SCALE = 1.0e-3  # applied to the 27 originals
SCALED_KEYS = ("alpha", "kappa", "planck", "rosseland")  # NOT sigma/sigma_raw
UNSCALED_KEYS = ("wn", "sigma", "sigma_raw", "T", "rho_kg_m3", "N", "dt_fs", "smooth_cm")

DATA_ROOT = Path("~/subneptune_data").expanduser()
FP_DIR = DATA_ROOT / "fp_opacity"
EXT_DIR = DATA_ROOT / "fp_ext"
ARCHIVE_DIR = DATA_ROOT / "fp_opacity_legacy_x1000"

_NA = 6.02214076e23  # mol^-1
_M_H2O = 18.015286  # g/mol


def number_density_cm3(rho_kg_m3: float) -> float:
    """Molecules per cm³ at ``rho_kg_m3`` (1 kg/m³ = 1e-3 g/cm³)."""
    return rho_kg_m3 * 1.0e-3 * _NA / _M_H2O


def alpha_over_n_sigma(d) -> float:
    """Median ``alpha / (n·sigma)``: 1000 for a legacy file, 1 for a correct one."""
    sigma = np.asarray(d["sigma"], float)
    alpha = np.asarray(d["alpha"], float)
    n = number_density_cm3(float(d["rho_kg_m3"]))
    m = sigma > 0
    return float(np.median(alpha[m] / sigma[m]) / n)


def verify_alpha_convention(paths: list[Path], expected: float, label: str) -> None:
    """Hard-fail unless every file follows ``alpha = expected · n · sigma`` (A50).

    This is the evidence the scaling decision rests on, so it is checked on the
    real files immediately before they are rewritten, not assumed.
    """
    for path in paths:
        with np.load(path) as d:
            got = alpha_over_n_sigma(d)
            kappa_check = float(
                np.median(np.asarray(d["kappa"], float)[np.asarray(d["sigma"], float) > 0]
                          / np.asarray(d["alpha"], float)[np.asarray(d["sigma"], float) > 0])
                * (float(d["rho_kg_m3"]) * 1.0e-3)
            )
        if not np.isclose(got, expected, rtol=1e-6):
            raise AssertionError(
                f"{path.name}: alpha/(n*sigma) = {got:.8g}, expected {expected:g} for "
                f"{label} files. The unit convention is not what A50 established; "
                f"stop and re-examine before rewriting any data."
            )
        if not np.isclose(kappa_check, 1.0, rtol=1e-6):
            raise AssertionError(
                f"{path.name}: kappa·rho/alpha = {kappa_check:.8g}, expected 1 "
                f"(kappa = alpha/rho_g_cm3 must hold in every file)"
            )


def _payload(d, scale: float) -> dict:
    """The rewritten file's arrays: ``SCALED_KEYS`` × ``scale``, the rest verbatim."""
    out = {k: np.asarray(d[k]) for k in UNSCALED_KEYS if k in d.files}
    for k in SCALED_KEYS:
        if k in d.files:
            out[k] = np.asarray(d[k], float) * scale
    for k in d.files:  # anything unanticipated is carried through untouched
        if k not in out:
            out[k] = np.asarray(d[k])
    out[UNITS_MARKER] = np.array(MARKER_VALUE)
    return out


def _write_npz(path: Path, payload: dict) -> None:
    """Write ``payload`` to ``path`` atomically (tmp file, then replace)."""
    tmp = path.with_suffix(".npz.tmp")
    with tmp.open("wb") as fh:
        np.savez(fh, **payload)
    tmp.replace(path)


def rewrite_legacy(path: Path, dry_run: bool = False) -> dict:
    """Divide one legacy file's α/κ/planck/rosseland by 1000 and mark it."""
    with np.load(path) as d:
        if UNITS_MARKER in d.files:
            raise RuntimeError(
                f"{path.name} already carries the {UNITS_MARKER!r} marker — refusing "
                f"to scale it a second time. Restore from {ARCHIVE_DIR} and re-run."
            )
        before = float(d["planck"])
        payload = _payload(d, LEGACY_SCALE)
    if not dry_run:
        _write_npz(path, payload)
    return {"name": path.name, "planck_before": before,
            "planck_after": before * LEGACY_SCALE}


def import_extension(src: Path, dst_dir: Path, dry_run: bool = False) -> dict:
    """Copy one extension file in unchanged, adding only the marker."""
    with np.load(src) as d:
        payload = _payload(d, 1.0)  # correct units already
        planck = float(d["planck"])
    if not dry_run:
        _write_npz(dst_dir / src.name, payload)
    return {"name": src.name, "planck_before": planck, "planck_after": planck}


def write_grid_json(fp_dir: Path, dry_run: bool = False) -> list[dict]:
    """Regenerate ``grid.json`` over every spectrum now in ``fp_dir``."""
    grid = []
    for path in sorted(fp_dir.glob("T*_opacity.npz")):
        with np.load(path) as d:
            t, rho = float(d["T"]), float(d["rho_kg_m3"])
        grid.append({"T": int(t) if t.is_integer() else t,
                     "rho": int(rho) if rho.is_integer() else rho})
    grid.sort(key=lambda e: (e["T"], e["rho"]))
    if not dry_run:
        (fp_dir / "grid.json").write_text(json.dumps(grid, indent=1) + "\n")
    return grid


def audit(fp_dir: Path) -> list[dict]:
    """Gate (c): every file loads, κ finite and positive; report wn spans."""
    rows = []
    for path in sorted(fp_dir.glob("T*_opacity.npz")):
        with np.load(path) as d:
            kappa = np.asarray(d["kappa"], float)
            wn = np.asarray(d["wn"], float)
            rows.append({
                "name": path.name, "T": float(d["T"]), "rho": float(d["rho_kg_m3"]),
                "wn_lo": float(wn.min()), "wn_hi": float(wn.max()),
                "kappa_lo": float(kappa.min()), "kappa_hi": float(kappa.max()),
                "n_nan": int((~np.isfinite(kappa)).sum()),
                "n_nonpos": int((kappa <= 0).sum()),
                "marked": UNITS_MARKER in d.files,
            })
    return rows


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--dry-run", action="store_true",
                    help="verify and report, write nothing")
    args = ap.parse_args(argv)

    legacy = sorted(FP_DIR.glob("T*_opacity.npz"))
    extension = sorted(EXT_DIR.glob("T*_opacity.npz"))
    print(f"fp_opacity/: {len(legacy)} legacy spectra")
    print(f"fp_ext/:     {len(extension)} extension spectra")
    if not legacy or not extension:
        print("nothing to do (one of the directories is empty)")
        return 1

    # --- the evidence the scaling rests on, checked on the real files (A50) ---
    print("\nverifying the alpha = k·n·sigma convention before touching anything ...")
    verify_alpha_convention(legacy, 1000.0, "legacy")
    verify_alpha_convention(extension, 1.0, "extension")
    print(f"  OK: all {len(legacy)} legacy files have alpha = 1000·n·sigma; "
          f"all {len(extension)} extension files have alpha = n·sigma")

    # --- 1. archive ---
    if ARCHIVE_DIR.exists():
        print(f"\narchive {ARCHIVE_DIR} already exists — leaving it untouched")
    elif args.dry_run:
        print(f"\n[dry-run] would copy {FP_DIR} -> {ARCHIVE_DIR}")
    else:
        shutil.copytree(FP_DIR, ARCHIVE_DIR)
        print(f"\narchived {FP_DIR} -> {ARCHIVE_DIR} "
              f"({len(list(ARCHIVE_DIR.glob('*')))} entries)")

    # --- 2. rewrite the 27 ---
    print(f"\nrewriting {len(legacy)} legacy spectra (x{LEGACY_SCALE:g} on "
          f"{', '.join(SCALED_KEYS)}; sigma/sigma_raw unchanged):")
    for path in legacy:
        r = rewrite_legacy(path, dry_run=args.dry_run)
        print(f"  {r['name']:<34s} planck {r['planck_before']:>12.4f} -> "
              f"{r['planck_after']:>10.4f} cm^2/g")

    # --- 3. import the 21 ---
    print(f"\nimporting {len(extension)} extension spectra unchanged:")
    for path in extension:
        r = import_extension(path, FP_DIR, dry_run=args.dry_run)
        print(f"  {r['name']:<34s} planck {r['planck_after']:>10.4f} cm^2/g")

    # --- 4. grid.json ---
    grid = write_grid_json(FP_DIR, dry_run=args.dry_run)
    t_vals = sorted({e["T"] for e in grid})
    r_vals = [e["rho"] for e in grid]
    print(f"\ngrid.json: {len(grid)} points; T = {t_vals}; "
          f"rho = {min(r_vals):g} - {max(r_vals):g} kg/m^3")

    if args.dry_run:
        print("\n[dry-run] no files written")
        return 0

    # --- gate (c) audit ---
    rows = audit(FP_DIR)
    print(f"\naudit of {len(rows)} files in {FP_DIR}:")
    print(f"  {'file':<34s} {'T':>6s} {'rho':>9s} {'wn lo':>8s} {'wn hi':>9s} "
          f"{'kappa lo':>11s} {'kappa hi':>11s} {'mark':>5s}")
    bad = 0
    for r in rows:
        bad += r["n_nan"] + r["n_nonpos"] + (0 if r["marked"] else 1)
        print(f"  {r['name']:<34s} {r['T']:>6.0f} {r['rho']:>9.4g} {r['wn_lo']:>8.4f} "
              f"{r['wn_hi']:>9.3f} {r['kappa_lo']:>11.5g} {r['kappa_hi']:>11.5g} "
              f"{'yes' if r['marked'] else 'NO':>5s}")
    wn_lo = max(r["wn_lo"] for r in rows)
    wn_hi = min(r["wn_hi"] for r in rows)
    print(f"\n  common wn span across all {len(rows)} files: "
          f"[{wn_lo:.6f}, {wn_hi:.6f}] cm^-1")
    print(f"  files with NaN/non-positive kappa or no marker: {bad}")
    return 0 if bad == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
