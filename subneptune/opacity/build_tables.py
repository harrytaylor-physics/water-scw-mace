"""Build gridded opacity tables from per-point spectra (SPEC §2, §3).

For a registered source: load its spectra, compute band-limited Planck and
Rosseland means over the common band (cm²/g), assert the factor-1000
plausibility check (SPEC §2), convert to SI, and wrap in a
:class:`GridOpacitySource`. Sources are cached by ``cache_key`` = (name, band,
rel_floor, source-specific knobs, data_root).
"""

from __future__ import annotations

if __name__ == "__main__":  # allow `python opacity/build_tables.py ...` as a CLI
    import sys as _sys
    from pathlib import Path as _Path

    _sys.path.insert(0, str(_Path(__file__).resolve().parents[1]))

import csv
from pathlib import Path

import numpy as np

from config import Config
from opacity.composite import CompositeOpacitySource
from opacity.loaders import find_spectra, load_spectrum
from opacity.means import planck_mean, rosseland_mean
from opacity.registry import REGISTRY, SourceSpec, get_spec
from opacity.source import GridOpacitySource, OpacitySource
from opacity.table import CM2_PER_G_TO_SI, OpacityTable

# SPEC §2: band-integrated kappa_P must land in ~1e1-1e4 cm²/g. Hard-fail bounds
# are one decade wider each side (avoids edge false-positives while still
# catching a factor-1000 error, which lands kappa_P at ~1e5-1e6 cm²/g).
_PLAUSIBLE_LO_CM2_G = 1.0e0
_PLAUSIBLE_HI_CM2_G = 1.0e5

_cache: dict[tuple, OpacitySource] = {}
_band_cache: dict[str, tuple[float, float]] = {}


def _source_variant(name: str, cfg: Config) -> tuple:
    """Config knobs, beyond (name, band, rel_floor), that change the built source.

    Keeping these in the cache key means a scan over the composite blend
    prescription (A35) rebuilds rather than silently reusing the first source.
    """
    op = cfg.opacity
    if name == "composite":
        return (
            op.composite_dilute_source,
            op.composite_dense_source,
            tuple(op.blend_window_kg_m3),
            op.blend_mode,
            op.fp_scale,
        )
    if name == "fp":
        return (op.fp_scale,)
    return ()


def cache_key(name: str, cfg: Config, band: tuple[float, float] | None) -> tuple:
    """Cache key for a built source. Band-agnostic injected sources use band=None."""
    if band is None:
        return (name, None, str(cfg.data_root_path))
    return (name, band, cfg.opacity.rel_floor, _source_variant(name, cfg),
            str(cfg.data_root_path))


def detect_common_band(cfg: Config, sources: list[str] | None = None) -> tuple[float, float]:
    """Largest wavenumber band the given sources share (A14, A21).

    = [max over spectra files of min(wn), min over files of max(wn)], verified
    per file. ``sources`` limits the intersection to those grid sources (band-
    agnostic sources like ``constant`` and undelivered sources are ignored);
    ``None`` uses every registered source with spectra on disk. Per-comparison
    banding (A21): a figure passes exactly the grid sources it plots, so a
    windowed dataset (e.g. POKAZATEL [50,4500]) only constrains comparisons it
    is part of.
    """
    names = sources if sources is not None else list(REGISTRY)
    key = (str(cfg.data_root_path), tuple(sorted(names)))
    if key in _band_cache:
        return _band_cache[key]
    mins: list[float] = []
    maxs: list[float] = []
    for name in names:
        if name not in REGISTRY:  # e.g. "constant" — band-agnostic, skip
            continue
        spec = REGISTRY[name]
        for path in sorted((cfg.data_root_path / spec.subdir).glob(spec.glob)):
            with np.load(path) as d:
                wn = np.asarray(d["wn"], dtype=float)
            mins.append(float(wn.min()))
            maxs.append(float(wn.max()))
    if not mins:
        raise FileNotFoundError(
            f"no opacity spectra found to detect a common band for sources={names}"
        )
    band = (max(mins), min(maxs))
    _band_cache[key] = band
    return band


def resolve_band(cfg: Config, sources: list[str] | None = None) -> tuple[float, float]:
    """Config band if pinned, else the auto-detected band for ``sources`` (A21).

    ``sources`` defaults to all registered sources with spectra. For a single
    source this yields its own full band; for a comparison, pass the plotted
    sources so every curve shares one exact band.
    """
    if cfg.opacity.band_cm is not None:
        return tuple(cfg.opacity.band_cm)  # type: ignore[return-value]
    return detect_common_band(cfg, sources)


def compute_point_means(
    name: str, cfg: Config, band_cm: tuple[float, float], rel_floor: float | None = None
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Per-point (T, ρ, κ_R, κ_P) means in cm²/g over ``band_cm``.

    ``rel_floor`` (default from config) applies the transparency-window floor
    before the Rosseland mean (A25). Reusable by the table builder and by the
    band/convention studies (docs/band_limitation.md, docs/opacity_crosscheck.md).
    """
    if rel_floor is None:
        rel_floor = cfg.opacity.rel_floor
    spec = get_spec(name)
    spectra_dir = cfg.data_root_path / spec.subdir
    paths = find_spectra(spectra_dir, spec)

    temps, rhos, kappa_r, kappa_p = [], [], [], []
    for path in paths:
        rec = load_spectrum(path, spec, band_cm, mtckd_ref_path=cfg.mtckd_ref_path)
        kappa_p.append(planck_mean(rec.wn_cm, rec.kappa_cm2_g, rec.temperature_k, band_cm, rel_floor))
        kappa_r.append(rosseland_mean(rec.wn_cm, rec.kappa_cm2_g, rec.temperature_k, band_cm, rel_floor))
        temps.append(rec.temperature_k)
        rhos.append(rec.rho_si)
    return (
        np.array(temps),
        np.array(rhos),
        np.array(kappa_r),
        np.array(kappa_p),
    )


def _assert_plausible(name: str, spec: SourceSpec, kappa_p_cm2_g: np.ndarray) -> None:
    """Hard-fail if band-integrated Planck means are implausible (SPEC §2)."""
    lo, hi = float(kappa_p_cm2_g.min()), float(kappa_p_cm2_g.max())
    if lo < _PLAUSIBLE_LO_CM2_G or hi > _PLAUSIBLE_HI_CM2_G:
        raise AssertionError(
            f"opacity source {name!r}: band-integrated kappa_P in "
            f"[{lo:.3g}, {hi:.3g}] cm^2/g is outside the plausible window "
            f"[{_PLAUSIBLE_LO_CM2_G:g}, {_PLAUSIBLE_HI_CM2_G:g}] cm^2/g. "
            f"Applied scale={spec.scale:g} ({spec.scale_reason}). "
            f"A ~1e3 offset indicates the factor-1000 scaling is wrong."
        )


def _line_only_name(name: str) -> str:
    """The bare line-list source a continuum-inclusive source is built on.

    Both read the same spectra directory; the continuum source differs only by
    its ``continuum`` field, so the line-only twin is the one registered entry
    with the same ``subdir`` and no continuum.
    """
    spec = get_spec(name)
    for other, cand in REGISTRY.items():
        if cand.subdir == spec.subdir and cand.continuum is None:
            return other
    raise KeyError(
        f"source {name!r} declares a continuum but no bare line-list twin is "
        f"registered for subdir {spec.subdir!r}"
    )


def build_composite(cfg: Config, band: tuple[float, float] | None = None) -> CompositeOpacitySource:
    """Build the composite source: line-list (dilute) blended with FP (dense).

    The two sub-sources share the band common to *them* (not to POKAZATEL or any
    other on-disk dataset), so the composite is not capped by an unrelated
    windowed source (A21).
    """
    op = cfg.opacity
    names = [op.composite_dilute_source, op.composite_dense_source]
    band = band or resolve_band(cfg, names)
    dilute = build_source(op.composite_dilute_source, cfg, band=band)
    dense = build_source(op.composite_dense_source, cfg, band=band)
    return CompositeOpacitySource(
        dilute=dilute, dense=dense, blend_rho_si=tuple(op.blend_window_kg_m3),
        mode=op.blend_mode,
    )


def _default_band(name: str, cfg: Config) -> tuple[float, float]:
    """A source's own band; raise the delivery-specific error if it has no spectra."""
    names = (
        [cfg.opacity.composite_dilute_source, cfg.opacity.composite_dense_source]
        if name == "composite"
        else [name]
    )
    try:
        return resolve_band(cfg, names)
    except FileNotFoundError:
        # Surface the clearer "not yet delivered" message per missing source.
        for n in names:
            find_spectra(cfg.data_root_path / get_spec(n).subdir, get_spec(n))
        raise  # pragma: no cover - find_spectra already raised


def build_source(name: str, cfg: Config, band: tuple[float, float] | None = None) -> OpacitySource:
    """Build a fresh opacity source for ``name`` (no caching).

    ``band`` fixes the Rosseland/Planck band; when ``None`` the source uses its
    *own* band (A21), so a pure-source run is independent of what else is on
    disk. Comparison contexts pass the shared band of the plotted sources.
    """
    if name == "composite":
        return build_composite(cfg, band=band)
    band_cm = band or _default_band(name, cfg)
    spec = get_spec(name)
    temps, rhos, kappa_r, kappa_p = compute_point_means(name, cfg, band_cm)
    # Plausibility is asserted on the *unperturbed* means: it guards the
    # factor-1000 data error (SPEC §2), not the deliberate fp_scale perturbation.
    # A continuum-inclusive source (A58) adds a physically-computed term that
    # legitimately carries kappa_P above the ceiling at high density, so the
    # guard is applied to its LINE-ONLY part — which is the delivered data the
    # guard exists to protect. The continuum itself is gated separately, against
    # AER's own binary (validation/test_mtckd_continuum.py).
    try:
        _assert_plausible(name, spec, kappa_p)
    except AssertionError:
        if spec.continuum is None:
            raise
        # Re-assert on the line-only twin: the delivered spectra are what the
        # guard protects, and an added continuum may legitimately carry the
        # total past the ceiling at high density.
        _assert_plausible(_line_only_name(name), spec,
                          compute_point_means(_line_only_name(name), cfg, band_cm)[3])
    # fp_scale (A37): a constant perturbation of the FP table, applied to both
    # means here so it is already in place when the composite blends (or when FP
    # is used pure). No effect on any other source.
    scale = cfg.opacity.fp_scale if name == "fp" else 1.0

    table = OpacityTable(
        temperature_k=temps,
        rho_si=rhos,
        kappa_rosseland_si=kappa_r * CM2_PER_G_TO_SI * scale,
        kappa_planck_si=kappa_p * CM2_PER_G_TO_SI * scale,
        band_cm=band_cm,
        dilute_flat=spec.dilute_flat,
    )
    return GridOpacitySource(name, table)


def get_source(name: str, cfg: Config, band: tuple[float, float] | None = None) -> OpacitySource:
    """Build or return a cached opacity source for ``name`` on ``band``.

    ``band`` None -> the source's own band (A21). Cached sources share a clamp
    counter; reset it before a structure integration via
    ``source.reset_clamp_counter()``.
    """
    # Band-agnostic injected sources (e.g. a constant-κ reference) are keyed with
    # band=None; they are placed in the cache directly and only looked up here.
    if name != "composite" and name not in REGISTRY:
        key = cache_key(name, cfg, None)
        if key in _cache:
            return _cache[key]
        raise KeyError(f"opacity source {name!r} is not registered or cached")
    band_cm = band or _default_band(name, cfg)
    # rel_floor changes the means; the variant carries the composite blend knobs.
    key = cache_key(name, cfg, band_cm)
    if key not in _cache:
        _cache[key] = build_source(name, cfg, band=band_cm)
    return _cache[key]


def clear_cache() -> None:
    _cache.clear()
    _band_cache.clear()


# --------------------------------------------------------------------------
# Band studies + CLI (A39): full-band vs common-band means for a line list
# --------------------------------------------------------------------------
def resolve_band_spec(cfg: Config, name: str, spec: str) -> tuple[float, float]:
    """Resolve a band spec for one source: ``common`` | ``full`` | ``lo,hi``.

    ``common`` is the band the science runs use (config ``opacity.band_cm`` when
    pinned, else the auto-detected common band). ``full`` is the source's own
    computed wavenumber range, straight from its spectra files — never wider than
    the data (SPEC §3's band-truncation quantification; A39).
    """
    if spec == "common":
        return resolve_band(cfg, [name])
    if spec == "full":
        return detect_common_band(cfg, [name])
    lo, hi = (float(x) for x in spec.split(","))
    return (lo, hi)


def write_point_means(
    name: str, cfg: Config, band_cm: tuple[float, float], tag: str,
    out_dir: Path | None = None, rel_floor: float | None = None,
) -> tuple[Path, Path]:
    """Compute and write per-grid-point means for one source/band (csv + npz)."""
    out_dir = Path(out_dir) if out_dir is not None else cfg.output_path
    out_dir.mkdir(parents=True, exist_ok=True)
    floor = cfg.opacity.rel_floor if rel_floor is None else rel_floor
    t, rho, kr, kp = compute_point_means(name, cfg, band_cm, rel_floor=floor)
    order = np.lexsort((rho, t))
    t, rho, kr, kp = t[order], rho[order], kr[order], kp[order]

    csv_path = out_dir / f"{name}_means_{tag}.csv"
    with csv_path.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["temperature_k", "rho_kg_m3", "kappa_R_cm2_g", "kappa_P_cm2_g",
                    "kappa_R_si", "kappa_P_si", "band_lo_cm", "band_hi_cm", "rel_floor"])
        for i in range(t.size):
            w.writerow([f"{t[i]:g}", f"{rho[i]:g}", f"{kr[i]:.6e}", f"{kp[i]:.6e}",
                        f"{kr[i] * CM2_PER_G_TO_SI:.6e}", f"{kp[i] * CM2_PER_G_TO_SI:.6e}",
                        f"{band_cm[0]:.12g}", f"{band_cm[1]:.12g}", f"{floor:g}"])
    npz_path = out_dir / f"{name}_means_{tag}.npz"
    np.savez_compressed(npz_path, temperature_k=t, rho_si=rho, kappa_R_cm2_g=kr,
                        kappa_P_cm2_g=kp, band_cm=np.array(band_cm), rel_floor=floor)
    return csv_path, npz_path


def _cli(argv=None) -> None:
    import argparse

    ap = argparse.ArgumentParser(
        description="Build gridded opacity means for one source on a chosen band."
    )
    ap.add_argument("config", type=Path, help="run config YAML")
    ap.add_argument("--source", default="hitemp", help="registered source name")
    ap.add_argument("--band", default="common",
                    help="common (science band) | full (the source's own computed "
                         "range) | lo,hi in cm^-1")
    ap.add_argument("--rel-floor", type=float, default=None,
                    help="override opacity.rel_floor for these means")
    ap.add_argument("--out-dir", type=Path, default=None)
    args = ap.parse_args(argv)

    from config import load_config

    cfg = load_config(args.config)
    band = resolve_band_spec(cfg, args.source, args.band)
    tag = args.band if args.band in ("common", "full") else f"{band[0]:g}-{band[1]:g}"
    csv_path, npz_path = write_point_means(args.source, cfg, band, tag,
                                           out_dir=args.out_dir, rel_floor=args.rel_floor)
    with np.load(npz_path) as d:
        kr, kp, t = d["kappa_R_cm2_g"], d["kappa_P_cm2_g"], d["temperature_k"]
    print(f"source={args.source}  band={args.band} = [{band[0]:.12g}, {band[1]:.12g}] cm^-1  "
          f"rel_floor={cfg.opacity.rel_floor if args.rel_floor is None else args.rel_floor:g}")
    print(f"{t.size} grid points; kappa_R {kr.min():.4g}-{kr.max():.4g} cm^2/g, "
          f"kappa_P {kp.min():.4g}-{kp.max():.4g} cm^2/g")
    print(f"wrote {csv_path}")
    print(f"wrote {npz_path}")


if __name__ == "__main__":
    _cli()
