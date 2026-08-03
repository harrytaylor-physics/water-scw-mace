"""Opacity source registry (SPEC §3).

Each source is described by a :class:`SourceSpec` holding its *intrinsic* data
properties: where its per-point spectra live (relative to ``data_root``), the
filename glob, and any on-load opacity scale factor. Adding a new source (e.g.
POKAZATEL once delivered) is a single entry here — no code changes elsewhere.

Locations for HITEMP/POKAZATEL are placeholders until their spectra are
delivered (AMENDMENTS A5); building them before then raises a clear error.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class SourceSpec:
    """Intrinsic properties of an opacity source's on-disk spectra."""

    name: str
    subdir: str  # spectra directory, relative to data_root
    glob: str  # filename pattern within subdir
    scale: float  # multiply kappa/alpha by this on load
    scale_reason: str
    # Below the grid's ρ floor, hold κ flat in ρ (dilute limit): pressure
    # broadening is negligible so mass opacity becomes ρ-independent (A26). True
    # for line lists (physically justified in the dilute regime); False for FP,
    # a condensed-phase method whose low-ρ extrapolation is not justified.
    dilute_flat: bool = False


REGISTRY: dict[str, SourceSpec] = {
    # FP carries the legacy factor-1000 error -> scale 1e-3 (AMENDMENTS A4).
    "fp": SourceSpec(
        name="fp",
        subdir="fp_opacity",
        glob="T*_opacity.npz",
        scale=1.0e-3,
        scale_reason="legacy factor-1000 error in FP files (AMENDMENTS A4)",
    ),
    # Line lists carry correct units -> no scaling. Filenames follow the
    # delivered POKAZATEL convention T{T}_rho{rho}.npz (no _opacity suffix);
    # generate_linelist_tables.py writes HITEMP the same way (AMENDMENTS A16).
    "hitemp": SourceSpec(
        name="hitemp",
        subdir="hitemp_opacity",
        glob="T*_rho*.npz",
        scale=1.0,
        scale_reason="line-list opacity already in correct units",
        dilute_flat=True,  # HITEMP delivered only for ρ≥0.1 (A26)
    ),
    "pokazatel": SourceSpec(
        name="pokazatel",
        subdir="pokazatel_opacity",
        glob="T*_rho*.npz",
        scale=1.0,
        scale_reason="line-list opacity already in correct units",
        dilute_flat=True,
    ),
}


def get_spec(name: str) -> SourceSpec:
    try:
        return REGISTRY[name]
    except KeyError:
        raise KeyError(
            f"unknown opacity source {name!r}; registered: {sorted(REGISTRY)}"
        ) from None
