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

from config.schema import FP_OPACITY_SUBDIR

# A87: which marker values belong with which first-principles directory. The
# registry's subdir and the files on disk have to agree, or spectra from one
# version load under another version's expectation -- which is how a "v3" ledger
# column came to carry v5 numbers. v3 carries two values because A74's four
# corner cells were stamped after A50's rewrite.
_EXPECTED_MARKERS: dict[str, tuple[str, ...]] = {
    "fp_opacity": ("cm2_per_g:A50", "cm2_per_g:A74"),
    "fp_opacity_v4": ("cm2_per_g:A80",),
    "fp_opacity_v5": ("cm2_per_g:A83",),
}


@dataclass(frozen=True)
class SourceSpec:
    """Intrinsic properties of an opacity source's on-disk spectra."""

    name: str
    subdir: str  # spectra directory, relative to data_root
    glob: str  # filename pattern within subdir
    scale: float  # multiply kappa/alpha by this on load
    scale_reason: str
    # npz key that must be present or the spectrum is REFUSED (A50). FP files
    # were rewritten in correct units by tools/fix_fp_units.py, which stamps this
    # marker; an un-stamped file is a legacy ×1000 one and must never load.
    units_marker: str | None = None
    # A87: the marker's expected *value*, not just its presence. The subdir and
    # the marker have to agree, or a file from one version is being read while
    # the registry points at another -- which is how a "v3" ledger column came to
    # carry v5 numbers. None = accept any value (line lists, which are versionless).
    units_value: tuple[str, ...] | None = None
    # Additive continuum model applied to every loaded spectrum, by name in
    # ``opacity.mtckd.CONTINUUM_MODELS``. None = bare line absorption. The
    # continuum is a function of the spectrum's own (T, ρ), so it is applied at
    # load time, before any band-limited mean is taken (A58).
    continuum: str | None = None
    # Below the grid's ρ floor, hold κ flat in ρ (dilute limit): pressure
    # broadening is negligible so mass opacity becomes ρ-independent (A26). True
    # for line lists (physically justified in the dilute regime); False for FP,
    # a condensed-phase method whose low-ρ extrapolation is not justified.
    dilute_flat: bool = False


REGISTRY: dict[str, SourceSpec] = {
    # FP units are now correct IN THE DATA (A50): tools/fix_fp_units.py divided
    # the 27 legacy spectra by 1000 at source, archived the originals in
    # fp_opacity_legacy_x1000/, and merged the 21 extension points (already in
    # correct units). No on-load scale survives; the marker is what enforces it.
    "fp": SourceSpec(
        name="fp",
        subdir=FP_OPACITY_SUBDIR,
        glob="T*_opacity.npz",
        scale=1.0,
        scale_reason="units fixed at source by tools/fix_fp_units.py (A50); spectra regenerated with the A81 estimator at sd = 10 cm^-1 (A83)",
        units_marker="units_fixed",
        units_value=_EXPECTED_MARKERS.get(FP_OPACITY_SUBDIR),

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
    # HITEMP line absorption PLUS the MT_CKD 4.3 water self-continuum evaluated
    # at each spectrum's own (T, ρ) — the second reference Fig. 2b reports, and
    # the one every planetary model actually uses (A58). Same spectra on disk as
    # ``hitemp``; the continuum is added on load. The continuum is extrapolated
    # far outside MT_CKD's fitted range in BOTH T and ρ — see opacity/mtckd.py.
    "hitemp_mtckd": SourceSpec(
        name="hitemp_mtckd",
        subdir="hitemp_opacity",
        glob="T*_rho*.npz",
        scale=1.0,
        scale_reason="line-list opacity already in correct units",
        continuum="mtckd_self",
        dilute_flat=True,
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
