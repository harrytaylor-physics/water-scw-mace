"""Where things live now that subneptune sits inside water-scw-mace (migration).

Every location outside this package is named here once, so a producer and the
gate that checks it cannot drift apart -- A87 found three producer groups writing
to hard-coded repository paths, and a gate reading one while a producer wrote
another. Each output directory keeps the environment override it already had.

    water-scw-mace/
      subneptune/                  ROOT  -- this package
      production_campaign/scw/     the ARC production pipeline, kept verbatim
      production_campaign/arc_scripts/  its Slurm scripts and the v4 training log
      data/fp_opacity_v5/          the 52 first-principles spectra (A83)
      data/release/                MESA-format kappa tables (A73)
      docs/figures/subneptune/     Figs 1-4 and Extended Data
      models/                      the MACE models and their error table
      benchmarks/reference/        reference.json (the 41/58 parity)
      docs/                        METHODS.md, LDA_CAMPAIGN.md

Raw inputs (the AQUA tables, MT_CKD, the line-list spectra) are not in the
repository (SPEC §1) and stay under ``paths.data_root`` in the config.
"""

from __future__ import annotations

import os
from pathlib import Path

ROOT = Path(__file__).resolve().parent
HOST = ROOT.parent


def _env(name: str, default: Path) -> Path:
    return Path(os.environ.get(name, default)).expanduser()


FIGS_DIR = _env("SUBNEPTUNE_FIGS_DIR", HOST / "docs" / "figures" / "subneptune")
RELEASE_DIR = _env("SUBNEPTUNE_RELEASE_DIR", HOST / "data" / "release")


# Spectra directories that are committed to this repository. Only these resolve
# here; every other registry source (line lists, the v3/v4 FP spectra the
# comparison tools read) stays under ``data_root``. A named list, not an
# existence check, so a missing copy fails instead of silently reading another.
COMMITTED_SPECTRA = ("fp_opacity_v5",)


def spectra_dir(data_root: Path, subdir: str) -> Path:
    """Where a registry source's spectra live; ``SUBNEPTUNE_FP_OPACITY_DIR``
    overrides the committed FP directory. The version marker in
    ``opacity/loaders.py`` checks whichever directory this names."""
    if subdir in COMMITTED_SPECTRA:
        return _env("SUBNEPTUNE_FP_OPACITY_DIR", HOST / "data" / subdir)
    return Path(data_root) / subdir


PRODUCTION_SCW = HOST / "production_campaign" / "scw"
ARC_SCRIPTS = HOST / "production_campaign" / "arc_scripts"
MLIP_ERRORS = HOST / "models" / "mlip_heldout_errors.json"
# the dipole model behind the 52-point v5 spectra (git-LFS)
DIPOLE_MODEL = HOST / "models" / "water_dipole_v3.model"
PARITY_REF = HOST / "benchmarks" / "reference" / "reference.json"
HOST_DOCS = HOST / "docs"


def rel(p: Path) -> str:
    """``p`` relative to the repository root, for provenance strings."""
    return Path(p).resolve().relative_to(HOST.resolve()).as_posix()
