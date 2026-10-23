"""AQUA table parsing, reshaping, and inversion pre-processing (SPEC §2, §4).

The AQUA ``*_rhot_*`` table is on a regular (ρ, T) grid (ρ log-spaced, fastest
axis T). This loader detects the header and grid shape robustly, reshapes every
column to (n_ρ, n_T), and precomputes the strictly-increasing per-column log P
needed to invert P(ρ,T) -> ρ(P,T) (AMENDMENTS A2, A8). Results are cached by
path so repeated construction is cheap.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import numpy as np

# Column order in aqua_eos_rhot_v1_0.dat (see table header).
_COLS = ("rho", "temp", "press", "ad_grad", "s", "u", "c", "mmw", "x_ion", "x_d", "phase")


@dataclass(frozen=True)
class AquaTables:
    """Parsed AQUA table on its native (n_ρ, n_T) grid, plus inversion aids."""

    logrho_grid: np.ndarray  # (n_ρ,) log10 ρ [ρ in kg m⁻³], strictly increasing
    logT_grid: np.ndarray  # (n_T,) log10 T [T in K], strictly increasing
    logP_strict: np.ndarray  # (n_ρ, n_T) log10 P, strictly incr along ρ per column
    press: np.ndarray  # (n_ρ, n_T) P [Pa], native (may have plateaus)
    ad_grad: np.ndarray  # (n_ρ, n_T) (dlnT/dlnP)_S
    entropy: np.ndarray  # (n_ρ, n_T) s [J kg⁻¹ K⁻¹]
    u: np.ndarray  # (n_ρ, n_T) internal energy [J kg⁻¹]
    sound: np.ndarray  # (n_ρ, n_T) bulk sound speed [m s⁻¹]
    phase: np.ndarray  # (n_ρ, n_T) phase id (int)

    @property
    def n_rho(self) -> int:
        return self.logrho_grid.size

    @property
    def n_T(self) -> int:
        return self.logT_grid.size


def _detect_header(path: Path) -> int:
    """Return the number of leading non-data lines."""
    with path.open() as fh:
        for i, line in enumerate(fh):
            tok = line.split()
            if len(tok) != len(_COLS):
                continue
            try:
                float(tok[0]), float(tok[1]), int(tok[-1])
                return i
            except ValueError:
                continue
    raise ValueError(f"no data rows found in {path}")


def _strictly_increasing(x: np.ndarray, tiny: float = 1e-12) -> np.ndarray:
    """Make a non-decreasing-ish array strictly increasing by nudging only the
    flat/decreasing entries (single-phase nodes are left unchanged; AMENDMENTS
    A8)."""
    x = np.maximum.accumulate(x)  # kill any tiny decreases
    dup = np.concatenate([[False], np.diff(x) <= 0.0])
    return x + np.cumsum(dup) * tiny


@lru_cache(maxsize=4)
def load_aqua_tables(path_str: str) -> AquaTables:
    """Load and pre-process an AQUA (ρ,T) table (cached by path)."""
    path = Path(path_str)
    header = _detect_header(path)
    data = np.genfromtxt(path, skip_header=header)
    if data.shape[1] != len(_COLS):
        raise ValueError(f"{path}: expected {len(_COLS)} columns, got {data.shape[1]}")

    rho_flat = data[:, 0]
    # T varies fastest -> n_T is the run length of the first ρ value.
    n_T = int(np.argmax(rho_flat != rho_flat[0]))
    if n_T == 0:
        raise ValueError(f"{path}: could not infer grid shape")
    n_rho, rem = divmod(data.shape[0], n_T)
    if rem != 0:
        raise ValueError(f"{path}: {data.shape[0]} rows not divisible by n_T={n_T}")

    col = {name: data[:, i].reshape(n_rho, n_T) for i, name in enumerate(_COLS)}
    rho_grid = col["rho"][:, 0]
    T_grid = col["temp"][0, :]
    if not (np.all(np.diff(rho_grid) > 0) and np.all(np.diff(T_grid) > 0)):
        raise ValueError(f"{path}: ρ/T axes are not both strictly increasing")

    logP_strict = np.column_stack(
        [_strictly_increasing(np.log10(col["press"][:, j])) for j in range(n_T)]
    )

    return AquaTables(
        logrho_grid=np.log10(rho_grid),
        logT_grid=np.log10(T_grid),
        logP_strict=logP_strict,
        press=col["press"],
        ad_grad=col["ad_grad"],
        entropy=col["s"],
        u=col["u"],
        sound=col["c"],
        phase=col["phase"].astype(int),
    )
