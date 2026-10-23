"""Gridded opacity table + interpolation with the extrapolation gate (SPEC §3).

The FP grid is scattered in (T, ρ), not a tensor-product grid (AMENDMENTS A6),
so "bilinear in (T, log ρ) on log κ" is realized as linear interpolation on a
Delaunay triangulation of (T, log₁₀ρ) over log₁₀κ. Extrapolation policy is a
gate, not a guess: outside the convex hull we clamp to the nearest grid value,
count every clamped call, and expose the clamped fraction so a run that clamps
too often is flagged unreliable in its output metadata.

Tables are stored in SI (κ in m² kg⁻¹, ρ in kg m⁻³); the cm²/g plausibility
assertion (SPEC §2) is applied earlier, at build time, in the native units.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator

# cm²/g -> m²/kg
CM2_PER_G_TO_SI = 0.1


@dataclass
class ClampCounter:
    """Counts opacity evaluations and how they were extrapolated.

    ``clamped_calls`` are genuine out-of-hull clamps (unreliable, gate the run);
    ``flat_rho_calls`` are dilute-limit ρ-floored calls (physically justified,
    reported separately, A26).
    """

    total_calls: int = 0
    clamped_calls: int = 0
    flat_rho_calls: int = 0

    def reset(self) -> None:
        self.total_calls = 0
        self.clamped_calls = 0
        self.flat_rho_calls = 0

    @property
    def fraction(self) -> float:
        return self.clamped_calls / self.total_calls if self.total_calls else 0.0

    @property
    def flat_rho_fraction(self) -> float:
        return self.flat_rho_calls / self.total_calls if self.total_calls else 0.0


@dataclass
class OpacityTable:
    """Scattered (T, ρ) opacity table with triangulated log-κ interpolation.

    Parameters
    ----------
    temperature_k, rho_si:
        Grid point coordinates (1-D, same length). ρ in kg m⁻³.
    kappa_rosseland_si, kappa_planck_si:
        Rosseland/Planck means at each grid point, in m² kg⁻¹.
    band_cm:
        Common band the means were integrated over (recorded for provenance).
    """

    temperature_k: np.ndarray
    rho_si: np.ndarray
    kappa_rosseland_si: np.ndarray
    kappa_planck_si: np.ndarray
    band_cm: tuple[float, float]
    # Dilute-limit flat-κ extrapolation below the grid ρ floor (A26): for line
    # lists, hold κ flat in ρ under ρ_min rather than clamping to nearest.
    dilute_flat: bool = False
    clamp: ClampCounter = field(default_factory=ClampCounter)

    def __post_init__(self) -> None:
        self.temperature_k = np.asarray(self.temperature_k, dtype=float)
        self.rho_si = np.asarray(self.rho_si, dtype=float)
        self.kappa_rosseland_si = np.asarray(self.kappa_rosseland_si, dtype=float)
        self.kappa_planck_si = np.asarray(self.kappa_planck_si, dtype=float)
        n = self.temperature_k.size
        if not (
            self.rho_si.size
            == self.kappa_rosseland_si.size
            == self.kappa_planck_si.size
            == n
        ):
            raise ValueError("OpacityTable arrays must all have the same length")
        if n < 3:
            raise ValueError("triangulated interpolation needs >= 3 grid points")
        if np.any(self.rho_si <= 0) or np.any(self.temperature_k <= 0):
            raise ValueError("grid T and rho must be positive")
        if np.any(self.kappa_rosseland_si <= 0) or np.any(self.kappa_planck_si <= 0):
            raise ValueError("grid opacities must be positive")

        # Interpolate log10 κ over (T, log10 ρ). rescale=True normalizes the
        # disparate T (~10³) and log10ρ (~2-3) scales before triangulation.
        pts = np.column_stack([self.temperature_k, np.log10(self.rho_si)])
        vals = np.column_stack(
            [np.log10(self.kappa_rosseland_si), np.log10(self.kappa_planck_si)]
        )
        self._lin = LinearNDInterpolator(pts, vals, rescale=True)
        self._near = NearestNDInterpolator(pts, vals, rescale=True)
        self._logrho_min = float(np.log10(self.rho_si.min()))

    @property
    def n_points(self) -> int:
        return self.temperature_k.size

    @property
    def rho_floor_si(self) -> float:
        """Grid ρ floor below which dilute-flat extrapolation applies (A26)."""
        return float(self.rho_si.min())

    def _eval_log10(self, temperature_k, rho_si) -> np.ndarray:
        """Return log10 κ (SI) for both means; shape (N, 2). Counts clamps.

        Dilute-flat (A26): for line lists, a query below the grid ρ floor has its
        ρ pinned to ρ_min (flat κ in ρ) and is counted as ``flat_rho`` — a
        physically-justified extrapolation, not an unreliable hull clamp.
        """
        t = np.atleast_1d(np.asarray(temperature_k, dtype=float))
        r = np.atleast_1d(np.asarray(rho_si, dtype=float))
        t, r = np.broadcast_arrays(t, r)
        logr = np.log10(r.ravel())

        n_flat = 0
        if self.dilute_flat:
            below = logr < self._logrho_min
            n_flat = int(below.sum())
            logr = np.where(below, self._logrho_min, logr)  # flat κ in ρ under floor

        pts = np.column_stack([t.ravel(), logr])
        out = self._lin(pts)  # (N, 2); NaN rows are outside the hull
        outside = np.isnan(out[:, 0])
        n_out = int(outside.sum())
        if n_out:
            out[outside] = self._near(pts[outside])

        self.clamp.total_calls += pts.shape[0]
        self.clamp.clamped_calls += n_out
        self.clamp.flat_rho_calls += n_flat
        return out

    def kappa_rosseland(self, temperature_k, rho_si):
        """Rosseland mean κ_R(T, ρ) in m² kg⁻¹ (scalar in -> float out)."""
        out = 10.0 ** self._eval_log10(temperature_k, rho_si)[:, 0]
        return float(out[0]) if np.isscalar(temperature_k) and np.isscalar(rho_si) else out

    def kappa_planck(self, temperature_k, rho_si):
        """Planck mean κ_P(T, ρ) in m² kg⁻¹ (scalar in -> float out)."""
        out = 10.0 ** self._eval_log10(temperature_k, rho_si)[:, 1]
        return float(out[0]) if np.isscalar(temperature_k) and np.isscalar(rho_si) else out

    def kappa_both(self, temperature_k, rho_si):
        """Both means in one evaluation: returns (κ_R, κ_P), each m² kg⁻¹."""
        out = 10.0 ** self._eval_log10(temperature_k, rho_si)
        if np.isscalar(temperature_k) and np.isscalar(rho_si):
            return float(out[0, 0]), float(out[0, 1])
        return out[:, 0], out[:, 1]

    def clamped_mask(self, temperature_k, rho_si) -> np.ndarray:
        """Boolean mask (per point) of queries outside the grid hull.

        Dilute-flat ρ-floored queries are NOT counted as clamped (they are
        reliable via the dilute limit, A26). Does not mutate the clamp counter
        (used by the composite source to attribute reliability to the used
        sub-source).
        """
        t = np.atleast_1d(np.asarray(temperature_k, dtype=float))
        r = np.atleast_1d(np.asarray(rho_si, dtype=float))
        t, r = np.broadcast_arrays(t, r)
        logr = np.log10(r.ravel())
        if self.dilute_flat:
            logr = np.maximum(logr, self._logrho_min)  # floored ρ -> in-hull, reliable
        pts = np.column_stack([t.ravel(), logr])
        return np.isnan(self._lin(pts)[:, 0])
