"""Opacity source interface (SPEC §3).

One interface, N sources. Every source exposes ``kappa_R(T, rho)`` and
``kappa_P(T, rho)`` in SI (m² kg⁻¹) and reports its extrapolation-clamp
statistics so the structure driver can flag unreliable runs. Concrete sources
are registered by name (see ``opacity.registry``); adding one requires no code
change outside a registry entry.
"""

from __future__ import annotations

from abc import ABC, abstractmethod

import numpy as np

from opacity.table import OpacityTable


class OpacitySource(ABC):
    """Abstract Rosseland/Planck opacity source in SI units."""

    name: str

    @abstractmethod
    def kappa_R(self, temperature_k, rho_si):
        """Rosseland mean opacity, m² kg⁻¹."""

    @abstractmethod
    def kappa_P(self, temperature_k, rho_si):
        """Planck mean opacity, m² kg⁻¹."""

    def kappa_both(self, temperature_k, rho_si):
        """Both means; default calls kappa_R/kappa_P (subclasses may optimize)."""
        return self.kappa_R(temperature_k, rho_si), self.kappa_P(temperature_k, rho_si)

    # -- extrapolation gate (SPEC §3) --
    @property
    @abstractmethod
    def clamp_fraction(self) -> float:
        """Fraction of evaluations that were clamped (outside the grid hull)."""

    @abstractmethod
    def reset_clamp_counter(self) -> None:
        """Reset the clamp counter (call before a structure integration)."""

    @property
    def flat_rho_fraction(self) -> float:
        """Fraction of evaluations dilute-flat ρ-floored (justified; A26). 0 by default."""
        return 0.0

    def clamped_mask(self, temperature_k, rho_si) -> np.ndarray:
        """Per-point mask of extrapolated queries. Default: none (no grid)."""
        shape = np.broadcast(np.asarray(temperature_k), np.asarray(rho_si)).shape
        return np.zeros(shape or (1,), dtype=bool)


class GridOpacitySource(OpacitySource):
    """Opacity source backed by a gridded :class:`OpacityTable`."""

    def __init__(self, name: str, table: OpacityTable) -> None:
        self.name = name
        self.table = table

    @property
    def band_cm(self) -> tuple[float, float]:
        return self.table.band_cm

    @property
    def n_points(self) -> int:
        return self.table.n_points

    def kappa_R(self, temperature_k, rho_si):
        return self.table.kappa_rosseland(temperature_k, rho_si)

    def kappa_P(self, temperature_k, rho_si):
        return self.table.kappa_planck(temperature_k, rho_si)

    def kappa_both(self, temperature_k, rho_si):
        return self.table.kappa_both(temperature_k, rho_si)

    def clamped_mask(self, temperature_k, rho_si):
        return self.table.clamped_mask(temperature_k, rho_si)

    @property
    def clamp_fraction(self) -> float:
        return self.table.clamp.fraction

    @property
    def flat_rho_fraction(self) -> float:
        return self.table.clamp.flat_rho_fraction

    @property
    def rho_floor_si(self) -> float:
        return self.table.rho_floor_si

    def reset_clamp_counter(self) -> None:
        self.table.clamp.reset()

    def __repr__(self) -> str:
        lo, hi = self.band_cm
        return (
            f"GridOpacitySource(name={self.name!r}, n_points={self.n_points}, "
            f"band_cm=[{lo:g}, {hi:g}])"
        )


class ConstantOpacitySource(OpacitySource):
    """Constant opacity, independent of (T, ρ). SI (m² kg⁻¹).

    Used for analytic limits — the G1 isothermal/Guillot-skin gate (κ → 0) and
    grey cross-checks — and never touches the extrapolation gate (no grid).
    """

    def __init__(self, kappa_rosseland_si: float, kappa_planck_si: float | None = None,
                 name: str = "constant") -> None:
        if kappa_rosseland_si <= 0:
            raise ValueError("kappa_rosseland_si must be positive")
        self.name = name
        self._kr = float(kappa_rosseland_si)
        self._kp = float(kappa_planck_si if kappa_planck_si is not None else kappa_rosseland_si)

    def kappa_R(self, temperature_k, rho_si):
        return np.full_like(np.asarray(temperature_k, float), self._kr) if not np.isscalar(
            temperature_k
        ) else self._kr

    def kappa_P(self, temperature_k, rho_si):
        return np.full_like(np.asarray(temperature_k, float), self._kp) if not np.isscalar(
            temperature_k
        ) else self._kp

    @property
    def clamp_fraction(self) -> float:
        return 0.0

    def reset_clamp_counter(self) -> None:
        pass
