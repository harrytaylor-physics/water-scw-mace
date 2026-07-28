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

    @property
    def clamp_fraction(self) -> float:
        return self.table.clamp.fraction

    def reset_clamp_counter(self) -> None:
        self.table.clamp.reset()

    def __repr__(self) -> str:
        lo, hi = self.band_cm
        return (
            f"GridOpacitySource(name={self.name!r}, n_points={self.n_points}, "
            f"band_cm=[{lo:g}, {hi:g}])"
        )
