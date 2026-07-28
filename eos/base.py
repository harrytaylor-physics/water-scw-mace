"""Equation-of-state interface (SPEC §4).

AQUA is the v1 envelope EOS, but the interface is written so a future custom
EOS (e.g. the group's MLIP-derived surface) drops in beside it: implement
:class:`EOS` and the structure code is unchanged.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass


@dataclass(frozen=True)
class EOSState:
    """Thermodynamic state at a point, SI units.

    rho             : density, kg m⁻³
    grad_ad         : adiabatic gradient (dlnT/dlnP)_S, dimensionless
    entropy         : specific entropy, J kg⁻¹ K⁻¹
    internal_energy : specific internal energy, J kg⁻¹ (evolution/G6)
    sound_speed     : bulk sound speed, m s⁻¹
    phase           : AQUA phase id (see the table header; e.g. 1 = vapor+liquid)
    """

    rho: float
    grad_ad: float
    entropy: float
    internal_energy: float
    sound_speed: float
    phase: int

    @property
    def two_phase(self) -> bool:
        """True in a coexistence/multi-phase region (phase id <= 2)."""
        return self.phase <= 2


class EOS(ABC):
    """Abstract EOS providing state as a function of (P, T), SI units."""

    name: str

    @abstractmethod
    def state(self, pressure_pa, temperature_k) -> EOSState:
        """Full thermodynamic state at (P, T)."""

    @abstractmethod
    def rho(self, pressure_pa, temperature_k):
        """Density ρ(P, T), kg m⁻³."""

    @abstractmethod
    def grad_ad(self, pressure_pa, temperature_k):
        """Adiabatic gradient (dlnT/dlnP)_S at (P, T), dimensionless."""

    @abstractmethod
    def entropy(self, pressure_pa, temperature_k):
        """Specific entropy s(P, T), J kg⁻¹ K⁻¹."""
