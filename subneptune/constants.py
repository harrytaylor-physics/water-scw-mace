"""Physical constants, defined once and reused everywhere (SPEC §1).

Sourced from ``scipy.constants`` (CODATA 2018) where available; astronomical
reference values use IAU 2015 nominal values. See ``docs/constants.md`` for
the full table, units, and citations. Do not re-hardcode any of these
elsewhere in the code.
"""

from __future__ import annotations

import scipy.constants as _sc

# --- CODATA 2018 via scipy.constants ---------------------------------------
G: float = _sc.G  # m^3 kg^-1 s^-2, Newtonian gravitation
SIGMA_SB: float = _sc.Stefan_Boltzmann  # W m^-2 K^-4, Stefan-Boltzmann
K_B: float = _sc.k  # J K^-1, Boltzmann
H_PLANCK: float = _sc.h  # J s, Planck
C_LIGHT: float = _sc.c  # m s^-1, speed of light
N_A: float = _sc.N_A  # mol^-1, Avogadro
M_U: float = _sc.m_u  # kg, atomic mass unit
BAR: float = _sc.bar  # Pa, one bar

# --- IAU 2015 nominal (not in scipy.constants) -----------------------------
GM_EARTH: float = 3.986004e14  # m^3 s^-2, nominal terrestrial mass parameter
M_EARTH: float = GM_EARTH / G  # kg, Earth mass (kept consistent with G)
R_EARTH: float = 6.3781e6  # m, nominal Earth equatorial radius

# --- H2O molecular data ----------------------------------------------------
MMW_H2O: float = 1.801527e-2  # kg mol^-1, matches AQUA `mmw` in dilute limit

__all__ = [
    "G",
    "SIGMA_SB",
    "K_B",
    "H_PLANCK",
    "C_LIGHT",
    "N_A",
    "M_U",
    "BAR",
    "GM_EARTH",
    "M_EARTH",
    "R_EARTH",
    "MMW_H2O",
]
