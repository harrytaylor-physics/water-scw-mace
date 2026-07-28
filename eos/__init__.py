"""EOS package: AQUA water EOS + rocky core (SPEC §4).

Public API::

    from eos import AquaEOS
    eos = AquaEOS.from_config(cfg)
    st = eos.state(P_pa, T_k)        # EOSState(rho, grad_ad, entropy, phase)
    rho = eos.rho(P_pa, T_k)         # kg/m^3

    from eos import core_radius_si
    R_c = core_radius_si(mass_mearth=5.0)   # metres
"""

from eos.aqua import AquaEOS
from eos.base import EOS, EOSState
from eos.core import (
    EARTH_LIKE_CMF,
    core_mass_si,
    core_radius_rearth,
    core_radius_si,
)
from eos.loader import AquaTables, load_aqua_tables

__all__ = [
    "EOS",
    "EOSState",
    "AquaEOS",
    "AquaTables",
    "load_aqua_tables",
    "core_radius_rearth",
    "core_radius_si",
    "core_mass_si",
    "EARTH_LIKE_CMF",
]
