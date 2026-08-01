"""Structure package: static interior models (SPEC §5).

Public API::

    from structure import build_structure, save_structure
    res = build_structure(cfg, t_eq_k=500, t_int_k=100, source_name="fp")
    print(res.radius_p_rearth, res.rcb_pressures_pa)
"""

from structure.boundary import (
    guillot_temperature,
    guillot_temperature4,
    skin_temperature,
    t_irr_from_t_eq,
)
from structure.integrator import (
    EnvelopeProfile,
    MatchState,
    integrate_interior,
    integrate_outer,
    luminosity_int,
)
from structure.model import (
    StructureResult,
    build_at_entropy,
    build_structure,
    save_structure,
)
from structure.shooting import ShootResult, shoot_radius

__all__ = [
    "guillot_temperature",
    "guillot_temperature4",
    "skin_temperature",
    "t_irr_from_t_eq",
    "MatchState",
    "EnvelopeProfile",
    "integrate_outer",
    "integrate_interior",
    "luminosity_int",
    "shoot_radius",
    "ShootResult",
    "build_structure",
    "build_at_entropy",
    "save_structure",
    "StructureResult",
]
