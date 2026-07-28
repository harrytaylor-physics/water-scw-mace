"""Opacity package: one interface, N sources (SPEC §3).

Public API::

    from opacity import get_source
    src = get_source("fp", cfg)
    kappa_R = src.kappa_R(T, rho)   # m^2/kg
    kappa_P = src.kappa_P(T, rho)   # m^2/kg
"""

from opacity.build_tables import build_source, clear_cache, compute_point_means, get_source
from opacity.registry import REGISTRY, SourceSpec, get_spec
from opacity.source import GridOpacitySource, OpacitySource
from opacity.table import CM2_PER_G_TO_SI, ClampCounter, OpacityTable

__all__ = [
    "REGISTRY",
    "SourceSpec",
    "get_spec",
    "OpacitySource",
    "GridOpacitySource",
    "OpacityTable",
    "ClampCounter",
    "CM2_PER_G_TO_SI",
    "build_source",
    "get_source",
    "compute_point_means",
    "clear_cache",
]
