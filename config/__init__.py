"""Run configuration package (SPEC §1)."""

from config.schema import (
    Boundary,
    Config,
    Core,
    Envelope,
    Interior,
    Irradiation,
    Numerics,
    Opacity,
    Paths,
    load_config,
    to_dict,
)

__all__ = [
    "Boundary",
    "Config",
    "Core",
    "Envelope",
    "Interior",
    "Irradiation",
    "Numerics",
    "Opacity",
    "Paths",
    "load_config",
    "to_dict",
]
