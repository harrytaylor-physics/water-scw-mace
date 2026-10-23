"""Run configuration package (SPEC §1)."""

from config.schema import (
    Anchor,
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
    "Anchor",
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
