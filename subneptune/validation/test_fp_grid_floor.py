"""The configured FP box corner must match the data (A54).

``opacity.fp_grid_t_floor_k`` / ``opacity.fp_grid_rho_floor_kg_m3`` define what
"inside the FP box" means for every scan, figure and diagnostic. They used to be
module constants (1000 K, 108.796) and the extension campaign moved them, which
is exactly the kind of drift that goes unnoticed. This gate ties them to the
built table's actual minima, so config cannot disagree with the data.
"""

from __future__ import annotations

import pytest

from figures._common import fp_floors
from opacity import get_source
from structure import FP_GRID_FLOOR_SI, FP_GRID_T_FLOOR_K


def test_config_floors_are_the_tables_actual_minima(cfg):
    table = get_source("fp", cfg).table
    t_floor, rho_floor = fp_floors(cfg)
    assert t_floor == pytest.approx(float(table.temperature_k.min()))
    assert rho_floor == pytest.approx(float(table.rho_si.min()))


def test_module_constants_mirror_the_schema_defaults(cfg):
    """The no-Config fallbacks must not drift from the configured values."""
    t_floor, rho_floor = fp_floors(cfg)
    assert FP_GRID_T_FLOOR_K == pytest.approx(t_floor)
    assert FP_GRID_FLOOR_SI == pytest.approx(rho_floor)


def test_the_blend_window_starts_at_or_above_the_fp_data_floor(cfg):
    """A54: below the FP grid's ρ floor the composite would hand κ to clamped FP.

    The blend's lower edge is where FP starts to carry weight, so putting it
    under the grid's own floor means the dense side is nearest-neighbour
    extrapolation by construction — which is why the proposed 3-10 window was
    rejected in favour of 10-30.
    """
    lo, _hi = cfg.opacity.blend_window_kg_m3
    _t, rho_floor = fp_floors(cfg)
    assert lo >= rho_floor, (
        f"blend window starts at rho={lo:g}, below the FP grid floor "
        f"{rho_floor:g}: FP would be clamped wherever it dominates")
