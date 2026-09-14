"""Evolution package: static-sequence thermal cooling (SPEC §6).

Public API::

    from evolution import cooling_track, cooling_time_ratio
    tr = cooling_track(cfg, "fp", t_eq_k=600)
    print(tr.radius_at_age(1.0), tr.energy_residual)  # R_p at 1 Gyr, G6 residual
"""

from evolution.history import HistoryComparison, compare_histories
from evolution.radiogenic import radiogenic_power_w, t_int_floor_k
from evolution.track import (
    SECONDS_PER_GYR,
    CoolingTrack,
    age_grid_problem,
    cooling_time_ratio,
    cooling_track,
)

__all__ = [
    "CoolingTrack",
    "cooling_track",
    "cooling_time_ratio",
    "age_grid_problem",
    "SECONDS_PER_GYR",
    "HistoryComparison",
    "compare_histories",
    "radiogenic_power_w",
    "t_int_floor_k",
]
