"""History-effect study: cooling under composite vs line-list opacity (A29).

The present-day structure was set by cooling under some opacity. If line lists
mislocate the RCB at depth during the hot early phase, the whole entropy
trajectory is biased. This compares two cooling tracks that share the initial
condition (same M_c, f_env, T_eq, same initial interior entropy) but use
different opacity, and reports how the histories diverge and whether they
reconverge.

Deliverables (per SPEC §6 evolution + owner request):
  - radius(age), RCB density(age), interior entropy difference vs time;
  - the FP-relevant epoch (when either track's radiative zone reaches ρ≥109);
  - peak radius divergence and whether it persists to Gyr;
  - the first 100 Myr specifically (atmospheric-escape sculpting window).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from config import Config
from eos import AquaEOS
from evolution.track import SECONDS_PER_GYR, CoolingTrack, cooling_track
from structure.model import FP_GRID_FLOOR_SI


@dataclass
class HistoryComparison:
    """Two cooling tracks on a common age grid + divergence diagnostics."""

    source_a: str
    source_b: str
    t_eq_k: float
    f_env: float
    age_gyr: np.ndarray  # common grid
    t_int_k: np.ndarray  # T_int on the age grid (≈ shared by both tracks)
    radius_a: np.ndarray
    radius_b: np.ndarray
    rcb_rho_a: np.ndarray  # deepest RCB density [kg/m³] vs age
    rcb_rho_b: np.ndarray
    entropy_a: np.ndarray
    entropy_b: np.ndarray

    @property
    def radius_frac_diff(self) -> np.ndarray:
        """(R_a − R_b)/R_b vs age."""
        return (self.radius_a - self.radius_b) / self.radius_b

    @property
    def peak_radius_frac_diff(self) -> float:
        return float(np.max(np.abs(self.radius_frac_diff)))

    @property
    def entropy_frac_diff(self) -> np.ndarray:
        return (self.entropy_a - self.entropy_b) / self.entropy_b

    def fp_relevant_mask(self) -> np.ndarray:
        """Ages where either track's RCB is at ρ ≥ FP floor (FP-relevant epoch)."""
        return (self.rcb_rho_a >= FP_GRID_FLOOR_SI) | (self.rcb_rho_b >= FP_GRID_FLOOR_SI)

    def fp_relevant_epoch_gyr(self) -> tuple[float, float] | None:
        """(age_start, age_end) over which the FP-relevant epoch spans, or None."""
        m = self.fp_relevant_mask()
        if not m.any():
            return None
        ages = self.age_gyr[m]
        return float(ages.min()), float(ages.max())

    def diff_at_age(self, age_gyr: float) -> float:
        """Radius fractional difference at a given age (interpolated)."""
        return float(np.interp(age_gyr, self.age_gyr, self.radius_frac_diff))

    def escape_window_extrema(self, window_gyr: float = 0.1) -> dict:
        """Largest |ΔR_p/R_p| and |Δs/s| inside the first ``window_gyr`` (A40).

        Atmospheric escape (not modelled in v1) sculpts planets early, so any
        composite-vs-line-list difference inside this window would be an
        irreversible imprint; the maxima are the numbers to quote. The window
        edge is included by interpolation so a sparse age grid cannot miss it.
        """
        mask = self.age_gyr <= window_gyr
        ages = np.concatenate([self.age_gyr[mask], [window_gyr]])
        dr = np.abs(np.concatenate([
            self.radius_frac_diff[mask], [self.diff_at_age(window_gyr)]]))
        ds = np.abs(np.concatenate([
            self.entropy_frac_diff[mask],
            [float(np.interp(window_gyr, self.age_gyr, self.entropy_frac_diff))]]))
        return {
            "window_gyr": float(window_gyr),
            "n_grid_points": int(mask.sum()),
            "age_min_gyr": float(self.age_gyr.min()),
            "max_abs_radius_frac": float(dr.max()),
            "age_at_max_radius_gyr": float(ages[int(dr.argmax())]),
            "max_abs_entropy_frac": float(ds.max()),
            "age_at_max_entropy_gyr": float(ages[int(ds.argmax())]),
        }

    def first_age_rcb_above(self, rho_si: float = FP_GRID_FLOOR_SI) -> dict:
        """First age at which each track's deepest RCB reaches ρ ≥ ``rho_si``.

        The crossing is interpolated in log age across the bracketing grid points;
        ``None`` means the track never gets there within its age range (A40).
        """

        def crossing(rcb_rho: np.ndarray) -> float | None:
            above = rcb_rho >= rho_si
            if not above.any():
                return None
            i = int(np.argmax(above))
            if i == 0:
                return float(self.age_gyr[0])
            r0, r1 = rcb_rho[i - 1], rcb_rho[i]
            a0, a1 = self.age_gyr[i - 1], self.age_gyr[i]
            if r1 == r0:
                return float(a1)
            w = (rho_si - r0) / (r1 - r0)
            return float(np.exp(np.log(a0) + w * (np.log(a1) - np.log(a0))))

        return {
            "rho_si": float(rho_si),
            self.source_a: crossing(self.rcb_rho_a),
            self.source_b: crossing(self.rcb_rho_b),
        }

    def radiogenic_floor_age_gyr(self, mass_core_mearth: float) -> float | None:
        """Age at which T_int reaches the radiogenic floor (A30): beyond this the
        planet is radiogenically supported and cools no further, truncating the
        FP-relevant epoch. None if the floor is never reached within the track.

        Solved self-consistently over the age grid: floor depends on R_p(age) and
        age; the crossing is where the (no-radiogenic) track T_int drops to it.
        """
        from constants import R_EARTH
        from evolution.radiogenic import t_int_floor_k

        floor = np.array([
            t_int_floor_k(mass_core_mearth, self.radius_b[i] * R_EARTH, self.age_gyr[i])
            for i in range(self.age_gyr.size)
        ])
        below = self.t_int_k <= floor  # T_int has cooled to/under the floor
        if not below.any():
            return None
        return float(self.age_gyr[np.argmax(below)])


def _deepest_rcb_rho(track: CoolingTrack, eos: AquaEOS, cfg: Config) -> np.ndarray:
    """Deepest RCB density at each track node (rebuild to read RCB densities)."""
    from dataclasses import replace

    from structure import build_structure

    c = replace(cfg, envelope=replace(cfg.envelope, f_env=track.f_env))
    out = []
    rp_guess = None
    for ti in track.t_int_k:
        res = build_structure(c, track.t_eq_k, float(ti), track.source, eos=eos,
                             n_sample=400, band=track.band_cm, rp_guess_rearth=rp_guess)
        out.append(max(res.rcb_densities_si) if res.rcb_densities_si else 0.0)
        rp_guess = res.radius_p_rearth
    return np.array(out)


def compare_histories(
    cfg: Config,
    source_a: str,
    source_b: str,
    t_eq_k: float,
    eos: AquaEOS | None = None,
    t_int_grid_k: np.ndarray | None = None,
    n_age: int = 60,
) -> HistoryComparison:
    """Compare cooling under two opacity sources on a common age grid.

    Both tracks share (M_c, f_env, T_eq) and the same T_int sequence, so they
    start from the same interior state; the opacity difference then drives the
    entropy/radius trajectories apart.
    """
    eos = eos or AquaEOS.from_config(cfg)
    band = None  # each track uses its source's default (config) band
    ta = cooling_track(cfg, source_a, t_eq_k, eos=eos, t_int_grid_k=t_int_grid_k, n_sample=400)
    tb = cooling_track(cfg, source_b, t_eq_k, eos=eos, t_int_grid_k=t_int_grid_k, n_sample=400)
    rho_a = _deepest_rcb_rho(ta, eos, cfg)
    rho_b = _deepest_rcb_rho(tb, eos, cfg)

    # Common age grid spanning both tracks (log-spaced; cooling covers decades).
    lo = max(ta.age_gyr[ta.age_gyr > 0].min(), tb.age_gyr[tb.age_gyr > 0].min())
    hi = min(ta.age_gyr.max(), tb.age_gyr.max())
    age = np.geomspace(lo, hi, n_age)

    def on_age(track, y):
        idx = np.argsort(track.age_gyr)
        return np.interp(age, track.age_gyr[idx], y[idx])

    return HistoryComparison(
        source_a=source_a, source_b=source_b, t_eq_k=t_eq_k, f_env=cfg.envelope.f_env,
        age_gyr=age, t_int_k=on_age(tb, tb.t_int_k),
        radius_a=on_age(ta, ta.radius_p_rearth), radius_b=on_age(tb, tb.radius_p_rearth),
        rcb_rho_a=on_age(ta, rho_a), rcb_rho_b=on_age(tb, rho_b),
        entropy_a=on_age(ta, ta.entropy_j_kg_k), entropy_b=on_age(tb, tb.entropy_j_kg_k),
    )
