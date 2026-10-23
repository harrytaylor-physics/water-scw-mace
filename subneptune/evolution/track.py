"""Static-sequence thermal evolution (SPEC §6).

Lopez & Fortney (2014); Ginzburg, Schlichting & Sari (2016). The envelope cools
along a sequence of static models parameterized by internal temperature T_int
(equivalently the deep-adiabat entropy). Each model supplies its total envelope
energy E and luminosity L = L_int; the age follows from dE/dt = −L:

    age(T_int) = ∫_{T_int}^{T_max} (dE/dT'_int) / L(T'_int) dT'_int,

since dE/dt = −L ⇒ dt = −dE/L and E increases with T_int (hotter = more energy).
Irradiation is held constant through T_eq; no escape in v1. Outputs R_p(t), L(t),
RCB(t) per source. Energy bookkeeping is gate G6.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.integrate import cumulative_simpson, simpson

# np.trapz removed in numpy 2.0.
try:
    from numpy import trapezoid as _trapz
except ImportError:  # numpy < 2.0
    from numpy import trapz as _trapz

from config import Config
from eos import AquaEOS
from structure import build_structure

SECONDS_PER_GYR = 3.15576e16


@dataclass
class CoolingTrack:
    """A cooling sequence for one source at fixed (M_c, f_env, T_eq)."""

    source: str
    t_eq_k: float
    f_env: float
    band_cm: tuple[float, float] | None

    t_int_k: np.ndarray  # ascending
    entropy_j_kg_k: np.ndarray
    radius_p_rearth: np.ndarray
    luminosity_w: np.ndarray
    energy_j: np.ndarray
    rcb_bar: np.ndarray
    age_gyr: np.ndarray  # age since the hot start (0 at T_int_max)

    energy_residual: float  # G6: |∫L dt − ΔE| / |ΔE|

    # Deepest RCB density at each node, read off the structures this track
    # already built (A56). Before, evolution.history rebuilt every structure a
    # second time just to recover this, doubling the cost of every history
    # comparison. None only for tracks constructed by hand in tests.
    rcb_rho_si: np.ndarray | None = None

    # -- interpolators (age is monotone increasing as T_int decreases) --
    def _age_sorted(self):
        idx = np.argsort(self.age_gyr)
        return idx

    def radius_at_age(self, age_gyr):
        idx = self._age_sorted()
        return np.interp(age_gyr, self.age_gyr[idx], self.radius_p_rearth[idx])

    def luminosity_at_age(self, age_gyr):
        idx = self._age_sorted()
        return np.interp(age_gyr, self.age_gyr[idx], self.luminosity_w[idx])

    def age_at_radius(self, radius_rearth):
        """Age when the planet has contracted to ``radius_rearth`` (R_p is
        monotone increasing in T_int, i.e. decreasing in age)."""
        idx = np.argsort(self.radius_p_rearth)
        return float(np.interp(radius_rearth, self.radius_p_rearth[idx], self.age_gyr[idx]))

    def f_env_at(self, radius_rearth, age_gyr) -> float:
        """Placeholder for the inferred-f_env shift (needs an f_env sweep;
        provided at the driver level, not on a single track)."""
        raise NotImplementedError("use evolution.inferred_f_env_shift over an f_env grid")


def cooling_track(
    cfg: Config,
    source_name: str,
    t_eq_k: float,
    eos: AquaEOS | None = None,
    t_int_grid_k: np.ndarray | None = None,
    band: tuple[float, float] | None = None,
    n_sample: int = 400,
) -> CoolingTrack:
    """Build the cooling track for one source at fixed (M_c, f_env, T_eq)."""
    eos = eos or AquaEOS.from_config(cfg)
    if t_int_grid_k is None:
        # Fine, log-spaced: the age integral is dominated by the cool, low-L end,
        # so resolution there controls the G6 energy residual. With the Simpson
        # quadrature n=60 gives ~2e-3, safely under the 1% gate.
        t_int_grid_k = np.geomspace(25.0, 500.0, 60)
    t_int_grid_k = np.sort(np.asarray(t_int_grid_k, dtype=float))

    # Warm-start each model from the previous R_p (T_int grid is ascending, so
    # R_p rises smoothly): successive shoots converge in a couple of iterations.
    rows = []
    rp_guess = None
    for ti in t_int_grid_k:
        res = build_structure(cfg, t_eq_k, float(ti), source_name, eos=eos,
                              n_sample=n_sample, band=band, rp_guess_rearth=rp_guess)
        rows.append(res)
        rp_guess = res.radius_p_rearth
    entropy = np.array([r.interior_entropy_j_kg_k for r in rows])
    radius = np.array([r.radius_p_rearth for r in rows])
    lum = np.array([r.luminosity_w for r in rows])
    energy = np.array([r.envelope_energy_j for r in rows])
    rcb = np.array([max(r.rcb_pressures_pa) / 1e5 if r.rcb_pressures_pa else np.nan for r in rows])
    rcb_rho = np.array([max(r.rcb_densities_si) if r.rcb_densities_si else 0.0
                        for r in rows])
    band_used = rows[0].band_cm

    # Cooling: dt = dE/L. Age from the hot end via a 4th-order Simpson cumulative
    # integral of 1/L against E (E monotone increasing with T_int). Simpson's
    # O(h⁴) accuracy lets a modest grid pass the G6 gate on the steep L(T_int).
    inv_l = 1.0 / lum
    time_from_bottom = cumulative_simpson(y=inv_l, x=energy, initial=0.0)  # s, 0 at coldest
    age_from_top = time_from_bottom[-1] - time_from_bottom  # 0 at hottest, grows as it cools
    age_gyr = age_from_top / SECONDS_PER_GYR

    # G6: independent Simpson quadrature of ∫L dt over the age array vs radiated ΔE.
    order = np.argsort(age_gyr)
    age_ord = age_gyr[order] * SECONDS_PER_GYR
    l_ord = lum[order]
    e_ord = energy[order]
    int_l_dt = simpson(y=l_ord, x=age_ord)
    delta_e = e_ord[0] - e_ord[-1]  # E(youngest) − E(oldest) = energy radiated
    energy_residual = abs(int_l_dt - delta_e) / abs(delta_e)

    return CoolingTrack(
        source=source_name, t_eq_k=t_eq_k, f_env=cfg.envelope.f_env, band_cm=band_used,
        t_int_k=t_int_grid_k, entropy_j_kg_k=entropy, radius_p_rearth=radius,
        luminosity_w=lum, energy_j=energy, rcb_bar=rcb, age_gyr=age_gyr,
        energy_residual=energy_residual, rcb_rho_si=rcb_rho,
    )


def age_grid_problem(track: CoolingTrack) -> str | None:
    """Message if a track's age integral came out unphysical, else None (A38).

    The cumulative-Simpson age integral overshoots when L(T_int) varies steeply
    across too few nodes, producing negative or non-monotone ages. Callers should
    re-run on a finer T_int grid rather than quietly interpolating nonsense.
    """
    order = np.argsort(track.t_int_k)[::-1]  # hottest first = youngest first
    age = np.asarray(track.age_gyr, float)[order]
    if np.any(age < 0):
        return (f"{track.source}: {int((age < 0).sum())} of {age.size} ages are negative "
                f"(min {age.min():.3f} Gyr) on a {age.size}-node T_int grid")
    if np.any(np.diff(age) < 0):
        return (f"{track.source}: age is not monotone in T_int on a {age.size}-node grid "
                f"({int((np.diff(age) < 0).sum())} decreasing steps)")
    return None


def cooling_time_ratio(track_a: CoolingTrack, track_b: CoolingTrack, radius_rearth: float) -> float:
    """Ratio of ages at which two sources have contracted to the same R_p
    (the definitive 'stay puffy longer' number: age_FP / age_LL)."""
    return track_a.age_at_radius(radius_rearth) / track_b.age_at_radius(radius_rearth)
