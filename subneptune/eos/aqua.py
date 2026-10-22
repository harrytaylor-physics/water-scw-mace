"""AQUA water EOS on (P, T) (SPEC §4).

Haldemann, Alibert, Mordasini & Benz (2020), A&A 643, A105. The native table
is (ρ, T)-indexed; ρ(P,T), grad_ad(P,T), s(P,T) are obtained by inverting
P(ρ,T) along T = const (AMENDMENTS A2) with the coexistence-plateau handling of
AMENDMENTS A8. Interpolation is bilinear in (log P, log T): per-column
interpolation against the strictly-increasing log P, then linear blend in
log T. log ρ and log P are interpolated in log space; s and grad_ad linearly.
"""

from __future__ import annotations

import numpy as np

from config import Config
from eos.base import EOS, EOSState
from eos.loader import AquaTables, load_aqua_tables


class AquaEOS(EOS):
    """AQUA H₂O equation of state, queried as a function of (P, T)."""

    name = "aqua"

    def __init__(self, table_path: str) -> None:
        self.table_path = str(table_path)
        self._t: AquaTables = load_aqua_tables(self.table_path)

    @classmethod
    def from_config(cls, cfg: Config) -> "AquaEOS":
        return cls(str(cfg.aqua_table_path))

    # -- core scalar interpolation --------------------------------------
    def _column_interp(self, log_p: float, j: int) -> tuple[float, float, float, float, float]:
        """(logρ, s, grad_ad, u, sound) at log P in T-column ``j``."""
        xp = self._t.logP_strict[:, j]
        log_rho = np.interp(log_p, xp, self._t.logrho_grid)
        s = np.interp(log_p, xp, self._t.entropy[:, j])
        grad = np.interp(log_p, xp, self._t.ad_grad[:, j])
        u = np.interp(log_p, xp, self._t.u[:, j])
        c = np.interp(log_p, xp, self._t.sound[:, j])
        return log_rho, s, grad, u, c

    def _state_scalar(self, pressure_pa: float, temperature_k: float) -> EOSState:
        log_p = np.log10(pressure_pa)
        log_t = np.log10(temperature_k)
        tg = self._t.logT_grid

        jt = int(np.clip(np.searchsorted(tg, log_t), 1, tg.size - 1))
        j0 = jt - 1
        wt = float(np.clip((log_t - tg[j0]) / (tg[jt] - tg[j0]), 0.0, 1.0))

        lr0, s0, g0, u0, c0 = self._column_interp(log_p, j0)
        lr1, s1, g1, u1, c1 = self._column_interp(log_p, jt)
        log_rho = (1.0 - wt) * lr0 + wt * lr1
        s = (1.0 - wt) * s0 + wt * s1
        grad = (1.0 - wt) * g0 + wt * g1
        u = (1.0 - wt) * u0 + wt * u1
        c = (1.0 - wt) * c0 + wt * c1

        jn = jt if wt > 0.5 else j0
        i = int(np.clip(np.searchsorted(self._t.logP_strict[:, jn], log_p), 0, self._t.n_rho - 1))
        phase = int(self._t.phase[i, jn])

        return EOSState(
            rho=float(10.0**log_rho),
            grad_ad=float(grad),
            entropy=float(s),
            internal_energy=float(u),
            sound_speed=float(c),
            phase=phase,
        )

    # -- public interface (scalar or array) -----------------------------
    def state(self, pressure_pa, temperature_k):
        """EOSState at (P,T). Scalar in -> EOSState; array in -> list in ravel order."""
        if np.isscalar(pressure_pa) and np.isscalar(temperature_k):
            return self._state_scalar(float(pressure_pa), float(temperature_k))
        p, t = np.broadcast_arrays(np.asarray(pressure_pa, float), np.asarray(temperature_k, float))
        return [self._state_scalar(float(pi), float(ti)) for pi, ti in zip(p.ravel(), t.ravel())]

    def _field(self, pressure_pa, temperature_k, attr: str):
        """Scalar or array-shaped access to one EOSState attribute."""
        if np.isscalar(pressure_pa) and np.isscalar(temperature_k):
            return getattr(self._state_scalar(float(pressure_pa), float(temperature_k)), attr)
        shape = np.broadcast(np.asarray(pressure_pa), np.asarray(temperature_k)).shape
        vals = [getattr(s, attr) for s in self.state(pressure_pa, temperature_k)]
        return np.array(vals).reshape(shape)

    def rho(self, pressure_pa, temperature_k):
        return self._field(pressure_pa, temperature_k, "rho")

    def grad_ad(self, pressure_pa, temperature_k):
        return self._field(pressure_pa, temperature_k, "grad_ad")

    def entropy(self, pressure_pa, temperature_k):
        return self._field(pressure_pa, temperature_k, "entropy")

    def internal_energy(self, pressure_pa, temperature_k):
        """Specific internal energy u(P,T), J kg⁻¹ (evolution/G6)."""
        return self._field(pressure_pa, temperature_k, "internal_energy")

    def sound_speed(self, pressure_pa, temperature_k):
        """Bulk sound speed c(P,T), m s⁻¹ (SPEC §4 sound-check helper)."""
        return self._field(pressure_pa, temperature_k, "sound_speed")

    # -- coverage --------------------------------------------------------
    @property
    def logP_range(self) -> tuple[float, float]:
        return float(self._t.logP_strict.min()), float(self._t.logP_strict.max())

    @property
    def T_range_k(self) -> tuple[float, float]:
        return float(10.0 ** self._t.logT_grid[0]), float(10.0 ** self._t.logT_grid[-1])
