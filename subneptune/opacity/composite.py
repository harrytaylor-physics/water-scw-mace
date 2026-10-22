"""Composite opacity: line-list (dilute) + FP (dense), log-blended (A11, A12).

FP is a condensed-phase method, valid at high density; line lists are the
physically correct opacity in the dilute regime. The composite uses the
line-list source below a density window and FP above it, blending log κ across
the window with a smooth (C¹) weight:

    w(ρ) = smoothstep( (log ρ − log ρ_lo) / (log ρ_hi − log ρ_lo) )
    log κ = (1 − w) · log κ_linelist + w · log κ_FP

Below ρ_lo: pure line-list; above ρ_hi: pure FP. The seam (geometric-mean
centre of the window) and its width are sensitivity-scan parameters; the
FP/line-list κ ratio at the seam is a quantified systematic (logged per run).
The default window is ρ = 100–300 kg/m³ so FP is used only above its real data
floor (~109 kg/m³); ρ = 30–100 is a documented sensitivity case (A15).

Two prescriptions (``mode``, A35): ``"linear"`` is the blend above (the default,
unchanged); ``"sharp"`` switches to FP discontinuously at ρ_lo with no blending,
i.e. w(ρ) = 1[ρ ≥ ρ_lo], for which the window may be zero-width. The composite's
surviving RCB sits *inside* the default window, so its location depends on the
prescription — scanned by ``figures/plot_s1_blend_sensitivity.py``.
"""

from __future__ import annotations

import numpy as np

from opacity.source import OpacitySource


def _smoothstep(x: np.ndarray) -> np.ndarray:
    x = np.clip(x, 0.0, 1.0)
    return x * x * (3.0 - 2.0 * x)


class CompositeOpacitySource(OpacitySource):
    """Density-blended line-list (dilute) + FP (dense) opacity."""

    def __init__(
        self,
        dilute: OpacitySource,
        dense: OpacitySource,
        blend_rho_si: tuple[float, float] = (100.0, 300.0),
        name: str = "composite",
        mode: str = "linear",
    ) -> None:
        if mode not in ("linear", "sharp"):
            raise ValueError(f"blend mode must be linear|sharp, got {mode!r}")
        lo, hi = blend_rho_si
        if not 0.0 < lo <= hi or (mode == "linear" and lo == hi):
            raise ValueError(
                f"blend window must satisfy 0 < lo < hi (0 < lo <= hi when "
                f"mode='sharp'), got {blend_rho_si} with mode={mode!r}"
            )
        self.name = name
        self.mode = mode
        self.dilute = dilute  # line-list (HITEMP/POKAZATEL)
        self.dense = dense  # FP
        self.blend_rho_si = (float(lo), float(hi))
        # Seam = where the FP weight passes 1/2: the window's geometric-mean
        # centre when blending, the switch density itself when sharp.
        self.seam_rho_si = float(lo) if mode == "sharp" else float(np.sqrt(lo * hi))
        # Both sub-sources are built on the same band (A21); expose it.
        self.band_cm = getattr(dense, "band_cm", getattr(dilute, "band_cm", None))
        self._log_lo = np.log10(lo)
        self._log_hi = np.log10(hi)
        self._total = 0
        self._clamped = 0
        self._flat_rho = 0

    def _weight(self, rho_si) -> np.ndarray:
        """FP weight w(ρ): 0 (pure line-list) → 1 (pure FP) across the window.

        ``mode='sharp'``: a step at ρ_lo (no blending), which also keeps a
        zero-width window well defined (A35).
        """
        logr = np.log10(np.asarray(rho_si, float))
        if self.mode == "sharp":
            return (logr >= self._log_lo).astype(float)
        return _smoothstep((logr - self._log_lo) / (self._log_hi - self._log_lo))

    def _blend(self, which: str, temperature_k, rho_si):
        scalar = np.isscalar(temperature_k) and np.isscalar(rho_si)
        t = np.atleast_1d(np.asarray(temperature_k, float))
        r = np.atleast_1d(np.asarray(rho_si, float))
        t, r = np.broadcast_arrays(t, r)
        w = self._weight(r)

        k_dilute = getattr(self.dilute, which)(t, r)
        k_dense = getattr(self.dense, which)(t, r)
        log_k = (1.0 - w) * np.log10(k_dilute) + w * np.log10(k_dense)
        out = 10.0**log_k

        # Reliability: attribute each point to its dominant source's clamp state.
        use_dense = w >= 0.5
        clamped = np.where(
            use_dense,
            self.dense.clamped_mask(t, r),
            self.dilute.clamped_mask(t, r),
        )
        self._total += clamped.size
        self._clamped += int(clamped.sum())
        # Dilute-flat (A26): dilute-dominated points below the line list's ρ floor.
        floor = getattr(self.dilute, "rho_floor_si", None)
        if floor is not None:
            self._flat_rho += int(np.sum((~use_dense) & (r.ravel() < floor)))

        return float(out[0]) if scalar else out

    def kappa_R(self, temperature_k, rho_si):
        return self._blend("kappa_R", temperature_k, rho_si)

    def kappa_P(self, temperature_k, rho_si):
        return self._blend("kappa_P", temperature_k, rho_si)

    def clamped_mask(self, temperature_k, rho_si):
        """Per-point clamp mask, attributed to the dominant sub-source (A26)."""
        t = np.atleast_1d(np.asarray(temperature_k, float))
        r = np.atleast_1d(np.asarray(rho_si, float))
        t, r = np.broadcast_arrays(t, r)
        use_dense = self._weight(r) >= 0.5
        return np.where(use_dense, self.dense.clamped_mask(t, r), self.dilute.clamped_mask(t, r))

    @property
    def clamp_fraction(self) -> float:
        return self._clamped / self._total if self._total else 0.0

    @property
    def flat_rho_fraction(self) -> float:
        return self._flat_rho / self._total if self._total else 0.0

    def reset_clamp_counter(self) -> None:
        self._total = 0
        self._clamped = 0
        self._flat_rho = 0
        self.dilute.reset_clamp_counter()
        self.dense.reset_clamp_counter()

    def seam_ratio(self, temperatures_k) -> dict:
        """FP / line-list κ ratio at the seam density, per temperature (A12).

        Returns the quantified systematic to log per run: how much FP and the
        line list disagree where the composite switches between them.
        """
        t = np.atleast_1d(np.asarray(temperatures_k, float))
        rho = np.full_like(t, self.seam_rho_si)
        kr_fp = np.atleast_1d(self.dense.kappa_R(t, rho))
        kr_ll = np.atleast_1d(self.dilute.kappa_R(t, rho))
        kp_fp = np.atleast_1d(self.dense.kappa_P(t, rho))
        kp_ll = np.atleast_1d(self.dilute.kappa_P(t, rho))
        return {
            "seam_rho_si": self.seam_rho_si,
            "temperatures_k": t.tolist(),
            "kappa_R_fp_over_ll": (kr_fp / kr_ll).tolist(),
            "kappa_P_fp_over_ll": (kp_fp / kp_ll).tolist(),
        }

    def __repr__(self) -> str:
        lo, hi = self.blend_rho_si
        return (
            f"CompositeOpacitySource(dilute={self.dilute.name!r}, dense={self.dense.name!r}, "
            f"blend_rho=[{lo:g}, {hi:g}] kg/m^3, mode={self.mode!r}, "
            f"seam={self.seam_rho_si:g})"
        )
