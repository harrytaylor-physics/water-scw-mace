"""
hitemp_opacity.py -- HITEMP H2O discrete-line kappa(nu) [cm^2/g] via RADIS, fixed reference broadening.

Reads the local RADIS HITEMP cache at $DATA/radisdb/hitemp (already built).
MUST be imported/run on a COMPUTE NODE -- loading HITEMP lines needs more memory
than the login node allows.

Fixed reference broadening: lines computed at a single low reference pressure P_REF_BAR
so we get the discrete-transition opacity, not GPa-broadened wings. Line strengths still
vary with T (physical). The intrinsic mass opacity kappa = abscoeff / rho_ref is independent
of P_REF in the linear (optically-thin-line) regime, since abscoeff ~ rho_ref.
"""
import os, numpy as np
import radis
radis.config["DEFAULT_DOWNLOAD_PATH"] = os.path.join(os.environ["DATA"], "radisdb")
from radis import calc_spectrum

NA = 6.02214076e23; M = 18.015; KB = 1.380649e-23
P_REF_BAR = 1.0

def hitemp_kappa(T, numin=10.0, numax=11000.0, wstep=0.01, p_ref_bar=P_REF_BAR):
    """Return (wn[cm^-1], kappa[cm^2/g]) for H2O at temperature T, fixed broadening p_ref_bar."""
    s = calc_spectrum(numin, numax, molecule="H2O", isotope="1,2,3",
                      Tgas=float(T), pressure=float(p_ref_bar), mole_fraction=1.0,
                      databank="hitemp", wstep=wstep, verbose=False)
    w, ac = s.get("abscoeff", wunit="cm-1")          # cm^-1 at (T, p_ref)
    p_pa = p_ref_bar * 1e5
    n_cm3 = (p_pa / (KB * T)) * 1e-6                  # H2O number density [cm^-3]
    rho_g_cm3 = n_cm3 * M / NA                        # mass density [g/cm^3]
    return np.asarray(w, float), np.asarray(ac, float) / rho_g_cm3   # cm^2/g


class HitempCalculator:
    """Load the HITEMP databank ONCE, then compute kappa(nu;T) per cell (no reload).
    Cuts a 27-cell run from ~1 h of repeated loads to ~15 min. Returns cm^2/g."""
    def __init__(self, numin=10.0, numax=11000.0, wstep=0.01, p_ref_bar=P_REF_BAR):
        from radis import SpectrumFactory
        self.p_ref_bar = p_ref_bar
        self.sf = SpectrumFactory(wavenum_min=numin, wavenum_max=numax,
                                  molecule="H2O", isotope="1,2,3",
                                  wstep=wstep, pressure=p_ref_bar, mole_fraction=1.0,
                                  verbose=False)
        self.sf.fetch_databank("hitemp")          # one-time load of the cached lines
    def kappa(self, T):
        s = self.sf.eq_spectrum(Tgas=float(T))
        w, ac = s.get("abscoeff", wunit="cm-1")
        n_cm3 = (self.p_ref_bar * 1e5 / (KB * T)) * 1e-6
        rho_g_cm3 = n_cm3 * M / NA
        return np.asarray(w, float), np.asarray(ac, float) / rho_g_cm3
