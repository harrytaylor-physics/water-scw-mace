"""
aqua_eos.py -- pressure P(T, rho) from the AQUA density-temperature table
(Haldemann et al. 2020; github.com/mnijh/AQUA).

On ARC, once only:
    git lfs install
    git clone https://github.com/mnijh/AQUA.git
    cd AQUA/Tables && unzip aqua_eos_rhot_v1_0.zip      # -> aqua_eos_rhot_v1_0.dat

Then:
    from aqua_eos import AquaEOS
    eos = AquaEOS("AQUA/Tables/aqua_eos_rhot_v1_0.dat")
    P_pa = eos.P(T_K=3000.0, rho_kgm3=108.8)            # Pascals

Table layout (from the AQUA Fortran module LOAD_TABLE_RhoT):
  - 21 header lines to skip
  - nRho=1501 (outer), nT=301 (inner); row order: for each rho, all T
  - columns: rho[kg/m3]  T[K]  P[Pa]  grad_ad  S  U  c_s  mu  x_ion  x_diss  phase
  - axes interpolated in log10(rho), log10(T)  (bilinear, matching AQUA)
"""
import numpy as np
from scipy.interpolate import RegularGridInterpolator

class AquaEOS:
    def __init__(self, dat_path, skip=21, nRho=1501, nT=301):
        raw = np.loadtxt(dat_path, skiprows=skip)
        if raw.shape[0] != nRho * nT:
            raise ValueError(f"AQUA rhoT table: expected {nRho*nT} rows, got {raw.shape[0]} "
                             f"(check skip={skip}, nRho={nRho}, nT={nT})")
        rho = raw[:, 0].reshape(nRho, nT)
        T   = raw[:, 1].reshape(nRho, nT)
        P   = raw[:, 2].reshape(nRho, nT)        # Pa
        self.logrho = np.log10(rho[:, 0])        # axis 0
        self.logT   = np.log10(T[0, :])          # axis 1
        self.logP   = np.log10(P)
        self._it = RegularGridInterpolator((self.logrho, self.logT), self.logP,
                                           bounds_error=False, fill_value=None)

    def P(self, T_K, rho_kgm3):
        T_K = np.atleast_1d(np.asarray(T_K, float))
        rho = np.atleast_1d(np.asarray(rho_kgm3, float))
        pts = np.column_stack([np.log10(rho), np.log10(T_K)])
        out = 10.0 ** self._it(pts)              # Pa
        return out[0] if out.size == 1 else out
