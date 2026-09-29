"""
mtckd_continuum.py -- MT_CKD water-vapour SELF continuum kappa(nu) [cm^2/g] for pure water.

Coefficients (MT_CKD 4.3, AER) are bundled in mtckd_ref.npz (extracted from
absco-ref_wv-mt-ckd.nc); no netCDF dependency needed. Implementation reproduces the
AER reference output to ~0.1% (validated against run_example/mt_ckd_h2o_output.nc).

Formula (per the AER mt_ckd_h2o_module.f90):
    sigma_self(nu) = self_absco_ref(nu) * (296/T)^self_texp(nu) * rho_rat * radterm(nu,T)
    radterm = nu*(1-e^-x)/(1+e^-x),  x = nu*1.4387752/T          [cm^2/molecule]
For PURE water we (a) drop the foreign term (no foreign gas) and (b) use the ACTUAL
water number-density ratio rho_rat = n_water/n_ref instead of the ideal-gas pressure
ratio the AER code uses -- correct for the dense, non-ideal regime. n_ref is the
reference number density (1013 mbar, 296 K). The self continuum then scales ~rho^2
in absorption (kappa ~ rho), as required for a collision-induced continuum.

NOTE: at our densities rho_rat ~ 1e2-1e3, far beyond MT_CKD's atmospheric calibration --
this is the *best empirical continuum extrapolated*, which is exactly what the breakdown
front measures MACE against.
"""
import numpy as np
NA=6.02214076e23; M=18.015; KB=1.380649e-23; RADCN2=1.4387752   # hc/k [cm K]

class MTCKD:
    def __init__(self, npz="mtckd_ref.npz"):
        d=np.load(npz)
        self.wn=d["wavenumbers"]; self.self_ref=d["self_absco_ref"]; self.texp=d["self_texp"]
        self.P_ref_mbar=float(d["ref_press_mbar"]); self.T_ref=float(d["ref_temp_K"])
        self.n_ref_cm3=(self.P_ref_mbar*100.0)/(KB*self.T_ref)*1e-6   # molecule/cm^3 at ref

    @staticmethod
    def _radterm(wn, T):
        wn=np.asarray(wn,float); x=wn*RADCN2/T
        return np.where(x<=0.01, 0.5*x*wn, wn*(1.0-np.exp(-x))/(1.0+np.exp(-x)))

    def kappa_self(self, wn, T, rho_kg_m3):
        """Pure-water self-continuum kappa(nu) [cm^2/g] on grid wn[cm^-1]."""
        wn=np.asarray(wn,float)
        n_w=rho_kg_m3*NA/(M*1e-3)*1e-6                    # actual water number density /cm^3
        rho_rat=n_w/self.n_ref_cm3
        sc=np.interp(wn,self.wn,self.self_ref); tx=np.interp(wn,self.wn,self.texp)
        sigma=sc*(self.T_ref/T)**tx*rho_rat*self._radterm(wn,T)   # cm^2/molecule
        return sigma*NA/M                                          # cm^2/g

    def _sigma_aer(self, wn, T, p_mbar, vmr):
        """AER convention (ideal-gas pressure rho_rat, with vmr) -- for validation only."""
        wn=np.asarray(wn,float); rho_rat=(p_mbar/self.P_ref_mbar)*(self.T_ref/T)
        sc=np.interp(wn,self.wn,self.self_ref); tx=np.interp(wn,self.wn,self.texp)
        return sc*(self.T_ref/T)**tx*vmr*rho_rat*self._radterm(wn,T)
