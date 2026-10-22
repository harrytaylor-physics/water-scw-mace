"""
opacity_means.py -- Planck and Rosseland mean opacities from a spectral kappa(nu).
Same definitions used for the MACE/AIMD pipeline, so all surfaces are consistent.
    wn   : wavenumber grid [cm^-1]   (ascending)
    kap  : kappa(nu)       [cm^2/g]
    T    : temperature     [K]
Constants cancel in the ratios; we keep them for clarity.
"""
import numpy as np
_trapz=getattr(__import__('numpy'),'trapz',None) or __import__('numpy').trapezoid
_H=6.62607015e-34; _C=2.99792458e10; _KB=1.380649e-23   # c in cm/s

def _x(wn,T): return _H*_C*np.asarray(wn,float)/(_KB*T)

def planck_mean(wn, kap, T):
    wn=np.asarray(wn,float); kap=np.asarray(kap,float)
    B=wn**3/np.expm1(_x(wn,T))
    return _trapz(kap*B,wn)/_trapz(B,wn)

def rosseland_mean(wn, kap, T):
    wn=np.asarray(wn,float); kap=np.maximum(np.asarray(kap,float),1e-30)
    x=_x(wn,T); w=wn**4*np.exp(x)/np.expm1(x)**2          # dB/dT shape
    return _trapz(w,wn)/_trapz(w/kap,wn)
