"""
exomolop_opacity.py -- POKAZATEL discrete-line opacity from a pre-computed ExoMolOP
cross-section file (fixed reference broadening => read at a low reference pressure).

ExoMolOP (Chubb et al. 2021) ships ready cross-sections sigma(nu; T, P). We use it at a
fixed reference pressure P_REF (lowest available, ~ discrete lines) and interpolate in T.

Supported here: TauREx HDF5 format (single .h5 with keys: 'bin_edges' or 'wno',
'xsecarr' [P,T,nu] in cm^2/molecule, 't', 'p'). petitRADTRANS folder format differs.

>>> VERIFY ON FIRST RUN <<<
  - Tell me which ExoMolOP format you downloaded (TauREx .h5 / petitRADTRANS / NEMESIS);
    this reader targets TauREx .h5. Confirm xsecarr units (cm^2/molecule) and axis order.
"""
import numpy as np, h5py
NA=6.02214076e23; M=18.015

def exomolop_kappa(h5_path, T, p_ref_bar=None):
    """Return (wn[cm^-1], kappa[cm^2/g]) at temperature T, fixed reference pressure."""
    with h5py.File(h5_path,"r") as f:
        keys=set(f.keys())
        wn = np.array(f["wno"]) if "wno" in keys else 0.5*(np.array(f["bin_edges"])[:-1]+np.array(f["bin_edges"])[1:])
        xs = np.array(f["xsecarr"])      # [P, T, nu] cm^2/molecule (CONFIRM)
        Tg = np.array(f["t"]); Pg = np.array(f["p"])
    ip = 0 if p_ref_bar is None else int(np.argmin(np.abs(Pg-p_ref_bar)))   # lowest P = discrete lines
    # interpolate in T at fixed pressure index
    sigma = np.array([np.interp(T, Tg, xs[ip,:,k]) for k in range(xs.shape[2])])
    return wn, sigma * NA / M            # cm^2/g
