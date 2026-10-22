"""opacity_pipeline.py -- validated absolute-opacity backend, source-agnostic."""
import numpy as np
from scipy.ndimage import gaussian_filter1d

_trapz = getattr(np, "trapz", None) or np.trapezoid

kB     = 1.380649e-23
h      = 6.62607015e-34
hbar   = 1.054571817e-34
c      = 2.99792458e8
c_cm_s = 2.99792458e10
eps0   = 8.8541878128e-12
debye_to_Cm = 3.33564e-30
NA     = 6.02214076e23
EANG_TO_DEBYE = 1.602176634e-29 / debye_to_Cm
DEFAULT_MOLAR_MASS_G_MOL = 18.01528


def autocorr_fft(x):
    n = len(x)
    nfft = 2 ** int(np.ceil(np.log2(2 * n)))
    fx = np.fft.rfft(x, n=nfft)
    ac = np.fft.irfft(fx * np.conjugate(fx), n=nfft)[:n]
    norm = np.arange(n, 0, -1, dtype=float)
    return ac / norm


def sigma_to_alpha_kappa(sigma_cm2, density_g_cm3, molar_mass_g_mol=DEFAULT_MOLAR_MASS_G_MOL):
    n_m3 = (density_g_cm3 * 1.0e6 / molar_mass_g_mol) * NA
    sigma_m2 = sigma_cm2 * 1.0e-4
    alpha_cm_inv = (n_m3 * sigma_m2) / 100.0
    return alpha_cm_inv, alpha_cm_inv / density_g_cm3


def planck_and_rosseland_from_kappa(wn_cm, kappa_cm2_g, T_K,
                                    ross_wn_min=100.0, ross_wn_max=5000.0, ross_floor_frac=1.0e-6):
    nu_hz = wn_cm * c_cm_s
    Bnu = (2.0 * h * nu_hz**3 / c**2) / np.maximum(np.exp(h * nu_hz / (kB * T_K)) - 1.0, 1.0e-300)
    xT = h * nu_hz / (kB * T_K)
    exp_xT = np.exp(np.clip(xT, None, 700))
    dB_dT = (2.0 * h * nu_hz**3 / c**2) * (exp_xT / np.maximum((exp_xT - 1.0) ** 2, 1.0e-300)) * (xT / T_K)
    kappa_planck = _trapz(kappa_cm2_g * Bnu, nu_hz) / _trapz(Bnu, nu_hz)
    rm = (wn_cm >= ross_wn_min) & (wn_cm <= ross_wn_max)
    nu_r = nu_hz[rm]; dB_r = dB_dT[rm]; k_r = kappa_cm2_g[rm]
    k_safe = np.maximum(k_r, ross_floor_frac * np.nanmax(k_r))
    kappa_ross = _trapz(dB_r, nu_r) / _trapz(dB_r / k_safe, nu_r)
    return kappa_planck, kappa_ross


def spectrum_from_dipoles(M, T_K, dt_fs, n_molecules, density_g_cm3,
                          units="debye", molar_mass_g_mol=DEFAULT_MOLAR_MASS_G_MOL,
                          max_wn_cm=6000.0, smooth_sigma_cm=20.0, wn_step_out=1.0,
                          ross_wn_min=100.0, ross_wn_max=5000.0, ross_floor_frac=1.0e-6, n_blocks=1):
    M = np.asarray(M, dtype=float)
    if units == "eAng":
        M = M * EANG_TO_DEBYE
    elif units != "debye":
        raise ValueError("units must be 'debye' or 'eAng'")
    beta = 1.0 / (kB * T_K)
    M_cm = M * debye_to_Cm
    n_frames = len(M_cm); dt_s = dt_fs * 1e-15

    def seg_smm(seg):
        L = len(seg); s = seg - seg.mean(0)
        Cmm = autocorr_fft(s[:, 0]) + autocorr_fft(s[:, 1]) + autocorr_fft(s[:, 2])
        return dt_s * np.real(np.fft.rfft(Cmm * np.hanning(L))), np.fft.rfftfreq(L, d=dt_s)

    if n_blocks <= 1:
        Smm, freq_hz = seg_smm(M_cm)
    else:
        L = n_frames // n_blocks; acc = None
        for i in range(n_blocks):
            smm, freq_hz = seg_smm(M_cm[i * L:(i + 1) * L])
            acc = smm if acc is None else acc + smm
        Smm = acc / n_blocks

    omega = 2.0 * np.pi * freq_hz; wn_cm = freq_hz / c_cm_s
    mask = (wn_cm > 0.0) & (wn_cm <= max_wn_cm)
    wn_cm = wn_cm[mask]; omega = omega[mask]; Smm = np.clip(Smm[mask], 0.0, None)
    x = beta * hbar * omega
    qcf = np.ones_like(x); mm = x > 1e-12; qcf[mm] = x[mm] / (1.0 - np.exp(-x[mm]))
    sigma_cm2_raw = (beta * omega**2 / (3.0 * eps0 * c * n_molecules)) * (Smm * qcf) * 1.0e4
    wn_out = np.arange(max(wn_cm.min(), 0.0), max_wn_cm + wn_step_out, wn_step_out)
    sig_raw = np.nan_to_num(np.interp(wn_out, wn_cm, sigma_cm2_raw, left=np.nan, right=np.nan), nan=0.0)
    sig_sm = gaussian_filter1d(sig_raw, sigma=smooth_sigma_cm / wn_step_out, mode="nearest")
    alpha, kappa = sigma_to_alpha_kappa(sig_sm, density_g_cm3, molar_mass_g_mol)
    kP, kR = planck_and_rosseland_from_kappa(wn_out, kappa, T_K, ross_wn_min, ross_wn_max, ross_floor_frac)
    return {"wn_cm-1": wn_out, "sigma_cm2_per_molecule_raw": sig_raw, "sigma_cm2_per_molecule": sig_sm,
            "alpha_cm-1": alpha, "kappa_cm2_g": kappa, "planck_mean_cm2_g": kP,
            "rosseland_mean_cm2_g": kR, "n_frames": n_frames}
