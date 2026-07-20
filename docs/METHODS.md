# Methods reference

## Labels (hybrid production)
CP2K, revPBE0-D3 with ADMM; CUTOFF 600 Ry, REL_CUTOFF 60, NGRIDS 5;
DZVP-MOLOPT-SR-GTH; GTH-PBE-q6/q1; OT (DIIS, FULL_SINGLE_INVERSE) below
3000 K; diagonalisation + Fermi-Dirac smearing at 3000 K; EPS_SCF 1e-6.
E0s: isolated spin-polarised atoms, matching settings.

## Training
train_v4.xyz: 1193 configurations — 800 bulk pools (650-2250 K),
344 RCB fill (2500-3000 K), 49 above 3000 K. Fine-tune from
mace-mpa-0-medium; SWA; energy/forces keys REF_energy/REF_forces.
Reactive extension >3000 K descoped (hybrid SCF non-convergence); the
potential is applied only at <= 3000 K.

## Production MD
64 molecules, NVT (CSVR-equivalent Langevin in ASE drivers), 0.5 fs,
60 ps per state point, 27-point grid (see deck slide 4). Pressure =
virial + N kB T / V (ideal-gas kinetic term).

## Dipoles and spectra
MLWF/Berry-phase dipole labels; MACE dipole model (water_dipole_v2)
evaluated per frame; ACF via FFT (unbiased normalisation), mean-removed;
Hanning window, maxlag 8000; nfft 2^17; quantum correction
alpha ~ w(1-exp(-hbar w / kB T)) I(w); white-noise floor from the
3900-5500 cm^-1 window; 20 cm^-1 Gaussian smoothing on output kappa.

## Opacity
sigma -> alpha -> kappa via number density (post-2026-06-29 code);
Planck mean over full band; Rosseland over 100-5000 cm^-1 with a 1e-6
relative floor.

## Emissivity
Self-consistent dielectric construction (k fixed by alpha; damped
iteration eps2=2nk, eps1=eps_inf+KK[eps2]); n_inf via Lorentz-Lorenz with
alpha_e = 1.47 A^3 (interim); normal-incidence Fresnel; e = 1 - R;
Planck-weighted e_P. Diagnostics: iteration convergence, Im sqrt(eps) vs k
mismatch (band truncation), KK cutoff drift, implied sigma_DC.
