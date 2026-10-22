# Methods reference

> ## Corrections (A90, 2026-10-01)
>
> This file was delivered with the ARC campaign. It is corrected rather than
> rewritten (subneptune A90), and what it said as delivered is quoted here so
> the delivery record survives.
> Each correction is against the delivered artefacts, which win over prose.
>
> | delivered text | corrected to | evidence |
> |---|---|---|
> | "60 ps per state point" | **30 ps** | every production dipole series (`subneptune/ed_data/arc_dipoles/*/dipoles.npy`, not committed; fingerprinted in `subneptune/ed_data/arc_dipoles_sha256.json`) is 60,000 frames at 0.5 fs |
> | "NVT (CSVR-equivalent Langevin in ASE drivers)" | **Langevin equilibration, then NVE production** | `production_campaign/scw/run_prod_md_05fs.py` line 48 runs Langevin for 20,000 steps, line 50 switches to `VelocityVerlet` at constant energy |
> | "27-point grid (see deck slide 4)" | **52 state points** | `data/fp_opacity_v5/` holds 52 `T*_opacity.npz`; subneptune treats them as ONE grid (`subneptune/SPEC.md`, and the `save_pdf` vocabulary gate) |
> | "64 molecules" | **64 or 216** | the `N` field of the 52 spectra takes both values |
>
> **The "27" is almost certainly the POKAZATEL coverage, not the MD grid** —
> `pokazatel_opacity/` holds exactly 27 spectra, and all 27 fall on FP state
> points (T = 1000-3000 K, rho = 108.796-2228.71). That identification is an
> inference from the counts, not something the pull states.
>
> **Not corrected, but superseded:** the "Dipoles and spectra" section below
> describes the v3 estimator — ACF windowing with a Hanning window at maxlag
> 8000, a white-noise floor taken from 3900-5500 cm^-1, and 20 cm^-1 smoothing.
> A81 showed that floor was spectral leakage rather than noise and must never be
> subtracted; A83 replaced the estimator with a Hann taper on the **time series**
> under `sum(w^2)` normalisation, the harmonic correction applied once, a
> continuous omega^2 prefactor and **10 cm^-1** smoothing. `subneptune/opacity/dipole_spectrum.py`
> is authoritative for v5.

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
64 or 216 molecules; Langevin equilibration (CSVR-equivalent in ASE
drivers) for 10 ps, then NVE production with VelocityVerlet; 0.5 fs,
30 ps per state point, 52-point grid. Pressure = virial + N kB T / V
(ideal-gas kinetic term). POKAZATEL spectra exist at 27 of the 52.
[Corrected A90 — see the block at the top of this file for the text as
delivered and the evidence for each change.]

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
