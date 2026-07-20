# Testing — outstanding, conditional fixes, and conditions of use

## Testing still to do
- Parity benchmark against a true held-out split (blocked on holdout file).
- Dipole-model tests: |M| distribution regression, ACF C(dt)/C(0) >= 0.997
  smoothness check as an automated test, spectrum band-position regression
  (libration/bend/stretch windows).
- EOS spot-checks at 2-3 more state points against IAPWS-95/AQUA.
- Emissivity: cross-check one point against a direct thin-film transfer
  calculation; KK band-truncation sensitivity beyond the 0.8-cutoff proxy.
- A "hot cell" stress test: short NVE at 3000 K / 2228 kg/m3, energy
  conservation and no dissociative runaway.

## Assuming the results are correct — technical fixes and checks
- Regenerate stale npz (TODO #1) and grep the manuscript for any absolute
  kappa quoted from pre-fix files (figures were corrected; text should be
  audited once).
- The _05fs production set is canonical; archive or clearly mark the
  non-05fs originals to prevent accidental use.
- alpha_e -> per-point DFT polarizabilities in the emissivity baseline.
- Record exact mace-torch/torch versions used for v4 (reproducibility).

## Assuming it is fine — how to use it usefully
- **Valid envelope:** T <= 3000 K, rho within the training pools
  (~100-2230 kg/m^3 along the production grid). Do not extrapolate above
  3000 K (reactive extension was descoped; hybrid labels did not converge
  there) or to ice phases (not trained).
- **Best uses:** RCB-region opacities for sub-Neptune cooling/structure
  models; interpolation of kappa_P / kappa_R between grid points (log-log
  in rho is smooth); IR emissivity for radiative boundary conditions;
  EOS points within the envelope.
- **Conditions of interest for new runs:** the low-density (100-300 kg/m3)
  branch where line lists and MD continuum overlap (crossover mapping);
  the 2500-3000 K dense branch for dissociation/conductivity onset; any
  TOI-270 d / K2-18 b adiabat traverses through the envelope.
- **Rule of thumb:** ratios from this pipeline are robust (scale factors
  cancel); absolute values require the post-fix data products.
