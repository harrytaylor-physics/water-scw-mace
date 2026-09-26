# First-principles steam opacity — release tables

Rosseland and Planck mean opacity of pure H2O from first-principles molecular
dynamics, on a regular (log R, log T) grid.

| | |
|---|---|
| Units | kappa in cm^2 g^-1, tabulated as log10; T in K; rho in kg m^-3 (MESA's log R uses g cm^-3) |
| Band | 100-4000 cm^-1, the widest band every source in the study shares |
| log R | 6.00 to 9.50, step 0.25 dex, R = rho / (T/10^6 K)^3 with rho in g cm^-3 |
| log T | 2.80 to 3.50, step 0.05 dex (631-3162 K) |
| Computed points | 48 states, T 650-3000 K, rho 10-2229 kg m^-3 |
| Coverage | grid 0.25 dex in log R, 0.05 dex in log T; 93 of 225 cells inside the 48-point hull |
| Outside the hull | `-99` in the MESA files, blank in the CSV — never extrapolated |

## Files

- `subneptune_kappa_R_mesa.txt`, `subneptune_kappa_P_mesa.txt` — the plain-text
  format MESA's `kap` module reads (`kap/private/load_kap.f90`,
  `form = 1`, `version = 37`). **Format, not a drop-in table**: a steam
  envelope sits at log R 6.0-9.5, where MESA's own tables run -8 to +1.
- `subneptune_kappa_grid.csv` — the same numbers with true blanks, plus T and rho.
- `subneptune_kappa_fit.json` — the two-regime power-law fit and this metadata.

Values are the study's own opacity table sampled onto the regular grid: linear
on the Delaunay triangulation in (T, log10 rho) over log10 kappa.

## Two caveats, one line each

- **Interaction-induced absorption.** These are first-principles molecular-dynamics opacities, so collision-induced and many-body absorption is present by construction and is not separable into a line list plus a pair continuum.
- **Dipole model.** The spectra come from a machine-learned dipole surface whose symmetric-stretch derivative is 2.8x the revPBE0 value at 0.957 A and 4.0x at 1.00 A, so the 3200-4000 cm^-1 stretch band is too strong and kappa_P is affected more than kappa_R.
