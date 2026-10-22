# fp_opacity_v5/ — first-principles H2O opacity spectra (v5)

52 absorption spectra of pure H2O from the production MD (revPBE0-D3 labels,
`models/water_scw_v4.model`), one per state point, on one grid. The dipoles
along the trajectories came from the MLWF-trained dipole model
`water_dipole_v3` (`production_campaign/arc_scripts/htc_opac*.sh` ->
`scw/mace_opacity.py`). Its binary is in this repository as
`subneptune/monomer_test/water_dipole_v3.model` (LFS); `models/` holds the
earlier `water_dipole_v2`. Regenerated from the production dipole series by
`subneptune/opacity/dipole_spectrum.py` via `subneptune/tools/v5_regenerate.py`
(subneptune A83). They supersede everything `production_campaign/scw/opacity_pipeline.py`
delivered; see [`production_campaign/README.md`](../../production_campaign/README.md).

## Grid

| | |
|---|---|
| state points | 52 |
| temperature | 650, 700, 800, 1000, 1500, 2000, 2250, 2500, 2750, 3000 K |
| density | 24 values, 10 to 2228.71 kg m^-3 |
| cells | 64 molecules at 39 points, 216 at 13 (field `N`) |
| MD | 0.5 fs step, 30 ps NVE production after 10 ps Langevin equilibration |
| list | `grid.json`, `[{"T": K, "rho": kg m^-3}, ...]` |

File names are `T{T}_rho{rho}_opacity.npz` (T in K, rho in kg m^-3).

## Contents of each `.npz`

| key | units | meaning |
|---|---|---|
| `wn` | cm^-1 | wavenumber grid, 6000 points, ~1 to 6000 cm^-1, 1 cm^-1 step |
| `sigma` | cm^2 per molecule | absorption cross-section, Gaussian-smoothed (sd = `smooth_cm`) |
| `sigma_raw` | cm^2 per molecule | the same, unsmoothed |
| `kappa` | cm^2 g^-1 | mass opacity |
| `alpha` | cm^-1 | absorption coefficient, `kappa * rho` |
| `planck`, `rosseland` | cm^2 g^-1 | mean opacities of this spectrum: Planck over the whole `wn` grid; Rosseland over 100-5000 cm^-1 with a floor of 1e-6 of the band maximum |
| `kappa_P_sigma`, `kappa_R_sigma` | cm^2 g^-1 | half the difference between the means of the two 15 ps halves of the record |
| `T`, `rho_kg_m3`, `N`, `dt_fs`, `smooth_cm` | K, kg m^-3, -, fs, cm^-1 | state point and processing |
| `units_fixed` | | **version marker**, `cm2_per_g:A83` |
| `version` | | `v5:tapered-periodogram:sd10` |
| `estimator` | | the estimator, in words |

The planetary model takes its own means over the 100-4000 cm^-1 band shared by
every source it compares (`subneptune/opacity/`), not the stored `planck` /
`rosseland`, which are the backend's own and kept for comparison with v3/v4.

`subneptune/opacity/loaders.py` refuses a spectrum whose `units_fixed` does not
match the value the registry expects for this directory (`cm2_per_g:A83`), so a
v3 or v4 file cannot be read as v5 by accident.

## Estimator

`sigma(w) = beta w^2 I2(w) / (6 eps0 c N)`, with `I2` the two-sided power
spectrum of the cell dipole after a Hann taper on the **time series**, normalised
by `sum(w^2)`; the harmonic quantum correction applied once; a continuous
`omega^2` prefactor; no clip at zero; 10 cm^-1 Gaussian smoothing.

Four checks; the first three are pinned by gates in `subneptune/validation/`:

| check | bound | measured | gate |
|---|---|---|---|
| synthetic classical harmonic oscillator vs the textbook band strength | 3% | 2.0% (1.0201) | `test_a79_normalisation.py`, `test_a80_v4_gates.py` |
| sum rule: band integral vs `beta <\|Mdot_band\|^2>` | 1% | 0.09% | `test_a83_v5_estimator.py` |
| leakage floor below the transparency-window minimum, all 52 cells | 100x | 476x at the worst, (3000 K, 108.796) | `test_a83_v5_estimator.py` |
| dilute-limit rotational band against the HITEMP line list | — | median FP/line list 1.043 (bend 1.483, stretch 2.589) | A83; the monomer test expects 0.97 (A79) |

## Also here

- `grid.json` — the 52 state points.
- `crossover_hitemp.json` — per-point HITEMP band ratios (via RADIS) that the
  planetary model's comparison figures read.

## Not here

The dipole series these were computed from are raw production output and are
not in the repository. Their SHA-256 fingerprints are in
`subneptune/ed_data/arc_dipoles_sha256.json` (67 series; 15 state points have an
older and a 0.5 fs run) — `tools/v5_regenerate.py` matches each spectrum to its
series by reproduction, not by name.
