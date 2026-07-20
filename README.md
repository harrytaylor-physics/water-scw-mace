# water-scw-mace

Machine-learned interatomic potential and opacity pipeline for warm dense
water (supercritical water, "SCW") in sub-Neptune interiors.

## What it is

`water_scw_v4` is a MACE potential fine-tuned from the `mace-mpa-0-medium`
foundation model on hybrid-functional (revPBE0-D3) labels computed in CP2K,
covering 650-3000 K and ~100-2230 kg/m^3. Alongside it: `water_dipole_v2`,
a MACE dipole model trained on MLWF/Berry-phase dipole labels, used to
compute IR absorption spectra from dipole autocorrelation along MD
trajectories, from which absolute opacities (Planck/Rosseland means) and
normal emissivities are derived.

## What it does

- MD of warm dense water at hybrid-DFT accuracy at classical-potential cost
- Equation-of-state points (pressure includes the ideal-gas kinetic term)
- IR absorption spectra and absolute opacity kappa(nu; T, rho)
- Opacity surface over the 27-point production grid (1000-3000 K)
- Normal spectral emissivity via a self-consistent dielectric construction

## Status

| Component | Status |
|---|---|
| water_scw_v4 potential | validated (parity, EOS, finite-size); production |
| water_dipole_v2 | production; behind the published opacity surface |
| Opacity surface (27 points) | complete; **absolute scale note below** |
| Emissivity (27 points) | complete (interim n_inf baseline) |
| LDA functional-sensitivity campaign | in progress (steps A-C running) |

**Known issue (fixed in code, stale in data):** a x1000 number-density unit
error in `opacity_pipeline.sigma_to_alpha_kappa` inflated all absolute
alpha/kappa in `opacity_v3.npz` files generated before 2026-06-29. The code
is fixed; published figures were corrected; the archived per-point npz and
`opacity_surface.npz` still carry the factor. Divide by 1000 or regenerate.
Ratios (MACE/AIMD, functional comparisons) were never affected.

## Layout

```
models/        trained .model files (git-lfs) + provenance
benchmarks/    rerunnable benchmarks with reference values
tests/         pytest suite (CPU-safe, skips if model absent)
convergence/   convergence evidence: cutoff, finite-size, timestep
studies/       satellite studies (e.g. dispersion/functional -> opacity)
docs/          results deck (deck.html), schematic, methods
```

## Quick start

```bash
pip install mace-torch ase numpy
python - << 'PY'
from mace.calculators import MACECalculator
from ase.io import read
calc = MACECalculator(model_paths="models/water_scw_v4.model", device="cuda")
atoms = read("benchmarks/reference/water64_sample.xyz")
atoms.calc = calc
print("E =", atoms.get_potential_energy(), "eV")
PY
```

For MD, use ASE Langevin/NVT at 0.5 fs; apply the potential only at
T <= 3000 K and within the trained density envelope (see TESTING.md,
"conditions of use").

## Methods (short)

Labels: CP2K revPBE0-D3 (ADMM), CUTOFF 600 Ry / REL_CUTOFF 60,
DZVP-MOLOPT-SR-GTH, GTH-PBE pseudos; OT below 3000 K, diagonalisation +
Fermi smearing at 3000 K. Training set train_v4.xyz: 1193 configurations
(800 bulk 650-2250 K, 344 RCB fill 2500-3000 K, 49 above 3000 K), E0s at
matching level. Fine-tune from mace-mpa-0-medium, SWA. Dipoles: MLWF/
Berry-phase labels, MACE dipole readout; spectra via dipole ACF with
Hanning window, quantum correction w(1-exp(-bhw)), 20 cm^-1 smoothing;
means over 100-5000 cm^-1. Full details: docs/METHODS.md and the deck.
