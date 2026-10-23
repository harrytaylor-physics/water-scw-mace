# water-scw-mace

Machine-learned interatomic potentials and an IR opacity/ emissivity
pipeline for warm dense water at
sub-Neptune interior conditions: 500-3000 K, ~100-2230 kg/m3.

**Why:** the radiative properties of warm dense H2O control cooling rates and
structural evolution in sub-Neptune interiors, but gas-phase line lists stop being valid at
high density. This repo holds the models and machinery to compute IR
opacities and emissivities from first-principles quality MD and MLWFs.


## What you can do with it

- MD of warm dense water at hybrid-DFT accuracy, with classical-potential cost
- EOS points (pressure includes the ideal-gas kinetic term)
- IR absorption spectra and absolute opacities kappa(nu; T, rho)
- Planck/ Rosseland mean opacities over the T-P grid
- Normal spectral emissivity via a self-consistent dielectric function construction

## Quick start

    git lfs install                  # required once; models are LFS objects
    git clone git@github.com:harrytaylor-physics/water-scw-mace.git
    pip install -r requirements.txt  # needs ase>=3.23 for the CSVR/Bussi protocol

Then see [models/README.md](models/README.md) for the three worked usage
patterns (evaluate / MD / EOS point). Apply the potentials only at
T <= 3000 K within the trained density envelope - full guidance in
[TESTING.md](TESTING.md).

## Map

| Where | What |
|---|---|
| [models/](models/) | trained models + mace architecture schematic + quickstart |
| [benchmarks/](benchmarks/) | benchmarks with pinned references (parity, EOS, timing) |
| [tests/](tests/) | pytest suite: invariances, dipole equivariance, golden values, unit guards |
| [convergence/](convergence/) | convergence protocols per stage + project reference results |
| [studies/](studies/) | satellite studies (functional/ dispersion sensitivity) |
| [docs/](docs/) | methods reference; results deck; [figures/subneptune/](docs/figures/subneptune/) — planetary-model figures |
| [subneptune/](subneptune/) | steam-envelope sub-Neptune interior and evolution model (see below) |
| [production_campaign/](production_campaign/) | the ARC production pipeline that made the 52-point grid, kept as run |
| [data/](data/) | the 52 v5 opacity spectra (LFS) and the MESA-format release tables |
| [MIGRATION.md](MIGRATION.md) · [CHANGELOG.md](CHANGELOG.md) | how subneptune came in; what changed since Jul 2026 |
| [TODO.md](TODO.md) | technical task list |
| [ROADMAP.md](ROADMAP.md) | understanding-oriented questions |
| [TESTING.md](TESTING.md) | outstanding tests + conditions of use |

## Status

| Component | Status |
|---|---|
| water_scw_v4 (revPBE0-D3 potential) | validated; production; apply <= 3000 K |
| water_dipole_v3 (dipole model) | production; behind the 52-point v5 spectra (`models/water_dipole_v3.model`) |
| water_dipole_v2 (dipole model) | the earlier dipole model, superseded by v3 |
| Opacity spectra, 52 points (v5) | complete; regenerated with the corrected estimator (Oct 2026, below) |
| Emissivity, 27 points | complete (interim n_inf baseline); computed before the v5 estimator, not regenerated |
| LDA sensitivity campaign | labels complete (1080/1080); water_lda_v1 training |
| Results deck | being rebuilt around the LDA results |

## Methods in one paragraph

Initial CP2K AIMD to generate appropriate configurations; CP2K labels on selected snapshots (revPBE0-D3 or any other functional; 600/60 Ry, DZVP-MOLOPT-SR-GTH),
MACE fine-tuned from mace-mpa-0-medium using CP2K labels, CSVR-equilibrated NVE production
MD, additional MLWF-trained dipole model, dipole-ACF
spectra with quantum correction factor, opacity means and normal emissivity
in post-processing module. Full details: [docs/METHODS.md](docs/METHODS.md).

## References

**MACE (architecture and foundation models)**
- Batatia et al., *MACE: Higher Order Equivariant Message Passing Neural Networks for Fast and Accurate Force Fields*, NeurIPS 35 (2022). arXiv:2206.07697
- Batatia et al., *A foundation model for atomistic materials chemistry*, arXiv:2401.00096 - the mace-mp/mace-mpa foundation-model line our fine-tunes start from
- Kovacs et al., *MACE-OFF23: Transferable Machine Learning Force Fields for Organic Molecules*, arXiv:2312.15211 - fine-tuning methodology

**Electronic structure (CP2K)**
- Kuhne et al., *CP2K: An electronic structure and molecular dynamics software package*, J. Chem. Phys. 152, 194103 (2020)
- VandeVondele & Hutter, *Gaussian basis sets for accurate calculations on molecular systems in gas and condensed phases* (MOLOPT), J. Chem. Phys. 127, 114105 (2007)
- Goedecker, Teter & Hutter, *Separable dual-space Gaussian pseudopotentials*, Phys. Rev. B 54, 1703 (1996)
- Guidon, Hutter & VandeVondele, *Auxiliary Density Matrix Methods for Hartree-Fock Exact Exchange* (ADMM), J. Chem. Theory Comput. 6, 2348 (2010)

**Functionals and dispersion**
- Zhang & Yang, *Comment on "Generalized Gradient Approximation Made Simple"* (revPBE), Phys. Rev. Lett. 80, 890 (1998)
- Adamo & Barone, *Toward reliable density functional methods without adjustable parameters: The PBE0 model*, J. Chem. Phys. 110, 6158 (1999)
- Grimme et al., *A consistent and accurate ab initio parametrization of density functional dispersion correction (DFT-D3)*, J. Chem. Phys. 132, 154104 (2010)
- Perdew & Zunger, *Self-interaction correction to density-functional approximations* (PADE/LDA parametrization), Phys. Rev. B 23, 5048 (1981); Ceperley & Alder, Phys. Rev. Lett. 45, 566 (1980)

**Thermostat and MD protocol**
- Bussi, Donadio & Parrinello, *Canonical sampling through velocity rescaling* (CSVR), J. Chem. Phys. 126, 014101 (2007)

**Reference EOS**
- Wagner & Pruss, *The IAPWS Formulation 1995 for the Thermodynamic Properties of Ordinary Water Substance* (IAPWS-95), J. Phys. Chem. Ref. Data 31, 387 (2002)

---
Maintainer: Harry Taylor (Oxford).

## LDA campaign additions (Jul 2026)

    lda_campaign/          LDA functional-sensitivity campaign:
      stepE_anchor/          on-trajectory validation, boundary measurement, dense-branch harvest
      stepG_dipoles/         LDA MLWF dipole chain, opacity + comparison figures
      eos_compare/           EOS evaluation (both functionals) + paper figure
      dispersion/            D3 on/off decomposition, dimer scans, MACE check
      emissivity_figures/    e_P surface + spectral waterfall plotters
      linelist/              self-healing ExoMol/BT2 runner
    docs/LDA_CAMPAIGN.md   campaign narrative: design, full gate record, results
    docs/AMENDMENTS_APPEND.md  ledger entries pending merge into AMENDMENTS.md

## Status (27 Jul 2026)

Hybrid (revPBE0-D3): 27-point opacity and emissivity surfaces complete; EOS
evaluated at all 27 points (IAPWS-anchored).
LDA: v1 potential trained and gated; validated envelope = bulk regime +
hot-dilute (boundary measured, density-driven: fails by ~730 kg/m3 at
3000 K, and by 1256 kg/m3 at 2250 K); 21-point production archive; three-point
full-band opacity/emissivity comparison vs hybrid delivered (T1000/rho130,
T2000/rho500, T~2807/rho109).
Headline results: EOS dP mean -47% (to -90% cold-dilute); LDA stretch
red-shift 150-200 cm-1; dispersion 1.7% of bulk forces, ~4:1
intermolecular. 28 Jul 2026.

## Planetary model and v5 opacity (Oct 2026)

Added in this update ([MIGRATION.md](MIGRATION.md) has the full mapping;
[CHANGELOG.md](CHANGELOG.md) the chronology):

    subneptune/              steam-envelope sub-Neptune interior + evolution model, with its history
    production_campaign/     ARC production pipeline (MD, labels, training, opacity), unmodified
    data/fp_opacity_v5/      52 first-principles absorption spectra, v5 estimator (git-LFS)
    data/release/            Rosseland/Planck mean tables (MESA kap format, CSV, fit JSON)
    docs/figures/subneptune/ Figs 1-4 and ten Extended Data items, PDF + 200-dpi PNG
    models/mlip_heldout_errors.json  error tables of the production potential

**The planetary model.** `subneptune` computes the interior structure and
thermal evolution of a sub-Neptune with a pure-H2O (steam) envelope.
- **Core:** inert rocky core (Zeng, Sasselov & Jacobsen 2016, Eq. 2).
- **Envelope:** AQUA equation of state.
- **Outer boundary:** Guillot (2010) irradiated boundary.
- **Opacity:** a composite of the first-principles spectra here, the HITEMP and
  POKAZATEL line lists, and the MT_CKD continuum.
- **Evolution:** static-sequence cooling.

Its physics choices are fixed in `subneptune/SPEC.md`. Every choice made without
a direct citation is logged in `subneptune/AMENDMENTS.md` (A1-A94).
`subneptune/validation/` is its test suite.

**The 52-point grid.** One grid of 52 state points, 650-3000 K and
10-2229 kg m^-3. It has 64 molecules per cell at 39 points and 216 at 13, with
0.5 fs steps and 30 ps NVE after 10 ps Langevin equilibration. The production
potential is `models/water_scw_v4.model`, the stage-one model.
Its held-out error is 57.8 meV/A and 20.8 meV/atom.

**The v5 spectral estimator.** The delivered spectra
(`production_campaign/scw/opacity_pipeline.py`) applied the harmonic quantum
correction twice, at line 82, and used a leaking, untapered periodogram. Both were
replaced, and all 52 spectra regenerated from the same dipole series
(subneptune A78-A83). The new estimator,
`subneptune/opacity/dipole_spectrum.py`, computes
`sigma = beta w^2 I2 / (6 eps0 c N)`:
- I2 is the two-sided power spectrum of the Hann-tapered dipole time series,
  normalised by `sum(w^2)`;
- one harmonic correction;
- a continuous omega^2 prefactor;
- 10 cm^-1 smoothing.

Four checks:

| check | result |
|---|---|
| harmonic-oscillator test: synthetic classical oscillator vs the textbook band strength | 2% (bound 3%) |
| sum rule: each band integral vs `beta <\|Mdot_band\|^2>` | 0.09% (bound 1%) |
| leakage floor below the transparency-window minimum | >= 476x at all 52 cells (bound 100x) |
| dilute-limit rotational band vs the HITEMP line list | median ratio 1.043 (+4%) |

Units, keys and the version marker are in [data/fp_opacity_v5/README.md](data/fp_opacity_v5/README.md).

**Reproducing.**

    git lfs install                        # spectra and models are LFS objects
    cd subneptune
    python -m pytest                       # the planetary model's gates
    bash tools/rebuild_producers.sh        # every product, in dependency order

- `rebuild_producers.sh` writes the figures to `docs/figures/subneptune/` (each
  as PDF + PNG), the release tables to `data/release/`, and the numbers ledger to
  `subneptune/outputs/manuscript_numbers.json`. The ledger is not committed.
- One figure: `python figures/fig1.py config/fiducial.yaml`. The ten Extended
  Data items: `python figures/ed_items.py config/fiducial.yaml` and
  `figures/ed_tables.py`.
- Release tables only: `python tools/release_tables.py config/fiducial.yaml`.
- The line lists, AQUA tables and MT_CKD are external inputs, read from
  `paths.data_root` in `subneptune/config/*.yaml` (default `~/subneptune_data`).
  Gates that need an absent input skip.
- The raw dipole series are not in the repository; SHA-256 manifests are in
  `subneptune/ed_data/`.
- **Expected `python -m pytest` results in `subneptune/`.** These assume the
  external inputs are at `paths.data_root`.
  - Fresh clone: 231 passed, 35 skipped, 4 xfailed. The raw dipole series are
    absent, and so is `subneptune/outputs/`, which the ledger gates read.
  - With the raw dipole series in place (`ed_data/arc_dipoles/`,
    `ed_data/arc_a78/T*/`) and after `bash tools/rebuild_producers.sh`:
    257 passed, 8 skipped, 5 xfailed.
