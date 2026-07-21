# water-scw-mace

Machine-learned interatomic potentials and an IR opacity/emissivity
pipeline for warm dense water ("scw" = supercritical water) at
sub-Neptune interior conditions: 1000-3000 K, ~100-2230 kg/m3.

**Why:** the radiative properties of warm dense H2O control cooling and
structure in sub-Neptune interiors, but line lists stop being valid at
high density. This repo holds the models and machinery to compute
opacities and emissivities from first-principles-quality MD instead.

![MACE architecture](models/architecture.svg)

## What you can do with it

- MD of warm dense water at hybrid-DFT accuracy, classical-potential cost
- EOS points (pressure includes the ideal-gas kinetic term)
- IR absorption spectra and absolute opacities kappa(nu; T, rho)
- Planck/Rosseland means over the 27-point production grid
- Normal spectral emissivity via a self-consistent dielectric construction

## Status

| Component | Status |
|---|---|
| water_scw_v4 (revPBE0-D3 potential) | validated; production; apply <= 3000 K |
| water_dipole_v2 (dipole model) | production; behind the published opacity surface |
| Opacity surface, 27 points | complete (absolute-scale note below) |
| Emissivity, 27 points | complete (interim n_inf baseline) |
| LDA sensitivity campaign | labels complete (1080/1080); water_lda_v1 training |
| Results deck | being rebuilt around the LDA results |

**Known issue (fixed in code, stale in archived data):** a x1000
number-density unit error inflated absolute alpha/kappa in npz files
generated before 2026-06-29. Code fixed; published figures corrected;
archived per-point npz and opacity_surface.npz still carry the factor -
divide by 1000 or regenerate. Ratios were never affected. A unit test
now traps this bug class (tests/test_units.py).

## Map

| Where | What |
|---|---|
| [models/](models/) | trained models + architecture schematic + quickstart |
| [benchmarks/](benchmarks/) | rerunnable benchmarks with pinned references (parity, EOS, timing) |
| [tests/](tests/) | pytest suite: invariances, dipole equivariance, golden values, unit guards |
| [convergence/](convergence/) | runnable convergence protocols per stage + project reference results |
| [studies/](studies/) | satellite studies (functional/dispersion sensitivity) |
| [docs/](docs/) | methods reference; results deck (in rebuild) |
| [TODO.md](TODO.md) | technical task list |
| [ROADMAP.md](ROADMAP.md) | understanding-oriented questions |
| [TESTING.md](TESTING.md) | outstanding tests + conditions of use |

## Quick start

    git lfs install                  # required once; models are LFS objects
    git clone git@github.com:harrytaylor-physics/water-scw-mace.git
    pip install -r requirements.txt  # needs ase>=3.23 for the CSVR/Bussi protocol

Then see [models/README.md](models/README.md) for the three worked usage
patterns (evaluate / MD / EOS point). Apply the potentials only at
T <= 3000 K within the trained density envelope - full guidance in
[TESTING.md](TESTING.md).

## Working on this repo

Pull before you edit; push when done (web edits count - pull them onto
the cluster before local work). Model files go through git-lfs, with a
provenance row in models/README.md for every addition. Pre-register
comparison metrics for new studies; log gate changes in the campaign
AMENDMENTS files.

## Methods in one paragraph

CP2K labels (revPBE0-D3/ADMM or LDA; 600/60 Ry, DZVP-MOLOPT-SR-GTH),
MACE fine-tuned from mace-mpa-0-medium, CSVR-equilibrated NVE production
MD (0.5 fs, 60 ps, 64 molecules), MLWF-trained dipole model, dipole-ACF
spectra with quantum correction, opacity means and Fresnel emissivity
downstream. Full details: [docs/METHODS.md](docs/METHODS.md).

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
Maintainer: Harry Taylor (Oxford). Paper in preparation; cite the repo
until then.
