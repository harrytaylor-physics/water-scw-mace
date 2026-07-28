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
- Planck/ Rosseland means over the 27-point T-P grid
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
| [docs/](docs/) | methods reference; results deck |
| [TODO.md](TODO.md) | technical task list |
| [ROADMAP.md](ROADMAP.md) | understanding-oriented questions |
| [TESTING.md](TESTING.md) | outstanding tests + conditions of use |

## Status

| Component | Status |
|---|---|
| water_scw_v4 (revPBE0-D3 potential) | validated; production; apply <= 3000 K |
| water_dipole_v2 (dipole model) | production; behind the published opacity surface |
| Opacity surface, 27 points | complete (absolute-scale note below) |
| Emissivity, 27 points | complete (interim n_inf baseline) |
| LDA sensitivity campaign | labels complete (1080/1080); water_lda_v1 training |
| Results deck | being rebuilt around the LDA results |

## Methods in one paragraph

Initial CP2K AIMD to generate appropriate configurations; CP2K labels on selected snapshots (revPBE0-D3 or LDA; 600/60 Ry, DZVP-MOLOPT-SR-GTH),
MACE fine-tuned from mace-mpa-0-medium using CP2K labels, CSVR-equilibrated NVE production
MD (0.5 fs, 60 ps, 64 molecules), additional MLWF-trained dipole model, dipole-ACF
spectra with quantum correction factor, opacity means and Fresnel emissivity
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
3000 K, by 1256 kg/m3 at 2250 K); 21-point production archive; three-point
full-band opacity/emissivity comparison vs hybrid delivered (T1000/rho130,
T2000/rho500, T~2807/rho109).
Headline results: EOS dP mean -47% (to -90% cold-dilute); LDA stretch
red-shift 150-200 cm-1; dispersion 1.7% of bulk forces, ~4:1
intermolecular. Presented at the AWE meeting, 28 Jul 2026.
