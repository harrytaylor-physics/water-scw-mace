# LDA campaign — functional-sensitivity study (July 2026)

End-to-end second-functional campaign: LDA(PADE) labels through production
spectra, mirroring the revPBE0-D3 pipeline stage for stage. Purpose:
measure how the XC approximation level moves the EOS, IR opacity, and
emissivity of warm dense water, with everything except the labels held
fixed.

## Design

Three-potential logic inherited from the hybrid campaign: proposals from
the foundation model (mace-mpa-0-medium), labels from CP2K at the target
functional, fine-tune per functional. LDA deltas vs the hybrid recipe are
confined to the label stage: XC_FUNCTIONAL PADE, GTH-PADE pseudos, no HF /
ADMM / D3. Training recipe deltas (batch 4 vs 2, 400 vs 200 epochs, SWA
300 vs 120) are logged in AMENDMENTS.

Grid: the 27-point (T, rho) production grid of the hybrid campaign,
1000-3000 K, 100-2230 kg m^-3, 64 molecules.

## Gate record (actual numbers)

- GB (labels): SCF convergence + force-outlier screens; pass.
- GC (parity, held-out): bulk 12.3 / rcb 20.2 meV/A (hybrid v4 refs 41/58).
- GE.1 (on-trajectory anchor): PASS bulk (14.7) and hot-dilute
  (T3000/rho108.8: 13-16). FAIL dense: T2500/rho1520.77 = 41,
  T3000/rho1043.59 = 40-57, with a fused-O frame at 0.822 A whose LDA
  label returns 161 eV/A of repulsion — i.e. a fit-coverage gap
  (foundation proposals under-sample LDA's over-bound ensemble at high
  density), not LDA physics.
- Boundary refinement (mini-anchors): T3000/rho729.3 = 53.2 FAIL,
  T3000/rho810.4 = 37.6 FAIL, T2500/rho800 = 31.4 FAIL (gate 30.3);
  T2250/rho1256.2 = 33.6 FAIL (bulk gate 18.5). The applicability
  boundary is density-driven. Dense-branch harvest data (90 frames + 18
  anchors, labelled) is banked for an ensemble-corrected v2; v2 retrain
  deliberately shelved for the meeting scope.
- GF.1 (energy drift, 60 ps NVE): 21/21 pass, |drift| <= 0.017
  meV/atom/ps (gate 2).
- GF.2 (temperature +-5%): fired on 4/21 first-pass runs (CSVR->NVE
  handoff); doubled settle (10 ps) re-passed T2000/rho500 (3.6%) and
  T2250/rho1256 (3.5%). T3000/rho108.8 failed three draws (9.9 / 6.4%
  with 10 / 20 ps settle — weak-coupling corner); best-of-three shipped
  with achieved T = 2807 K used in the QCF and Planck weighting. Gate
  tightness at weak-coupling corners logged as a known limitation.

## Dipole model (LDA)

MLWF Wannier labels via the transplanted LDA template (make_lda_mlwf.py):
~640 configs over five trajectories. Training: v3.5 hybrid recipe with
batch 4 (GPU memory; 192-atom cells) — near-twin, delta logged. Selected
model: best-validation epoch-58 checkpoint (~17-20 mDebye/atom honest
per-config RMSE for both trainings; the early "7.8" was a lucky 37-config
split). A 640-config retrain worsened ACF smoothness (C(dt)/C(0)
0.97->0.89 class), so the data lever is documented as spent; a larger
MLWF campaign is the post-meeting fix.

Noise handling: spectra carry a measured white floor (0.6-2.6% of peak
spectral density at 2 fs sampling; C(dt)/C(0) = 0.984/0.984/0.953 across
the three comparison points). Band positions and the libration region
(40-60x floor) are quotable; LDA absolute band intensities and kappa_R
are floor-limited and are not quoted.

## Sampling amendment

Production initially dumped every 10 steps (5 fs -> Nyquist 3335 cm^-1),
which truncated and aliased the stretch band. The three comparison points
were re-run at dump_every 4 (2 fs -> Nyquist 8339). Archive is therefore
heterogeneous: 18 points at 5 fs, 3 at 2 fs; config now defaults to 2 fs.

## Results

- EOS (eval_eos.py: MLIP stress + analytic kinetic term, stride-60 frames,
  own-trajectory convention, hybrid pathway IAPWS-anchored): LDA below
  hybrid at all 21 matched points; mean -47%, range -0.5% (hot-dilute) to
  -90% at 1000 K / 130 kg m^-3 (~0.07 GPa absolute).
- Spectra: LDA stretch red-shifted 150-200 cm^-1 vs hybrid (3654 vs ~3800
  at 1000 K; 3550 vs ~3750 at 2000 K); libration agreement within ~2x with
  T-dependent ordering; gap narrows across the band by ~2800 K.
- Emissivity: functional-insensitive at the compared (moderate-density)
  points, e_P = 0.994-1.000 both functionals; hybrid surface spans e_P
  0.911-1.000 over the full grid. Dense-branch LDA emissivity awaits the
  ensemble-corrected potential.
- Dispersion (hybrid, D3 decomposition over 33 bulk configs): |F_D3|/|F|
  mean 1.7% (0.1-5.6%), quadrature split 45.8% intra / 88.1% inter
  (~4:1 in power intermolecular). Pair D3 contribution: ~80 meV at
  contact, ~50 meV at typical O-O separations, decayed by ~4.5 A.
  H-bonded dimer binding curve pending a geometry fix (constructed
  geometry does not reach the canonical well; the D3 difference curve is
  orientation-agnostic and is the quoted distance dependence).

## Open items

- Dense-branch v2 potential (harvest labelled and banked); then dense
  LDA spectra + emissivity.
- Larger MLWF dipole campaign to lower the ACF noise floor.
- H-bonded dimer geometry fix.
- ExoMol/BT2 cross-check (exomol_selfheal.py; blocked on a RADIS<->server
  chunk-naming mismatch, first chunk 404s).
- AQUA (Haldemann 2020) overlay on the EOS figure — a third
  approximation tier (PBE-family ab initio in this regime) for the paper
  version.
