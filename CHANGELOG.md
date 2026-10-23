# CHANGELOG

Work since the last commit on `main` (`27193f4`, 2026-07-28), oldest first.
Amendment numbers (A*n*) refer to [`subneptune/AMENDMENTS.md`](subneptune/AMENDMENTS.md),
which holds the full reasoning, numbers and evidence for each.

## 2026-07-29 to 2026-08-04 — the planetary model

- `subneptune` built as a separate code: a pure-H2O steam-envelope sub-Neptune
  interior and evolution model.
  - **Core:** inert rocky core, Zeng, Sasselov & Jacobsen 2016 Eq. 2.
  - **Envelope:** AQUA H2O EOS, with the rho(P,T) inversion through the
    vapour-liquid dome (A2, A8).
  - **Outer boundary:** Guillot (2010) two-stream (A10).
  - **Opacity:** a composite opacity, blending the first-principles (FP) spectra
    and the HITEMP and POKAZATEL line lists in log space (A11-A16).
  - **Evolution:** static-sequence evolution (A24).
- Opacity conventions:
  - FP factor-1000 units issue identified (A4);
  - common band 100-4000 cm^-1 with a Rosseland floor (A25);
  - dilute-limit flat-kappa extrapolation (A26).

## 2026-09-13 to 2026-09-20 — systematics and the enlarged grid

- Blend prescription made a config knob and a stated systematic. Sensitivity
  to an FP scale of ±20-25% (A35-A40).
- Domain scans and coverage, 108 cells. Every RCB but two is within 8% of
  flipping (A41-A49).
- **FP factor-1000 error fixed in the data** (`tools/fix_fp_units.py`); the
  originals were archived. 21 extension points merged: **grid 27 -> 48 points**
  (A50). Crossover, blend window, scan and cooling rebuilt on 48 points
  (A51-A56). Structure builds that could hang now fail cleanly (A57).

## 2026-09-23 to 2026-09-26 — continuum, line lists, Extended Data, release tables

- MT_CKD 4.3 water self-continuum added as an opacity source (A58, A63).
  Broadening-convention sensitivity (A59). HITEMP rows at 2250 and 2750 K (A60).
- The dilute FP/line-list excess is a frequency-dependent tilt pointing at the
  dipole derivative; monomer test bundle (A61, A64, A65).
- The FP points are treated as **one grid** in every output and figure. This is
  enforced by a vocabulary gate in `save_pdf` and in both Extended Data ledgers
  (A62, A67).
- Numbers ledger `tools/manuscript_numbers.py` -> `outputs/manuscript_numbers.json`
  (A66). Extended Data items (A67). LDA functional sensitivity: band positions
  only (A68). ARC pull tracked (A69).
- rho^2 scaling test (A70). **Line-wing cut-off**: the largest line-list
  systematic found; kappa_R varies by up to 21.8x across 25-500 cm^-1 (A71).
  AQUA entropy is the uncorrected one (A72).
- **Release tables** in MESA's `kap` plain-text format, plus CSV and fit JSON
  (A73).

## 2026-09-29 — corner cells

- Four cells at (700 K, 10/30/60 kg m^-3) and (800 K, 30 kg m^-3) close the
  coverage gap the boundaries needed: **grid 48 -> 52 points**,
  650-3000 K, 10-2229 kg m^-3 (A74).

## 2026-09-30 — the spectral estimator (A75-A83)

- A75-A77: the dilute "noise floor" and a grid-wide sum-rule violation
  investigated. No quantum-correction calibration fixes it.
- **A78:** `production_campaign/scw/opacity_pipeline.py` line 82 applies the
  harmonic quantum correction **twice**. The dipole model and the MD are
  exonerated: `beta <|Mdot|^2>` is flat in T.
- **A79:** the one-sided normalisation is correct (factor 1), pinned by a
  synthetic classical oscillator.
- **A80:** the fixed estimator built and gated (v4).
- **A81:** the "floor" is spectral leakage from the untapered periodogram. A Hann
  taper on the time series moves the window margin from 1.8x to 1243x.
- **A82:** the smoothing width limit; 10 cm^-1 chosen.
- **A83:** all 52 spectra regenerated (**v5**) and the whole downstream rerun.
  Each spectrum carries a half-record uncertainty. kappa_R(v5)/kappa_R(v3) ranges
  0.089-1.687 (median 0.807). Band-level physics unmoved.
- **A84:** two artefacts in A83's own numbers corrected (a config mix-up and a
  stale CSV).

## 2026-09-30 to 2026-10-01 — ledger guards, the definitive comparison, figures

- **A85: ledger staleness guards.** The ledger refuses inputs older than the
  spectra or than the outputs it depends on. `tools/rebuild_producers.sh` rebuilds
  every producer in dependency order. Producers write atomically.
- **A86:** a structure solve stalling 0.3 K from water's critical point; a single
  build-level retry at looser tolerance, flagged.
- **A87:** registry version markers per spectra directory, checked by the
  loader. Figure and release output directories made redirectable, so the v3
  vs v5 comparison can run in isolation.
- **A88-A89:** Extended Data rebuilt and folded to ten items (eight figures, two
  tables). The window-validation ratios are reproduced from the file that exists.
  The definitive v3 vs v5 comparison, Extended Data included: 38 of 213
  numeric keys changed. The rho^2 claim survives v5 (+3.73 sigma at 2500 K,
  +3.86 at 3000 K).
- **A90:** MLIP error quantities separated: held-out split error, training-set
  parity by regime, NVE drift rate. `docs/METHODS.md` corrected on four counts:
  30 ps, Langevin then NVE, 52 points, 64 or 216 molecules. Pooled-sigma
  homogeneity tested.
- **A91-A93:** Fig. 1b's "notch" is RADIS's line-wing truncation. The panel now
  shows both cut-offs and quotes medians. Three figure labels changed, with no
  measurement moving.

## 2026-10-22 — production-model error table; migration into this repository

- **A94:** the production potential is the stage-one model (epoch 109,
  `models/water_scw_v4.model`, sha256 `dc4a8977...`): held-out
  **20.8 meV/atom, 57.8 meV/A**. The error A90 quoted for it came from the
  stage-two (SWA) table of the same training log; that model was not used for
  production. Verified by evaluation (23.28 meV/atom over
  `train_v4.xyz.gz`). The 41/58 parity in `benchmarks/reference/reference.json`
  is reproduced (41.0/58.3) by the same model. Error table:
  `models/mlip_heldout_errors.json`.
- `benchmarks/run_benchmarks.py`: the parity benchmark read a `REF_forces` key the
  reference file does not store, and copied each frame before reading its
  calculator, so it skipped every frame. Fixed; it reproduces 41.0 / 58.3 meV/A.
- `rebuild_producers.sh` now also runs `tools/aqua_entropy_check.py` and
  `tools/aguichine2025_compare.py`, and the ledger's staleness guard checks
  their outputs. Both run cooling tracks on the composite opacity, so they
  depend on the FP spectra. They had been exempt as "independent of the FP
  spectra", and their 26 Sep outputs were the v3, 48-point values. Regenerated
  on v5; the old values are reproduced exactly by the 26 Sep code on the 26 Sep
  inputs.
- `subneptune` imported with its history (filtered; see
  [MIGRATION.md](MIGRATION.md)). Production pipeline in `production_campaign/`.
  v5 spectra and release tables in `data/`. Figures in `docs/figures/subneptune/`.
  `docs/METHODS.md` carries the A90 corrections.
