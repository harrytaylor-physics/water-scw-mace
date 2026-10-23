# AMENDMENTS

Every design decision made without a direct citation, or that resolves a
gap/ambiguity in `SPEC.md`, is logged here with a one-line rationale
(SPEC §10). Dated; newest first.

## 2026-10-23 (four stored outputs were still the v3, 48-point values)

- **A95 — `aqua_entropy_check.json`, `aguichine2025_compare.json` and
  `crossover_v3_{before,after}_a60.json` were built from the 48 v3 spectra and
  had never been rebuilt on v5. All four are now v5 values, and the staleness
  guard checks them.**

  This amendment exists only in the water-scw-mace copy of subneptune
  (`subneptune/`, branch `integrate-subneptune-v5`); `~/subneptune` stops at A94.

  **Why they survived.** A85's audit exempted all four from the staleness guard
  because each had "zero first-principles references". That was a text search of
  the producers, and the dependence is not in their text: the AQUA check and the
  Aguichine tracks run cooling tracks on the composite opacity, and the A60 tool
  takes its ratios from `crossover_v3_rows`. Of A85's five exemptions only
  `broadening_sensitivity.json` (line lists only; no input newer than it) was
  right.

  **(1) The AQUA check and the Aguichine 2025 tracks.**

  *Attribution.* Neither file carries a commit stamp. They were written on
  26 Sep, 13:22-13:29 (AQUA) and 13:32-13:39 (Aguichine): after `162624b`
  (12:26) and before `024e7f9` (16:11), the commit that added both producers.
  - **Producers, `024e7f9` to now:** the only change is the atomic write (A86).
  - **Modules they import:** the migration's paths (`layout.py`), the atomic
    write, A86's one-time retry in `structure/model.py` (it runs only after a
    solve has already failed), the v5 spectra directory (A83) and A87's version
    markers. None changes the calculation.
  - **Inputs on 26 Sep:** the v3 spectra (48 files, dated 20 Sep) with the
    48-point `grid.json` (kept in `~/subneptune/outputs/pre52/`). `fiducial.yaml`
    is unchanged, and every other input predates the files.
  - **Reproduced:** `024e7f9`, run in a worktree on a data root rebuilt to the
    26 Sep state, gives the Aguichine file byte for byte and the AQUA file
    identical except for two path strings. The current code on v5 gives
    byte-identical output on two runs.

  *AQUA check* (M_c = 10 M_E, water fraction 0.5, T_eq = 500 K, composite
  opacity). Old = v3, 48 points; new = v5.

  | quantity | definition | old | new |
  |---|---|---|---|
  | entropy offset | corrected minus AQUA entropy along the fiducial track, 0.001-10 Gyr | 1.646-3.194% | 1.646-3.219% |
  | density check | interpolated against the corrected table, max | 1.18e-4 | 1.18e-4 |
  | adiabatic gradient | median / max relative difference | 3.5-3.9e-4 / 3.94e-2 | 3.5-3.9e-4 / 3.39e-2 |
  | age change via energy | age integral redone with the corrected internal energy, same luminosities; 1 / 10 Gyr | -3.65 / -6.48% | -3.65 / **-6.02%** |
  | age shift at fixed entropy, 10 Gyr | entropy offset / (ds/d age) at the 10 Gyr node; order of magnitude only (A72) | 8.76 Gyr | **7.45 Gyr** |
  | radiated energy | total over the 44-node track down to T_int = 12 K | -2.75% | **-2.40%** |
  | supercritical fraction | | 49.5% | 49.5% |

  The age shift moves because the slope at 10 Gyr steepens (T_int 32.34 ->
  33.63 K, ds/d age -35.3 -> -41.7 J/kg/K per Gyr). The radiated energy is
  dominated by the track far past the age of the Universe: its coldest node sits
  at 387 Gyr on v3 and 1093 Gyr on v5. The same check on `hitemp`, which uses no
  first-principles spectra, gives -3.66% at 1 Gyr, -4.46% at 10 Gyr, -2.10%
  radiated and a 4.56 Gyr shift. So the composite-minus-`hitemp` difference
  narrows on v5 at 10 Gyr (2.0 -> 1.6 points) and in radiated energy (0.65 ->
  0.30 points), and is unchanged at 1 Gyr.

  *Aguichine 2025* (residual = our radius / theirs at 1 µbar, minus 1; 12 cases
  x 3 ages):

  | quantity | old | new |
  |---|---|---|
  | median | 4.05% | **4.00%** |
  | worst | 5.91% | 5.91% |
  | at 1 / 5 / 10 Gyr | 4.56 / 3.61 / 3.19% | 4.56 / 3.58 / 3.17% |
  | ours larger | 36 of 36 (min +2.11%) | 36 of 36 (min +2.05%) |
  | their R_20mbar / R_b | 1.012-1.017 | unchanged (their data) |
  | opacity source, case 1, `hitemp`/composite - 1, at 1 / 5 / 10 Gyr | +0.007 / -0.303 / -0.624% | +0.005 / -0.233 / -0.477% |

  The `hitemp` residuals, +2.90 / +2.05 / +1.63%, are unchanged and reproduce
  A72's.

  **(2) The A60 before/after crossover files.** Written 24 Sep 14:08, with 48
  points: the v3 grid. Their FP opacities at the four affected points are v3's
  (`~/subneptune/outputs_v3/crossover_v3.json`: Planck ratio 1.407, 4.155,
  1.430, 11.448). Regenerated on v5: 52 points. At the 48 shared points every
  line-list column is identical, and only the FP opacities and the four ratios
  move. The FP values match `outputs/crossover_v3.json`, and two runs are
  byte-identical.

  The ledger key `opacity.hitemp_grid_rows_2250_2750.before_after`, ratio
  [before, after]:

  | T [K] | rho [kg/m^3] | FP/HITEMP Planck | FP/HITEMP Rosseland | FP/(H+MT) Planck | FP/(H+MT) Rosseland |
  |---|---|---|---|---|---|
  | 2250 | 175.5 | [1.407, 1.407] -> [2.129, 2.129] | [1.413, 1.370] -> [2.578, 2.500] | [1.334, 1.327] -> [2.018, 2.008] | [1.399, 1.356] -> [2.552, 2.475] |
  | 2250 | 1256 | [4.110, 4.155] -> [6.380, 6.449] | [14.08, 11.03] -> [23.92, 18.74] | [2.887, 3.198] -> [4.481, 4.963] | [13.42, 10.40] -> [22.80, 17.67] |
  | 2750 | 222.7 | [1.430, 1.430] -> [2.478, 2.477] | [1.359, 1.335] -> [2.659, 2.612] | [1.369, 1.374] -> [2.372, 2.381] | [1.343, 1.319] -> [2.627, 2.580] |
  | 2750 | 1841 | [11.31, 11.45] -> [20.26, 20.52] | [34.74, 30.22] -> [61.54, 53.53] | [9.272, 9.821] -> [16.62, 17.60] | [32.37, 27.88] -> [57.36, 49.40] |

  **What A60 measured does not move.** Before and after share the FP numerator,
  so after/before depends on the line list alone. It is identical to two
  decimals on both vintages: -21.66% in Rosseland at 2250 K / 1256 kg/m^3,
  -13.01% at 2750 K / 1841 kg/m^3, and no point away from the added temperatures
  moves (largest 0.00%). What moved is the absolute ratios, by the v3 -> v5
  change at these hot points.

  **Ledger.** 764 entries before and after. **Exactly nine values change**:
  - `aqua_entropy.age_frac_change_from_energy`;
  - `aqua_entropy.age_shift_at_fixed_entropy_gyr`;
  - `aqua_entropy.delta_s_frac_range`;
  - `aqua_entropy.radiated_energy_frac_change`;
  - `g5.aguichine2025.by_age`, `.median_abs_rel`, `.median_signed_rel`,
    `.worst_abs_rel`;
  - `opacity.hitemp_grid_rows_2250_2750.before_after`.

  Every other difference is a `commit` stamp. The ED merge was rerun after the
  ledger, as `rebuild_producers.sh` orders it; its 355 `ed.*` entries are
  unchanged. **No figure or Extended Data item reads these outputs or carries
  one of these keys** (ED Fig. 7 reads `g5_residuals.json`, the Aguichine 2021
  gate on `hitemp`; Figs. 3 and 4 measure the 2250/2750 K sensitivity from the
  current tables). So none was regenerated.

  **Changed.**
  - `tools/rebuild_producers.sh` runs all three producers: the A60 tool as
    step 7, after the crossover; AQUA and Aguichine as steps 13 and 14.
  - `tools/manuscript_numbers.py`: `_STALE_EXEMPT` now holds only
    `broadening_sensitivity.json`.

  **Any manuscript text quoting** an age change of 3.7-6.5% (now 3.7-6.0%), an
  8.76 Gyr shift, -2.75% radiated energy, a 4.05% Aguichine median, or the v3
  absolute ratios at 2250/2750 K **is now out of step with the ledger.**

## 2026-10-22 (the MLIP error table belongs to the production model)

- **A94 — A90 quoted the wrong model's error. The production potential's
  held-out split error is 57.8 meV/A and 20.8 meV/atom, not 44.4 and 7.0.**

  The training log `train_v4_7924966.out` is now in the repository
  (`ed_data/arc_scripts/logs/`, with the ARC submission scripts). It prints TWO
  "Error-table on TRAIN and VALID" blocks:

  | model | lines | valid_Default E / F | train_Default E / F |
  |---|---|---|---|
  | stage one, epoch 109 (`water_scw_v4_run-1_epoch-109.pt`) | 120-133 | **20.8 meV/atom / 57.8 meV/A** | 23.5 / 56.9 |
  | stage two, SWA, epoch 186 (`..._epoch-186_swa.pt`) | 135-152 | 7.0 / 44.4 | 6.7 / 32.7 |

  A90 transcribed the stage-two table. Production MD loaded `water_scw_v4.model`,
  ARC sha256 dc4a8977..., which differs from `water_scw_v4_stagetwo.model`
  (380e20a8...) -- H. Taylor, from the ARC checksums. The water-scw-mace copy of
  `models/water_scw_v4.model` has sha256
  `dc4a8977d44e58cf959901b39a594abd3daf69fb09e5622bbe7e8d3e74b7f06c` (its LFS oid,
  recomputed locally: identical), so the file that ships is the one that ran.
  The log's own naming agrees: stage one was compiled to
  `water_scw_v4_compiled.model`, the name water-scw-mace ships.

  **Checked by evaluation, not only by checksum.** The dc4a8977 model on all
  1193 frames of `train_v4.xyz.gz` (CPU, fp64): energy RMSE **23.28 meV/atom**
  -- the stage-one tables mixed 90/10 predict 23.2; stage two would give ~7 --
  and force RMSE by regime **41.0 (bulk) / 58.3 (RCB) meV/A**, which reproduces
  `reference.json`'s 41 / 58. **So the 41/58 parity was computed with the
  production stage-one model too**; nothing in the pull was computed with the
  SWA model except the log's second table. Most of the energy error is a
  uniform offset (mean -18.9 meV/atom), invisible to forces and to the MD.
  Recorded in `mlip_heldout_errors.json` under `reproduction_a94`.

  Incidentally, `benchmarks/run_benchmarks.py` in water-scw-mace reads
  `a.arrays["REF_forces"]`, but `train_v4.xyz.gz` stores `forces`, which ASE
  moves into a SinglePointCalculator: as written its parity benchmark skips
  every frame. Not fixed here (it is not this repository's code).

  **Changed.** `ed_data/mlip_heldout_errors.json`: `errors` is now the stage-one
  table, with a `production_model` block (stage, checkpoint, sha256, log lines);
  the stage-two table is kept under `stage_two_not_used_for_production`.
  `source_log` is the repository path and `source_in_repository` is true.
  `tools/pipeline_facts.py` and Fig. 1a's caption follow. **Ledger: 764 entries
  before and after; exactly four numeric changes** --
  `pipeline.mlip_heldout_force_mev_per_a` 44.4 -> 57.8,
  `pipeline.mlip_heldout_energy_mev_per_atom` 7.0 -> 20.8,
  `pipeline.mlip_train_force_mev_per_a` 32.7 -> 56.9,
  `pipeline.mlip_train_energy_mev_per_atom` 6.7 -> 23.5 -- plus the Fig. 1a
  caption, facts and notes under `figures.fig1.a.*`. Every other difference is
  the `commit` stamp. Gate: `validation/test_a94_mlip_errors.py` re-reads both
  tables from the log and fails if the JSON quotes stage two as production.

  **What this does to the A88-A90 reading.** The held-out force error (57.8) now
  sits at the training-set RCB parity (58), not below it; the held-out split is
  no longer "better than parity" in any sense. The three-quantity distinction of
  A90 (split error, parity by regime, NVE drift rate) is unchanged. **Any text
  outside this repository quoting 44.4 meV/A or 7 meV/atom for the production
  potential is now out of step with the ledger.**

## 2026-10-01 (three figure labels; the colour-bar extension checked)

- **A93 — label changes only. No measurement moves: the figure ledger has zero
  numeric differences against A92.**

  **(1) Fig. 1a**, last stage box: "dipole autocorrelation, mass opacity" ->
  "dipole autocorrelation, κ(ν)". **(2) Fig. 3b** y-axis: "Factor" -> **"Ratio or
  threshold"**, because the axis carries two different quantities -- the measured
  ratio and the threshold it has to beat -- and "Factor" named neither. The
  `nums.axes` label was changed with it, so the ledger and the axis agree; they
  would otherwise have disagreed, since the ledger's `y_label` comes from the
  `nums.axes` argument rather than from `set_ylabel`.

  **(3) Fig. 2a's colour-bar extension, checked rather than assumed.** The bar's
  limits are taken from the **model** surface on the display grid, 2.782510 to
  7583.703531 cm^2/g, while the panel draws the Clough-Tocher display surface
  (A91), so the two can disagree at the edges -- which is what `extend="both"` was
  quietly covering.

  | end | outside the limit? | by how much |
  |---|---|---|
  | maximum | **no** | largest computed point 7583.703530966697 vs vmax 7583.703530966707 — **-5.8e-16 dex**, floating-point identical |
  | minimum | **yes** | coldest, thinnest computed point 2.740488 is **0.0066 dex** below vmin; the displayed surface reaches **1.8629 cm^2/g**, **0.174 dex** below, over **1195 of 182583** in-hull cells (0.65%) |

  So the lower arrow is earned and the upper one is not: `extend="min"`. The
  minimum undershoot is not an error -- the display grid does not sample the
  corner point exactly, and a C1 interpolant may dip below the data between
  nodes -- but it is real, so the arrow stays. Recorded as
  `figures.fig2.a.colourbar_extend` and `colourbar_extend_reason`.

  **Ledger diff against A92**: two keys added (both the provenance above), one
  text value changed (`figures.fig3.b.y_label`), `_n_entries` 758 -> 760 as the
  count of those two. **Zero numeric measurements changed.** Figure gates: 42
  passed.

## 2026-10-01 (Fig. 1b states the cut-off dependence instead of hiding it)

- **A92 — Fig. 1b now plots the line list at both wing cut-offs and quotes medians
  over a window. The single-channel 2550 cm^-1 marker is gone: A91 showed it sat
  on the pedestal of the line at 2543.07, so what it measured was the cut-off.**

  A91 established that the comparison does not converge in the truncation and
  then left the panel quoting one number from one cap with the sensitivity buried
  in the ledger. That was the wrong way round. The panel now shows the line list
  at **50 cm^-1** (the delivered convention, which every line-list table in the
  manuscript is built at) and at **500 cm^-1** (the widest A71 scanned, dashed and
  lighter), shades 2400-2600 cm^-1, and draws a median bar across that window for
  each of the three curves.

  | curve | median over 2400-2600 cm^-1 | first principles / line list |
  |---|---|---|
  | first principles (v5) | **370.5** | — |
  | HITEMP, 50 cm^-1 cut-off | **13.17** | **x28.1** |
  | HITEMP, 500 cm^-1 cut-off | **30.26** | **x12.2** |

  Every one is computed from the three spectra when the figure is built, not
  carried over from A91's text. The two ratios differ by more than a factor of
  two, and that *is* the statement: a reader sees the cut-off dependence on the
  face of the figure rather than meeting a single "x35" that happens to be where
  a pedestal edge fell.

  **Ledger.** `figures.fig1.b.marked_wavenumber_cm`, `marked_kappa_*`,
  `marked_ratio`, `caption_correction` and `truncation_sensitivity` are **removed**;
  `median_window_cm`, `median_first_principles_cm2_g`, `median_line_list_cm2_g`,
  `median_ratio_fp_over_line_list`, `why_two_cut_offs`,
  `single_channel_marker_removed` and `producer_for_the_500_spectrum` are **added**;
  `line_list_truncation_cm` goes 50.0 -> [50.0, 500.0]. **Nothing else in the
  ledger moves** — the diff is confined to `figures.fig1.b.*`, because no spectrum
  or table was regenerated.

  The panel raises `SystemExit` naming the exact
  `tools/broadening_sensitivity.py` invocation if either line-list spectrum is
  absent, so the second cut-off cannot silently go missing on another machine.

## 2026-10-01 (the notch is the wing cut-off; six figure changes)

- **A91 — Fig. 1b's "notch" is RADIS's line-wing truncation, not line structure,
  and raising the cap does not converge the comparison: the marked ratio is 35 at
  50 cm^-1, 27 at 200 and 12 at 500. The figure therefore keeps the delivered
  convention and names it. Five other figure changes, and one product I corrupted
  and restored.**

  **(1) The notch.** At (2000 K, 1200 kg/m^3) the spectrum has flat plateaus with
  sharp rectangular edges. The edges pair at **exactly 100.01 cm^-1 = 2 x 50**:
  lines at 2473.76 and 2543.07 switch on at nu_0 - 50 and off at nu_0 + 50. That
  is the `truncation` cap, and `tools/broadening_sensitivity.py`'s own docstring
  predicts it -- RADIS renormalises each truncated profile to unit area, so where
  the true HWHM greatly exceeds the cap (here O(1e3) cm^-1 against 50) a line
  becomes a ~100 cm^-1-wide rectangle carrying its full strength. **The marked
  2550 cm^-1 channel sits on the pedestal of the 2543.07 line.** Artefact, not
  physics.

  Regenerated at 200 and 500 cm^-1 (A71's range; `--truncation`, new):

  | cap | kappa(2550) | FP/LL at 2550 | band median 2400-2600 | FP/LL median | steps >2% | LL kappa_R |
  |---|---|---|---|---|---|---|
  | **50** (delivered) | 9.5213 | **35.0** | 13.1662 | 28.1 | 8 | 50.64 |
  | 200 | 12.1568 | **27.4** | 13.3448 | 27.8 | 1 | 63.48 |
  | 500 | 28.9124 | **11.5** | 30.2576 | 12.2 | 0 | 116.27 |

  The pedestals do disappear (8 steps -> 1 -> 0) but **nothing converges**: the
  ratio moves by a factor of 3 and kappa_R more than doubles. A band median is not
  the escape either -- it looks stable from 50 to 200 and then jumps 2.3x at 500.
  Monotone and physical: the wider the retained wing, the more far-wing absorption
  from the bend and stretch bands fills the transparency window.

  **So Fig. 1b keeps 50 cm^-1** -- every line-list table in the manuscript is built
  at that cap, and changing this panel alone would make it inconsistent with them
  -- and now **records the cap and the sensitivity** in its ledger entry, with the
  caption naming it. This is A71's systematic ("any 'N times the line list' claim
  has to name the cut-off") landing on a specific figure. **No downstream number
  changes**, because nothing downstream was regenerated: the delivered spectra and
  tables are untouched.

  **A product I corrupted, and the guard that stops it recurring.** The 200 cm^-1
  run silently damaged `outputs/broadening_sensitivity.json`. `--out-dir`
  redirects the spectra but the summary path was `out_dir.parent`, so a
  single-point, single-convention run at a non-default cap replaced one row of the
  six-point file -- whose header still said `truncation_cm: 50.0` -- and dropped
  that row's `ref1bar` and `aqua` conventions. Restored exactly from the intact
  50 cm^-1 spectra (`kappa_R_ideal` back to 50.6427, matching the ledger to every
  digit). A run that is not the full default set now writes its summary beside its
  own spectra and is stamped `partial_run`; the banner and the JSON header record
  the cap **actually used** rather than the module constant, which they did not.

  **(2) Fig. 1a** drops "drift 0.017". It is a stability check on the integrator,
  not a validation of the opacity, and it was the one caption number a reader
  could not place. Still in the ledger as `pipeline.nve_drift_mev_per_atom_per_ps`.

  **(3) Fig. 3b** gains the bounded dispute as a second marker per boundary --
  ratio / 1.400 at rho >= 60, / 3.039 below. Both are ledger quantities, not
  choices made in the figure: 1.400 is `window_validation.ratio.window`, measured
  against Wannier-centre dipole labels, and 3.039 is A88's largest v5 dilute
  stretch-band excess, at (2000 K, 30 kg m^-3). Nine of ten boundaries stay above
  their flip threshold once bounded; the tenth (1008 K, 44 kg m^-3) does not, and
  now shows it. **Fig. 3a** gets a fixed 1-1e4 kg/m^3 axis and reads "first
  principles" where it said "composite".

  **(4) Fig. 2b** is Rosseland only, against HITEMP and HITEMP + MT_CKD; it
  carried three references x two means = six overlapping series, and kappa_R is
  what the structure model integrates. POKAZATEL and the Planck ratios move to a
  new **ED Fig. 2c**, which also forced ED Fig. 2 into two rows -- three panels
  with colour bars do not fit across one journal width, and b's colour-bar label
  ran into c's y-label. The b label aligns with a (dy 1.30 -> 1.06).

  **Fig. 2a's faceting was real.** The model interpolates linearly on a Delaunay
  triangulation, so the facets *are* the model, not a rendering artefact. The
  surface is therefore re-interpolated **for display only** with a C1
  Clough-Tocher scheme on the same 52 points, keeping the model's own hull mask,
  and `display_vs_model_max_dex` records how far the displayed surface departs
  from the one the model uses. The stepped lower edge was genuinely only the
  `pcolormesh` cell size (NGRID 320 -> 480, `contourf`).

  **(5) Fig. 4**: colour-bar ticks at 2, 5, 10, 20, 50; "rho below the grid" moves
  from orange to blue, because vermillion and orange are adjacent hues and read as
  one colour at print size; and the zero-count category leaves the legend, since a
  category with no members is not a key to anything. **(6) Fig. 1b**'s first band
  reads "rotation/libration" -- it is "rotation" in the manuscript and "libration"
  in `opacity/bands.py`, and this is where a reader meets it first.

  **Two faults of my own in the figures, both caught by looking at them.**
  `contourf` puts a tick on every contour level, so Fig. 2a's colour bar printed
  `1.00634 x 10^3` and nothing else; and ED Fig. 2 at three panels across one width
  had overlapping labels. Both fixed before commit.

  validation: 251 passed, 8 skipped, 5 xfailed.

## 2026-10-01 (the pooling assumption tested; the held-out errors found)

- **A90 — the pooled sigma survives a test on the diagnostic that matters, so
  A70's convention stands; the held-out MLIP errors A89 called absent exist in a
  log that was never delivered; METHODS.md corrected against the artefacts.**

  **(1) The pooling assumption, tested rather than asserted.** `_scatter_pooled`
  pools residuals over every isotherm, which is what buys the ~23 dof A70 quoted
  its errors on; the single-pair alternative carries one. `tools/rho_scaling.py`
  now fits each isotherm with the same quadratic in ln rho, takes the residuals,
  and runs Levene (robust, median-centred) and Bartlett (normal-theory) on their
  variances, writing `scatter_homogeneity`.

  | diagnostic | Levene p | Bartlett p | variance ratio | pooling |
  |---|---|---|---|---|
  | **window_2550** | **0.597** | **0.378** | 10.2x | **not rejected** |
  | window_2550_pm50_mean | 0.023 | 0.0014 | 70.7x | **rejected** |
  | window_1000_1200 | 0.340 | 0.015 | 60.3x | not rejected (tests disagree) |
  | rosseland | 0.054 | 0.027 | 21.6x | marginal |

  **On A70's own diagnostic the assumption holds**, so the pooled error bar — and
  with it the rho^2 claim — is defensible where the manuscript makes it. It fails
  on the +-50 cm^-1 variant and is marginal for kappa_R, so the pooled sigma
  should not be carried to those without saying so. Seven isotherms, 1-5 dof each
  (n - 3 per isotherm); Bartlett is the more sensitive and the less robust, and
  where the two disagree the groups are small enough that Levene is the one to
  read. The A89 fork is **narrowed, not closed**: the test licenses the pooled
  convention for `window_2550`, it does not make the single-pair estimate wrong.

  **The significance at 2550 cm^-1 under v5, both ways:**

  | | slope | pooled sigma | single-pair sigma |
  |---|---|---|---|
  | 2500 K | 2.638 | **+3.73** | +1.88 |
  | 3000 K | 2.434 | **+3.86** | +1.95 |

  **(2) The held-out errors exist; A89's finding was right about the pull and
  wrong about the work.** They come from the training job's own stdout,
  `train_v4_7924966.out`, which was never part of the delivery: `valid_Default`
  **7.0 meV/atom and 44.4 meV/A**, `train_Default` 6.7 and 32.7. Recorded in
  `ed_data/mlip_heldout_errors.json` — beside the pull, not inside it, following
  `arc_dipoles_sha256.json`, and carrying `source_in_repository: false` so the
  chain is visible. Four `pipeline.mlip_*` ledger entries cite it. **These do not
  replace the 41/58 meV/A parity**; that is training-set parity BY REGIME, which
  `reference.json` itself calls optimistic, and the 58 at RCB conditions remains
  the number for the regime the structural claim lives in. Three different error
  quantities, all three now in the ledger, and Fig. 1a shows the held-out pair
  with the regime parity beneath it.

  **(3) METHODS.md corrected against the artefacts, with the delivered text
  kept.** It is inside the ARC pull, so it is corrected rather than rewritten and
  a dated block at the top quotes what it said: 60 ps -> **30 ps** (every dipole
  series is 60,000 frames at 0.5 fs); NVT -> **Langevin equilibration then NVE**
  (`run_prod_md_05fs.py` lines 48 and 50); 27-point grid -> **52 state points**;
  and a fourth the brief did not list, 64 molecules -> **64 or 216** (the `N`
  field takes both). **The "27" is almost certainly POKAZATEL coverage**:
  `pokazatel_opacity/` holds exactly 27 spectra and all 27 fall on FP state
  points, T = 1000-3000 K and rho = 108.796-2228.71. That identification is an
  inference from the counts, not something the pull states. The file's "Dipoles
  and spectra" section is flagged as superseded rather than rewritten — it
  describes the v3 estimator, including the white-noise floor A81 showed was
  leakage and must never be subtracted.

  **(4) `config/candidate.yaml`: M_c = 5 M_earth, f_env = 0.05, T_eq = 700 K.**
  There is no `wmf` key anywhere in `config/`; `f_env` is defined in
  `config/schema.py` as `M_env / M_p`, and because the envelope is pure H2O on a
  `zeng2016_rock` core the water mass fraction **is** `f_env`. **A86 and A89 do
  not conflict**: A86's "f_env 0.02-0.10" is `F_ENV_GRID = (0.02, 0.05, 0.10)` in
  `figures/plot_s6_candidate_evolution.py`, the scan S6 runs **around** the
  config value; A89's `w0 = 0.05` is the config value itself, the scan's centre.
  One consequence worth stating: w0 = 0.05 sits exactly on the bottom of the
  tested set {0.05, 0.20, 0.50} and the candidate's shift is downward at every
  age, so "outside the tested range" there is structural, not a pathology. To
  bracket it, the grid needs points below 0.05.

## 2026-10-01 (A88 finished; the remaining manuscript numbers, from v5)

- **A89 — A88's parts 4-7 are done and every number below is rebuilt from its
  producer against v5, both columns from producers. The headline: **A70's
  "faster than rho^2" survives v5** at +3.73 sigma (2500 K) and +3.86 sigma
  (3000 K) on A70's own pooled-sigma convention — 3000 K strengthens — but the
  producer's own `sigma_above_2` field uses a different sigma and reports +1.9;
  which to quote is a statistics fork for Harry, below. Two validation numbers in
  the A88 brief are not in this repository. And a gate has been asserting one
  thing instead of four since A83.**

  **(4) Water-fraction inversion** — `tools/wmf_inversion.py`, new. At fixed core
  mass and age, `dwmf = w0 - w_LL = -(R_comp - R_LL)/(dR/dwmf)`; negative means the
  line list would have overstated the water fraction. `np.gradient` and a plain
  secant agree to five decimals, and the per-wmf radii are stored so the gradient
  can be audited.

  | case | 1 Gyr | 5 Gyr | 10 Gyr |
  |---|---|---|---|
  | water world M_c = 10, w0 = 0.50 | no shift | -0.00004 (-0.007%) | **-0.00125 (-0.251%)** |
  | water world M_c = 15, w0 = 0.50 | no shift | no shift | **-0.00080 (-0.159%)** |
  | candidate M_c = 5, w0 = 0.05 | -0.00041 (-0.814%) | -0.00277 (-5.534%) | **-0.00399 (-7.977%)** |

  "No shift" is `|dwmf| < 1e-5`, the solver-noise floor: dR between the tracks is
  ~1e-6 R_E where they agree and dR/dwmf is ~2, so nothing smaller is real. **The
  candidate needs wmf < 0.05 at all three ages**, below the bottom of the tested
  set {0.05, 0.20, 0.50}; reported, not extrapolated. One grid point
  (M_c = 15, wmf = 0.30) does not build under the integration guard. These
  supersede the one-off numbers quoted mid-session (-0.00284, -0.00191 at 10 Gyr),
  whose grid is not reconstructible.

  **(5) Fig. 1a — four stages, not six.** The six-box version drew the dipole
  autocorrelation as a stage of its own, giving it the weight of the whole DFT
  campaign. Each box now carries the number that validates it, and those numbers
  come from `tools/pipeline_facts.py`, which is also where the ledger's
  `pipeline.*` entries come from — one source, so the schematic cannot quote a
  number the ledger does not carry.

  **Two numbers in the A88 brief are not in this repository.** It gave a held-out
  force error of 44 meV/A and 7 meV/atom. Neither string appears anywhere in the
  pull. `benchmarks/reference/reference.json` carries **41 meV/A (bulk) and
  58 meV/A (RCB)** and calls them TRAINING-SET parity, "optimistic vs true
  held-out error", frames above 3000 K excluded; the only held-out parity in the
  pull is 12.3 / 20.2 meV/A and belongs to the **LDA** model, not this one. The NVE
  number is not an energy error but a **drift rate**, 0.017 meV/atom/ps worst over
  21 runs against a 2 meV/atom/ps gate. The figure states what the files state.

  **(6) Extended Data is ten items.** ED Fig. 9's seven scaling panels are now ED
  Fig. 4c-i: the dilute band ratios and the rho^2 scaling are one question asked at
  two densities. 183 x 236 mm, 118 recorded numbers, nine panels in the ledger. Two
  stale counts fixed in passing — ED Fig. 4i's legend said "computed points (48)"
  and `ed_index.py` called ED Table 1 "the 48 state points"; the grid has been 52
  since Tasks 42-44.

  **(7) Both columns rebuilt from producers. 736 keys each, 731 shared, 213
  numeric, 38 changed (17.8%).**

  | group | numeric | changed |
  |---|---|---|
  | agreement | 6 | **6** |
  | candidate_fiducial | 18 | **18** |
  | evolution | 3 | **3** |
  | s6_closure | 3 | **2** |
  | ed | 50 | 6 |
  | figures | 13 | 3 |
  | broadening | 56 | 0 |
  | window_validation | 15 | 0 |
  | g5 | 15 | 0 |
  | wing_cutoff | 13 | 0 |
  | pipeline | 6 | 0 |
  | boundaries, opacity, planets, aqua_entropy, direct | 28 | 0 |

  `pipeline.*` unchanged is the right answer: those describe the MD campaign, not
  the spectra. The S9 case set is **not** the same on both sides — v3 built
  `M_c5_wmf0.50` and v5 built `M_c15_wmf0.20` — so five keys each way are reported
  as present-in-one-column, not as unchanged. `tools/column_diff.py`, new, keeps
  that distinction, which is the distinction A86 lost. Two scripts failed under v3
  only (`plot_f3_cooling_ratio.py`, `plot_sensitivity.py`); neither writes ledger
  keys, so the comparison is unaffected.

  **Item 3 — the rho^2 excess survives v5, and there are two sigmas in the file.**
  `outputs/rho_scaling.json` carries both: `scatter.per_diagnostic` is a
  **single-pair** estimate (1 dof, 0.5436 at 2550 cm^-1 under v5) and
  `scatter_pooled.per_diagnostic` is pooled over the isotherms (**23 dof**, 0.2743).
  The stored `standard_error` and `sigma_above_2` use the pair one; **A70's text
  says "errors from the pooled sigma", and ED Fig. 4g/h rescale by pooled/pair
  before drawing.** Running today's producer against the **v3** spectra and
  applying that rescaling returns A70 exactly — 2.759 +- 0.201 = **+3.78 sigma**
  at 2500 K and 2.423 +- 0.132 = **+3.20 sigma** at 3000 K, against A70's
  2.76 +- 0.20 = 3.7 and 2.42 +- 0.13 = 3.1. **Nothing changed between A70 and
  now.** Under v5, on that same convention:

  | 2550 cm^-1 | slope | pooled sigma (A70, ED 4g/h) | pair sigma (JSON field) |
  |---|---|---|---|
  | v3, 2500 K | 2.759 | **+3.78** | +2.08 |
  | v5, 2500 K | 2.638 | **+3.73** | +1.88 |
  | v3, 3000 K | 2.423 | **+3.20** | +1.76 |
  | v5, 3000 K | 2.434 | **+3.86** | +1.95 |

  **FORK, not resolved here.** Pooled has the dof but assumes the scatter is
  common across isotherms; the pair estimate is local but rests on one degree of
  freedom. A70 chose pooled and the figure follows it, so the manuscript is
  self-consistent as it stands and the claim holds. Flagged for Harry rather than
  decided. Every other `sigma_above_2` quoted in this entry — including the
  +-50 cm^-1 variant's +3.39 — is the **pair** figure.

  Where v5 moves the physics is the cold end, and upward: 650 K goes
  1.578 -> 2.238, 800 K 1.719 -> 2.153, 1000 K 1.675 -> 1.917, so the cold
  isotherms now sit at rho^2 where v3 had them sub-rho^2. That part is spectral.
  **Two corrections of my own this session: I first reported this group as
  "+0.25 sigma", and then as a 3.7 -> 1.9 collapse caused by an error model. Both
  were wrong; the first read the wrong statistic and the second read the wrong
  field.**

  **Item 6 — the first 100 Myr is identical, and the reason is the composite's
  definition.** Radius and interior entropy agree to |dR| <= 1e-6 R_E and
  |dS| <= 0.023 J/kg/K at 10, 30, 50 and 100 Myr for both water worlds and the
  candidate. Not a clamped interpolation this time: the composite blends first
  principles in only above `blend_window_kg_m3 = [10, 30]`, and over that span the
  radiative-convective boundary sits at rho = 0.017-3.97 kg/m^3, **below the
  window**, so the composite evaluates to the pure line list and the tracks are
  identical by construction. They first part company at 4909 Myr (M_c = 10),
  5604 Myr (M_c = 15) and 497 Myr (candidate) — in each case as rho_rcb enters the
  window. `tools/early_evolution.py`, new.

  **A gate has been asserting one thing instead of four since A83.** A83 restated
  `test_agreement_with_the_delivered_analysis`'s Planck pin for v5 and left a bare
  `return` in front of the three assertions after it, which made them unreachable
  rather than restated — the comment said "Restated, not loosened" and the code did
  neither. All three would have failed: measured under v5 the deviations from the
  v3-era RADIS reference are **82.5% (on grid), 102.1% (T-interpolated) and 96.2%
  (density-clamped)** in kappa_R. Restored and restated on A83's own logic, that
  the delivered analysis describes the v3 spectra so the useful pin is that the two
  now differ substantially. This is **not** the density-clamp cost, a different
  comparison that stays pinned at <= 6.5% by
  `test_the_density_clamp_costs_at_most_six_percent`.

  **Gate changed.** `test_the_ledger_never_names_the_delivery_history` now reads
  prose tokens rather than raw substrings, as the Extended Data gate beside it
  already did: Fig. 1a records `.../docs/LDA_CAMPAIGN.md` as the source of its
  drift number, and naming the file a number came from is the opposite of hiding
  its provenance.

  **Still open.** The S6 blend-window claim (the sliver spans rho = 10.2-49.7) and
  the band-diagnostic monotonicity xfails are untouched. `METHODS.md` contradicts
  the delivered data on two counts — it says 60 ps per state point and a 27-point
  grid, where the dipole series give 30 ps and the spectra directory holds 52
  cells; `tools/pipeline_facts.py` records the data and says so per fact.

## 2026-09-30 (the window validation reduced; the bound measured; T_int against age)

- **A88, parts 1-3 — ED Fig. 5e now reads the file that exists and reproduces all
  five band ratios to 0.0004, which unblocks the Extended Data group A87 had to
  exclude. The "conservative bound" turns out not to exist in this repository at
  all, which is why A83-A87 kept reporting it unchanged. Replaced with two
  measured factors, 9 of 10 boundaries survive. Parts 4-7 are not done and are
  listed at the end.**

  **(1) ED Fig. 5e and the window numbers.** The panel expected pre-reduced
  `band_centre_cm` and `A_pred_over_A_label` arrays; the delivered file carries
  the spectra (`nu`, `A_label`, `A_pred`, `A_resid`) and both dipole series
  (`M_label`, `M_pred`, 501 frames at 2 fs = 1 ps). The reduction now happens in
  `figures/ed_items.py`, and `figures/ed_items.py` runs to completion for the
  first time since the file changed format — all nine ED figures build.

  **The reduction is a plain channel sum, not a trapezoid**, and that matters:

  | band (cm^-1) | reproduced | reported | residual/label | reported |
  |---|---|---|---|---|
  | libration 300-1000 | **1.0234** | 1.023 | **0.0042** | 0.0042 |
  | bend 1400-1800 | **1.3589** | 1.359 | **0.0332** | 0.0332 |
  | window 1900-3100 | **1.3998** | 1.400 | **0.1709** | 0.1709 |
  | stretch 3200-4000 | **1.3464** | 1.346 | **0.0850** | 0.0850 |
  | overtone 4200-8000 | **0.9599** | 0.960 | **0.6680** | 0.668 |

  All five within **0.0004**. With trapezoids instead, the first four move by
  0.001-0.023 and the overtone by **0.086** — the grid is uniform at 33.29 cm^-1,
  so the two differ only in the endpoint half-weights, but that is enough where a
  band edge lands on a steep part of the spectrum. The producing script used sums.

  The two-halves values come from the dipole series through the same estimator:
  they agree with the spectral arrays to 0.002-0.007 in the four physical bands
  (libration 1.0254, bend 1.3657, window 1.3769, stretch 1.3464) and diverge in
  the overtone (0.8655), which is the band whose residual is 67% of its label.
  Their spread gives the uncertainty, and it is now in the ledger:
  **window 1.400 +- 0.021**, libration +- 0.018, bend +- 0.035, stretch +- 0.005,
  overtone +- 0.192. Sixteen `window_validation.*` entries in all.

  **(2) The conservative bound, and why it kept reporting unchanged.** There is
  **no conservative-bound entry in the ledger and no 4.4 anywhere in the
  repository.** Every key matching "bound" is a boundary quantity — dispute
  ratios, flip thresholds, enclosed-mass changes. The 4.4 lives in the manuscript
  draft, which is not in this tree.

  So the honest answer to why A83 through A87 reported it unchanged is: **they did
  not measure it. My own grouping regex matched zero scalars and I printed
  "unchanged" for an empty group** — a reporting fault, repeated four times. An
  empty group should have read "no entries", and A88 is where that is corrected.

  **The two measured factors.** At rho >= 60 kg m^-3 the window-validation ratio,
  **1.400**; at rho < 60 the largest v5 dilute stretch-band excess over the line
  list, which is **3.039** at (2000 K, 30 kg m^-3) — the fourteen dilute points
  span 2.244 to 3.039 with a median of 2.642. Applied to all ten boundaries on
  computed data (in box, unclamped):

  | T (K) | rho | dispute | factor | bounded | threshold | exceeds? |
  |---|---|---|---|---|---|---|
  | 1008.0 | 43.8 | 2.55 | 3.039 | **0.84** | 1.19 | **no** |
  | 732.2 | 54.1 | 6.29 | 3.039 | 2.07 | 1.12 | yes |
  | 739.1 | 54.6 | 6.18 | 3.039 | 2.03 | 1.03 | yes |
  | 735.0 | 61.3 | 6.77 | 1.400 | 4.84 | 1.05 | yes |
  | 735.4 | 81.9 | 8.97 | 1.400 | 6.41 | 1.08 | yes |
  | 744.8 | 90.1 | 9.55 | 1.400 | 6.82 | 1.05 | yes |
  | 741.0 | 96.3 | 10.58 | 1.400 | 7.56 | 1.01 | yes |
  | 734.1 | 107 | 12.01 | 1.400 | 8.58 | 1.09 | yes |
  | 726.1 | 110 | 12.83 | 1.400 | 9.17 | 1.21 | yes |
  | 736.5 | 124 | 12.95 | 1.400 | 9.25 | 1.03 | yes |

  **Nine of ten survive.** The bounded dispute ratio runs **0.84 to 9.25**, and
  the smallest bounded dispute, 0.84, is the one that fails: (1008 K, 43.8), the
  hottest and least dense of the ten, where the dilute factor of 3.039 is nearly
  as large as the dispute itself. Every boundary at rho >= 60 clears its threshold
  by a factor of 4.6 to 9.2 even after dividing by 1.400.

  **(3) Internal temperature against age.** T_int interpolated on each cooling
  track, line list against composite (v5):

  *Water-world fiducial, M_c = 10, wmf = 0.50, T_eq = 500 K —
  `config/fiducial.yaml`:*

  | age (Gyr) | 0.1 | 0.5 | 1 | 2 | 5 | 10 |
  |---|---|---|---|---|---|---|
  | line list | 125.22 | 94.29 | 77.16 | 57.27 | 44.68 | 37.80 |
  | composite | 125.22 | 94.29 | 77.16 | 57.27 | **43.88** | **33.61** |
  | difference | 0.00 | 0.00 | 0.00 | 0.00 | **-0.80 K, -1.78%** | **-4.19 K, -11.09%** |

  *Water world, M_c = 15, wmf = 0.50, T_eq = 500 K — `config/fiducial.yaml`:*

  | age (Gyr) | 0.1 | 0.5 | 1 | 2 | 5 | 10 |
  |---|---|---|---|---|---|---|
  | line list | 134.57 | 102.49 | 84.47 | 63.26 | 48.73 | 41.25 |
  | composite | 134.57 | 102.49 | 84.47 | 63.26 | **48.31** | **37.57** |
  | difference | 0.00 | 0.00 | 0.00 | 0.00 | **-0.42 K, -0.86%** | **-3.68 K, -8.93%** |

  *Candidate, M_c = 5, f_env = 0.05, T_eq = 700 K — `config/candidate.yaml`:*

  | age (Gyr) | 0.1 | 0.5 | 1 | 2 | 5 | 10 |
  |---|---|---|---|---|---|---|
  | line list | 52.07 | 34.23 | 29.07 | 24.62 | 19.76 | 16.88 |
  | composite | 52.07 | 33.75 | 27.19 | 22.35 | 17.50 | 14.94 |
  | difference | 0.00 | **-1.40%** | **-6.46%** | **-9.23%** | **-11.44%** | **-11.50%** |

  The two water worlds are **indistinguishable out to 2 Gyr** and separate only
  late, reaching 8.9-11.1% by 10 Gyr; the candidate separates from 0.5 Gyr and
  saturates near 11.5%. The composite is cooler at fixed age throughout, which is
  the sign a more opaque envelope gives: it retains heat longer, so at any age it
  has cooled less far along its own track and sits at a *lower* T_int than the
  line-list model at the same age.

  **One caveat, and it changed the numbers.** The candidate's track on a T_int
  grid floored at 20 K ends at **4.71 Gyr** (line list) and 2.99 Gyr (composite),
  so 5 and 10 Gyr fell outside it and `np.interp` silently clamped both to the
  floor — reporting 20.00 K and a difference of exactly zero. The grid floor is
  lowered to 12 K, which extends the tracks to 32.5 and 21.6 Gyr and gives the
  real values above. **The water-world tracks span 242 and 280 Gyr and never
  needed it.**

  **Parts 4 to 7 are not done.** Part 4 (water-fraction inversion at fixed radius
  and age) needs a two-dimensional solve with the core mass tracking
  `M_p (1 - wmf)`, part 5 is a new vector schematic for Fig. 1a, part 6
  restructures the Extended Data to ten items, and part 7's final table depends
  on 4 to 6 being in place. None is blocked; they are simply not built yet, and
  reporting three parts measured rather than seven part-measured is the honest
  split. The ED rebuild in both columns that part 1 asks for also waits on part 7,
  since it is the same pair of runs.

## 2026-09-30 (A86's table was still v5 against v5 in every figure group)

- **A87 — The direct version check caught it immediately: the spectra differ by
  1.13x at 2550 cm^-1 and the ledger reported the **same number in both
  columns**. Three producer groups write outside `cfg.output_path` and could not
  be isolated, so A86's v3 column carried v5 figure numbers. Fixed, rerun, and
  the honest count is **32 changed, not 23** — including the window value at
  2550 cm^-1, which A86 reported as unchanged.**

  **(1) The check, and what it found.** Straight from the npz, common band:

  | T | rho | kappa(2550) v3 | v5 | kappa_R v3 | v5 | kappa_P v3 | v5 |
  |---|---|---|---|---|---|---|---|
  | 2000 | 1200 | 296.1 | 333.6 (x1.126) | 498.3 | 788.6 (x1.58) | 915.9 | 1296 (x1.42) |
  | 700 | 10 | 14.68 | 1.405 (**x0.096**) | 28.65 | 2.74 (x0.096) | 423.1 | 486.9 (x1.15) |
  | 1000 | 10 | 5.877 | 6.282 (x1.07) | 11.67 | 8.663 (x0.74) | 379.1 | 445.5 (x1.18) |
  | 3000 | 2228.71 | 4426 | 8246 (x1.86) | 4094 | 7584 (x1.85) | 4183 | 8052 (x1.93) |

  Against that, A86's ledger gave
  `figures.fig1.b.marked_kappa_first_principles_cm2_g` = **333.5783 in both
  columns**. The v3 column was the v5 number. After the fix it reads
  **296.139 -> 333.578, x1.1264** — the npz ratio exactly.

  **(2) How each producer finds its spectra, and where it writes.** Every
  producer locates spectra through the **registry**, via `cfg.fp_opacity_path`
  and `FP_OPACITY_SUBDIR`; none hard-codes a spectra path. The fault was on the
  *output* side, and there were three:

  | | reads spectra | writes to | isolated before A87? |
  |---|---|---|---|
  | scan, evolution, crossover, band diagnostic, text checks, scaling, wing cut-off | registry | `cfg.output_path` | yes |
  | all figures and ED items/tables | registry | **`figs/`**, hard-coded in `figures/_manuscript.py` | **no** |
  | release tables and fit | registry | **`data/release/`**, hard-coded | **no** (A86 caught this one) |
  | the ledger's own `figures.*` block | — | — | **no**: read `ROOT / "figs"` directly |

  `FIGS_DIR` (and `ED_DIR`, which derives from it), `RELEASE_DIR`, and the
  ledger's figure-root now all honour `SUBNEPTUNE_FIGS_DIR` /
  `SUBNEPTUNE_RELEASE_DIR`, defaulting to the repo paths so ordinary runs are
  unchanged. The v3 driver sets all three.

  **And the version guard the brief asks for.** The registry entry now carries
  the marker *value* expected for the directory it points at —
  `fp_opacity` -> `cm2_per_g:A50` or `:A74` (A74's four corner cells were stamped
  later), `fp_opacity_v4` -> `:A80`, `fp_opacity_v5` -> `:A83` — and
  `opacity/loaders.py` **raises** when a file's marker is not among them.
  Verified three ways: v5 loads, v3 loads once the subdir is flipped (the
  expectation follows it), and a deliberate mismatch refuses.

  **(3) The corrected table.** Both columns rebuilt from every producer,
  including the figures, into isolated roots. **697 entries, 141 scalars
  comparable, 32 changed (22.7%).**

  | group | scalars | changed |
  |---|---|---|
  | window value / factor at 2550 cm^-1 | 4 | **2** |
  | boundary disputes / thresholds / shifts | 5 | **1** |
  | candidate fiducial (S6) | 18 | 18 |
  | water worlds + old fiducial (S9) | 3 | 3 |
  | S6 closure | 3 | 2 |
  | crossover agreement | 6 | 6 |
  | wing cut-off ratios | 14 | **0** |
  | other figures | 7 | **0** |
  | dilute ranges, scaling indices, conservative bound, release/fit | — | **0** |

  **What A86 got wrong**, now that the figure groups are honestly compared:

  | number | v3 rebuilt | v5 | A86 said |
  |---|---|---|---|
  | `fig1.b` kappa at 2550 cm^-1 | **296.139** | **333.578** (x1.126) | unchanged |
  | `fig1.b` marked ratio | **31.103** | **35.035** (x1.126) | unchanged |
  | `fig3.a` boundary shift (dex) | **-0.68995** | **-0.572077** (x0.829) | unchanged |

  Everything A86 reported about the candidate fiducial, the water worlds and the
  closure stands unchanged: the candidate's 18 numbers all fall to 0.59-0.76x at
  T_int = 30 K, the old fiducial's cooling ratio 6.750 -> 4.052, and **the two
  wmf = 0.50 water worlds still move only 0.4% and 0.1%.** The six
  `agreement.*` entries move from sub-percent to 60-100%, which is the expected
  consequence of comparing v5 against the owner's v3-era RADIS analysis (A83).

  **One group is still not comparable, and it is not the isolation's fault.**
  `figures/ed_items.py` **cannot run to completion**: ED5 panel (e) reads
  `outputs/window_residual.npz` expecting `band_centre_cm` and
  `A_pred_over_A_label`, and the file on disk carries `nu`, `A_label`, `A_pred`,
  `A_resid`, `M_label`, `M_pred`, `N`, `dt`. It fails identically under v3 and
  v5, so it is a pre-existing breakage, not a version effect. The consequence is
  that **neither column can freshly build the ED figures**: the 341 `ed.*`
  entries in the v5 column are inherited from an earlier run whose provenance
  cannot be established from here, and the v3 column has only the three the
  tables supply. **The Extended Data group is therefore excluded from the table
  rather than reported as unchanged**, which is what A86 did with it. A88 should
  regenerate `window_residual.npz` from its producer, or fix ED5 to read the
  format that exists.

  **(4) The band tilt, re-pinned to the monomer test.** The gate asked for
  `stretch/bend` over the line list **> 3**, which is what the v3 spectra
  happened to give rather than what the monomer test predicts: the dipole
  model's squared derivative ratios put it between about **1.2 and 3.0**. v5
  gives **1.588 to 2.283** over the six dilute points — inside that window at
  every one, where v3's values were above it. The gate now pins 1.2 to 3.0 and
  passes; the tilt is the claim and it survives, its magnitude never was. The two
  *monotonicity* assertions are a different claim and stay xfailed.

  **The rho = 49.7 layer changes no quoted planetary number.** It is the
  composite's in-box radiative sliver in the candidate structural case, T =
  854.7-967.9 K, rho = 10.2-49.7, **9.38e-4 of the envelope mass**. Exactly three
  ledger entries derive from in-box radiative quantities, and all three are the
  closure numbers already recorded (`s6_closure.*`). The 22 quoted planetary
  numbers come from `plot_s6_candidate_evolution.py` and
  `plot_s9_case_evolution.py`, which are integrated structure and evolution
  results computed from profiles that **already contain** the layer wherever it
  sits. Nothing needs correcting: what moved is the closure *statement*, which is
  recorded as a number, not any planetary quantity.

## 2026-09-30 (the stuck case was the critical point; the definitive table, at last)

- **A86 — The case that blocked A85 stalls on water's critical point: an equation
  of state that is smooth pointwise but violently stiff, at a tolerance no step
  can meet. It is an integrator limitation, the solution is well determined, and
  one bounded retry fixes it. With that, both ledger columns were rebuilt from
  their producers and **23 of 197 comparable numbers changed — every headline
  group is untouched.**

  **(1) The stuck case.** M_c = 15, wmf = 0.20, T_eq = 500 K. First, it is not
  the case that fails but **one node on its cooling track**: of the 44 T_int
  nodes from 12 to 450 K, **exactly one fails — T_int = 12.00 K, the cold
  endpoint — and only under the composite.** The line list converges at all 44.
  A85 reported the case as failing at T_int = 30 K; it converges there in 344
  right-hand-side evaluations.

  Where it stops, from an instrumented right-hand side:

  | | value |
  |---|---|
  | T | **647.391 K** |
  | P | **22.343 MPa** |
  | rho | 405.8 kg m^-3 |
  | r | 1.630e7 m |
  | m | 1.120e26 kg |
  | grad_ad | 0.0867 |

  Water's critical point is **T_c = 647.096 K, P_c = 22.064 MPa**. The stall sits
  **0.30 K (0.046%) and 0.28 MPa (1.3%) away from it.** What triggers the A57
  guard is not divergence but **step collapse**: consecutive evaluations move the
  mass coordinate by 8.6e17 kg out of 1.1e26, a relative step of **8e-9**, with P
  creeping 7 Pa and T 2e-5 K at a time, until the 200,000-evaluation budget runs
  out.

  **What it is not.** AQUA reports `two_phase = False` at the stall and on both
  sides of it, so it is not a phase boundary; `grad_ad` varies smoothly through,
  so it is not a layer flipping between radiative and convective; the state is
  interior to the table, not on an edge. **What it is:** the critical region's
  compressibility. Sweeping P at the stall temperature, `dln(rho)/dln(P)` runs
  18.66, 3.73, 2.73, 2.37, 1.80 across 21 -> 25 MPa; sweeping T at the stall
  pressure, rho goes 506 -> 155 kg m^-3 between 640 and 660 K. The equation of
  state is single-valued and continuous but its derivatives are enormous, and at
  `rtol = 1e-8` no step is accepted.

  **It is the tolerance, and nothing else.** Fails under **LSODA, BDF and Radau**
  alike, and at **1x, 2x, 5x and 10x** the evaluation budget with identical ~23 s
  wall time each. But at `rtol = 1e-7, 1e-6, 1e-5` it converges in 1-2 s, and the
  answer is stable: R_p = **2.5818, 2.5821, 2.5820 R_earth** — a spread of 1e-4
  R_earth, **0.01%**, over three decades of tolerance.

  **Plain statement: an integrator limitation in the critical region, not physics
  and not a broken equation of state.** So a fix is warranted, and the case is
  neither dropped nor swapped.

  **The fix, and why it sits where it does.** `build_structure` retries **once**
  at ten times looser `rtol`/`atol`, and only after the budget has actually been
  exhausted. Loosening globally would move every number in the repository;
  retrying only the *interior* integration does **not** work, because the outer
  Guillot solve shares the tolerance and a looser outer solve gives a different
  match state — the interior-only retry failed where the global loosening
  succeeded, which is how the retry ended up at build level. Cases that converge
  today never reach the path and are bit-identical.
  `StructureResult.loosened_tolerance` records a model that needed it, so one is
  never silently mixed with the rest. Verified: the stuck node now converges to
  R_p = 2.5818 with the flag set; a node that always worked returns
  2.626818 with the flag clear.

  **(2) Atomic writes.** `atomicio.py` provides `write_text` and an `open_write`
  context manager that writes to a temporary file beside the target and
  `os.replace`s it on clean exit only; on an exception the temporary is removed
  and the previous good file is untouched. Applied to the **four streaming CSV
  producers** — the scan and the three evolution figures, which write row by row
  over minutes and are where A85's partial file came from — and to **18 tool JSON
  writers**, so there is one rule rather than two. A failed producer can no
  longer leave a plausible-looking partial that reads as current.

  **(3) The S6 closure gate, restated on its claim.** The 20% is retired. It was
  never derived: commit `0ed64e7` measured 377 ppm against the line list's 3025
  ppm and put 20% above that 12.5% as margin. The gate now pins the claim the
  figure actually supports — the composite's in-box radiative mass is **strictly
  positive** (a shrinkage, not an elimination) and **strictly below the line
  list's** (there is a shrinkage at all) — and the measured ratio is recorded in
  the ledger as `s6_closure.composite_over_hitemp`, to be quoted rather than met.

  One thing fell out of this that is **not** the same claim. The test also
  asserted that whatever survives sits below the blend window's top, rho = 30,
  where the composite is the line list by construction. Under v5 the sliver spans
  **rho = 10.2 to 49.7**, so part of it lies where the composite is fully first
  principles. That is a separate statement about *where* the residue sits rather
  than *how big* it is, so it has been split into its own test and left as an
  xfail naming A86, not folded into the restated gate.

  **(4) The labelling.** `water_world_fiducial.*` is renamed
  **`candidate_fiducial.*`**, and the case string corrected from "M_c = 10
  M_earth, water mass fraction 0.50, T_eq = 500 K" — an S9 case — to what S6
  actually runs: **M_c = 5, f_env 0.02/0.05/0.10, T_eq = 700 K**. A
  `candidate_fiducial.producer` entry now names the script and config. Where each
  planetary quantity comes from:

  | ledger prefix | producer | case | config |
  |---|---|---|---|
  | `candidate_fiducial.*` | `figures/plot_s6_candidate_evolution.py` | candidate, M_c = 5, f_env 0.02-0.10, T_eq = 700 K | **`config/candidate.yaml`** |
  | `evolution.M_c10_wmf0.50_*`, `evolution.M_c15_wmf0.50_*` | `figures/plot_s9_case_evolution.py` | the wmf = 0.50 water worlds, T_eq = 500 K, T_int = 30 K | `config/fiducial.yaml` |
  | `evolution.old_fiducial_*` | same | M_c = 5, wmf = 0.05, T_eq = 500 K, T_int = 20 K | `config/fiducial.yaml` |
  | `s6_closure.*` | `figures/plot_s4_domain_scan.py --extended` | the candidate cell of the 108-cell scan | `config/fiducial.yaml` |

  There is **no M_c = 5, wmf = 0.50 water world**: the scan's three strongest
  cells are M_c = 15 and 10 at wmf 0.50 and M_c = 15 at wmf 0.20 (the once-stuck
  one), and the only M_c = 5 case S9 carries is the old fiducial at wmf = 0.05.

  **A fault in my own rebuild driver, found by counting.** The first rebuild
  produced a 356-entry ledger where the stale snapshot had 693. `ed_index.py`
  **merges** its numbers into `outputs/manuscript_numbers.json` under `ed.*`, so
  the ledger has to be built *before* it; my driver ran the ledger last and
  overwrote 341 Extended Data entries. Both drivers now run `ed_tables` ->
  ledger -> `ed_index`, with the reason in a comment, and both columns come to
  **697 entries including all 341 `ed.*`**.

  **(5) The definitive table.** Both columns rebuilt from their producers — v5
  live, and v3 by pointing `FP_OPACITY_SUBDIR` at the v3 spectra with
  `output_dir` redirected to `outputs_v3/`, so the v3 side is a genuine v3 value
  and not, as in A83, whatever file happened to be on disk.

  **One leak in that isolation, caught by the gates.** `release_tables.py` writes
  to `data/release/` through `Path(__file__).parents[1]`, a hard-coded repo path
  rather than `cfg.output_path`, so the v3 pass overwrote the tracked v5 release
  products and `test_release_table_is_the_models_own_table` and
  `test_release_fit_matches_its_own_residual` failed. Rebuilding them under v5
  restored both. The v3 driver no longer runs that producer, with the reason
  recorded at the call: the ledger carries no release scalar, so the v3 column
  loses nothing by skipping it. **197 scalars are
  present and numeric in both. 23 changed, 11.7%.**

  **Every headline group is unchanged**, and this is the result:

  | group | scalars | changed |
  |---|---|---|
  | crossover / Planck-Rosseland ratio ranges | — | **0** |
  | window value and factor at 2550 cm^-1 | 7 | **0** |
  | dilute ranges | 1 | **0** |
  | rho^2 scaling indices | 9 | **0** |
  | wing cut-off ratios | 16 | **0** |
  | boundary categories, disputes, thresholds, shifts | 7 | **0** |
  | conservative bound | — | **0** |
  | Extended Data | 40 | **0** |
  | release tables and fit | — | **0** |
  | candidate fiducial (S6) | 18 | **18** |
  | water worlds + old fiducial (S9) | 3 | **3** |
  | S6 closure | 3 | **2** |

  The candidate fiducial, all of it in the same direction — v5 **reduces** every
  difference at T_int = 30 K and leaves T_int = 40 K alone:

  | number | v3 rebuilt | v5 | ratio |
  |---|---|---|---|
  | tau_KH ratio, f_env 0.02, T_int 30 K | **2.2304** | **1.4917** | 0.669 |
  | tau_KH ratio, f_env 0.05, T_int 30 K | 2.0023 | 1.4309 | 0.715 |
  | tau_KH ratio, f_env 0.10, T_int 30 K | 1.8114 | 1.3723 | 0.758 |
  | tau_KH ratio, T_int 40 K (all f_env) | 1.0020 | 1.0005-1.0020 | ~1.000 |
  | delta R / R at 5 Gyr | 0.0055488 | 0.0034695 | 0.625 |
  | delta R / R at 7 Gyr | 0.0066264 | 0.0039328 | 0.594 |
  | delta R / R at 10 Gyr | 0.0076126 | 0.0047743 | 0.627 |
  | delta S / S at 5 Gyr | 0.054678 | 0.035386 | 0.647 |
  | delta S / S at 10 Gyr | 0.086292 | 0.057329 | 0.664 |

  The water worlds and the old fiducial:

  | number | v3 rebuilt | v5 | ratio |
  |---|---|---|---|
  | old fiducial cooling-time ratio | **6.7500** | **4.0520** | 0.600 |
  | M_c = 10, wmf 0.50 cooling-time ratio | 1.00993 | 1.00591 | 0.996 |
  | M_c = 15, wmf 0.50 cooling-time ratio | 1.00227 | 1.00106 | 0.999 |

  **The two water worlds barely move** — a 0.4% and a 0.1% change on a ratio
  within 1% of unity — so the A56 reading that their small differences are the
  real signal survives v5 intact. The old fiducial's cooling-time ratio falls by
  40%, and the candidate's differences by a third.

  The closure: `s6_closure.composite_over_hitemp` **0.3654 -> 0.4770** with the
  line-list side unchanged, so the composite's in-box radiative mass grows by
  1.31x. (These are the scan's candidate cell; the S6 test builds its own models
  at its own parameters and reads 12.5% -> 31% on the same quantity. Both say the
  same thing and the ledger's is the one to quote.)

  **A83's 34-entry diff should be discarded.** Its v3 side was stale and its v5
  side incomplete; the honest count is these 23, and the difference between the
  two lists is what staleness was contributing.

## 2026-09-30 (the ledger now refuses stale inputs; the definitive table is blocked on one case)

- **A85 — The audit found one materially stale producer, and the ledger now
  refuses to run on a stored file older than the spectra it derives from.
  Rebuilding that producer changed which cells S9 selects, and one of the new
  ones **does not integrate under v5** — the same failure at 1x, 2x, 5x and 10x
  the evaluation budget, so it is stuck rather than slow. The v5 ledger and
  therefore the definitive v3-versus-v5 table cannot be completed until that
  case is resolved. Everything that does converge is reported below. No
  manuscript number and no gate changed.**
  `tools/manuscript_numbers.py` (freshness guard), `tools/rebuild_producers.sh`,
  `tools/rebuild_v3_column.sh`.

  **(1) The audit.** The ledger reads **15** stored files and rebuilds none of
  them. That is the mechanism behind A84's finding, and it had a second victim.

  | stored file | status |
  |---|---|
  | `S4_domain_scan_extended.csv` | **materially stale** — written 09-29 10:24, before the v5 switch |
  | `S6_candidate_evolution.csv` | was stale from commit 4834c22 (A84); rebuilt |
  | `crossover_hitemp_v2.json`, `crossover_v3.json`, `band_diagnostic.json`, `text_checks.json`, `rho_scaling.json`, `wing_cutoff.json`, `S8_rcb_flip.csv`, `S9_case_evolution.csv` | fresh — rebuilt during A83's propagation |
  | `broadening_sensitivity.json`, `aqua_entropy_check.json`, `aguichine2025_compare.json`, `crossover_v3_before_a60.json`, `crossover_v3_after_a60.json` | old but **exempt**: each has **zero** first-principles references, so the v5 switch cannot move them |

  **The stale one matters most.** `plot_s4_domain_scan.py` writes the extended
  108-cell scan only when passed `--extended`, which costs about five minutes and
  which A83's propagation never passed — it ran the script bare. So the entire
  boundary group, categories included, was still the v3-era scan while
  everything around it had moved to v5.

  **The guard.** Every `_csv`/`_json` read now compares the stored file's mtime
  against the newest spectrum under `cfg.fp_opacity_path` and the HITEMP
  directory, and **raises** if the file is older. It is an mtime check and not a
  dependency graph — it catches exactly the failure that bit us, a spectra switch
  without a producer rerun, and nothing subtler. The five FP-independent files
  are exempt by name, each with its reason in the source. Verified: the ledger
  refuses on the stale scan and names it.

  `tools/rebuild_producers.sh` runs all eleven producers in dependency order
  before the ledger, with the S6 config trap commented in place so the A83/A84
  mistake cannot recur silently.

  **The first version of the guard was not enough, and finding that out is part
  of the result.** Checking each stored file against the *spectra* catches a
  spectra switch, but some producers depend on **another producer's output**: S9
  and S8 take their case list from the extended scan, so rebuilding the scan
  invalidates them while the spectra have not moved. The first guard passed a
  ledger that mixed boundary numbers from the new scan with evolution numbers
  computed against the old one. A small `_OUTPUT_DEPS` map now covers that, and
  it was checked against the source rather than guessed: `plot_s8_rcb_flip.py`
  and `plot_s9_case_evolution.py` reference the scan, `plot_s6_candidate_evolution.py`
  does **not** (it takes a fixed f_env grid from `config/candidate.yaml`), so S6
  is deliberately not listed — an entry for it produced a spurious refusal and
  was removed.

  **A hole remains that mtime cannot close.** A producer that fails partway
  leaves a *partial* output with a *new* mtime, which looks fresher than its
  inputs and passes every check here. That happened during this amendment:
  diagnosing the S9 failure rewrote `S9_case_evolution.csv` with the **3 cases
  that converge** where it previously held 4, and the guard cannot tell. The
  producers should write atomically — temp file then rename — so a failed run
  leaves no output at all. That is for A86; until then the case count is stated
  wherever S9 is quoted.

  **(2) Blocked, and by a real failure rather than a stale file.** Rebuilding the
  scan changed its ranking, and S9 takes "the 3 strongest scan cells". The third
  is now **M_c = 15, wmf = 0.20, T_eq = 500 K, T_int = 30 K**, which was not in
  the stale scan's top three. It **fails the A57 integration guard**:

  | evaluation budget | result |
  |---|---|
  | 200,000 (default) | fails after 29 s |
  | 400,000 | fails after 23 s |
  | 1,000,000 | fails after 23 s |
  | 2,000,000 | fails after 23 s |

  Identical wall time at ten times the budget means the interior integration is
  **not making progress**, which is precisely what the guard exists to surface.
  It is not a tuning problem. So `plot_s9_case_evolution.py` cannot complete, the
  `evolution.*` ledger group cannot be rebuilt under v5, and the full ledger — and
  with it the definitive v3-versus-v5 table the brief asks for — **is not
  produced here.** Excluding the case would silently change the case set and make
  the comparison the same kind of apples-to-oranges A84 was cleaning up, so it is
  reported instead.

  Note what this is *not*: a regression. With the stale scan S9 picked a
  different third case and ran. The failure appeared because the scan was
  un-staled, so it belongs to v5 and has simply never been exercised before.

  Steps 1-3 of the v5 rebuild (extended scan, corner boundaries, S8) completed;
  the extended scan itself converged at all 108 cells. The v3 column driver is
  written and untested, since it would stop at the same step. `S9_case_evolution.csv`
  currently holds **3 of its 4 cases** — the two wmf = 0.50 water worlds and the
  old fiducial — with the M_c = 15, wmf = 0.20 case absent because it does not
  integrate. Any `evolution.*` number read from it is a three-case number.

  **(3) The planetary numbers, under v5, each with its config named.** The two
  water worlds converge and reproduce the stored values exactly, so these are
  genuine v5 numbers.

  `figures/plot_s9_case_evolution.py`, **`config/fiducial.yaml`** (dR and ds are
  composite minus line list, fractional; 44/88 nodes bracket the age-grid
  sensitivity):

  | case | dR 1 Gyr | dR 5 Gyr | dR 10 Gyr | ds 5 Gyr |
  |---|---|---|---|---|
  | M_c=10 wmf=0.50 T_eq=500 T_int=30 | -5.9e-08 / -1.1e-07 | 3.20e-05 / 4.00e-05 | 8.37e-04 / 9.17e-04 | 4.07e-04 / 4.59e-04 |
  | M_c=15 wmf=0.50 T_eq=500 T_int=30 | -1.1e-07 / 5.5e-08 | -3.51e-06 / 9.76e-06 | 5.38e-04 / 5.71e-04 | 2.42e-05 / 1.43e-04 |
  | old fiducial M_c=5 wmf=0.05 T_eq=500 T_int=20 | 3.08e-03 / 3.16e-03 | 9.80e-03 / 9.73e-03 | 1.24e-02 / 1.23e-02 | 1.67e-01 / 1.62e-01 |

  **There is no M_c = 5, wmf = 0.50 water world to report**: the scan's three
  strongest cells are M_c = 15 and 10 at wmf 0.50 and M_c = 15 at wmf 0.20, and
  the only M_c = 5 case S9 carries is the old fiducial at wmf = 0.05. The brief's
  "M_c = 5 and 15 water worlds" exists only as the M_c = 15 one.

  `figures/plot_s6_candidate_evolution.py`, **`config/candidate.yaml`**
  (M_c = 5, f_env 0.02/0.05/0.10, T_eq = **700 K** — not 500):

  | quantity | 5 Gyr | 7 Gyr | 10 Gyr |
  |---|---|---|---|
  | delta R (R_earth) | 0.006102 | 0.006904 | 0.008363 |
  | delta R / R | 0.003470 | 0.003933 | 0.004774 |
  | delta S (J/kg/K) | 343.34 | 413.56 | 531.14 |
  | delta S / S | 0.035386 | 0.043494 | 0.057329 |

  | tau_KH ratio | f_env 0.02 | 0.05 | 0.10 |
  |---|---|---|---|
  | T_int = 30 K | 1.4917 | 1.4309 | 1.3723 |
  | T_int = 40 K | 1.0020 | 1.0005 | 1.0020 |

  **A third mislabelling, and it is in the ledger.** These S6 rows are stored
  under the key prefix `water_world_fiducial.*`, and the ledger asserts the case
  is "M_c = 10 M_earth, water mass fraction 0.50, T_eq = 500 K". That is an **S9**
  case. S6 runs the candidate at M_c = 5, f_env 0.02-0.10, T_eq = 700 K. The
  label and the numbers under it have never matched. Every quantity the
  manuscript takes from `water_world_fiducial.*` is a *candidate* number, not a
  water-world one.

  **(4) The S6 closure, diagnosed.** What it compares: the envelope mass fraction
  in layers that are both **radiative** and **inside the first-principles box**
  (T >= 650 K, rho >= 10 kg m^-3), composite against line list. The claim is that
  handing kappa to the first principles inside its own box removes the radiative
  layers the line list leaves there.

  **The 20% was chosen, not derived.** Commit `0ed64e7`, which set the blend
  window to 10-30 and wrote the gate, records the measurement in its own message:
  *"The candidate's closure becomes partial: 3025 ppm in-box under HITEMP, 377 ppm
  under the composite, at the box's dilute edge."* That is **377/3025 = 12.5%**,
  and 20% was put above it as margin. Under v5 the line-list side is unchanged at
  3025 ppm and the composite side is **938 ppm**, so the ratio is **31%** — the
  composite's in-box radiative mass has grown 2.5x and consumed the margin.

  **Which of the three it is: the physics, read through a chosen threshold.** The
  gate's own docstring already says the closure is partial because the box floor
  (A55) reaches down to rho = 10, into the range where the composite is the line
  list *by construction* — the blend window is 10-30 with its seam at rho = 17.3.
  A84 measured what changed there: kappa_R first-principles over line list at the
  seam fell from 5.74-9.67 to 2.32-2.82. A less opaque composite at the dilute
  edge keeps more mass radiative, which is exactly what 377 -> 938 ppm is. So it
  is **not** a method error and **not** a gate error in the sense of being wrongly
  written; it is the physics of v5 moving a quantity past a bound that was set
  with 1.6x of margin over a 12.5% measurement. The qualitative claim the gate
  supports — shrinkage, not elimination — still holds at 31%. Whether 20% should
  become a looser number or the claim should be restated is Harry's call, and the
  xfail stays until he makes it.

  **What A86 needs.** Decide the `M_c = 15, wmf = 0.20` case: whether its
  non-convergence is a physical pathology (a two-phase envelope, a hull-clamp
  cycle) or an integrator limitation. Until then no v5 ledger and no v3
  comparison column is defensible, and the 34-entry diff A83 quoted should not be
  quoted at all — its v3 side is stale and its v5 side is incomplete.

## 2026-09-30 (checking v5: the 6.07 was mine, and the v3 baseline was stale)

- **A84 — The tau_KH jump A83 reported as 1.00 -> 6.07 was an artefact of a
  bracket I widened, and the v3 number it was compared against was stale by
  eleven commits. Measured properly, same code and same config, the change is
  **2.230 -> 1.492** — the ratio *falls* by a third. The solve machinery is
  sound, the composite has no discontinuity at any seam, and the remaining
  headline claim (kappa_R above the line list everywhere) holds. No gate and no
  manuscript number was changed here; the six xfails stay.**
  Read-only checks; `figures/plot_s6_candidate_evolution.py` reverted to its
  default entropy bracket and its outputs regenerated with the config its own
  docstring specifies.

  **(1) FP over the line list, all 52 points under v5.**

  | quantity | min | median | max | points below 1 |
  |---|---|---|---|---|
  | kappa_R bare | 1.573 | 5.987 | 96.96 | **none** |
  | kappa_R + MT_CKD | 1.570 | 3.522 | 71.23 | **none** |
  | kappa_P bare | 1.246 | 2.158 | 34.11 | **none** |
  | kappa_P + MT_CKD | 0.213 | 1.760 | 30.24 | **13** |

  The Rosseland claim is unqualified: **the first-principles opacity exceeds the
  line list at every one of the 52 points**, by a median of 6.0 bare and 3.5 with
  the continuum on the line-list side. The Planck mean does not: with MT_CKD
  included, 13 points fall below 1, all of them cold and dense — (650, 500) at
  **0.213**, then (650, 300) 0.302, (650, 150) 0.390, (800, 500) 0.447,
  (650, 100) 0.539, (800, 300) 0.570, (800, 150) 0.693, (700, 60) 0.695,
  (1000, 500) 0.817, (1000, 300) 0.877, (800, 100) 0.884, (1000, 129.744) 0.982
  and (700, 30) 0.992. Where the continuum dominates the Planck mean it exceeds
  what the first-principles spectra produce, by up to a factor of five.

  **(2) The tau_KH trace, and what A83 got wrong.** Two errors, both mine.

  *The bracket.* A83 widened `build_at_entropy`'s T_int floor from 20 K to 12 K
  because the solve failed. It failed because the script was run with
  `config/fiducial.yaml`, whose first irradiation point is T_eq = 500 K, where
  its own docstring specifies `config/candidate.yaml` at T_eq = 700 K. With the
  right config the solve lands at **T_int = 27.15 K**, comfortably inside the
  original bracket, and the widened bracket is what let the 500 K case run away
  to 6.07. The bracket is reverted and the outputs regenerated.

  *The baseline.* The v3 number, 1.001989, carries commit **4834c22** in the
  ledger, and `manuscript_numbers.py` reads `outputs/S6_candidate_evolution.csv`
  from disk without regenerating it. That CSV had not been rebuilt since
  4834c22 — before the corner cells (A74) and everything after. **It was not a
  v3 number at all.** Running the same script against the v3 spectra now gives
  **2.230**, not 1.002. The 1.002 is, to four figures, what the identity control
  in (3) returns, which is consistent with the composite having been
  indistinguishable from the line list at that case when the CSV was written.

  **The sound comparison, same code and same config, v3 spectra against v5:**

  | | v3 | v5 |
  |---|---|---|
  | composite T_int (K) | 24.55 | **27.15** |
  | tau_KH ratio | **2.230** | **1.492** |
  | kappa_R FP/LL at the seam | 5.74-9.67 | **2.32-2.82** |
  | kappa_P FP/LL at the seam | 1.32-1.34 | 1.49-1.55 |
  | max radiative rho | 20.4 | 22.9 |

  and the mechanism is monotone: v5's first-principles opacity is a smaller
  multiple of the line list **at the seam** (rho = 17.3, where the blend lives)
  than v3's was, so the composite is less opaque there, needs a *higher* T_int to
  hold the reference entropy, and its cooling time is correspondingly closer to
  the line list's. My earlier confusion came from comparing the FP/LL ratio at
  (700, 10), which the envelope never visits; what matters is the seam.

  **No discontinuity anywhere.** Both versions give a **single** radiative-
  convective transition over 799 layers with no reversal; 407 layers sit below
  the rho = 10 grid floor and read the A26 flat clamp, identically in both; the
  clamp fraction is 0.5% with **zero** clamped radiative layers and
  `clamp_reliable` true; and the seam ratios march monotonically across the
  10-30 kg m^-3 blend in both. The composite is behaving.

  **(3) Controls, all passing.**

  | control | result |
  |---|---|
  | the line list matched against itself | T_int 29.986 K, ratio **1.0019** (identity, to solver tolerance) |
  | the same, with kappa scaled x2 and x5 on both sides | 1.0019, 1.0017 — scale-invariant, as it must be |
  | line list x1, x2 against the unscaled reference | T_int 27.15 -> 22.82, ratio **1.492 -> 2.986** — monotone |
  | line list x5 | entropy unreachable in the default bracket: a genuine bracket limit at extreme contrast, the same failure A83 met and mis-diagnosed |

  **(4) Plain statement.**

  *The 6.07 was not physical.* It was the widened bracket acting on a case the
  script was never meant to run. The physical number is **1.492**, against a
  properly recomputed v3 value of **2.230** — so v5 moves this quantity **down**
  by a third, not up by a factor of six. The solve is sound: it reproduces the
  identity to 2e-3 and responds monotonically to a scaled opacity.

  *The radius change is real but smaller than A83 said, and A83 mislabelled it.*
  The keys `delta_radius_frac.5.0` and the rest are indexed by **age in Gyr, not
  core mass**; A83's table called the first of them "M_c = 5" and that is wrong.
  With the bracket reverted and the right config, the fractional radius
  difference at 5 Gyr is **3.985e-05 -> 0.00347**, not the 0.0108 A83 reported,
  and the entropy difference **8.2e-4 -> 0.0354**. Both still rest on the same
  stale 4834c22 baseline, so the honest reading is that the v5 values stand on
  their own and the "v3" column should not be quoted until it is recomputed.

  *The S6 closure failure is real and is not any of the above.* It reproduces
  unchanged under the corrected config: the composite's in-box radiative mass is
  9.38e-04 against the line list's 3.03e-03, i.e. **31%** where the gate requires
  under 20%. Nothing about the bracket, the config or the baseline touches it. It
  is a genuine consequence of v5 being less opaque at the seam, and it remains
  for physics review.

  **What this says about the other A83 comparisons.** Every number whose v3 value
  came from a CSV that `manuscript_numbers.py` reads rather than rebuilds is
  suspect in the same way. A85 should regenerate the v3 column by checking out
  the v3 spectra and rerunning the *producers*, not the ledger, before any of the
  34 changed entries is quoted.

## 2026-09-30 (v5 is live: the estimator propagated, and three physics claims moved)

- **A83 — All 52 spectra regenerated with the A81 estimator at sd = 10 cm^-1,
  the registry switched, and the whole downstream list rerun. The band-level and
  scaling results barely move; the structure and evolution results move a great
  deal. Three of the manuscript's physics claims changed materially and **no gate
  encoding a physics claim has been re-pinned** — those are Harry's calls.**
  `opacity/dipole_spectrum.py`, `tools/v5_regenerate.py` ->
  `outputs/v5_regenerate.json`, two gates in
  `validation/test_a83_v5_estimator.py`.

  *Corrected by A84 in three places. The tau_KH ratio of 6.07 was an artefact of
  the entropy bracket widened below, run on a case the script was not meant to
  take; the real change is 2.230 -> 1.492, downward. The v3 column for every S6
  quantity came from a CSV last rebuilt at commit 4834c22 and is not a v3 number.
  And the `delta_radius`/`delta_entropy` keys are indexed by **age in Gyr**, not
  core mass, so the row below labelled "M_c = 5" is at 5 Gyr. The bracket has
  been reverted and those outputs regenerated. See A84.*

  **(1) v5.** Estimator settled: Hann-tapered periodogram, `sum(w^2)`
  normalisation, the harmonic correction once, **continuous `omega^2`**. A82's
  discrete-derivative form is a finite-difference diagnostic and is not used;
  it would move kappa_R by **-0.71% to -1.26%**, and that difference is recorded
  here so the choice is on the record rather than implicit. Smoothing is
  **sd = 10 cm^-1 at every cell**, replacing v3's 20 (and 15 at four points).

  Each v5 file now carries its own error bars, `kappa_R_sigma` and
  `kappa_P_sigma`, from the two 15 ps halves with `sigma_30 = |a - b| / 2`:

  | | min | median | max |
  |---|---|---|---|
  | sigma(kappa_R)/kappa_R | 0.1% | **5.0%** | **22.5%** |
  | sigma(kappa_P)/kappa_P | 0.00% | 1.68% | 8.32% |

  **A82's framing of this was wrong and is corrected here.** A82 called it a
  dilute-corner problem; it is not. The least certain cells are
  (2500, 197.586) at **22.5%**, (2750, 222.666) at 18.9% and (2000, 60) at
  18.7% — **mid-density, 2000-2750 K**, where the windows are deep and the
  spectrum is busiest. At several of these the full-record kappa_R lies outside
  the range of its own two halves, which a harmonic mean can do because it is
  noise-biased low, but it means the single-pair sigma is itself only known to
  about a factor of two.

  **v5 against v3, per cell:** kappa_R **0.089 to 1.687** (median 0.807),
  kappa_P **0.960 to 1.677** (median 1.139). The cold cells fall hard — (700, 10)
  from 28.57 to 2.55, a factor of **11** — and the hot ones rise to 1.69x. All
  three v3 defects filled the transparency windows (the spare `Q_harmonic`, the
  leakage, and sd = 20), and kappa_R is a harmonic mean set by the minima, so
  removing all three lowers it most where the windows are deepest.

  **The band-level physics is unmoved**, as it must be, since band integrals are
  insensitive to smoothing: the dilute band medians are **1.043** (rotation),
  **1.483** (bend) and **2.589** (stretch), against A81's v4 values of 1.044,
  1.515 and 2.658 and A79's prediction of 1.028, 1.539 and 2.659.

  **(2) Propagation.** Registry switched via a new `FP_OPACITY_SUBDIR` in
  `config/schema.py`, so the whole pipeline follows one symbol. The owner's
  independent `crossover_hitemp.json` is a source-independent input and carries
  over; the v3-derived crossover *output* is recomputed. Everything on the list
  ran: band diagnostic, crossover bare and continuum-inclusive, rho^2 scaling,
  wing cut-off, dilute band table, release tables and fit, all 29 figures, the
  ED items and tables, and `manuscript_numbers.json`.

  One intervention was needed. `figures/plot_s6_candidate_evolution.py` matches
  the composite to the line list's interior entropy, and under v5 the composite
  reaches that entropy at **T_int = 19.1-19.8 K**, below the solver's default
  20 K bracket floor, so the solve failed outright. The floor was widened to
  **12 K** — the smallest value that both brackets the root and integrates, since
  at 10 K the interior integration stops converging against A57's guard. The
  solutions are interior roots well away from the new floor, not boundary
  artefacts.

  **Of 693 manuscript entries, 34 changed.** They are not spread evenly:

  | group | changed / total |
  |---|---|
  | dilute band numbers | **0 / 7** |
  | rho^2 scaling indices | **0 / 12** |
  | wing cut-off ratios | **0 / 21** |
  | conservative bound | **0 / 4** |
  | boundary numbers | 2 / 29 |
  | window at 2550 cm^-1 | 2 / 8 |
  | **fiducial / evolution** | **19 / 24** |

  The headline movers:

  | number | v3 | v5 |
  |---|---|---|
  | `fig1.b` kappa at 2550 cm^-1 (cm^2/g) | 296.1 | **333.6** |
  | `fig1.b` marked ratio | 31.10 | **35.04** |
  | `fig3.a` boundary shift (dex) | -0.685 | **-0.572** |
  | `fig3.b` number of boundaries | 8 | **10** |
  | water world tau_KH ratio, T_int = 30 K | 1.002 | **6.074** |
  | water world tau_KH ratio, T_int = 40 K | 0.998 | **1.234** |
  | water world delta R/R, M_c = 5 | 4.0e-5 | **0.0108** |
  | water world delta S/S, M_c = 5 | 8.2e-4 | **0.201** |
  | old fiducial cooling-time ratio | 6.750 | **4.052** |

  **(3) The fork, and why A83 stops at it.** Three claims moved materially.

  1. **The water-world fiducial.** The tau_KH ratio at T_int = 30 K goes from
     1.00 to **6.07**, and the arithmetic is transparent: the entropy-matched
     composite sits at T_int = 19.1 K against the line list's 30 K, and
     `tau_KH ~ T_int^-4` gives `(30/19.1)^4 = 6.1`. So the number is a
     consequence of the matched solve moving, not of a bug — **but the direction
     needs checking, because v3's first-principles kappa_R was 42x the line
     list at (700, 10) against v5's 3.8x, so v3's composite should have needed
     an even lower T_int and yet it solved at ~30 K.** Either v3's solve was
     wrong or the composite's blend does something I have not traced. This is
     the single thing A84 should settle first.
  2. **The S6 closure claim.** The composite's in-box radiative mass is
     9.38e-04 against the line list's 3.03e-03 — 31%, where the gate required
     under 20%. The composite no longer closes the radiative zone as tightly.
  3. **The dilute band tilt (A61/A65).** `stretch/bend` for the first principles
     over the line list falls from above 3 to **1.61**, and the ratio is no
     longer monotone in band centre at (1000, 10). The tilt survives in
     direction but not in magnitude.

  Also for review: the cold-cell kappa_R excess over the line list is **3.5 to
  8.4x** bare and **2.0 to 5.4x** with MT_CKD on the line-list side, so the
  continuum accounts for roughly a third to a half of it and the window claim
  should be quoted continuum-inclusive.

  **Gate disposition.** Mechanical changes made; physics claims left failing.

  | gate | action |
  |---|---|
  | `test_a80_v4_gates::test_nothing_has_been_propagated` | **inverted** to `test_v5_is_the_live_source` |
  | `test_sum_rule::..._still_violate_it` | **inverted** to `..._now_satisfy_the_sum_rule`, \|p\| < 0.3 |
  | `test_s7_dilute_floor` (3 tests) | **retired** by module skip: A81 showed the white-noise-in-M model wrong in kind. File kept for the trail |
  | `test_crossover_v2::test_agreement_with_the_delivered_analysis` | **restated**: the delivered analysis is v3-era, so the pin is now that v5 differs by > 10% |
  | `test_opacity::test_known_grid_point_means` | **re-pinned** 489.8 -> 565.6 |
  | `test_fp_units` baseline kappa_R | **re-pinned** 35.409 -> 17.498 |
  | `test_s6_candidate`, `test_band_diagnostic` (3), `test_crossover_v2::test_spectral_contrast`, `test_scaling_and_release::test_where_the_pair_estimate_is_not_conservative` | **left failing as xfail**, each naming A83 and awaiting physics review |

  **ED Table 1 gains two uncertainty columns**, `kappa_P_sigma_cm2_g` and
  `kappa_R_sigma_cm2_g`, at all 52 rows. One care: the stored sigma is absolute
  on the pipeline's own 100-5000 cm^-1 Rosseland mean, while the table's kappa_R
  is band-limited to 100-4000, and the two differ by about 10%. The **fractional**
  error is what transports between band definitions -- it is a property of the
  record length, not of the integration limits -- so the table takes the fraction
  from the file and applies it to its own mean, which reproduces v5's
  0.1-22.5% (median 5.0%) exactly. Quoting the absolute sigma against the
  band-limited mean would have been wrong by that 10%.

  **The new gate** (`validation/test_a83_v5_estimator.py`) pins what the whole
  A78-A81 hunt bought: the sum rule holds to **1%** (measured 0.09%, against a
  taper-consistent reference), and the leakage floor stays **100x** below the
  transparency-window minimum at every one of the 52 cells (measured 476x at the
  worst, (3000 K, 108.796)).

## 2026-09-30 (the smoothing width cannot be chosen from a 30 ps record)

- **A82 — No width satisfies both of the brief's criteria, and the reason is not
  the width. sd = 10 cm^-1 meets them at five of the six cells; the sixth,
  (1000 K, 10 kg m^-3), carries a 4.4% smoothing bias against an **8.0%
  statistical uncertainty**, so the 2% criterion is being applied below the
  noise. The 30 ps records cannot determine kappa_R to 2% anywhere in the dilute
  corner. Reporting and stopping before propagation, as the brief directs.**
  `tools/smoothing_width.py` -> `outputs/smoothing_width.json`. Nothing
  propagated; the registry is still on v3.

  **(1) The study.** Six cells, seven widths, A81's estimator throughout
  (Hann-tapered periodogram, `sum(w^2)` normalisation, one harmonic correction)
  with the **discrete-derivative prefactor** `[2 sin(w dt/2)/dt]^2` the brief
  specifies. Four quantities per cell and width:

  | quantity | how |
  |---|---|
  | as delivered | kappa_R of the smoothed tapered periodogram |
  | noise-corrected | x `nu/(nu-2)`, since kappa_R integrates `1/kappa` and `E[1/S_hat] = (1/S) nu/(nu-2)` |
  | statistical uncertainty | from the two 15 ps halves |
  | smoothing bias | against a datum at sd = 2 |

  with `nu = 6 . 2 s sqrt(pi) / (1.5 dnu_fft)` — 6 for three Cartesian
  components of a chi^2_2 periodogram, the Gaussian kernel's effective bin
  count, and 1.5 for the Hann taper's noise-equivalent bandwidth. `nu` is the
  same for every cell (1.0159 at sd = 10) because every record is 60002 frames
  at 0.5 fs. The analytic `nu` is used only for the bias it predicts; every
  error bar is the empirical half-record one.

  **Two corrections to how those were first computed.** The half-record spread
  is `sqrt(2) sigma_15 / kappa` and `sigma_30 = sigma_15/sqrt(2)`, so the
  uncertainty on the full record is the spread **divided by two**, not the
  spread. And the Welch-8 reference, though nearly noise-free, has **5.0 cm^-1**
  intrinsic resolution against the full record's 1.11, so its sd = 2 value is
  already smoothed to an effective ~5.4 cm^-1: **it cannot certify a 2% bias for
  any width near 5, whatever the data say.** The verdict below therefore uses
  the *full-record noise-corrected* value at sd = 2 as the bias datum, which is
  properly resolved at the cost of being noisy; Welch is kept as a cross-check
  and gives the same ordering.

  **(2) The verdict.** Worst case over the six cells:

  | sd (cm^-1) | worst step | worst bias | worst statistical | cells passing |
  |---|---|---|---|---|
  | 5 | 2.2% | 2.2% | 8.3% | 5/6 |
  | 7.5 | 1.2% | 3.4% | 8.1% | 5/6 |
  | **10** | **1.0%** | **4.4%** | **8.0%** | **5/6** |
  | 15 | 1.7% | 6.2% | 7.7% | 5/6 |
  | 20 | 1.8% | 8.1% | 7.5% | 3/6 |
  | 30 | 4.5% | 13.0% | 7.5% | 2/6 |
  | 40 | 5.8% | 19.5% | 7.5% | 2/6 |

  **No width passes at all six cells, so A82 stops.** But the failure is one
  cell, not a trend: at sd = 10 the bias is 0.7-1.5% at (700, 60), (1500, 300),
  (2000, 1200), (2500, 500) and (3000, 2228.71), and **4.4% at (1000, 10)**
  alone — the grid's most dilute corner, where the transparency windows are
  deepest and the harmonic mean is most exposed to them. Widths of 5 to 15 all
  pass five of six; only at 20 and above does the bias spoil the denser cells
  too.

  **The binding constraint is the record length, not the width.** The
  statistical uncertainty on kappa_R from a 30 ps record is **1.1% to 7.9%,
  median 6.2%**, and 8.0% at the failing cell — three to four times the
  tolerance being tested. Reaching 2% needs **9.5x the record: ~285 ps at the
  median cell and ~463 ps at the worst**, against the 30 ps in hand. No choice
  of smoothing width can recover what the record does not contain.

  **The error model closes, which is the positive result.** At (1000, 10) the
  noise-corrected kappa_R rises **19.5%** from sd = 2 to sd = 40, and the
  independently measured Welch smoothing bias over the same range is **18.4%**.
  The two agree to 1.1%, so the whole width dependence of kappa_R is smoothing
  bias plus a well-modelled chi^2 noise term and nothing else. **The 30 ps
  record's spectral *shape* is sound; only its *length* is limiting.** That is a
  different and much better position than A80 or A81 were in.

  **(3) Not done.** Propagation, v5, the registry switch and the whole downstream
  list are untouched, per the brief's instruction to stop if (2) fails. The
  registry still points at v3 and the A80 gate that asserts it still passes.

  **What A83 should do, and the choice is Harry's.**
  (a) Extend the dilute trajectories to ~300-500 ps — the only route to the 2%
  the criteria ask for, and an ARC job of about ten times the original cost for
  the dilute cells alone (the dense cells are already at 1-4%).
  (b) Adopt **sd = 10 cm^-1** and propagate with the systematic stated: +4.4%
  smoothing bias worst case, 8% statistical at the dilute corner, both quoted
  beside A71's wing cut-off. sd = 10 is the widest width that is stationary to
  1.0% and passes five of six cells, and it is half the delivered 20.

  **One inconsistency to settle either way.** The brief specified the
  discrete-derivative prefactor; A81's estimator and the v4 files use the
  continuous `w^2`. The two differ by **-0.71% to -1.26% on kappa_R** across the
  six cells, always negative. It is small, but v5 must pick one and every
  downstream number must come from the same choice.

## 2026-09-30 (there is no noise floor: it was spectral leakage)

- **A81 — The floor that held A80 is not dipole-model noise and not
  discontinuities. It is the untapered periodogram's own leakage: the record's
  end-to-end mismatch acts as a step, leaks as omega^-2, and the omega^2 of the
  absorption formula turns that into a *flat floor in sigma* straight across the
  transparency windows. A Hann taper on the time series removes it. The margin
  goes from 1.8x to **1243x** and all three A80 gates pass.**
  `tools/floor_identify.py` -> `outputs/floor_identify.json`; the estimator in
  `opacity/dipole_spectrum.py` now tapers; `ed_data/arc_dipoles_sha256.json` is
  the tracked fingerprint of the bundle. **Nothing propagated.**

  **(1) Not discontinuities.** Across all 52 cells and **3,109,989** consecutive
  frame steps there is **not one** outlier above 10 sigma. The largest |dM|
  anywhere is **0.199** of the box length, against the ~0.4-0.8 L a wrapped
  atomic charge would give — at the dilute cells that would be a 450-sigma
  event, and there are none. The magnitudes do not cluster at multiples of
  q x L; the largest step in each cell sits at a median of **5.53** MAD-sigma
  (range 4.20-7.59), against the ~5.1 expected for the maximum of 60000 draws
  from a chi distribution with three degrees of freedom — a slightly heavier
  tail than Gaussian, worth noting, but nothing like a wrap. The dipole model
  is not being evaluated across a periodic image, so there is nothing to unwrap.

  **(2) The tail exponent, and the trap in measuring it.** Fitting
  `sigma ~ nu^m` above 6000 cm^-1 gave a slope that *rose* with frequency —
  -0.21, +0.26, +0.67, +1.36 across the four bands to Nyquist at (700, 30) —
  which is not a power law at all. The cause was the prefactor: writing
  `sigma ~ omega^2 S` uses the **continuous** omega^2 on a series sampled at
  0.5 fs, where the right factor is `[2 sin(omega dt/2)/dt]^2`. Under the
  discrete factor the same bands give **-0.30, +0.01, +0.05, +0.04** — flat.
  Across all 52 cells the near-Nyquist discrete exponent has median **-0.04**.

  Flat in the discrete sigma means `S ~ omega^-2`, which is the brief's *step*
  signature — and (1) says there are no large steps, so it had to be many small
  ones. But dM is **smooth**: its lag-1 autocorrelation is 0.91-0.98, with only
  0.2-4% of its power in the flat tail. Something was producing an omega^-2
  tail without any step in the data.

  **The step is the record boundary.** A discrete Fourier transform treats the
  record as periodic, so the mismatch between the first and last frame *is* a
  step discontinuity, once per record, and it leaks as omega^-2 — exactly the
  observed shape. Two tests confirm it and exclude the alternatives:

  - **It scales with the boundary mismatch, which is the proof rather than the
    inference.** Leakage from a step of size `Delta = M[last] - M[first]` goes as
    `Delta^2 / n`. Across all 52 cells the measured untapered floor tracks that
    prediction with a log-log correlation of **r = 0.973** and a slope of
    **0.896** against the exact 1.0 — one mechanism, over four orders of
    magnitude in floor level.
  - **It falls with record length.** Halving and quartering the series raises
    the apparent floor (60002 frames against 15000: **0.06 to 0.69**). Genuine
    additive model noise is a bias in the expected spectrum and would not move;
    leakage from a single boundary mismatch dilutes as the record lengthens.
  - **It is not the thermostat.** Production is **NVE VelocityVerlet**
    (`run_prod_md_05fs.py` line 50); Langevin appears only in equilibration at
    lines 47-48. And it is not float32: the implied step is 4000-25000x the
    float32 ulp of the stored dipoles.

  **The repair is a taper, and it is exact.** Applying a Hann window to the time
  series and normalising by `sum(w^2)` — the normalisation that keeps Parseval
  exact for a tapered record — drops the high-frequency floor by **three to five
  orders of magnitude** (ratios 1e-5 to 0.08) while leaving every band integral
  unchanged:

  | gate | A80 (untapered) | A81 (tapered) |
  |---|---|---|
  | (i) oscillator vs textbook band strength | 0.9939 | 1.0201 |
  | (ii) band integral vs `beta<Mdot^2>` | 0.14% | **0.09%** |
  | (iii) floor margin at the window minimum | **1.8x** | **1243x** |

  Gate (i) moves the wrong way, from 0.6% to 2% off, and the reason is the
  taper's ~1.5x noise-equivalent bandwidth: the synthetic test's line is 27 cm^-1
  wide, comparable to the widened resolution, so the trapezoid over it integrates
  slightly differently. On the real spectra, whose features are far broader than
  the resolution, the effect is smaller — which gate (ii) confirms at 0.09%.

  Gate (ii) needs one care: the reference `beta <|Mdot_band|^2>` must be taken on
  the *same tapered* record. Compared against the untapered series it shows a
  spurious 3-7%, which is what the A80 tool reported on first being handed a
  tapered spectrum; taper-consistent, it is 0.09%.

  **(3) What subtraction would have done, and why it was the wrong repair.**
  Before the taper was found, the measured flat floor was subtracted as a
  constant (its level the lowest 2000-cm^-1 binned median between 6000 cm^-1 and
  Nyquist, uncertainty 0.3-4.1%). It removed essentially the whole window: the
  margin went from 0.8x to **0.0x**, kappa_R fell by a median of **62%** at
  rho <= 60 and 35% at 60-300 while kappa_P moved 1%, and 25 of 52 cells lost
  more than 20% of kappa_R. That is the signature of subtracting something that
  is not additive noise but an artefact of the estimate itself. The record is
  kept here because it is the measurement that showed the subtraction route was
  wrong.

  **A75 and A52 are both superseded.** Both modelled the dipole model's residual
  as **white noise in M**, which gives `sigma ~ nu^2` and a floor that must be
  scaled down from the 3900-5500 cm^-1 estimation window to 2550. Neither part
  of that survives: the residual is a smooth function of configuration, so it
  distorts the bands rather than raising a floor, and the flat tail that was
  taken for it was the estimator's leakage. **There is no white-noise floor on
  these spectra.** A52's original worry, A75's reconstruction of it and A80's
  blocker are one artefact seen three times.

  **(4) Plain statement.** The cause is spectral leakage from the untapered
  periodogram, not the dipole model, not the dynamics and not periodic-image
  wrapping. Corrected with a Hann taper, the floor margin at the worst cells is
  **1243x** at the window minimum and **2390x** at 2550 cm^-1, against the 10x
  required — three orders of magnitude of headroom, at every one of the 14
  dilute points. **Propagation can proceed on this criterion**, and A81 does not
  propagate because it was told not to: the registry still points at v3, no
  cache is cleared and nothing downstream has been rerun. A gate asserts it.

  **The window depth is real, not the taper's main lobe.** A Hann taper trades
  -13 dB sidelobes for -31 dB but widens the main lobe, which could fill a deep
  window from its edges. Against a Welch estimate (8 segments, 50% overlap —
  coarser resolution, independent noise statistics) the tapered sigma at
  2550 cm^-1 agrees to **6-12%** at (700, 30), (1000, 10) and (2000, 30). The
  windows are as deep as they look.

  **One thing to carry into A82, and it is not small.** With the leakage gone
  the windows are genuinely deep, so kappa_R is now sensitive to the smoothing
  width: sd 10 against sd 40 moves it by **13.2%**, where untapered v4 gave 4.3%
  and kappa_P is unmoved at 0.3%. Worse, it does not converge. At (700, 30):

  | sd (cm^-1) | 5 | 10 | 15 | 20 | 30 | 40 |
  |---|---|---|---|---|---|---|
  | kappa_R | 3.405 | 3.494 | 3.561 | 3.629 | 3.797 | 4.026 |

  kappa_R rises monotonically and is **still moving by 2.6% between sd 5 and 10**,
  so there is no plateau to sit on and the delivered sd = 20 is a choice rather
  than a converged value — it is 6.6% above sd 5 and 11% below sd 40. Which
  width is right is a physics question (the true linewidth supportable by a 30 ps
  record), not a taste one. **A82 must settle it before the release tables are
  rebuilt**, and until then it belongs beside A71's wing cut-off as a stated
  systematic of the same size.

  **The bundle's fingerprint.** `ed_data/arc_dipoles/` is gitignored as raw data
  (A80), so `ed_data/arc_dipoles_sha256.json` carries a tracked SHA-256 per cell
  with byte counts — 67 series — so the inputs behind A80 and A81 can be
  verified without the 41 MB being in git.

## 2026-09-30 (the estimator replaced; propagation held at the noise floor)

- **A80 — The fixed estimator is built, gated and written as v4. Two of the
  three pre-propagation gates pass. The third fails, and it is the one that
  matters: with the estimator corrected, the dilute noise floor is **1.8x**
  below the deepest point of the transparency window at (1000 K, 10 kg m^-3),
  not the 10x required — 55% of what is there is dipole-model noise. Nothing has
  been propagated. The registry still points at v3.**
  `opacity/dipole_spectrum.py` (the estimator), `tools/v4_regenerate.py` ->
  `outputs/v4_regenerate.json`, five gates in
  `validation/test_a80_v4_gates.py`. No manuscript changes.

  **Manifest, before anything else.** All **52 delivered cells matched**, none
  missing. Fifteen tags have two directories — an older ~30002-frame run and an
  `_05fs` 60002-frame one — so matching is done by *reproduction* rather than by
  name: each candidate is pushed through the unfixed production backend and the
  one that reproduces the delivered `sigma_raw` wins. All 52 resolve to the
  `_05fs` series, median error **1.4e-8**, worst **7.9e-5**. The 15 bare
  directories are unused.

  **(1) The estimator.** `opacity/dipole_spectrum.py` implements

      sigma(w) = beta w^2 Itilde_2(w) / (6 eps0 c N)

  with `Itilde_2` the **two-sided** power spectrum of the mean-removed cell
  dipole, estimated **window-free** by the periodogram of the full series, the
  harmonic correction entering **once** (as the cancellation A79 derived), and
  **no zero clip** — the clip is dropped because a periodogram is non-negative
  by construction, so it was provably inert here, not because it was forgotten.
  Everything downstream of the estimator is kept as delivered: the 1 cm^-1 grid,
  the Gaussian smoothing (sd 20, and 15 at the four 1500 K points — confirmed
  from the files: exactly four carry 15, at rho = 150, 300, 800, 1200), the
  number density and the mass-opacity conversion. So a v3/v4 comparison isolates
  the estimator and nothing else. It reproduces A79's independent implementation
  to **3e-12**.

  **Smoothing sensitivity, kappa_R and kappa_P at sd 10 and 40:**

  | T | rho | sd used | kappa_R 10 / used / 40 | kappa_P 10 / used / 40 |
  |---|---|---|---|---|
  | 700 | 30 | 20 | 15.24 / 15.42 / 15.93 | 455.9 / 456.1 / 456.9 |
  | 1000 | 100 | 20 | 29.20 / 29.48 / 30.16 | 437.3 / 437.4 / 437.7 |
  | 1500 | 300 | 15 | 57.07 / 57.41 / 58.80 | 450.0 / 450.0 / 450.2 |
  | 2000 | 1200 | 20 | 558.0 / 560.9 / 564.1 | 906.9 / 906.8 / 906.5 |
  | 2500 | 500 | 20 | 170.5 / 172.3 / 174.5 | 423.7 / 423.6 / 423.5 |
  | 3000 | 2228.71 | 20 | 6709 / 6742 / 6758 | 6603 / 6596 / 6584 |

  **kappa_R moves by at most 4.3% and kappa_P by 0.3% across a factor of four in
  smoothing width.** The delivered choice is not a sensitive one, and this is a
  marked improvement on v3, where the clip made the windows noise-dominated and
  therefore smoothing-dependent.

  **(2) The gates.**

  | gate | result | |
  |---|---|---|
  | (i) A79 oscillator through the v4 estimator | 0.9939 against 1 | **PASS** |
  | (ii) band integral = `beta <|Mdot_band|^2>` at the three rho = 30 cells | worst deviation **0.14%** against 2% | **PASS** |
  | (iii) A75 dilute floor >= 10x below the window | **1.8x** | **FAIL** |

  For (ii) the dipole derivative is taken **spectrally**, by multiplying by
  `i w` in Fourier space, not by central differences: the central-difference
  estimator carries a known `[sin(w dt)/(w dt)]^2` bias worth 5% at the top of
  the stretch band, which is a property of the test rather than of the estimator
  under test. With it, every band at every cell agrees to 0.1%.

  **(iii) is the blocker, and the reason is not what A75 implies.**

  | T | rho | margin at the window minimum | noise fraction there |
  |---|---|---|---|
  | 1000 | 10 | **1.8x** | **55%** |
  | 700 | 60 | 2.2x | 46% |
  | 700 | 30 | 3.2x | 31% |
  | 1000 | 30 | 6.0x | 17% |

  A white residual carrying the lag test's uncorrelated fraction
  `r^2 = 7.26e-5` of the dipole power has a flat two-sided spectrum
  `Itilde = r^2 <|M|^2> dt`, which is pushed through the *same* estimator as the
  signal. At 2550 cm^-1 — where A75 quoted its margins — the range is
  **5.3x to 108.9x**; at the deepest point of the 1900-3100 window it is
  **1.8x to 83.7x**. Both fail the 10x bound at the cold, dilute end.

  **A75's margin was not wrong; it was not well posed, and A80 supersedes it.**
  A75 compared a clean analytic floor against the delivered v3 kappa. But v3's
  window opacity was itself largely artefact: the clip at zero rectified noise
  upward while the symmetric window suppressed the signal, which is exactly why
  A79 found that no scalar could repair kappa_R. **v4 is the first version in
  which the window opacity is a clean estimate, and therefore the first in which
  the floor comparison means anything.** A75 did not establish a large margin;
  it established nothing.

  **The floor quoted here is a lower bound.** Against the spectrum's own plateau
  at 4500-5800 cm^-1, above every water band, the observed level runs 1.3x to
  22x the lag-test prediction, and the excess grows with temperature — real hot
  bands live up there at 1500 K and above. But at the cold dilute cells there is
  almost no real absorption there, so their excess is noise the lag test did not
  see: it was measured at (800, 61) and (1000, 130), and the rho = 10 cells are
  noisier. Using each cell's plateau-implied `r^2` instead, the worst margin
  falls from 1.8x to **1.4x**. Both numbers are reported and neither is averaged
  into the other.

  **(3) Not propagated.** The brief gates propagation on (2), and (2) fails, so
  `tools/v4_regenerate.py --write` puts the 52 v4 spectra in their own directory
  (`fp_opacity_v4/`, units marker `cm2_per_g:A80`, `version = v4:fixed-estimator`)
  for inspection, and **the registry is untouched, no cache is cleared and
  nothing downstream has been rerun**. A gate asserts that
  (`test_nothing_has_been_propagated`), so v4 cannot be wired in while the floor
  gate fails without something turning red. v3's `crossover_hitemp.json` is
  v3-derived and would be silently wrong against v4, so it is carried into the
  v4 directory only under the name `crossover_hitemp_v3_DO_NOT_USE.json`.

  **What v3 -> v4 does, which is the table that can be given.** The two defects
  partly cancel — the spare `Q_harmonic` made v3 too large, the window and clip
  too small — so the *magnitude* barely moves. The **temperature slope** does:

  | T (K) | n | median kappa_R v4/v3 | median kappa_P v4/v3 |
  |---|---|---|---|
  | 650 | 4 | **0.746** | 1.088 |
  | 700 | 3 | 0.862 | 1.093 |
  | 800 | 5 | 0.741 | 1.017 |
  | 1000 | 8 | 0.821 | 1.024 |
  | 1500 | 8 | 0.992 | 1.146 |
  | 2000 | 8 | 1.157 | 1.383 |
  | 2250 | 2 | 1.363 | 1.436 |
  | 2500 | 5 | 1.472 | 1.520 |
  | 2750 | 2 | 1.529 | 1.564 |
  | 3000 | 7 | **1.638** | 1.608 |

  Over the 52 cells kappa_R moves by a median of 1.081 (range 0.608-1.695) and
  kappa_P by 1.158 (0.989-1.655) — but **the temperature dependence of kappa_R
  changes by a factor of 2.2 across the grid**, from 0.75x at 650 K to 1.64x at
  3000 K. Every downstream number inherits that, which is why propagation is
  worth doing properly rather than quickly. The manuscript-number table the
  brief asks for is what propagation would produce; it is not given here because
  nothing was propagated.

  **The dilute band medians under v4**, FP/HITEMP at the 13 dilute points plus
  (800, 100): rotation **1.044** [0.980, 1.139] against 0.97; bend **1.515**
  [1.306, 1.717] against 1.25-1.32; stretch **2.658** [2.207, 2.963] against
  1.7-3.9. These match A79's prediction from the measured 3.5155 (1.028 / 1.539
  / 2.659) to within 2%, which is an independent confirmation that the estimator
  replacement does what A79 said a correct one would.

  **Why this is not a precautionary halt.** The boundaries are not in the dense
  interior. Of the 108 radiative-convective boundaries, **89 sit at
  rho <= 60 kg m^-3** (42 of them on computed data), the median boundary is at
  **rho = 12.9**, and **47 sit below the grid's rho floor of 10** and therefore
  read the dilute edge cells through the A26 flat-in-rho clamp. The corner where
  the floor gate fails is the corner the boundary physics lives in. Propagating
  now would build the boundary result on window opacities that are 17-55%
  dipole-model noise at their coldest, most dilute points.

  **The two options, and the choice is Harry's.**
  (a) Propagate v4 anyway and flag the dilute window opacities as noise-limited,
  the way A71 flagged the wing cut-off — defensible only if the boundary result
  is shown to be insensitive to the dilute edge, which the 47-of-108 count above
  argues against.
  (b) Reduce the floor first — longer trajectories, or a dipole model with a
  smaller uncorrelated residual — then propagate. Both are ARC jobs. This is the
  recommendation.

  **Where the dipole bundle lives.** `ed_data/arc_dipoles/` is 41 MB of raw
  trajectory product across 67 cells — the raw input SPEC §1 keeps out of the
  repo, unlike A78's 6.8 MB of scripts and spectra, which is provenance. It is
  gitignored and stays on disk beside `~/subneptune_data`; the A80 gates read it
  the way every opacity gate reads the data root, and skip cleanly when it is
  absent. Say so if you would rather it were tracked under git-LFS.

  A81 should carry whichever of (a) and (b) is chosen. The
  `opacity_pipeline.py` fix belongs upstream with the reruns and is still not
  applied to the provenance copy.

## 2026-09-30 (the normalisation is right; the estimator is not, and A80 is blocked on it)

- **A79 — The production prefactor is exactly correct: the derived
  normalisation factor is 1. The factor of ~3.4 A78 left open is not a
  normalisation at all but the line-66 window and the line-79 clip, worth a
  measured 3.52. With it, the dilute bands land on the monomer test's
  expectation. The N audit is clean. But A80's recipe cannot proceed as
  written: a scalar on `sigma_raw` leaves a **40% temperature swing in kappa_R**
  along the rho = 30 row, the quantity the manuscript rests on.**
  `tools/normalisation_check.py` -> `outputs/normalisation_check.json`.
  Existing data only; no manuscript changes.

  **(a) The derivation.** Conventions, stated once: SI throughout; `M` is the
  **total cell** dipole in C m and `N` the molecule count;
  `C(t) = <M(0).M(t)>` is the **trace**, summed over the three Cartesian
  components and not averaged, so the orientational average of the field
  projection supplies the explicit 1/3;
  `Itilde_2(w) = int_{-inf}^{inf} dt e^{iwt} C(t)` is the **two-sided**
  transform, equal to `2 Itilde_1(w)` since C is real and even; the refractive
  index is 1; and sigma is per molecule, so the cell result is divided by N.

  The quantum expression is
  `alpha(w) = w (1 - e^-x) Itilde_qm(w) / (3 eps0 hbar c V)` with
  `x = beta hbar w`. The harmonic correction is `Itilde_qm = Q_h Itilde_cl`,
  `Q_h = x/(1 - e^-x)`, so the `(1 - e^-x)` **cancels**, leaving

  | form | expression |
  |---|---|
  | two-sided | `sigma(w) = beta w^2 Itilde_cl,2(w) / (6 eps0 c N)` |
  | one-sided | `sigma(w) = beta w^2 Itilde_cl,1(w) / (3 eps0 c N)` |

  The 6-versus-3 is exactly the one-/two-sided convention and nothing else.
  Production line 82 writes the **one-sided** form, and its `Smm` at line 66 is
  `dt Re rfft(Cmm . w)` — a one-sided cosine transform. **The prefactor
  `beta w^2 / (3 eps0 c N)` and the `1e4` (m^2 -> cm^2) are both correct, term
  for term. The derived factor by which production differs is 1.**
  `spectra.py` is not in this chain at all: it computes
  `alpha = w (1 - e^-x) I` for diagnostics and writes `spectrum.npz`; nothing
  downstream of the opacity grid reads it.

  **This is not taken on trust from a remembered reference form.** The
  one-/two-sided factor of two is exactly where these derivations go wrong, so
  it is settled numerically: a thermostatted classical harmonic oscillator of
  known `dmu/dq` and reduced mass, integrated and compared against the textbook
  classical band strength `A = (dmu/dq)^2 / (12 eps0 c^2 m)` (McQuarrie; the
  same constant behind the km/mol intensities quantum-chemistry codes print).
  The two-sided/6 form reproduces it to **0.994**; the two-sided/3 form gives
  **1.988**. Equipartition in the test integrator holds to 0.1%.

  **So the ~3.4 A78 could not name is not a physical constant.** It is the
  spectral estimator: the symmetric `np.hanning(L)` on a one-sided ACF (line 66)
  followed by the clip at zero (line 79). Measured against the window-free
  periodogram estimate of `Itilde_2`, which satisfies Parseval exactly and needs
  no lag window, the shortfall is

  | cell | N | rho | rotation | bend | stretch |
  |---|---|---|---|---|---|
  | T700_rho30 | 216 | 30 | 3.487 | 3.525 | 3.699 |
  | T1000_rho30 | 216 | 30 | 3.574 | 3.618 | 3.508 |
  | T2000_rho30 | 216 | 30 | 3.390 | 3.473 | 3.516 |
  | T1000_rho300 | 64 | 300 | 3.363 | 3.594 | 3.470 |
  | T2000_rho300 | 64 | 300 | 3.581 | 3.505 | 3.643 |

  median **3.5155**, range 3.363-3.699, **10% spread over two densities, three
  temperatures and three bands**. The brief asked for a derived factor and said
  not to fit one; the honest report is that the derived factor is 1 and this one
  is **measured**, so it is quoted as a measurement with its scatter and not as
  a constant.

  **(b) The dilute-limit check, and it agrees with (a).** Dividing out the spare
  `Q_harmonic` and applying 3.5155, FP/HITEMP over the 13 dilute points:

  | band | median | range | monomer-test expectation | ratio |
  |---|---|---|---|---|
  | rotation 50-1100 | **1.028** | [0.946, 1.182] | 0.97 | 1.06 |
  | bend 1000-2400 | **1.539** | [1.274, 1.667] | 1.25-1.32 | 1.20 |
  | stretch 2800-4000 | **2.659** | [2.187, 3.039] | 1.7-3.9 | 0.95 |

  The rotation band, whose expectation is the tightest (the model's mu^2 ratio,
  0.974), lands within **6%**. The stretch lands inside its range. The bend is
  **20% above** the top of its 1.25-1.32 window. That residual is real rather
  than estimator noise: truncating the trajectories to 15 ps and 7.5 ps leaves
  the per-band factors at 3.14-4.08 with no drift and no systematic spread
  between bands, so the flatness of 3.49/3.52/3.52 is structural and not a
  property of this record length. 20% is four times the estimator's own 5%
  scatter — but it is also well inside the uncertainty of the expectation
  itself, which is a squared derivative ratio read at two bond lengths, so it
  is recorded and not interpreted. A76's factor-2.6 rotation deficit and
  A78's uniform 0.3 are both fully explained: they were the estimator, not the
  dipole model. **(a) and (b) agree, so the brief's stop condition is not
  triggered.**

  **(c) The N audit is clean; A78's flag is withdrawn.** After both corrections,
  the rotation-band integral per molecule across all 52 files rises smoothly
  with density from 4.7e-17 to 2.1e-16 cm^2 cm^-1 — a factor of 4.5, and in the
  direction expected, since the condensed-phase dipole grows with density. The
  13 N=216 files span 4.87-6.07e-17 (max/min **1.25**) and the 39 N=64 files
  4.72-21.3e-17; where the two groups meet in density the ratio of medians is
  **0.934**, against the **3.375** (or 0.296) a misassigned N would produce.
  **No cell is flagged.** The step would be unmissable and there is none, at any
  density. Three of the 13 N=216 cells have local run records (the A78 bundle);
  all three record N=216, and their mean cell dipole implies 1.87 D per
  molecule at N=216 against 3.43 D at N=64. The two box sizes are two
  legitimate production batches.

  **A80 is blocked, on different grounds from the brief's stop condition.**
  A80 would rescale the delivered `sigma_raw` by a constant. Tested against the
  honest recomputation from the dipole series, on the five cells where that is
  possible:

  | cell | kappa_R scalar/recomputed | kappa_P scalar/recomputed |
  |---|---|---|
  | T700_rho30 | **0.761** | 0.988 |
  | T1000_rho30 | **0.870** | 0.985 |
  | T2000_rho30 | 1.061 | 0.998 |
  | T1000_rho300 | 1.041 | 0.996 |
  | T2000_rho300 | 1.043 | 0.980 |

  **The kappa_R error runs with temperature: 0.761 at 700 K to 1.061 at 2000 K
  along the rho = 30 row — a 40% swing.** That is larger than the factor-2.1
  bend trend that started this whole chain at A76, and it would have gone
  straight into every cooling track. kappa_P survives to 2%; kappa_R does not. A
  scalar restores band integrals, which is what kappa_P is; it cannot restore a
  lineshape, which is what kappa_R is. kappa_R is set by the transparency
  windows, where the true sigma is smallest and the clip at zero rectifies noise
  upward, so the relative bias is largest exactly where kappa_R has its weight.
  Since kappa_R is what sets the radiative-convective boundaries, the cooling
  tracks and most of the manuscript's headline numbers, regenerating 52 spectra
  this way would propagate a 25% temperature-dependent error into all of them,
  and it would be invisible afterwards because every downstream number would
  move together and look self-consistent.

  **What A80 needs instead.** The correct operation is not a rescaling but a
  recomputation: `sigma = beta w^2 Itilde_2 / (6 eps0 c N)` with `Itilde_2`
  estimated window-free, per cell, from its own dipole series. That needs each cell's own
  dipole series. This machine has five usable (plus the two `*_BAD64` series,
  which are different runs at N = 256 and are not usable), so **47 are
  requested**: `$DATA/scw_prod/<tag>/dipoles.npy` for

  > T650_rho100, T650_rho150, T650_rho300, T650_rho500, T700_rho10, T700_rho60,
  > T800_rho30, T800_rho100, T800_rho150, T800_rho300, T800_rho500, T1000_rho10,
  > T1000_rho60, T1000_rho100, T1000_rho129.744, T1000_rho150, T1000_rho500,
  > T1500_rho10, T1500_rho30, T1500_rho60, T1500_rho150, T1500_rho300,
  > T1500_rho471.343, T1500_rho800, T1500_rho1200, T2000_rho10, T2000_rho60,
  > T2000_rho150, T2000_rho500, T2000_rho800, T2000_rho1200, T2250_rho175.545,
  > T2250_rho1256.23, T2500_rho197.586, T2500_rho300, T2500_rho500,
  > T2500_rho800, T2500_rho1520.77, T2750_rho222.666, T2750_rho1841.02,
  > T3000_rho108.796, T3000_rho729.344, T3000_rho810.382, T3000_rho1043.59,
  > T3000_rho1343.91, T3000_rho1730.66, T3000_rho2228.71

  The local directories spell these with a trailing `_05fs`; either resolves,
  since `mace_opacity.py` parses `T(\d+(?:\.\d+)?)` and `rho([\d.]+)` out of the
  basename. **Nothing else is needed**: `N` and `dt_fs` are already carried by
  each delivered `.npz` (all 52 are `dt_fs = 0.5`; 39 are `N = 64` and 13
  `N = 216`). With those 47 files A80 becomes a clean rerun of the fixed
  pipeline and needs **no factor at all**. A79 therefore stops here and
  regenerates nothing.

  `opacity_pipeline.py` is a copy of an ARC file inside the provenance pull, so
  it is left byte-identical; the fix belongs upstream, where the reruns happen,
  and A80 should carry it. Nothing in this amendment applies the 3.5155 to any
  delivered product: it is a measurement used to test (b) and to decide (a)
  against (b), not a correction.

  **Gate.** `validation/test_a79_normalisation.py` pins the one-/two-sided
  factor of two with the synthetic oscillator, since the whole conclusion —
  that the production prefactor is correct — turns on it and no amount of
  re-reading the code settles it.

## 2026-09-30 (the production formula, read: the harmonic correction is applied twice)

- **A78 — `opacity_pipeline.py` multiplies `beta omega^2 Smm`, which is already
  the harmonic-corrected result, by `Q_harmonic` a second time. The delivered
  spectra carry the harmonic quantum correction twice. That spare factor is the
  whole of the temperature trend A76 found and A77 measured grid-wide: at the
  source, `beta <|Mdot|^2>` obeys the classical sum rule exactly.**
  `tools/dipole_sum_rule.py` -> `outputs/dipole_sum_rule.json`; four gates in
  `validation/test_a78_dipole_sum_rule.py`. Existing data only.

  **(1) The formula, line by line.** `ed_data/arc_a78/scw/opacity_pipeline.py`,
  `spectrum_from_dipoles`, imports only `numpy` and
  `scipy.ndimage.gaussian_filter1d` — there is no third module to chase.
  `spectra.py` imports `torch`, `ase.io` and `mace` and does nothing but
  evaluate the dipole model frame by frame and cache `dipoles.npy`; it contains
  no opacity arithmetic. The chain, with the line numbers:

  | line | operation |
  |---|---|
  | 54-58 | `M` in e.A -> Debye (`EANG_TO_DEBYE`), mean removed at line 64 |
  | 60 | `M_cm = M * 3.33564e-30` (Debye -> C m) |
  | 65 | `Cmm = sum_k autocorr_fft(M_k)`, **unbiased** (`ac / (n - lag)`) |
  | 66 | `Smm = dt_s * Re rfft(Cmm * np.hanning(L))` |
  | 77-79 | keep `0 < wn <= 6000`; `Smm` clipped at 0 |
  | 80 | `x = beta hbar omega` |
  | 81 | `qcf = x / (1 - exp(-x))`  — **this is `Q_harmonic`** |
  | 82 | `sigma_raw = (beta omega^2 / (3 eps0 c N)) * Smm * qcf * 1e4` |
  | 83-84 | linear interpolation onto a 1 cm^-1 grid, NaN -> 0 |
  | 85 | `gaussian_filter1d(sigma, sd = smooth_cm)` |
  | 86 | `alpha = n sigma`, `n = rho/M_w N_A`; `kappa = alpha / rho_cgs` |

  **It is neither of A77's candidates.** A77 carried form A
  (`sigma ~ w (1-e^-x) I_cl`) and form B (`sigma ~ beta w^2 I_cl`). Production is
  **form B times a second `Q_harmonic`**. The algebra is three lines: the
  standard relation is `alpha ~ w (1 - e^-x) I_qm`; the harmonic correction is
  `I_qm = Q_harmonic I_cl` with `Q_harmonic = x/(1 - e^-x)`; so
  `alpha ~ w (1-e^-x) . x/(1-e^-x) . I_cl = beta hbar w^2 I_cl` — **the
  `(1-e^-x)` has already cancelled against `Q_harmonic`'s denominator.** Line 82
  starts from that cancelled result and applies `Q_harmonic` again.

  Read, not paraphrased: running the function as it stands on each supplied
  `dipoles.npy` reproduces the delivered `sigma_raw` to **8e-8 relative** (float
  round-off on float32 input), and the three rho = 30 cells are bit-identical to
  the delivered grid files, so these are the production cells.

  **A second defect in the same function, line 66.** `np.hanning(L)` is a
  *symmetric* window and `Cmm` is a *one-sided* autocorrelation, so the zero-lag
  value — the largest and best-determined point — is multiplied by **zero**, and
  the peak weight lands at lag L/2 where the unbiased estimator `ac/(n-lag)` is
  at its noisiest. `spectra.py` line 46 does the same job correctly, with
  `np.hanning(2*maxlag)[maxlag:]`. This distorts the lineshape and the overall
  scale; it does **not** introduce a temperature trend, which (3) below
  establishes independently. Its size cannot be quoted without also settling the
  `maxlag` convention (`spectra.py` truncates at 8000 lags, the pipeline does
  not), so A79 should settle both together.

  **(2) The sum rule at its source, with no formula in between.** From each
  `dipoles.npy`: `Mdot` by central differences at the stored 0.5 fs, an FFT mask
  per band, then `beta <|Mdot_band|^2>` per molecule. This is the quantity the
  classical sum rule fixes, since `beta <|Mdot|^2> = beta (kT/m) sum (dmu/dQ)^2`.
  The three rho = 30 cells, 60002 frames, N = 216:

  | band | 700 K | 1000 K | 2000 K | index p |
  |---|---|---|---|---|
  | rotation 50-1100 | 4.915e-12 | 5.050e-12 | 5.587e-12 | **+0.125** |
  | bend 1000-2400 | 1.663e-12 | 1.728e-12 | 1.578e-12 | **-0.061** |
  | stretch 2800-4000 | 1.729e-12 | 1.813e-12 | 2.114e-12 | **+0.196** |

  **The dipole series satisfies the sum rule.** Against that, the delivered band
  integrals give +0.049, **-0.741** and **-0.737**. Removing line 81-82's `qcf`
  and nothing else recovers +0.158, **-0.049** and **+0.229** — which track
  `beta <|Mdot|^2>` to within 0.03 in every band, at a ratio constant to
  **3.5-5.4%** across the three temperatures. (That constant absorbs the line-66
  window and the unbiased-estimator normalisation; only its constancy is the
  claim, and its constancy is what shows the window carries no temperature
  dependence.) The trend attributable to the spare factor is **-0.108**,
  **-0.692** and **-0.966**. Taking `Q_harmonic` on its own, weighted by each
  band's own (un-doubled) intensity, gives **-0.108, -0.692 and -0.966** — the
  same numbers to three decimals. The spare factor does not mostly account for
  the trend; it *is* the trend.

  **(3) Plain statement.** **No** — the temperature trend is *not* present in
  `beta <|Mdot|^2>`; that quantity is flat to |p| <= 0.2 in all three bands.
  **The defect is in the formula, not the dipole series.** So the brief's
  fallback does not apply and no equipartition check is needed: the dynamics and
  the dipole model are exonerated for this finding, and the ARC velocities are
  not required. What is required is one line of arithmetic in
  `opacity_pipeline.py`.

  **Corrections to A75, A76 and A77.**

  - **A75 is wrong where it says the delivered spectra carry `Q = 1` on the
    lineshape.** That was read from `spectrum_convergence.py`, which is not the
    production code. Production carries `Q_harmonic` **twice**. A75's noise-floor
    result does not depend on this and stands.
  - **A76's conclusion was right for the wrong reason.** It found that none of
    the three corrections improves the match and that each steepens the trend.
    Of course: the spectra already carried one too many. The operation needed is
    to **divide** by `Q_harmonic`, which A76 never tried because A75 had told it
    the factor was absent. Redone with the formula corrected, A76's calibration
    comes out flat:

    | band | corrected ratio / dipole factor | range | dln/dT per 1000 K | rho_s |
    |---|---|---|---|---|
    | rotation | 0.302 | [0.275, 0.346] | **-0.000** | +0.08 |
    | bend | 0.341 - 0.323 | [0.288, 0.382] | **-0.015** | +0.01 |
    | stretch | 0.365 - 0.196 | [0.296, 0.423] | +0.191 | +0.82 |

    The bend's temperature trend goes from -0.553 per 1000 K to **-0.015**, and
    its Spearman rho from -0.89 to **+0.01**. That the trend vanishes under a fix
    derived independently, from reading the source, is the strongest evidence
    that the fix is the right one.
  - **A77's grid-wide violation is entirely this.** Its measured indices stand;
    its open A/B question is closed against both. Its statement that the trend
    "enters at or above the ACF -> cross-section step" was correct, and the
    answer is the step itself.

  **What the fix does not repair, and it is now the open question.** With the
  spare factor removed, the first-principles bands sit at a **uniform ~0.3 of the
  line list** once the monomer dipole factor is divided out — 0.30 in rotation,
  0.34 in bend, 0.37 in stretch — flat in temperature and consistent across
  bands. A76's rotation deficit therefore deepens from 2.6x to **3.3x** rather
  than closing. A frequency-independent factor of three is the signature of a
  normalisation, not of a dipole surface. The line-66 window sets an absolute
  constant and is the obvious suspect, but it cannot simply be swapped: the
  half-window alternatives change the band integrals by 4.5x, 30.5x and 82.5x
  — strongly frequency-dependent, nothing like a uniform 3.3x — and they fail
  the constancy test above, so they are not a correction. **A79 should derive
  the absolute normalisation of `beta omega^2 Smm / (3 eps0 c N)` from first
  principles and compare it against the measured 3.0e-6, the one number that
  would settle both the window and the factor of three.** Do not treat the 0.3
  as a dipole-model result until then.

  **Two provenance notes, neither affecting the above.** (i) The supplied
  `T2000_rho300_05fs/opacity_v3.npz` carries the **legacy x1000** error in
  `alpha` and `kappa` while its `sigma` is correct and identical to the delivered
  file; the internal check `alpha = n sigma` catches it, and **all 52 delivered
  files pass that check**, so the paper's grid is clean and only the bundled cell
  is stale. (ii) The delivered grid is **mixed in box size**: 39 files carry
  N = 64 and 13 carry N = 216, and `mace_opacity.py` takes N from the command
  line with a default of 64. For the three rho = 30 cells N = 216 is right — the
  mean cell dipole 5.267 e.A implies 1.87 D per molecule, in the range expected
  for near-isolated water (gas phase 1.85 D), whereas N = 64 would imply 3.43 D,
  which no water molecule has; this discriminates N but says nothing about how
  well the model extrapolates — but the two
  `*_BAD64` directories show the pipeline has been run with a wrong N before
  (they carry N = 256 and different dipole series, and were not used here).
  Harry should confirm the box size of the 13 N = 216 cells against their run
  records; `sigma ~ 1/N`, so a wrong N is a factor of 3.4 in kappa.

  **Gates.** `validation/test_a78_dipole_sum_rule.py` adds four: the formula as
  read reproduces the delivered `sigma_raw` (<= 1e-5; measured 1.7e-7); the
  dipole series satisfies the sum rule (|p| < 0.3; measured <= 0.196); the band
  integrals are recovered from the dipole series at a temperature-independent
  ratio (spread < 10%; measured <= 5.4%); and the delivered bend integral still
  carries the spare correction (index < -0.5; measured -0.741), so an upstream
  fix cannot pass unnoticed.

  These gates read the bundle, so `ed_data/arc_a78/` is tracked in full, and
  `.gitignore`'s global `*.npz` is lifted inside `ed_data/` exactly as A69
  lifted the PNG rule there: the seven `opacity_v3.npz` are gate inputs and
  provenance, not run outputs. 6.8 MB, largest file 720 kB, so no git-LFS.

## 2026-09-30 (tracing the sum-rule violation: grid-wide, and half of A76's index was a formula assumption)

- **A77 — The violation is real, it is grid-wide rather than dilute-only, and it
  survives under either candidate production formula. But A76's specific index
  came from a script whose equivalence to the production one is unverified, and
  the production script is not on this machine. Half of A76's magnitude rests on
  that unverified equivalence, not on a measurement.**
  `tools/sum_rule_trace.py` -> `outputs/sum_rule_trace.json`; two gates in
  `validation/test_sum_rule.py`. Existing data only; no manuscript changes.

  *Superseded in part by A78, which read the production script: it is neither
  candidate form, but form B times a second `Q_harmonic`. The measured indices
  here stand; the open A/B question is closed.*

  **(1) The production formula cannot be read here, and A76 cited a convergence
  script for it.** The production chain is
  `spectra.py` -> `dipoles.npy` -> `~/scw/opacity_pipeline.py ::
  spectrum_from_dipoles` -> `*_opacity.npz`, named in
  `ed_data/arc/water-scw-mace/lda_campaign/stepG_dipoles/make_opacity_lda.py`
  lines 2-4 and 17. **Neither `~/scw/opacity_pipeline.py` nor `spectra.py` is in
  the ARC pull or anywhere on this machine** (`~/scw` does not exist), and
  nothing else in the pull references them. A76 took its formula from
  `convergence/04_spectra/spectrum_convergence.py`, where line 28 states form A
  unconditionally — so A76 **read** its formula rather than assuming it. What it
  took on trust was that script's own docstring, which describes it as working
  from a saved dipole series with "same conventions as the production
  `spectra.py`": a claim of matching conventions, not of being the production
  code, and one nothing in the pull substantiates.

  **They demonstrably differ.** The convergence script FFTs onto a `nfft = 2^17`
  `rfft` grid (0.509 cm^-1 spacing at `dt` = 0.5 fs) and applies no smoothing;
  the delivered files are on a uniform **1.000 cm^-1** grid of 6000 points from
  1.1118 to 6000.1 cm^-1 and carry a 20 cm^-1 smoothing. So the answer to "do
  the two agree line by line" is **no** at the grid and the smoothing, and
  undetermined at the formula.

  **What is recoverable is recoverable exactly.** Every production step still
  visible in the delivered arrays reproduces on all 52 files:

  | step | verified to |
  |---|---|
  | `sigma_raw = max(... - noise floor, 0)` | 12.0-21.6% of channels clipped, **zero negatives** |
  | `sigma = gaussian_filter1d(sigma_raw, sd = 20 cm^-1)` | **exactly** (median rel. error 0.0) |
  | `alpha = n sigma`, `n = rho / M N_A`, M = **18.015280 g/mol** | 4.4e-16 |
  | `kappa = alpha / rho_cgs` | 4.4e-16 |

  Note the smoothing width is a **standard deviation**, not a FWHM (assuming
  FWHM misses by 7.7%) — so the effective FWHM is 47 cm^-1, broader than
  METHODS.md's "20 cm^-1 Gaussian smoothing" reads. The one unreadable step is
  the ACF -> cross-section map, and that is exactly where the temperature and
  quantum factors live.

  **The two candidates for it differ by exactly `Q_harmonic = x/(1-e^-x)`:**

  | | form | classical quantity |
  |---|---|---|
  | **A** | `sigma ~ w (1-e^-x) I_cl` (the convergence script, line 28; A76 used this) | `kappa x Q_harmonic` |
  | **B** | `sigma ~ beta w^2 I_cl` (the textbook classical form) | `kappa` itself |

  Both reduce to `w^2 I / T` as x -> 0 and differ only at large x, where I_cl is
  not available — so **no test on the delivered arrays can separate them.** The
  stretch band is where they differ most (x = 6.6-8.2 at 700 K, so A and B stand
  about 7x apart) and it would discriminate if the first-principles stretch
  intensity were otherwise right. It is not: the delivered stretch ratio to the
  line list is already 5.2 at (700 K, 10 kg m^-3) before any correction, outside
  the monomer test's 2.1-3.9 dipole range, so that comparison is dominated by
  the very violation this amendment is tracing. The script is needed.

  **(2) The violation is grid-wide, and survives either way.** Index p of
  `int (classical quantity) dnu ~ T^p`, which must be 0, over every density row
  with three or more temperatures:

  | rho | nT | rot A | B | bend A | B | stretch A | B | full A | B |
  |---|---|---|---|---|---|---|---|---|---|
  | 10 | 4 | -0.24 | -0.11 | -1.52 | -0.83 | -1.77 | -0.81 | -1.53 | -0.59 |
  | 30 | 5 | -0.11 | +0.02 | -1.47 | -0.76 | -1.74 | -0.77 | -1.46 | -0.50 |
  | 60 | 4 | -0.37 | -0.23 | -1.34 | -0.63 | -1.69 | -0.72 | -1.46 | -0.55 |
  | 100 | 3 | -0.42 | -0.18 | -1.53 | -0.68 | -2.41 | -1.39 | -2.15 | -0.95 |
  | 150 | 5 | -0.20 | -0.03 | -1.44 | -0.72 | -1.90 | -0.92 | -1.61 | -0.61 |
  | 300 | 6 | -0.36 | -0.17 | -1.38 | -0.71 | -1.99 | -1.04 | -1.67 | -0.70 |
  | 500 | 5 | -0.43 | -0.23 | -1.29 | -0.63 | -1.98 | -1.04 | -1.65 | -0.71 |
  | 800 | 3 | -0.10 | +0.02 | -0.80 | -0.31 | -1.66 | -0.82 | -1.21 | -0.42 |

  **Every row fails, from rho = 10 to 800 kg m^-3 — this is not a dilute-corner
  problem**, which is what A76 left open and is the main new fact. The line list,
  measured identically over eight rows from rho = 0.1 to 1000 and up to 13
  temperatures, gives **|p| <= 0.063 in every row and every window**. The control
  is clean, so the measurement is sound and the defect is on the
  first-principles side.

  **A76's index was formula-dependent and this halves it.** A76 reported T^-1.5
  and a 4.4x fall, which is column A. Over the six rows that span a factor of two
  or more in temperature (see the gates note below; the rho = 100 and rho = 800
  rows are too narrow to carry an index), column A gives full-band **-1.46 to
  -1.67, median -1.57**, and column B gives **-0.50 to -0.71, median -0.60** —
  still a violation in every row, still far outside the |p| <= 0.063 the line
  list achieves, but roughly **half the magnitude A76 claimed**. A76's
  qualitative conclusion stands; its number should be quoted as "-0.6 to -1.6
  full-band depending on the production formula" until the script is read.
  Weak evidence favours **B**: under B the rotation band is flat to +0.02/-0.23
  across all eight rows, while under A it sits at -0.10 to -0.43. The rotation
  band is the most nearly classical here and should obey the rule best. This is
  suggestive, not proof — the rotational sum rule is the loosest of the three.

  **Two mechanisms that could have faked it, both excluded.** (i) *Intensity
  leaving the window*: extending A76's 100-4000 cm^-1 range to the full
  100-6000 cm^-1 delivered grid changes p by at most 0.19, and not
  systematically towards zero (rho = 30: full -1.46, wide -1.56), so the
  strength is lost, not displaced above 4000 cm^-1. (ii) *The zero clipping*: `sigma_raw` is clipped at zero and
  ~1 channel in 6 is clipped, which rectifies noise positively — but the clipped
  fraction is **12.0-21.6% with no temperature trend** in any band or row
  (rho = 30 bend: 17.5, 19.8, 17.2, 16.0, 12.2% at 700-2000 K), so it cannot
  produce a monotone factor-of-several trend.

  **(3) and (4) Dynamics and dipole-model checks: stopped, as the brief
  directs.** No trajectories and no dipole time series are on this machine — the
  ARC pull carries scripts, figures and logs only, and `~/subneptune_data` holds
  processed opacity tables. Needed:

  - `~/scw/opacity_pipeline.py` — `spectrum_from_dipoles`, the ACF -> sigma step
    (this is the step-1 blocker, and it decides A vs B)
  - `~/scw/spectra.py` — writes `dipoles.npy`; the ACF convention
  - `$DATA/scw_prod/T700_rho30/dipoles.npy`, `.../T1000_rho30/`,
    `.../T2000_rho30/` — plus the matching trajectories with **velocities**
    (64 molecules, 0.5 fs, 60 ps) for the equipartition and dipole-regression
    checks

  rho = 30 kg m^-3 is the dilute row with the most temperatures (700, 800, 1000,
  1500, 2000 K), so it is the cheapest row spanning the trend. Directory names
  follow `make_opacity_lda.py`'s `--trajdir` convention, `T{T}_rho{rho}`.

  **(5) Plain statement — where the T^-p behaviour first appears.** As far upstream as this
  machine reaches, which is the delivered `sigma_raw` array. Everything below
  that point is verified exact and temperature-free: the clipping has no
  temperature trend, the smoothing is a fixed 20 cm^-1 standard deviation, and
  `alpha = n sigma` and `kappa = alpha / rho` introduce no temperature factor at
  all at fixed density. **So the trend enters at or above the ACF ->
  cross-section step, and this machine cannot say which of the three candidates
  it is** — the dipole model's temperature response, the trajectory dynamics
  (a sub-equipartition vibrational amplitude), or the formula in
  `opacity_pipeline.py` itself. Those are exactly steps (3), (4) and (1), and all
  three are blocked on the files listed above. What can be said without them:
  the defect is **not** in post-processing, it is **not** confined to the dilute
  corner, it is **not** an artefact of the integration window or the clipping,
  and it is **not** present in the line list processed through the same
  measurement.

  **Gates.** `validation/test_sum_rule.py` adds two. The line list must satisfy
  the rule (|p| < 0.1 in every row; measured <= 0.063) — the control, which if it
  ever fails means the measurement broke, not the spectra. And the
  first-principles bend must keep failing it (every qualifying row p < -0.40 and the
  median < -0.55; measured -0.628 to -0.825, median -0.715), so that a silent upstream fix
  cannot invalidate A76/A77 unnoticed. The second gate is pinned to the
  **delivered** kappa with no correction applied — the quantity the files
  actually contain — which keeps it independent of the unresolved A/B question.

  A row qualifies only if its temperatures span a factor of two or more, since a
  log-log index over less than that is not a measurement of a trend. Two rows
  fail, and both are three hot points: rho = 100 (650-1000 K, 1.54x) and
  rho = 800 (1500-2500 K, 1.67x). The second shows why the rule is needed —
  dropping its coldest point moves its index from -0.31 to **+0.07**, so it
  carries no information either way, and it is the row that made the raw range
  in the table above look wider than the measurement supports. Across the six
  qualifying rows no single leave-one-out fit rises above -0.535.

  Also fixed here: `tools/floor_check.py` cited the convergence script by a path
  relative to the ARC pull rather than to the repo. A78 should read
  `opacity_pipeline.py` and settle A vs B before anything else.

## 2026-09-30 (calibrating the quantum correction: it is not the free parameter)

- **A76 — None of the three corrections works, and none can. The temperature
  trend the calibration was meant to null is not in the lineshape: the
  first-principles vibrational band strengths fall as T^-1.5 where the classical
  sum rule they should obey requires T^0.**
  `tools/qcf_calibration.py` -> `outputs/qcf_calibration.json`. Existing data
  only; no new runs.

  *Superseded in part by A78: the production spectra carry `Q_harmonic` twice,
  so the operation needed was a division, not any of the three multiplications
  tested here. The conclusion -- that no correction helps -- stands; the reason
  given for it does not.*

  A75 established that the delivered spectra carry `alpha = w (1 - exp(-x)) I_cl`
  — the detailed-balance factor with **Q = 1 on the lineshape** — and that the
  three standard corrections raise kappa_R by x1.7-2.0 (standard), x2.2-4.8
  (harmonic) and x2.5-10.2 (Schofield). Which one belongs there is not decidable
  from the spectra alone. It is, in principle, decidable against the line list,
  because HITEMP is quantum by construction and the dipole model's own
  band-intensity error is known independently from the monomer test. So:
  `int kappa(nu) dnu` over each band for both sources at the 13 points with
  rho <= 60 plus (800, 100); FP/HITEMP under each of the four forms; then divide
  out the monomer test's intensity factor. What is left should be 1, flat in T.

  **The rotation and bend bands are the diagnostic.** Their intensity factors are
  tight (0.974; 1.25-1.32, the squares of the measured bend derivative ratios
  1.12 at r_OH = 0.957 A and 1.15 at 1.00 A). The stretch factor spans 2.1-3.9
  — 0.27 dex on its own — and cannot discriminate; it is reported and not used.

  **Corrected ratio = (FP/HITEMP) / intensity factor**, median [min, max] over
  the 14 points, with the log-space temperature slope, the hottest/coldest ratio
  and the Spearman rho:

  | band | correction | median | [min, max] | dln(r)/dT /1000 K | hot/cold | rho_s |
  |---|---|---|---|---|---|---|
  | rotation | delivered | **0.388** | [0.327, 0.441] | **-0.10** | 0.93 | -0.45 |
  | | standard | 0.481 | [0.384, 0.553] | -0.19 | 0.83 | -0.68 |
  | | harmonic | 0.494 | [0.389, 0.569] | -0.21 | 0.81 | -0.71 |
  | | Schofield | 0.502 | [0.391, 0.580] | -0.22 | 0.79 | -0.75 |
  | bend | delivered | **0.86 / 0.81** | [0.560, 1.198] | **-0.55** | 0.48 | -0.89 |
  | | standard | 1.55 / 1.47 | [0.858, 2.309] | -0.74 | 0.38 | -0.95 |
  | | harmonic | 2.18 / 2.06 | [0.946, 4.103] | -1.09 | 0.24 | -0.95 |
  | | Schofield | 2.69 / 2.55 | [0.997, 6.266] | -1.36 | 0.17 | -0.95 |
  | stretch | delivered | 1.75 / 0.94 | [1.102, 2.742] | -0.58 | 0.44 | -0.90 |
  | | standard | 3.48 / 1.88 | [2.057, 5.477] | -0.63 | 0.41 | -0.90 |
  | | harmonic | 9.53 / 5.13 | [3.141, 19.42] | -1.35 | 0.16 | -0.92 |
  | | Schofield | 26.5 / 14.3 | [4.152, 124.8] | -2.51 | 0.03 | -0.94 |

  (Where two medians are given they are the two ends of the factor range; the
  bracketed range is at the low end. Per-point, per-band, per-correction ratios
  are in the JSON.)

  **Every correction makes the bend worse on both criteria, and this is forced,
  not accidental.** All three Q >= 1, and all three fall towards 1 as T rises at
  fixed nu. The delivered bend ratio already falls with T. So every correction
  must steepen the trend, and the data do exactly that: -0.55 -> -0.74 -> -1.09
  -> -1.36 per 1000 K. **No choice of Q can flatten a ratio that is already
  falling.** The brief's question presupposes one of them works; none does.

  **The trend is in temperature, not density, and is not a window artefact.**
  At fixed T the bend ratio is flat in density — d ln(ratio)/d ln(rho) = +0.014,
  range -0.017 to +0.101 over the four isotherms — while at fixed rho = 30 the
  temperature slope is -0.548, indistinguishable from the pooled -0.553. The
  first-principles bend sits 5-26 cm^-1 above the line list's, so a fixed
  1400-1800 window could in principle catch a drifting fraction of it; widening
  it makes the trend *steeper*, not flatter (-0.553 at 1400-1800, -0.607 at
  1200-2100, -0.641 at 1000-2400).

  **What it is instead: the band-strength sum rule.** `Q_harmonic x (1 - e^-x)
  = x` exactly, so the harmonic-corrected form is `alpha ~ beta w^2 I_cl` — the
  plain classical absorption formula. For a classical trajectory
  `int beta w^2 I_cl dw ~ beta <|mu_dot|^2> = beta (kT/m) sum (dmu/dQ)^2` is
  **independent of temperature**, and the line list's net band strength is
  T-independent by the same sum rule. Fitting `int kappa_harmonic dnu ~ T^p` over
  windows wide enough to hold each whole band:

  | window (cm^-1) | FP rho=10 | rho=30 | rho=60 | HITEMP |
  |---|---|---|---|---|
  | rotation 50-1100 | -0.24 | -0.11 | -0.37 | **-0.00** |
  | bend 1000-2400 | **-1.52** | **-1.47** | **-1.34** | +0.06 |
  | stretch 2800-4000 | **-1.77** | **-1.74** | **-1.69** | -0.04 |
  | full 100-4000 | **-1.53** | **-1.46** | **-1.46** | +0.02 |

  HITEMP obeys the sum rule to **2.7% over the full band** across 700-2000 K, as
  it must. The first-principles side **falls by 4.4x** over the same range
  (normalised to 700 K at rho = 30: 1.000, 0.822, 0.565, 0.298, 0.226 against
  HITEMP's 1.000, 1.004, 1.010, 1.019, 1.027). The loss is not redistribution
  between bands — the full-band index matches the bend and stretch indices — and
  it is density-independent, so it is not a normalisation error, which would hit
  every band equally. The identification is exact rather than assumed: `ac()` in
  the generation script returns the **unbiased** estimator `a / (n - lag)` of the
  mean-removed dipole ACF, with no C(0) normalisation and no temperature factor,
  so the delivered `alpha` is exactly `w (1 - e^-x) FFT[C]` and multiplying by
  Q_harmonic gives `(hbar/kT) w^2 FFT[C]` and nothing else.

  The obvious objection does not rescue it. Equipartition fixes <Qdot^2>, not
  <|grad mu|^2>: on an anharmonic dipole surface the sum rule is
  `beta <|mu_dot|^2> = <|grad mu|^2>` averaged over the configurations actually
  sampled, which need not be strictly T-independent. But the monomer test says
  grad mu **grows** with bond length (sym 2.84 at 0.957 A to 4.01 at 1.00 A), so
  hotter sampling should push p *positive*. The observed -1.5 is anomalous under
  either reading, and the band dependence settles it in any case: a global
  prefactor error would give the same p in every band, and the rotation band's
  p ~ -0.2 against the vibrational -1.5 rules that out. **The rotation band, which is near-classical here
  (intensity-weighted <x> = 0.64 at 700 K falling to 0.35 at 2000 K), roughly
  obeys the rule; the bend and stretch, where x runs 1.0-8.2, do not.** The defect is upstream of the lineshape and no QCF can reach it.

  **So "delivered is closest" is arithmetic, not physics.** Q_harmonic is larger
  at low T (3.42 at 1600 cm^-1, 700 K; 1.68 at 2000 K), and dividing it back out
  cancels part of a temperature dependence the spectra should not have in the
  first place. That is the whole of the delivered form's advantage. Forcing the
  bend ratio to 1 point by point would need an intensity factor running from
  **1.50 at 700 K to 0.71 at 2000 K** (at rho = 30) — the dipole model would have
  to go from 50% too strong to 29% too weak in the bend — which a fixed monomer dipole
  surface cannot do. The one T-dependence the monomer test does supply goes the
  wrong way: hotter states sample longer bonds, so the bend factor should rise
  from 1.25 towards 1.32, making the corrected ratio *smaller* at high T. It is
  also 5% against a factor of 2.1.

  **The rotation band says something separate and cleaner.** Both sides obey the
  sum rule there, so the comparison is well posed — and the first-principles
  rotation band is a **flat factor of 2.6 too weak** after the dipole factor
  (0.388, range 0.327-0.441, hot/cold 0.93, and +0.013 per 1000 K at fixed
  rho = 30, i.e. no trend at all). No correction closes it: the largest,
  Schofield, reaches 0.50. This is an intensity deficit in the model's permanent
  dipole, not a lineshape question.

  **Two consistency notes.** (i) Neither side carries MT_CKD here — the loader
  returns the stored kappa — so the line list is lines only while the
  first-principles spectrum, being bulk MD, in principle contains the continuum.
  That asymmetry flatters the first-principles side. Adding the self-continuum to
  the rotation denominator is worth 0.5-69% of the line integral (largest cold
  and dense) and deepens the deficit from 0.32-0.43 to **0.24-0.39**, so the
  factor of 2.6 is conservative. (ii) These band integrals are immune to A71's
  wing cut-off, the largest line-list systematic: RADIS renormalises every line
  profile to unit area (`lineshape /= area`, unconditional), so across the whole
  25-500 cm^-1 truncation range kappa_P moves by **0.5-9.0%** (0.5% at the
  nearest state, 1000 K and 30 kg m^-3) where kappa_R moves by up to **21.8x**.

  **Plain statement.** Of the four forms, the delivered one — Q = 1 on the
  lineshape, as shipped — brings the corrected bend ratio closest to 1 (median
  0.81-0.86 across the factor range, 0.08 dex from 1 at its midpoint) and is also the flattest in temperature. But none of the
  three corrections removes the temperature trend, and none can: all three are
  >= 1 and fall with T, while the ratio they act on already falls with T, so each
  one steepens it. The correct answer to the question as posed is that no
  correction achieves it, because the trend is not a lineshape effect. The
  first-principles bend and stretch band strengths violate the classical sum rule
  by a factor of 4.4 across 700-2000 K while the line list holds it to 2.7%; the
  delivered form looks best only because 1/Q_harmonic happens to cancel part of
  that spurious temperature dependence. **The quantum correction is not the free
  parameter that reconciles the two sources.** Until the sum-rule violation is
  traced upstream in the dipole time series, no lineshape correction can be
  calibrated from this comparison, and the delivered Q = 1 should stand — as the
  form that requires the fewest unearned assumptions, not as one that has been
  shown correct. A76 opens that as the question A77 should answer.

## 2026-09-30 (the dilute noise floor, measured rather than assumed)

- **A75 — There is a floor and it is three orders of magnitude below the
  spectra. A52's worry does not survive contact with the lag-test data.**
  `tools/floor_check.py` -> `outputs/floor_check.json`. Existing data only.

  A52 built a floor from an amplitude ratio the campaign quoted — 5.75 against
  12.4 mDebye/atom, so r^2 = 0.215 of the dipole power — spread white over the
  0.5 fs Nyquist bandwidth, and found that the resulting floor **exceeded** the
  delivered kappa in the transparency window. A52 recorded that as a
  contradiction, because an additive white residual cannot exceed the signal it
  is part of, and made no absolute reliability claim from it.

  **The contradiction was in the input, not the construction.**
  `ed_data/arc/scw/lagtest_v2.json` carries the residual's own autocorrelation,
  and the residual is mostly **not** white: ACF at the first lag is 0.9365 at
  (800 K, 61 kg m^-3) and 0.9552 at (1000, 130), so only **4.5-6.4%** of its
  variance is uncorrelated. Taking the uncorrelated part,
  r^2 = residual_variance x (1 - ACF_1) / signal_variance, the two cold cells
  give **7.26e-5** and **6.79e-5** — agreeing to 7%, so this does not rest on
  one cell. (The two 2250 K cells give 3.0e-3 and 4.3e-2; they are hot and dense
  and are not used for the dilute corner.)

  **Both differences from A52 are needed to explain it, and the amplitude ratio
  is the larger.** Through the same production normalisation:

  | r^2 | floor in the window (cm^2/g) | kappa(2550)/floor |
  |---|---|---|
  | A52, 0.215 | 143 - 492 | 0.02 - 0.13 |
  | A52's ratio x the uncorrelated fraction, 0.0116 | 7.8 - 26.6 | 0.35 - 2.33 |
  | lag test, 7.0e-5 | **0.046 - 0.157** | **58 - 386** |

  Correcting only the uncorrelated fraction leaves the floor at the kappa level
  at 9 of the 17 points. The lag test's own total-residual ratio at the cold
  cells is 1.3e-3 against A52's 0.215, a factor of **162**; the uncorrelated
  fraction is a further **18.5**. A52's 5.75/12.4 mDebye/**atom** and the lag
  test's resid/signal are not stated to be the same quantity, so both are
  reported rather than reconciled.

  **The verdict, over the 13 dilute points and the nine vertices of the ten
  strongest boundaries' triangles** (17 states):
  kappa(2550) sits **58-386x** above the floor, the 1000-1200 cm^-1 gap
  **440-8400x**, and kappa_R **185-458x**. **No point is within a factor of two
  of the floor; the nearest is 58x.** The Rosseland relative floor
  (1e-6 x max = 0.0013-0.007 cm^2/g) is two to three orders below the true window
  minimum and never binds. The minimum of kappa(nu) falls at 2593-2772 cm^-1 at
  every point — the floor of the transparency window, resolved on 200 samples,
  not an artefact of the 20 cm^-1 smoothing.

  **What the delivered spectra carry, read from the code.**
  `ed_data/arc/water-scw-mace/convergence/04_spectra/spectrum_convergence.py`
  line 28 is
  `alpha = w * (1 - np.exp(-HBAR*w/(KB*T))) * I` with `I` the classical dipole
  autocorrelation spectrum: the detailed-balance factor and **no quantum
  correction on the lineshape**. The ARC prose calls that factor itself "the
  quantum correction", which is where the name "harmonic" comes from in
  discussion; the code is the authority and it is Q = 1.

  Applying the three standard lineshape corrections raises kappa_R by x1.7-2.0
  (standard), x2.2-4.8 (harmonic) and x2.5-10.2 (Schofield), and the line list is
  quantum by construction and gets none — so the first-principles/HITEMP
  Rosseland ratio can only **grow**: 1.08-37.8 delivered, to 1.9-73, 2.4-156 and
  2.7-273. The choice cannot make the excess go away. What the spread measures is
  how much the cold spectra depend on an unverifiable choice, and it is largest
  exactly where the strongest claims sit. A76 calibrates it.

  **Caveat.** The uncorrelated fraction is a small difference of large numbers
  (1 - 0.955) on 200 pairs per lag and is quoted to one figure. The floor uses
  each spectrum's own total dipole power as its denominator, so it is a per-point
  quantity, not a constant per density.

## 2026-09-29 (four corner cells: the grid is 52 points)

- **A74 — Four cells at (700, 10/30/60) and (800, 30) close the gap the boundaries
  were falling into. The first-principles grid is 52 points.**
  `tools/corner_cells.py`, `tools/corner_boundaries.py`.

  **Checked before the marker, not after.** `units_marker` is what stops a legacy
  factor-1000 spectrum loading (A50), so an unstamped file is refused by the
  loader by design; the checks therefore read the npz directly and only then is
  `units_fixed = cm2_per_g:A74` written. On the common band:

  | cell | kappa_P | kappa_R | vs 1000 K at the same rho | the 48-point model there |
  |---|---|---|---|---|
  | (700, 10) | 423.1 | 28.65 | x1.12 | 379.1, clamped to (1000, 10) |
  | (700, 30) | 404.7 | 18.04 | x1.05 | 384.0, clamped to (1000, 30) |
  | (700, 60) | 404.7 | 19.30 | x0.97 | 468.9, clamped to (650, 100) |
  | (800, 30) | 400.2 | 22.90 | x1.04 | 384.0, clamped to (1000, 30) |

  and (800, 30) against (700, 30) is x0.989 in kappa_P. Every cell is within 12%
  of its 1000 K neighbour against a refusal threshold of a factor of 10, so
  nothing resembling the A50 units error is present. Three further checks passed
  that the brief did not ask for and that a factor-of-10 test would not catch:
  `alpha / rho` reproduces `kappa` exactly, `kappa / sigma` equals the production
  files' constant, and each spectrum spans the band. **All four states were
  clamped before**, and the embedded `planck`/`rosseland` scalars were not used —
  they are full-range (1.1-6000 cm^-1), not band-limited (A4).

  **The four new crossover rows** (line list read off its own table, as all 52
  rows are):

  | state | P/bare | R/bare | P/+MT_CKD | R/+MT_CKD | MT_CKD share |
  |---|---|---|---|---|---|
  | (700, 10) | 1.247 | 37.84 | 1.113 | 28.75 | 0.240 |
  | (700, 30) | 1.192 | 21.47 | 0.876 | 11.89 | 0.446 |
  | (700, 60) | 1.191 | 22.09 | 0.664 | 8.91 | 0.597 |
  | (800, 30) | 1.287 | 14.41 | 1.049 | 12.00 | 0.168 |

  Three of the four sit on the 90-point line-list grid, so the line list was also
  taken **directly** at those states as a check: table and direct agree to
  **1.0000** at all three, and to 0.9925 in kappa_R at (700, 60), the one that
  needed a new spectrum (`hitemp_extra/T700_rho60.npz`). The 90-point line-list
  grid is unchanged. Two of these rows have a **sub-unity** Planck ratio against
  the continuum-inclusive reference — 0.876 and 0.664 — and the continuum carries
  up to 60% of the combined Rosseland mean there.

  **All six boundaries that were falling between computed points are now inside**,
  and the category counts go **43 -> 49** on computed data, 18 and 41 unchanged,
  and **6 -> 0** between computed points. Profiled by the same code path the ten
  strongest use, all six still exceed their flip threshold — but **the clamp was
  wrong in both directions**:

  | case | T, rho | FP/LL before | after | threshold |
  |---|---|---|---|---|
  | M_c=5 wmf 0.20 T_eq 500 T_int 30 | 739.1, 54.6 | 32.61 | **18.70** | 1.034 |
  | M_c=5 wmf 0.50 T_eq 500 T_int 30 | 732.2, 54.1 | 34.09 | **19.26** | 1.122 |
  | M_c=10 wmf 0.05 T_eq 500 T_int 50 | 710.7, 11.0 | 14.26 | **34.13** | 1.100 |
  | M_c=10 wmf 0.50 T_eq 500 T_int 50 | 730.2, 11.5 | 12.52 | **28.35** | 1.063 |
  | M_c=15 wmf 0.05 T_eq 500 T_int 50 | 712.5, 13.1 | 13.86 | **31.30** | 1.136 |
  | M_c=15 wmf 0.20 T_eq 500 T_int 50 | 729.5, 13.1 | 12.42 | **28.44** | 1.044 |

  At rho ~ 54 the clamp reached a colder, denser point and **overestimated** the
  first-principles opacity by 1.8x; at rho ~ 11-13 it reached (1000, 10) and
  **underestimated** it by 2.3x. Neither error was visible as anything but a flag.

  **All ten of the strongest boundaries changed triangle**, two of them from
  clamped to interpolated, and the dispute ratios fell across the board — the
  maximum from 54.27 to **47.08**, and the two formerly clamped from 34.09/32.61
  to 19.26/18.70. All ten still exceed their thresholds, so ED Table 2 is now
  **10 interpolated, 0 clamped**, ratios 6.15-47.08.

  **What the merge does NOT do: move the boundaries.** Their positions come from
  the line-list opacity, which is untouched. The stored scan CSV differs from
  today's by up to 1.2 K in the boundary temperature and 1.3e-4 R_earth in radius,
  which looked at first like the merge. It is not: re-running the scan with the
  four cells **moved back out** reproduces the same offsets against the stored
  CSV, and two runs at the same state are bit-identical. The offsets predate this
  task; the shooting solver is warm-start sensitive at 1e-7 in R_p, which moves
  which sample of a 1500-point profile is the deepest radiative one. Recorded
  because it was measured, not because it matters here: a 1.2 K jitter cannot
  produce a 6-boundary category change.

  **The water-world fiducial's triangle changed**, from (650,100)-(800,100)-(1000,60)
  to (650,100)-(700,60)-(800,100), so E2 was rerun. The composite-against-line-list
  differences it reports move by about a tenth of themselves: the radius difference
  at 10 Gyr from 0.00121 to **0.00136** (+12.6%), the entropy difference from
  0.01346 to **0.01510** (+12.2%), the cooling-time ratio from 1.00875 to
  **1.00993**, and the track end from 387 to 417 Gyr. In absolute terms the
  fiducial is still a 0.14% radius effect; what moved is the size of a small
  number, not the conclusion drawn from it.

  **Density scaling on the two new isotherms.** 700 K now has three densities
  (10/30/60), all below 100, so it has no dense segment to fit; its local kappa_R
  indices are **0.58** (10 -> 30) and **1.10** (30 -> 60) — the dilute limit, as
  the line list shows there. 800 K gains rho = 30 below its 100-500 run, and its
  indices are 1.48 (30 -> 100), 1.97, 1.52, 1.72, with the fitted dense slope
  unchanged at 1.67 +- 0.15 because the fit starts at 100. **The 30 -> 100 pair
  crosses the N = 216 to N = 64 boundary** and so measures that change as well as
  the density; the per-isotherm homogeneity flag only looked at rho >= 100 and
  would have missed it, so an explicit per-pair flag was added. It shows the same
  problem at 1000 K (60 -> 100), 1500 K (three pairs) and 2000 K (60 -> 150).

  **Three gates re-pinned, all narrowings the new cells caused.** The dilute row
  count 9 -> 13; `floor_rosseland_window_share` from > 0.8 to > 0.70, because
  (700, 10) is the only one of the thirteen below 0.8 at 0.724; the cold
  Rosseland ratio floor from > 20 to > 12, because the new cold **dilute** rows
  run 14.4-37.8 where every cold row used to be dense; and the FP/HITEMP
  window-contrast separation from >= 15 to >= 10, because the new cells have the
  deepest first-principles windows anywhere cold (contrast 14.8-22.4 against
  4.8-11.2 for the cold-dense rows) and the separation narrows to 11.2 at
  (800, 30). Each claim survives with less margin; each bound is moved to what
  the data give rather than left where it was. Three more in the dilute-floor
  gates for the same reason: the cell count 9 -> 13, the Rosseland window share
  from > 0.8 to > 0.70 ((700, 10) at 0.724), and the fraction of the band cores
  above the floor from > 0.5 to > 0.35 ((700, 30) at 0.394). The contrast that
  carries A52's claim — under 5% of the *window* above the floor against 39-75% of
  the *cores* — is untouched. And two in the in-box clamp accounting, both
  improvements: one cell (M_c = 5, wmf 0.05, T_eq = 500, T_int = 30) now clamps on
  **0.2%** of its in-box samples and clears the 1% limit, so "all 108 violate the
  unbounded corner gate" becomes 107 of 108, and the full-box violation count goes
  28 -> 27. The ceilings have not moved; a hole below them has been filled, and
  this is the first cell to notice.

  **The release fit's dilute regime changes, the dense one does not.** With four
  more points and a temperature range that now starts at 700 K rather than 1000,
  kappa_0 goes 28.9 -> **31.6**, a from +0.256 to **+0.134**, b from +1.537 to
  **+1.013**, and the RMS from 0.093 to **0.124** dex over 13 points. The dense
  regime is untouched (39 points, identical coefficients). The MESA-format grid
  gains coverage: **100 of 225 cells** inside the hull, up from 93.

## 2026-09-26 (rho^2 scaling; wing cut-off; the AQUA entropy typo; release tables)

- **A70 — The absorption is faster than rho^2 in the window, and only there.**
  `tools/rho_scaling.py` -> `outputs/rho_scaling.json`, ED Fig. 9.
  Every continuum model on offer is a pair quantity, so alpha ~ rho^2 and
  alpha/rho^2 = kappa/rho is flat. Fitting d ln alpha / d ln rho over rho >= 100
  along each isotherm with at least three densities:

  | source | 2550 cm^-1 window | rotation-bend gap | kappa_R |
  |---|---|---|---|
  | HITEMP | 0.64-1.01 | 1.07-1.10 | 0.99-1.04 |
  | HITEMP + MT_CKD | 0.65-1.99 | 1.07-1.95 | 1.02-1.97 |
  | first principles | **1.58-2.76** | 1.53-2.06 | 1.53-2.22 |

  The two line-list rows are the control and they behave: bare HITEMP is rho^1
  everywhere (isolated lines), and MT_CKD's continuum is rho^2 at 500 K and falls
  to rho^1 by 1250 K, which is where it stops dominating. Against that, the
  first-principles absorption at 2550 cm^-1 reaches **2.76 +- 0.20 at 2500 K** and
  2.42 +- 0.13 at 3000 K — **3.7 and 3.1 sigma above 2** (errors from the pooled
  sigma). Nothing else comes close: kappa_R's largest excess is 2.22 +- 0.12 at
  2500 K, 1.9 sigma, and the rotation-bend gap never exceeds 2 by more than 0.5
  sigma. The excess is specific to the transparency window, which is where a
  continuum shows itself and where a line list has least to say.

  **Scatter, two ways.** The brief's estimator is the 729.344/810.382 kg m^-3 pair
  at 3000 K, the closest-spaced pair on the grid, which gives
  sigma_ln_alpha = 0.33 for kappa_R and 0.59 at 2550 cm^-1 — on one degree of
  freedom. Pooling the residuals of a quadratic in ln rho over every isotherm with
  four or more densities gives 0.19 and 0.33 on 23. The pair is the more
  pessimistic of the two for kappa_R and for the 2550 cm^-1 window — the two the
  claims rest on — and is used for the per-pair breakdown; the pooled value is
  used for the fitted slopes, and both are recorded. It is **not** conservative
  for the rotation-bend gap (0.118 against a pooled 0.141), which is harmless
  because no breakdown density is quoted for that diagnostic, and is pinned by a
  gate so the exception stays visible.

  **Per-isotherm breakdown density** (first adjacent pair whose index exceeds 2 by
  more than the pair-derived noise at that spacing), kappa_R: 650, 800, 1000 and
  2000 K, not within the grid; 1500 K at 1200 kg m^-3; 2500 K at 1521; 3000 K at
  1044. At 2550 cm^-1: 1500 and 2000 K at 1200, 2500 K at 1521, 3000 K at 1044.

  **Two run parameters are not uniform over the 48 points**, which matters for a
  window measurement: N = 216 molecules at rho <= 60 on the 1000, 1500 and 2000 K
  isotherms and 64 elsewhere, and the spectral smoothing is 15 cm^-1 at four of
  the 1500 K points and 20 cm^-1 at the rest. Restricting the fits to rho >= 100
  removes the N split; the smoothing split survives on the **1500 K** isotherm
  alone, which is flagged in the output. Neither of the significant results
  (2500 and 3000 K) sits on a mixed isotherm.

  **Lambda depends on density and not at all on temperature.** With
  lambda_mfp = 1/(sqrt(2) n sigma), sigma = pi d^2 and d = 2.8 A,
  Lambda = 1/(sqrt(2) pi d^3 n): **Lambda = 1 at rho = 307 kg m^-3** and
  Lambda = 0.1 at 3067. They are vertical lines on the (T, rho) plane, drawn as
  such. Converted through AQUA, the published tables stop well short of them:
  Freedman et al. (2014) at 300 bar reaches rho = 68 kg m^-3 at 1000 K and 32 at
  2000 K (Lambda 4.5 and 9.4); Kempton et al. (2023) at 1000 bar reaches 254 and
  108 (Lambda 1.2 and 2.8). Interpolating each limit linearly in T between those
  two isotherms and holding it flat outside them, **41 of the 48 computed points
  and 18 of the 108 boundaries are denser than the Freedman limit**, and 31 and 6
  respectively are denser than the Kempton limit — so the states this study
  computes are mostly outside the documented range of the tables the field uses,
  while the boundaries themselves mostly are not.

- **A71 — The line-wing cut-off is worth more than the broadening convention.**
  `tools/wing_cutoff.py` -> `outputs/wing_cutoff.json`, ED Fig. 3d.
  RADIS truncates each line profile at `truncation` cm^-1 and renormalises the
  remainder to unit area, so strength is conserved and shape is not. Varying only
  that cut-off — 25, 50 (delivered), 100, 300, 500 cm^-1, everything else at the
  delivered convention — the 50 cm^-1 row reproduces the delivered kappa_R
  **exactly** at four of the six points and to 1e-9 at a fifth, which is what
  makes the rest of the scan a measurement rather than a re-derivation. The
  sixth agrees to 1.0e-3, for the reason below.

  Two efficiencies worth recording, one of which is not quite free. The spectral
  range is cut to [50, 4600] cm^-1, which is 3.4x faster; with
  `neighbour_lines = 0` a line above 4000 + truncation cannot reach the band, so
  no *line* is lost for truncation <= 500. It is not exact all the same: RADIS
  builds its lineshape-interpolation grid from the set of loaded lines, so
  dropping 4600-11000 cm^-1 perturbs every profile a little. Measured rather than
  assumed — the reduced range reproduces the delivered kappa_R exactly at four of
  the six points, to 1e-9 at a fifth, and to **1.0e-3** at (3000 K, 2228.7), the
  state whose Lorentzian widths span the widest range. Recomputed there over the
  delivered [50, 11000] cm^-1 the agreement is exact. 0.1% is two orders below the
  effect being measured, so the scan stands, but the reduced range is an
  approximation with a number on it, not a free one.

  And the cost turns out to be **flat in the truncation** with RADIS's default
  line-shape optimisation, not the quadratic its documentation warns of — 170 s
  per spectrum at 25 and at 500 cm^-1 alike.

  **Renormalisation cannot be switched off.** `lineshape /= area` is unconditional
  at five places in `radis/lbl/broadening.py` (1526, 1592, 1770, 1892, 1989) in
  0.16.3, and no `SpectrumFactory` argument or `radis.config` key exposes it. The
  brief's "where RADIS allows" resolves to nowhere, and the library is not patched
  to find out.

  **The numbers, and they are large.** kappa_R rises monotonically with the
  cut-off at every point. Across 25 to 500 cm^-1:

  | state | kappa_R (cm^2/g) | spread | kappa_P spread |
  |---|---|---|---|
  | 650 K, 300 | 0.545 - 11.88 | **21.8x** | 1.09x |
  | 1000 K, 129.7 | 4.02 - 18.06 | 4.50x | 1.02x |
  | 1500 K, 471.3 | 20.26 - 67.42 | 3.33x | 1.01x |
  | 2000 K, 1200 | 47.94 - 116.3 | 2.43x | 1.007x |
  | 1000 K, 30 | 3.97 - 8.79 | 2.21x | 1.005x |
  | 3000 K, 2228.7 | 102.4 - 152.2 | 1.49x | 1.02x |

  kappa_P barely moves — 0.5% to 9%, and 9% only at the coldest point — because an arithmetic mean does not care where
  inside the band the redistributed area lands. kappa_R is a harmonic mean set by
  the transparency windows, and filling them is exactly what a truncated,
  renormalised profile does. The largest spread is at the **cold, dense** corner,
  where the Lorentzian half-width most exceeds the cut-off.

  **This is the largest line-list systematic found so far, and it eats into the
  headline.** The first-principles / line-list Rosseland ratio at these six states
  spans **1.3 to 180** across the cut-off range: at (1500 K, 471.3) it falls from
  4.32 at 25 cm^-1 to **1.30** at 500, which is agreement, not a factor of four.
  At (3000 K, 2228.7) it stays 26.9-40 and at (650 K, 300) it stays 8.25-180, so
  the excess does not vanish everywhere — but the line-list reference is not a
  fixed quantity to be exceeded by a definite factor, and any statement of the
  form "the first-principles opacity is N times the line list" has to name the
  cut-off. The delivered 50 cm^-1 is RADIS's default, not a physical choice.

- **A72 — AQUA's entropy is the uncorrected one, and the cost is 3% in the label,
  3-6% in the age.** `tools/aqua_entropy_check.py` -> `outputs/aqua_entropy_check.json`.

  The table on disk is `aqua_eos_rhot_v1_0.dat`, SHA-256 `c78c4523…`. The upstream
  repository offers only the v1_0 tables — one release, tagged 2020-09-02, last
  touched 2020-12-30 — so it predates the August 2024 cut-off Aguichine et al.
  (2025) footnote 3 names, and carries the uncorrected Mazevet et al. (2019)
  entropy: a sign error on the first two terms of F_T in their Eq. (13), and a
  reference entropy S0 of 4.9 rather than 9.8 n_at k_B.

  **The corrected values are read, not reimplemented.** PALEOS (Zenodo
  10.5281/zenodo.22776069) re-grids AQUA with the sign fix applied analytically
  and states that rho and P are untouched. Measured along the water-world
  fiducial's own profiles, at five ages:

  | age | s (AQUA) | s (corrected) | shift | E shift | age shift at fixed s |
  |---|---|---|---|---|---|
  | 0.001 Gyr | 13942 | 14171 | +1.65% | +0.58% | +0.001 Gyr |
  | 0.01 | 13382 | 13616 | +1.75% | +0.54% | +0.009 |
  | 0.1 | 12669 | 12911 | +1.91% | +0.47% | +0.051 |
  | 1 | 11002 | 11273 | +2.46% | +0.35% | +0.30 |
  | 10 | 9681 | 9990 | +3.19% | +0.25% | +8.8 |

  **Half the envelope is in the affected region**: 49.5% of the profile points on
  this track fall in AQUA's supercritical/superionic block, the part taken from
  Mazevet et al. (2019), at every one of the five ages.

  Two independent checks passed: the interpolated density agrees to 1.2e-4, and
  the adiabatic gradient differs by a **median 3.9e-4** (the maximum, 3.9e-2, is
  one point at a phase boundary) — so the correction really does leave the
  derivatives alone, as stated.

  **What it costs depends on which route you take.** The cooling age integrates
  dt = dE/L with E the envelope internal *energy* and the adiabat set by grad_ad,
  so ages, radii and luminosities never read the entropy column; entropy is a
  label. Through the energy the total radiated energy falls **2.75%** and ages at
  1-10 Gyr shorten by **3.7-6.5%**. Through the label — which is what
  `structure.build_at_entropy` and the fixed-entropy cooling-time ratio use — a
  3% entropy shift is worth **+8.8 Gyr at 10 Gyr**, because ds/d(age) has almost
  flattened by then — that one is Delta s divided by a derivative going to zero,
  so it is an order-of-magnitude statement, not a three-figure one. Nothing is
  substituted: the model still reads the table it has, and this is what that
  costs.

  **G5's reference is now Aguichine et al. (2025)**, matched age for age at their
  1 microbar radius with their Earth-like core (cmf 0.325). Twelve cases across
  three equilibrium temperatures, two water fractions and two masses:
  **median |dR/R| 0.0405**, worst 0.0591, improving with age (0.0456 at 1 Gyr,
  0.0319 at 10). Against the 2021 relation's 0.0844 that is half the residual,
  and **all 36 are positive** where all eight of the 2021 residuals were negative
  — so the two disagree about the sign, not just the size. The 2021 comparison is
  kept as the second row.

  Two systematics measured rather than assumed: their own M-star and G-star tracks
  differ by **0.86%** at 5 Gyr, a fifth of the residual, so the choice of the G
  rows is not what is being measured; and their R_20mbar/R_1000bar is 1.013 for a
  20 M_E water world, so the radius convention is worth about 1%, not tens of
  percent. The comparison reads their 1 microbar radius, which is the level this
  study's anchor already used.

  **The two rows use different opacity sources, and it does not matter.** The 2021
  row runs `hitemp` (full radiative-zone coverage; the radius there is
  EOS-dominated), the 2025 row the production `composite`. Measured at the first
  case, the source moves the radius by +0.007% at 1 Gyr, -0.30% at 5 and -0.62% at
  10, so the 2025 residual would read +2.90 / +2.05 / +1.63% under `hitemp`
  against +2.90 / +2.36 / +2.27% under `composite` — still positive, still about
  half the 2021 residual. The disagreement in sign between the two anchors is not
  the opacity source.

  Two limits stated rather than discovered later: the gate reads
  `outputs/aguichine2025_compare.json` and the reference tracks live outside the
  repository, so it **skips on a fresh clone** until both are built; and this
  model's track starts at T_int = 200 K, about 0.01 Gyr after a hot start, which
  together with their 1 Myr formation offset is negligible at the 1-10 Gyr the
  comparison uses.

- **A73 — Release tables in MESA's format, which is not the same as a MESA table.**
  `tools/release_tables.py` -> `data/release/`. The format is taken from the
  reader (`kap/private/load_kap.f90`, `Read_Kap_X_Table`), with
  `form = kap_table_fixed_metal_form = 1` and `version >= min_version = 37`, and
  validated by a round trip through a strict parser of that spec (agreement
  4.6e-5 dex, which is half a unit in the written last place).

  MESA's R is rho/T6^3 with **rho in g cm^-3**, so a steam envelope at
  10-2229 kg m^-3 and 650-3000 K sits at **log R 6.0-9.5**, where MESA's own
  tables run -8 to +1. No MESA table covers this regime — which is the reason to
  release one — but the files are a format, not a drop-in replacement, and say so
  in their own header. 93 of the 225 cells lie inside the 48-point hull; the rest
  carry -99 (the OPAL no-data value; MESA's parser requires a number in every
  column) and are blank in the companion CSV. Values are the model's own
  interpolant — linear on the Delaunay triangulation in (T, log10 rho) over
  log10 kappa — sampled onto the regular grid, so the release reproduces the
  paper's own table to 4.6e-5 dex rather than being a second fit to it.

  The two-regime fit kappa_R = kappa_0 (rho/100)^a (T/1000)^b:

  | regime | kappa_0 (cm^2/g) | a | b | RMS | n | validity |
  |---|---|---|---|---|---|---|
  | dense (rho >= 100) | 32.5 | +1.023 | +0.559 | 0.192 dex | 39 | 100-2229 kg m^-3, 650-3000 K |
  | dilute (rho < 100) | 28.9 | +0.256 | +1.537 | 0.093 dex | 9 | 10-60 kg m^-3, 1000-2000 K |

  The dense exponent a = 1.02 is the same statement as A70's kappa_R slope near 2
  (alpha = rho kappa), from an independent fit. **The dense fit is poor**: 0.192
  dex is a factor of 1.56, and a single power law cannot follow the curvature
  ED Fig. 9c shows across a 22x density range. It is published as a convenience
  with its residual attached, not as a substitute for the table.

  `data/release/` is exempted from the SPEC §1 "data is never committed" rule:
  these are outputs of this code, not raw inputs. The exemption needed
  `/data/*` rather than `/data/` in `.gitignore`, because git does not descend
  into an excluded directory and the negation was unreachable.

## 2026-09-26 (housekeeping: the pull, the model binary, one ledger key)

- **A69 — The ARC pull and the dipole-model binary now live in the repository;
  one ledger key renamed.** Three decisions, all the owner's:

  `ed_data/` is tracked, all **92** files. The `*.png` rule is lifted inside it
  (`!ed_data/**/*.png`), because 41 of those files are ARC's own figures and they
  are the provenance for panels this repository cannot rebuild — ED Fig. 1 and the
  functional-sensitivity numbers read from this tree, so a clone without it cannot
  reproduce them. 7.1 MB.

  `monomer_test/water_dipole_v3.model` (14 MB) is tracked through **git-LFS**
  (`.gitattributes`: `monomer_test/*.model filter=lfs`). `monomer_test/eval_model.py`
  needs it and nothing else in the repository can regenerate it. There is no remote
  yet, so the pointer is local until one is added — a remote that does not serve LFS
  will reject the push.

  `opacity.grid_extension.*` is renamed **`opacity.hitemp_grid_rows_2250_2750.*`**,
  and `tools/grid_extension_effect.py` to
  `tools/hitemp_grid_rows_2250_2750_effect.py`. The old name read as the
  original/extension distinction the manuscript does not make (A62); it never meant
  that — it is the two **line-list** rows of A60, not the 48 first-principles
  points — but a name that has to be explained is the wrong name.

  Also corrected in A67: the pull holds **41** PNGs, not 33. 33 is its script count.

## 2026-09-25 (Extended Data; the Fig. 4 categories; what the ARC pull does and does not contain)

- **A67 — Extended Data: 8 figures and 2 tables, built from what the repository
  and the ARC pull actually hold.** `figures/ed_items.py`, `figures/ed_tables.py`,
  `figures/ed_index.py` → `figs/ed/`. Same `save_pdf` guard as the main figures, so
  the vocabulary rule, the 7 pt floor and the exact width are enforced here too —
  it caught four 6.9 pt labels of mine in ED Fig. 5 before they reached a PDF.

  **What the ARC pull contains is mostly code, not data.** 92 files: scripts,
  READMEs, campaign documentation and 41 PNGs, and exactly **four** data files —
  `scw/lagtest_v2.json`, `scw/md_pressure.csv`, `scw/opacity_surface_means.csv`
  and `benchmarks/reference/reference.json`. That decides which ED panels exist.
  Two notes on the CSVs, neither of which is used here: `opacity_surface_means.csv`
  carries **legacy ×1000 units** and full-band rather than common-band means (A50),
  and `md_pressure.csv` is MD-vs-AQUA pressure per state point, which is not the
  finite-size check ED Fig. 1c asks for.

  **Three disagreements between the text and the data, all in ED Fig. 1.**

  1. The text describes the 41 / 58 meV Å⁻¹ force parity as against **held-out**
     labels. Both `benchmarks/README.md` and `reference.json` state the opposite:
     these are *training-set* parity by regime, "optimistic vs true held-out
     error", with T > 3000 K frames excluded. The only held-out parity in the pull
     is **12.3 / 20.2** meV Å⁻¹ and it belongs to the **LDA** model, not this one.
     The values 41/58 are right; the description of them is not.
  2. The ~49 and ~376 meV Å⁻¹ of ED Fig. 1b appear **nowhere** in the pull. Panel
     b is a labelled placeholder.
  3. The finite-size sweep (ED Fig. 1c) has no output in the pull — only
     `finite_size.py` and the README's verdict ("flat over a 32× volume range,
     64 → 2048; 6400 OOM'd on 44 GB GPUs"). Panel c is a labelled placeholder.

  ED Fig. 1d is a **bound, not a distribution**: the pull records |drift| ≤ 0.017
  meV atom⁻¹ ps⁻¹ over 21 NVE runs of the LDA campaign against a gate of 2, not a
  value per state point.

  **ED Fig. 5 settles Task 26's question, and the answer is yes.** The monomer test
  has been run (`monomer_test/pyscf_monomer.json`, revPBE0 = 0.25 HF + 0.75 revPBE
  exchange + PBE correlation, three bases). Model / revPBE0 dipole-derivative ratio
  at aug-cc-pVTZ:

  | mode | at r_OH = 0.957 Å | at 1.00 Å |
  |---|---|---|
  | symmetric stretch | **2.84** | **4.01** |
  | antisymmetric stretch | 1.24 | 1.72 |
  | bend | 1.12 | 1.15 |

  The permanent dipole is fine — model 1.822 D against revPBE0 1.847 and
  experiment 1.855 — and the **bend** derivative is fine. The **symmetric stretch**
  derivative is nearly 3× too large at equilibrium and **4× too large at 1.00 Å**,
  the length the bonds sample at 1000 K. Band intensity goes as |dμ/dQ|², and the
  intensity proxy the campaign reports moves from **1.69 to 2.97** between the two
  references. That is the A61/A65 stretch excess, and it is anharmonic — worse
  where the spectrum lives than at equilibrium, which is exactly why the stretched
  reference was included in the geometry set (A64).

  **ED Fig. 7 qualifies a number the text quotes.** The anchor residuals are
  **8 of 8 negative**: the median *signed* residual is −0.0844, not ±0.0844. The
  model is systematically 4-18% smaller than the published relation, so
  "median |rel| 0.084" describes a one-sided bias, not scatter about the anchor.
  The worst is 0.179 and the 1-100 bar matching-pressure spread is 0.35%, both as
  the text has them.

  **ED Fig. 6b shows fewer boundaries than asked, for a stated reason.** The flip
  threshold ∇_ad/∇_rad needs both gradients at the point, and the scan CSV stores
  neither; only the ten profiled boundaries have it. The dispute ratio is available
  for all 43 on computed data. Panel b therefore shows the ten. With the brief's
  conservative divisor of 4.4 applied, the ratio still exceeds the flip threshold
  at **all ten** (bound 1.41-12.33 against thresholds 1.010-1.208).

  **Main Fig. 4 now carries the same four categories** as ED Fig. 6a — filled where
  the opacity table's own triangulation contains the boundary, open where the model
  clamps, with the clamp cause distinguished, and the counts in the legend:
  **43 on computed data, 18 T below the grid, 41 ρ below the grid, 6 between
  computed points.**

  **The vocabulary guard now covers the ledgers, not just the figures.**
  `save_pdf` scans the drawn text; it cannot see the JSON written beside it, and
  every ED number is merged into `outputs/manuscript_numbers.json`, so a note in
  the ledger reaches the read-through exactly as a figure label would. Four notes
  had gone round the guard that way. They are reworded, and
  `test_the_extended_data_ledgers_never_name_the_delivery_history` now gates both
  ED ledgers, exempting only tokens that are file paths — a ledger has to be able
  to say a number came from `docs/LDA_CAMPAIGN.md` without that counting as the
  manuscript adopting the word.

  **Five layout defects found by rendering each figure rather than reading its
  code.** The 7 pt floor and the width are checked on the file; overlap is not, and
  it is invisible in a thumbnail. ED Fig. 1d and ED Fig. 4a had the bold panel
  letter printed on top of a long y-label; ED Fig. 6a had it inside the legend;
  ED Fig. 5 had the two displacement axis labels running into each other, a y-label
  overflowing its panel, six two-line tick labels colliding, an annotation clipped
  at the panel edge, and no legend at all for its four sources and three modes.
  ED Fig. 5a,b now also share one y-scale, because the point of the pair is that
  the slopes grow with r_OH and separate autoscales hid it.

  **ED Fig. 7b drew a line that implied the opposite of the finding.** It marked
  both +8.4% and -8.4%, which reads as scatter about the anchor. Every residual is
  negative, so it now marks the signed median once, labelled "median -8.4% (8/8
  below)".

  **ED Table 1 and 2 keep trailing zeros at three significant figures** (1.30, not
  1.3), round rather than pad between 10^3 and 10^5 (4180, not 4182.78 and not
  `4.18e+03`), and use a mantissa and a power beyond that, because `%.3g` on a
  pressure of 2.83e7 prints all eight digits, which is not three significant
  figures. ED Table 2's LaTeX provenance column now names the clamp cause with the
  same words as ED Fig. 6a and main Fig. 4 ("clamped, between pts"), not just
  "clamped".

  **The file-level figure gates now cover the ED items.** `test_each_figure_is_a
  _vector_pdf` and `test_each_figure_is_exactly_the_column_width` ran over the four
  main figures only. Extending them caught six ED PDFs written without
  `set_manuscript_style()` — Type 3 text, default metrics — from an interactive
  rebuild of mine; they are rebuilt through the CLI. All eight ED pages measure
  183.00 mm.

- **A68 — Functional sensitivity: band positions only; the means cannot be
  quoted.** From `docs/LDA_CAMPAIGN.md`: the stretch band is red-shifted under LDA
  by **146 cm⁻¹ at 1000 K (3654 vs ~3800)** and **200 cm⁻¹ at 2000 K (3550 vs
  ~3750)**. The third campaign point, 3000 K / 108.8 kg m⁻³, has **no stretch-peak
  pair in the pull**; its thermostat handoff failed three draws and the shipped run
  achieved **2807 K**, used in the quantum correction and the Planck weighting.

  **κ_P and κ_R under each functional are not reported, for two independent
  reasons.** No spectra are in the pull — `lda_campaign` holds only scripts. And
  the campaign documentation states that LDA absolute band intensities and κ_R are
  **floor-limited and are not quoted**: the measured white floor is 0.6-2.6% of
  peak spectral density, and only band positions and the libration region (40-60×
  floor) are quotable. Computing a mean from a floor-limited spectrum would produce
  a number the campaign itself declines to quote. The files that would be needed
  are named in `outputs/manuscript_numbers.json` under `ed.functional_sensitivity`.

## 2026-09-24 (isolated-molecule dipole test; manuscript ledger; text checks)

- **A64 — The dilute band tilt points at the dipole derivative, and the isolated
  molecule is the test. `monomer_test/` is the ARC-ready bundle.** A61 found the
  first-principles absorption at 0.35-0.50 of the line list's in the rotational
  band, 0.7-1.1 in the bend and 2.4-3.9 in the stretch, with a stretch/bend
  balance 3.3-3.6x the line list's — **and the same at 1000, 1500 and 2000 K**. A
  ratio that is mode-dependent but temperature-independent is not a population
  effect; it is what a wrong dμ/dQ looks like. `water_dipole_v3` was trained
  almost entirely on dense configurations, so it may return a liquid-like dμ/dr
  when a bond stretches in the gas, and band intensity goes as |dμ/dQ|², so a
  stretch derivative too large by ~1.8x would produce the observed excess alone.

  **13 geometries, two references.** Equilibrium (0.9572 Å, 104.52°), then ±
  pairs for symmetric stretch (±0.01 Å), asymmetric stretch (±0.01 Å opposed) and
  bend (±1.0°), about equilibrium *and* about **r_OH = 1.00 Å** — the length the
  bonds sample at 1000 K. The second reference is not decoration: an error that is
  anharmonic (right at equilibrium, wrong where the spectrum lives) shows only
  there. Orientation is identical across all thirteen, so dipole *vectors*
  compare, not just magnitudes.

  **Both the model and CP2K are on ARC, so the deliverable is a bundle**, not a
  result. What was nonetheless verified from here:

  * all 13 geometries build, and `config_from_atoms` → `AtomicData.from_config` →
    `DataLoader` → `batch.to_dict()` runs for every one of them
    (`eval_model.py --dry-run`, which needs no checkpoint). **Only the forward
    pass and the `"dipole"` key are untested locally** — `torch`, `mace` and `ase`
    all import here, but `water_dipole_v3.model` does not exist on this machine.
  * `compare.py`'s arithmetic was checked against synthetic inputs carrying a
    known **1.8x** stretch bias: it recovered 1.800 on every stretch derivative,
    1.000 on bend, and an intensity-proxy ratio of exactly **3.240 = 1.8²**.

  Two safety properties are deliberate. `harvest_dft.py` obtains `dip()` by
  *reading* `build_dipole_xyz_v3.py` and executing it **only up to `out=[]`**,
  never importing it — everything after that marker is the harvest loop that
  rewrites the training set. And `make_cp2k_inputs.py` edits the MLWF template in
  exactly one place, the cell, with the same `re.sub(r"ABC\s+[\d. ]+", ...)` the
  label builder uses; the template is not reproduced or guessed, because a
  settings difference would surface as a model/DFT disagreement.

  The sbatch array wrapper follows the resources the brief specifies but is **not**
  copied from `htc_dipole_label3.sh`, which was not readable from here; the README
  says so and says to prefer that file's partition, account and `srun` form.

- **A65 — Two text numbers checked: the window value is right, and the band tilt
  is indifferent to state.** `tools/text_checks.py`.

  **(a)** At (2000 K, 1200 kg/m³), 2550 cm⁻¹: first principles **296.1** cm²/g
  against **9.521** for the line list at identical conditions — **×31.1**. The
  text's "9.5 vs ~300, ×31" is confirmed. The 1 bar-broadened line list gives
  **3.013** there, which is where the draft caption's "4" came from and is not the
  comparison the text should make (×98 against it). The first-principles spectrum
  is the delivered one, smoothed to 20 cm⁻¹ on a 1 cm⁻¹ grid, exactly as plotted.

  **(b)** The line list carries no ρ = 60 spectrum — its grid is 0.1, 1, 10, 30,
  100, 300, 1000 — so three were computed by `tools/extra_linelist_points.py`,
  which inherits its conventions from `generate_linelist_tables` **by import, not
  by copy**. They live in a separate directory and are deliberately **not** part of
  the production grid: a density column would change the triangulation and every
  ratio, figure and gate downstream, for a diagnostic that needs spectra rather
  than a table.

  | band | min | max | mean | spread | monotone in ρ at every T? |
  |---|---|---|---|---|---|
  | libration | 0.348 | 0.510 | 0.411 | 1.47× | no |
  | bend | 0.706 | 1.122 | 0.869 | 1.59× | no |
  | window | 1.060 | 2.891 | 1.691 | 2.73× | no |
  | stretch | 2.376 | 4.425 | 3.174 | 1.86× | no |

  **The finding survives the extension and sharpens.** The ordering
  libration < bend < window < stretch holds at all nine states; no band's ratio is
  monotone in density at any temperature; each band scatters by only **1.5-2.7×**
  about its own value over a six-fold density and two-fold temperature range,
  while the bands stay **5.7-9.7× apart from one another** throughout. A tilt fixed
  in frequency and indifferent to the state is what a wrong dipole derivative
  gives; a population or collision effect would track T or ρ. The continuum changes
  none of it — its largest share of any band here is **23.0%**, at (1000 K, 60).

- **A66 — `outputs/manuscript_numbers.json`: one file, 303 entries, nothing
  missing.** Every number the manuscript quotes from this repository, each with
  the script that produced it, the file it is stored in, and the commit that last
  touched that script — so a number traces to the code that made it rather than to
  the day it was read. Sections: opacity ratios by group and over all 48 points,
  the two direct computations, the six-point broadening sensitivity, continuum
  shares, agreement with the independent crossover, the band diagnostic, G5
  residuals, the 108-planet boundary ranges and the 10 analysed boundaries, the
  water-world fiducial, evolution differences and escape maxima, radiogenic
  floors, and every value in `figs/fig_numbers.json`.

  **Two gaps were closed rather than recorded.** The scan CSV has no boundary
  pressure column, but the boundary pressure *is* P(T_rcb, ρ_rcb) and AQUA is
  tabulated on (ρ, T) with P as an output (A2) — a lookup, not a fit, so it is
  computed. And the **reachable bounds** existed only in A56's prose; they are now
  columns that `plot_s9_case_evolution.py` writes, because the floor crossing age
  comes from the track and belongs with the track.

  **The reachable bound needed two corrections before it was right**, both caught
  before it reached the ledger. It is the difference **at** the last age the planet
  is still cooling, which is the earliest of: the radiogenic floor crossing, the
  end of the track, and the **10 Gyr reporting horizon**. Taking the *maximum* over
  that window instead gives 2.1% radius for the M_c = 15 water world against the
  0.08% it actually ends with — the difference is not monotonic in age. And
  omitting the horizon lets "end of track" mean **549 Gyr**, since the T_int grid
  reaches 12 K. With both fixed it reproduces A56 exactly: the water worlds at
  **0.078 / 0.121 / 0.237%** in radius and **0.96 / 1.35 / 2.15%** in entropy, and
  the old fiducial at **0.079% / 0.705%** by its 0.386 Gyr floor against
  1.65% / 35% at 10 Gyr — the inversion that picked the recommended case.

## 2026-09-24 (one grid; band attribution; HITEMP rows at 2250/2750 K; manuscript figures)

- **A60 — HITEMP grid extended to T = 2250 and 2750 K, and the new rows are
  windowed to ν ≥ 50 cm⁻¹ because the 0-50 cm⁻¹ line file is not cached.** A59
  showed the two 20-30% deviations from the delivered crossover analysis sit at
  exactly the two temperatures absent from HITEMP's grid, so those points were
  both ρ-clamped *and* T-interpolated. Adding the rows removes the T
  interpolation. `generate_linelist_tables.py` now carries T = 2250 and 2750 in
  `T_GRID_K` (13 temperatures), its `P_AQUA_PA` table is regenerated for the wider
  grid (the 11 original columns reproduce the committed values to 4.7e-5), and it
  gained `--temperatures`, `--densities` and `--wavenum-min` so a partial rebuild
  is one command. Two RADIS-version fixes came with it, both from A59: the slotted
  `SpectrumFactory` (pressure is a per-call argument now) and the registration
  abort on a fully-cached offline load.

  **The wavenumber floor is a deviation and is stated as one.** The delivered rows
  were computed over [10, 11000] cm⁻¹. The local HITEMP cache has no 0-50 cm⁻¹
  file and fetching one needs a HITRAN login, so the new rows are computed over
  **[50, 11000]**. With `truncation = 50` a line below 50 cm⁻¹ cannot reach past
  100 cm⁻¹, so the science band [100, 4000] is untouched — measured, not assumed:
  ∫κ dν over [10, 4000] and [50, 4000] differ by 0.06% at (2000 K, 300), the
  [100, 4000] integral is identical, and A59's [50, 11000] rerun reproduced a
  delivered [10, 11000] spectrum to 1e-6 on both means. What it *does* change is
  the auto-detected band: `detect_common_band(["fp", "hitemp"])` moves from
  [10, 6000] to [50, 6000]. Nothing in the science path uses it (the band is pinned
  to [100, 4000]), and the affected gate is updated with this reason rather than
  loosened. Re-running with `--wavenum-min 10` on a machine that has the file
  restores exact parity.

  **A scattered-grid property worth knowing.** With one new row present
  (2250 K, 0.1 kg/m³) and thirteen still to come, the FP/HITEMP Rosseland ratio at
  (650 K, 500) moved by 0.4%; with four present it was back to its original value.
  Adding a point to a Delaunay triangulation can change triangles well away from
  it, so a *partially* rebuilt grid gives numbers that belong to no grid at all —
  and intermediate states do not converge towards the answer, they are just
  different. Every output and figure here was regenerated only after all 14 rows
  landed, and `tools/hitemp_grid_rows_2250_2750_effect.py` refuses to run on a partial
  set.

  **THE RESULT.** All 14 spectra ran, none failed (~2.7 h). The before/after is
  captured by moving the added spectra aside and back, so both states come from
  the same code on the same day:

  | T [K] | ρ | FP/HITEMP Planck | FP/HITEMP Rosseland | FP/(H+MT) Planck | FP/(H+MT) Rosseland |
  |---|---|---|---|---|---|
  | 2250 | 175.5 | 1.407 → 1.407 | 1.413 → **1.370** | 1.334 → 1.327 | 1.399 → 1.356 |
  | 2250 | 1256 | 4.110 → 4.155 | 14.081 → **11.032** | 2.887 → 3.198 | 13.422 → 10.402 |
  | 2750 | 222.7 | 1.430 → 1.430 | 1.359 → **1.335** | 1.369 → 1.374 | 1.343 → 1.319 |
  | 2750 | 1841 | 11.306 → 11.448 | 34.737 → **30.216** | 9.272 → 9.821 | 32.374 → 27.881 |

  The Rosseland ratios fall by 1.8-21.7%; the bare Planck ratios move by at most
  1.3%. **Nothing else on the grid moved at all** — the largest change at any point
  not at an added temperature is 0.00%, so on the *complete* grid the effect is
  exactly local, and the grid-wide ranges are unchanged (bare Planck 1.290-17.721,
  bare Rosseland 1.077-213.432, continuum-inclusive Planck 0.212-15.708,
  continuum-inclusive Rosseland 1.076-38.449).

  **This closes the A53/A59 puzzle completely.** Deviation from the owner's
  independent RADIS crossover, by how κ_LL was obtained:

  | | n | Rosseland, before | Rosseland, after |
  |---|---|---|---|
  | both coordinates on the line-list grid | 3 | 0.24% | 0.24% |
  | T interpolated | 15 | 2.98% | **0.77%** |
  | ρ clamped | 9 | **29.65%** | **5.89%** |

  A53 attributed the 30% to the density clamp; A59 measured the clamp at ≤6% and
  identified the outliers as the two T-interpolated temperatures instead; removing
  that interpolation removes the outliers. The worst residual, **5.89% at
  (3000 K, 2229)**, is exactly the clamp error measured directly against a
  line-by-line run at the true density (A59). Every remaining disagreement with the
  independent analysis is now accounted for by a measured effect.

- **A63 — ARC's MT_CKD reference is bit-for-bit the file A58 built, so A58's one
  open question is closed.** A58 could not reach ARC (`scp` needs an interactive
  password) and built the self-continuum from AER's release 4.3 instead — the
  upstream primary source. The one thing that left open was whether ARC's
  `mtckd_ref.npz` holds the same coefficients. It has since been copied across and
  `tools/mtckd_ref_diff.py` compares them at three levels:

  | | result |
  |---|---|
  | keys and shapes | names differ (`wavenumbers`/`ref_press_mbar`/`ref_temp_K` vs ours), content maps one to one; 2003 points either way |
  | `wavenumber_cm`, `self_absco_ref`, `self_texp` | **max &#124;rel diff&#124; = 0** on all three |
  | `ref_press_mb`, `ref_temp_k` | 1013 mb, 296 K in both |
  | derived C(ν,T)·R(ν,T) at 296/650/1000/2000 K on [100, 4000] cm⁻¹ | **max &#124;rel diff&#124; = 0** at every temperature |

  **Identical, not merely close.** Nothing in Tasks 19-23 changes, and A58's model
  rests on exactly the coefficients ARC uses. ARC's file carries no version or
  provenance string; ours carries AER's title and the sha256 of the netCDF.

  **One honest limit on the version identification.** The *self* arrays alone
  cannot separate MT_CKD 4.0 from 4.1-4.3: the self continuum was last revised in
  4.0, and 4.1-4.3 changed only the **foreign** term, which neither file carries.
  So "4.3" rests on our file's provenance, not on the numbers. For a pure-steam
  envelope the distinction is empty anyway — the foreign term carries
  (1 − x_H2O) = 0.

- **A61 — The dilute FP/line-list excess is a frequency-dependent tilt, not an
  intensity scale.** At ρ = 10-60 the FP Planck mean exceeds bare HITEMP by
  1.32-1.53× while the Rosseland ratio falls to 1.08, and a band-integrated mean
  cannot say where that sits. `opacity/bands.py` integrates
  A = ∫ ν(1 − e^{−hcν/k_BT}) κ(ν) dν over the libration (300-1000), bend
  (1400-1800), window (1900-3100) and stretch (3200-4000) cm⁻¹ bands, each source
  on its own delivered grid but with the integrand interpolated onto the exact band
  edges so the limits are identical (FP is on 1 cm⁻¹, HITEMP on 0.01; a mask-only
  integral would use different limits for each).

  | T [K] | ρ | libration | bend | window | stretch |
  |---|---|---|---|---|---|
  | 1000 | 10 | 0.501 | 1.075 | 2.510 | 3.587 |
  | 1000 | 30 | 0.405 | 1.122 | 2.154 | 3.935 |
  | 1500 | 10 | 0.372 | 0.708 | 1.413 | 3.074 |
  | 1500 | 30 | 0.399 | 0.809 | 1.387 | 2.804 |
  | 2000 | 10 | 0.348 | 0.721 | 1.060 | 2.376 |
  | 2000 | 30 | 0.445 | 0.706 | 1.136 | 2.524 |

  **The ratio rises monotonically with band centre at every one of the six points
  and crosses unity between the bend and the window.** None of the three outcomes
  the brief anticipated fits: it is not near-uniform (so not an overall dipole-
  surface intensity scale), and it is not an excess confined to the stretch on an
  otherwise matching spectrum — there is a **deficit of 2-3× in the libration
  band** as well. The scale-free statement is the same: FP's own stretch/bend
  balance is 9.4-15.7 against the line list's 2.8-4.3, a factor **3.3-3.6**, which
  no overall intensity factor can produce. Attribution therefore points at whatever
  acts differently on high- and low-frequency modes — the classical dynamics plus
  the harmonic quantum correction — rather than at the dipole surface. The window
  column (1.06-2.51) is the collision-induced contribution and is expected to
  exceed unity whatever the answer, since a line list has almost nothing there.

- **A62 — The 48 points are one grid, in every output and figure.** The manuscript
  makes no distinction between how or when a point was computed, so neither do the
  outputs. `crossover_v3_rows` no longer emits `is_extension_point` or
  `in_delivered_set`, and `GROUPS` is now defined by *where a point sits*:
  `dilute` (ρ ≤ 60), `cold_dense` (T ≤ 1000 K and ρ ≤ 500), `dense` (the rest).
  Counts move from 27/12/9 to **26/13/9** — (1000 K, 129.744) is cold-dense under a
  physical definition and was only ever grouped elsewhere because it happened to be
  in the delivered set. The delivered-JSON cross-check survives untouched: it is a
  check on *our pipeline*, not a label on the grid.

  `figures/_manuscript.py` holds the journal-fixed style (183 mm width, ≥ 7 pt
  sans at final size, Okabe-Ito for categorical marks, vector PDF with live text)
  and a ledger that writes `figs/fig_numbers.json` — every number a reader could
  take off a panel, so a caption can be checked against the figure rather than
  against memory. `save_pdf` refuses to write a figure whose text contains
  "original", "extension" or "campaign", so the vocabulary rule is enforced at the
  point of writing rather than in review. The two continuous colour scales use
  viridis: Okabe-Ito is categorical and cannot carry a continuous variable.

## 2026-09-23 (line list + empirical continuum; broadening convention)

- **A58 — MT_CKD 4.3 water SELF-continuum added as a fourth opacity source,
  `hitemp_mtckd`, implemented from AER's own source and gated against AER's own
  binary.** Fig. 2b's ratios are FP against *bare* HITEMP — line absorption, no
  continuum. Every planetary radiative-transfer model adds an empirical continuum
  on top of its line list, so both references are now computed.

  **Provenance, and a deviation from the brief.** The brief said to read
  `~/scw/mtckd_continuum.py` and `mtckd_ref.npz` from ARC and reimplement what
  they do. **Neither file could be copied: `scp` to `htc-login.arc.ox.ac.uk`
  requires an interactive password and this session cannot supply one.** Rather
  than block, the model was implemented from the *upstream* source, which is
  better provenance than a reimplementation of a reimplementation: AER's MT_CKD
  release 4.3 at `~/software/MT_CKD` — `src/mt_ckd_h2o_module.f90` for the
  formulation and `data/absco-ref_wv-mt-ckd.nc` (sha256 `69944eb8…87ddb3`) for
  the coefficients. Every question the brief wanted the ARC read to settle is
  settled here from the primary source. **What remains unverified is only whether
  ARC's script agrees with AER's**; `mtckd_v4.3_self_ref.npz` is deliberately
  *not* named `mtckd_ref.npz`, so Harry's later `scp` cannot overwrite it and a
  one-line diff is possible.

  **Which version, which terms.** MT_CKD **4.3**; **self only, with the radiation
  term**; no foreign term. Self-only is not a modelling choice — AER's foreign
  coefficient carries a factor (1 − x_H2O) and the envelope is pure steam, so the
  foreign continuum is identically zero by composition. Implemented exactly as
  the Fortran computes it (Mlawer et al. 2023 Eqs. 1-4):

      C(ν,T)   = C_ref(ν) · (296/T)^n(ν)                       Eq. 2
      α_cont   = n_H2O · C(ν,T) · (n_H2O/n_ref) · R(ν,T)       Eqs. 1, 3
      κ_cont   = α_cont / ρ                                    [cm²/g]
      R(ν,T)   = ν·tanh(hcν/2k_BT)                             Eq. 4 (LBLRTM RADFN)

  with n_ref the number density of AER's reference state (1013 mb, 296 K;
  ρ_ref = 0.7418 kg/m³). AER's Fortran scales by `(p/p_ref)·(T_ref/T)`, which *is*
  n/n_ref for an ideal gas; steam at envelope densities is not ideal, so the
  number-density ratio is formed directly and pressure is never used — which is
  also what the brief specified. **α_cont is quadratic in density, so the mass
  opacity κ_cont is linear in ρ**, unlike line absorption, whose mass opacity is
  ρ-independent at fixed broadening. That is the whole reason it matters here.
  AER's XINT 4-point interpolation off the 10 cm⁻¹ reference grid is reproduced
  rather than replaced by a linear one; the choice is worth up to 11% pointwise
  and <1% on the band-limited means (gated both ways).

  **The gate is against AER's compiled binary, not against ourselves.** At the
  brief's stated gate state (296 K, ρ ≈ 0.7418 kg/m³ = 1 bar-equivalent) the
  density scaling is unity by construction, so that comparison alone tests
  nothing. Two gates against AER's `cntnm` v4.3 do the real work, both over the
  whole 100-4000 cm⁻¹ band, against a 1% requirement:

  | gate | what it exercises | agreement |
  |---|---|---|
  | `WATER.COEF` at **296, 650, 1000, 2000 K** | the reference table, (296/T)^n(ν), the radiation term | **4.8e-5** |
  | `CNTNM.OPTDPT` at **(1000 K, 129.744)** and **(650 K, 300 kg/m³)** | all of the above *plus the ρ² density scaling*, at 175× and 404× ρ_ref | **8e-4** |

  The second is the one that matters, because the density extrapolation is the
  part of this model that is genuinely being pushed. Its 8e-4 residual is not
  noise: it is positive at ~97% of wavenumbers with a median of 3.1e-4, which is
  exactly cntnm's own column constant (`WTOT = ALOSMT·(P/1013)·(273/T)·L`, with
  `ALOSMT = 2.68675e19` and 273 rather than 273.15) against CODATA. That
  constant, and the fact that the *density ratio* cancels identically between the
  two implementations, are pinned by their own gates rather than absorbed into a
  tolerance. Six reference outputs are committed under `validation/data/mtckd/`
  with the exact recipes (`cntnm` reads its state from **stdin** and *ignores* the
  `INPUT` file its own example ships — a trap worth recording).

  **The continuum is extrapolated, hard, in both variables, and this is stated
  wherever the numbers are.** MT_CKD is fitted to terrestrial-atmosphere
  measurements: T ≈ 250-350 K, densities at or below the reference. It is used at
  **T = 500-3000 K** (a power law in T fitted over ~100 K carried over ~2700 K,
  with exponents n(ν) = 1.43-10.3) and **ρ = 0.1-1000 kg/m³, up to ~1350× the
  reference density**. Nothing about that is endorsed by AER. It is done because
  it is what a planetary model does, and the point is to show what that costs.

  **Implementation.** `opacity/mtckd.py` (the model), `tools/extract_mtckd_ref.py`
  (one-off netCDF → npz, so netCDF4 stays out of the runtime dependencies, SPEC
  §1), a new `SourceSpec.continuum` field, and one registry entry pointing at the
  *same* `hitemp_opacity/` spectra — adding the source needed no change to the
  table builder, the interpolator or the structure driver (SPEC §3). The continuum
  is added to κ(ν) **at load time, before any mean is taken**, so the harmonic
  Rosseland mean sees one spectrum rather than a sum of two means. The SPEC §2
  plausibility ceiling (κ_P ≤ 1e5 cm²/g) guards the delivered-data ×1000 error, not
  a physically computed additive term; a continuum source that trips it is
  re-asserted on its line-only twin instead. It does not currently trip
  (max κ_P = 1.74e4 cm²/g, at 500 K / 1000 kg/m³).

  **THE RESULT (Task 20; `outputs/crossover_v3.json`, `figures/plot_p2_crossover.py`).**
  The 48 points split into the three groups the campaigns delivered — the 27
  original, the 12 cold-dense (T ≤ 1000 K, ρ = 100-500), the 9 dilute (ρ ≤ 60):

  | group | n | FP/HITEMP Planck | FP/HITEMP Rosseland | FP/(H+MT) Planck | FP/(H+MT) Rosseland | MT_CKD share of κ_R |
  |---|---|---|---|---|---|---|
  | original   | 27 | 1.38-17.7 | 1.36-41.9 | 0.92-15.7 | 1.34-38.4 | 0.006-0.083 |
  | cold-dense | 12 | 1.29-2.27 | **7.80-213** | **0.21-1.10** | **7.33-21.8** | **0.059-0.944** |
  | dilute     |  9 | 1.32-1.53 | 1.08-7.99 | 1.21-1.44 | 1.08-7.65 | 0.001-0.043 |

  Three things follow, and all three belong in the manuscript.

  **(1) The continuum changes nothing away from the cold-dense corner.** At the 27
  original points it carries 0.6-8.3% of κ_R, and at the 9 dilute points
  0.1-4.3% — it is quadratic in density and those cells are thin. The published
  ratios are unaffected there.

  **(2) At the 12 cold-dense points the continuum IS the line-list opacity.** Its
  share of κ_R runs to **94.4%** (650 K, 500 kg/m³): the line list's own absorption
  is then a *minority* of what a planetary model would use. The FP/line-list
  Rosseland ratio collapses from **7.8-213× to 7.3-21.8×** — A53's headline
  "24-213× at the cold points" survives only against the bare line list. The sign
  of the conclusion does not change (FP is still the more opaque), but the size of
  the disagreement does, by up to an order of magnitude, and the reader must be
  told which reference a quoted number is against.

  **(3) Eleven of the 48 points have FP/(HITEMP+MT_CKD) < 1, all on the Planck
  mean** — the extrapolated empirical continuum is *more* absorbing than the
  first-principles calculation there. They are the 8 cold-dense points at
  T = 650/800 K (ratios 0.21-0.83), the 2 at (1000 K, 300/500) and (1000 K, 129.7),
  the last being an *original* grid point. No point falls below 1 on the Rosseland
  mean. This has to be stated: it is the one place where "FP is more opaque" fails,
  and it fails against a continuum model used ~674× above its reference density
  and ~350 K above its reference temperature. Whether that inversion is physical or
  an artefact of extrapolating MT_CKD cannot be settled here — the honest reading
  is that at T ≤ 1000 K the *empirical continuum's extrapolation uncertainty is
  comparable to the FP/line-list difference the manuscript is about*.

  **THE CAVEAT THAT HAS TO TRAVEL WITH ALL OF IT — the decomposition is not well
  posed on most of this grid.** MT_CKD is not "extra absorption" in the abstract.
  It is *defined* as the residual beyond a line-by-line calculation done to
  LBLRTM's convention: every line truncated at **25 cm⁻¹** from centre with the
  Lorentz value at that bound (the "plinth") subtracted. That split only means
  anything while the Lorentz half-width is small compared with 25 cm⁻¹. Taking
  HITEMP's own broadening parameters — intensity-weighted over 2.8×10⁶ H2O lines
  across 150-4150 cm⁻¹, γ_self = 0.455 cm⁻¹/atm at 296 K with exponent 0.529 —
  and evaluating γ = γ_self·(P/1 atm)·(296/T)^n at each cell's ideal-gas pressure:

  | group | n | γ [cm⁻¹] | inside the 25 cm⁻¹ bound |
  |---|---|---|---|
  | original   | 27 | 141 - 4070 | **0 / 27** |
  | cold-dense | 12 | 88.9 - 544 | **0 / 12** |
  | dilute     |  9 | 10.9 - 90.5 | 3 / 9 |

  **45 of the 48 points have lines broader than the truncation the continuum is
  defined against**, by one to two orders of magnitude. Two consequences:

  * The delivered HITEMP spectra are RADIS at `truncation = 50` cm⁻¹, **with the
    profile renormalised to unit area** (A59) and **no plinth subtraction**. So
    the two conventions disagree about the 25-50 cm⁻¹ region twice over: MT_CKD
    assumes it has been cut and its plinth removed, while RADIS retains it and
    additionally folds the whole >50 cm⁻¹ wing back into it. The mismatch is real
    and it is not small, but **its sign is not determined** — renormalisation
    conserves strength, so this is a redistribution, not an excess, and whether
    it raises or lowers κ_R depends on where the area lands relative to the
    windows. An earlier draft of this amendment asserted it biased the combined
    opacity high (and so the FP ratios low); that claim is withdrawn.
  * More fundamentally, at these pressures the band is a single merged profile and
    "which part is line and which is continuum" has no operational meaning. The
    right reading of the continuum-inclusive column is **"what a planetary model
    would compute if it followed standard practice here"**, not "the physically
    correct opacity". Every output carries `lorentz_hwhm_cm` and
    `decomposition_well_posed` per point so this cannot be quoted without it.

  This is a caveat on the *reference*, not on FP, and it does not touch the bare
  ratios. But it is the first thing a careful reader will ask, and it is gated
  (`test_the_line_continuum_decomposition_is_not_well_posed_on_most_of_the_grid`).

  Interaction with the open A52 question: none numerically — the continuum is added
  to the HITEMP side, and A52 concerns the FP side's noise floor — but the two
  compound. The cold Rosseland ratio is simultaneously (a) possibly measuring FP's
  unresolved windows (A52/A53) and (b) sensitive at the order-of-magnitude level to
  whether the line-list reference carries an extrapolated continuum. Neither is
  resolved in this repository.

  The nine ρ > 1000 points keep the existing nearest-neighbour clamp (A53) rather
  than a hybrid κ_line(clamped) + κ_cont(true ρ), which would be a third
  convention and is not additive in the harmonic mean. The continuum's
  contribution is therefore *under*-counted there by **ρ/1000** — linear, not
  quadratic: the absorption coefficient α_cont goes as ρ², but the tabulated
  quantity is the mass opacity κ = α/ρ, which goes as ρ. That is up to **2.23×**
  at (3000 K, 2229), and the factor is reported per point as
  `continuum_underscale`.

- **A59 — Broadening-convention sensitivity: the brief's premise is wrong, the
  delivered HITEMP spectra are NOT at a fixed 1 bar reference, and the real
  question is ideal-gas P versus AQUA P.** Task 21 was written to quantify "the
  effect of broadening at the cell's actual pressure" against "the 1 bar-reference
  values" the bare HITEMP spectra were said to use. They do not use one.

  **The evidence, three independent lines.** (a) Every delivered file carries
  `pressure_mode = 'ideal'` and a `p_bar` equal to ρRT/M at its own density —
  0.4615 bar at (1000 K, 0.1 kg/m³), 4615 bar at (1000 K, 1000). (b) κ(ν) at fixed
  T differs by up to **143×** between ρ = 0.1 and ρ = 1000, which is impossible
  under a ρ-independent broadening: κ = α/ρ would then be identical at every
  density, since α ∝ n. (c) Re-running RADIS locally at the ideal-gas pressure
  reproduces the delivered `T1000_rho30.npz` to **1e-6 on both means** and 3.5e-5
  monochromatically across 100-4000 cm⁻¹. **A17 already recorded this for the
  delivered POKAZATEL set** — the `HitempCalculator` *snippet* used a 1 bar
  reference, the production run did not — and the HITEMP delivery followed the
  same convention.

  So the brief's comparison is not available, and the meaningful one is **ideal-gas
  P (what the tables were built at) versus AQUA P(T, ρ) (what the structure model's
  own EOS says)**. All three conventions are computed, so the brief's literal
  question is answered as well:

  | convention | P at the six points [bar] | what it is |
  |---|---|---|
  | `ref1bar` | 1 | the convention the brief assumed |
  | `ideal` | 138 - 30 900 | P = ρRT/M, the delivered tables' convention |
  | `aqua` | 136 - 486 000 | the structure model's EOS |

  AQUA/ideal spans **0.253×** at (650 K, 300) — inside the vapour-liquid region,
  where the EOS pins P far below ideal — to **15.7×** at (3000 K, 2229). Every
  point sits in AQUA phase 5 except (1000 K, 30), which is phase 3.

  **Mass opacity is convention-independent by construction.** RADIS returns α at
  its own ideal number density n = P_broad/(k_BT), and the pipeline forms
  κ = α/(n·M/N_A) ≡ σ·N_A/M. The broadening pressure therefore moves *only* the
  Voigt widths. That is what makes this a clean one-variable experiment.

  **The Lorentzian half-width cap the brief asks to be stated: RADIS's
  `truncation`, 50 cm⁻¹** (its default, and what the delivered spectra used). The
  nominal self-broadened HWHM at these pressures runs far past it — water is
  ~0.46 cm⁻¹/atm at 296 K, so at ρ = 2229 kg/m³ it is O(10³) cm⁻¹.

  **CORRECTION (same day, before any of this was quoted anywhere): the cap does
  not discard absorption, it redistributes it.** An earlier draft of this
  amendment said the high-pressure spectra under-count and that the ratios are
  lower bounds. That is wrong. RADIS renormalises each truncated profile to unit
  area (`radis/lbl/broadening.py::_calc_lorentzian_lineshape`, `lineshape /=
  area`), so line strength is conserved exactly. Measured on the delivered
  spectra, which settles it independently of reading the code: at T = 1000 K,
  ∫κ(ν)dν over the delivered range is **constant to 1 part in 10⁴** from
  ρ = 0.1 to 1000 kg/m³ — a 10⁴× range in broadening pressure (0.46 → 4615 bar).

  What the truncation actually does is pile the area that physically belongs
  beyond ±50 cm⁻¹ back into the window. Where HWHM ≫ 50 cm⁻¹ a line becomes, in
  effect, a ~100 cm⁻¹-wide rectangle carrying its full strength — not a
  Lorentzian. **The sign of that error on κ_R is not determined a priori**: it
  depends on where the redistributed area lands relative to the transparency
  windows, and κ_R is a harmonic mean dominated by exactly those windows. So the
  ratios here are **not bounds in either direction**. They carry an unquantified
  line-shape systematic wherever HWHM > 50 cm⁻¹ — which, by A58's half-width
  table, is **45 of the 48 FP grid points**, and includes every point in the
  delivered HITEMP tables at ρ ≳ 10. Gated by
  `test_line_strength_is_conserved_across_broadening_pressure`.

  Any point RADIS refuses outright has its refusal text recorded in the output,
  since that text is itself the answer to "what half-width does RADIS accept".

  **Run environment (local, not ARC).** RADIS 0.16.3 and a 7.6 GB HITEMP cache are
  both on the Mac, so the brief's ARC fallback was not needed;
  `tools/broadening_sensitivity.sbatch` is written anyway and is the route if the
  run has to move. Two environment facts are worth recording because they will bite
  again:

  * **RADIS ≥ 0.16 made `SpectrumFactory` slotted.** `sf.pressure = ...` — what
    `~/scw/hitemp_opacity.py` and `generate_linelist_tables.py` do between points —
    now raises `AttributeError`. The pressure is a per-call `eq_spectrum(pressure=)`
    argument. The new form works on the older API too, so no version branch.
  * **`fetch_hitemp` aborts on a fully-cached offline load.** It calls
    `HITEMPDatabaseManager.register()` whenever the databank is absent from
    `~/radis.json`, and that method formats `self.wmin`/`self.wmax`, which are only
    set when files are actually downloaded (`TypeError: unsupported format string
    passed to NoneType`). Registration is pure bookkeeping, so the driver skips it
    when nothing was downloaded rather than writing an entry into the user's home
    config.
  * The local cache lacks HITEMP's 0-50 cm⁻¹ file, so lines are loaded from
    **50 cm⁻¹** rather than the delivered run's 10. With truncation = 50 a line
    below 50 cm⁻¹ cannot reach past 100 cm⁻¹, so the science band [100, 4000] is
    unaffected — which the 1e-6 reproduction above confirms empirically.
    `--wavenum-min 10` restores the delivered convention exactly where the file
    exists (the sbatch wrapper passes it).

  `tools/broadening_sensitivity.py` (driver, resumable, per-point failure isolation,
  JSON rewritten after every point), `figures/plot_ed1_broadening.py` (the Extended
  Data panel: κ_R ratio against the 1 bar reference, and the residual AQUA/ideal
  systematic actually present in the delivered tables).

  **RESULTS — COMPLETE. All 18 spectra ran, none refused.** ~2.5 h wall on the Mac.
  `outputs/broadening_sensitivity.json` + `outputs/ED1_broadening.png`.

  | T [K] | ρ [kg/m³] | P ideal / AQUA [bar] | κ_R 1 bar | κ_R ideal | κ_R AQUA | ideal/1 bar | AQUA/ideal | κ_P 1 bar | κ_P ideal |
  |---|---|---|---|---|---|---|---|---|---|
  | 1000 | 30 | 138.46 / 135.88 | 2.051 | 4.432 | 4.427 | **2.161** | 0.999 | 277.4 | 277.5 |
  | 1000 | 129.7 | 598.60 / 553.96 | 2.051 | 4.660 | 4.652 | **2.272** | 0.998 | 277.4 | 277.8 |
  | 650 | 300 | 899.97 / 227.67 | 0.7715 | 0.6561 | 0.6313 | **0.850** | **0.962** | 357.7 | 358.8 |
  | 1500 | 471.3 | 3262.7 / 4275.7 | 12.54 | 22.07 | 22.32 | **1.760** | **1.011** | 251.8 | 253.5 |
  | 2000 | 1200 | 11077 / 51080 | 34.23 | 50.64 | **112.3** | **1.479** | **2.218** | 245.8 | 247.2 |
  | 3000 | 2228.7 | 30858 / 485710 | 86.35 | 103.7 | **98.06** | **1.201** | **0.945** | 237.3 | 229.4 |

  Over the six points: κ_R(ideal)/κ_R(1 bar) = **0.850-2.272**, κ_R(AQUA)/κ_R(1 bar)
  = **0.818-3.281**, κ_R(AQUA)/κ_R(ideal) = **0.945-2.218**.

  **The ρ-clamp error is now MEASURED, and A53's "up to 30%" was mis-attributed.**
  Two of these points sit above HITEMP's ρ = 1000 ceiling, so the table holds κ at
  the ceiling node while the RADIS run computes the true density directly. That is
  the first direct measurement of the clamp:

  | T [K] | ρ | ρ/ρ_ceiling | κ_R clamped | κ_R true | clamp error | deviation from the delivered analysis |
  |---|---|---|---|---|---|---|
  | 2000 | 1200 | 1.20 | 49.82 | 50.64 | **−1.62%** | +1.88% |
  | 3000 | 2228.7 | 2.23 | 97.59 | 103.7 | **−5.89%** | +5.89% |

  The measured clamp error reproduces the deviation from the owner's independent
  RADIS crossover to **better than 0.3 percentage points at both points**, so A53's
  mechanism is confirmed — where it applies. Across the nine clamped points
  (ρ/ρ_ceiling = 1.04-2.23) the clamp is therefore a **≤6% effect on κ_R**, which
  is negligible beside the FP/HITEMP ratios of 7-42× those points carry.

  **What it does not explain.** A53 reported "ρ clamped → up to 30% on the
  Rosseland mean" and attributed all of it to the clamp. The deviations are not
  monotone in ρ/ρ_ceiling at all: 1.04 → +0.07%, 1.20 → +3.19% and +1.88%,
  **1.26 → +29.65%**, 1.34 → +0.19%, 1.52 → +2.56%, 1.73 → +2.00%,
  **1.84 → +19.91%**, 2.23 → +5.89%. The two outliers are at T = 2250 and 2750 K —
  **the only two temperatures off HITEMP's own T grid**, so those points are
  simultaneously ρ-clamped *and* T-interpolated. The clamp is worth a few percent
  there, like everywhere else; the remaining ~25% is something about the
  T = 2250/2750 K interpolation, or about the delivered analysis at those states,
  and this repository cannot separate the two. **A53's sentence should read "the
  clamp costs ≤6%, measured; the two 20-30% outliers are the T-interpolated points
  and are unexplained".**

  **A resolution caveat, and it falls on the 1 bar column — not the high-pressure
  one.** RADIS raised exactly one accuracy warning in the whole run, during the
  **T = 2000 K, 1 bar** spectrum: *"Some lines are too narrow (FWHM ~ 0.03 cm⁻¹)
  for the current spectral grid (wstep = 0.01)"*. At 1 bar the narrowest lines are
  resolved by only ~3 grid points. This is the reverse of the worry the brief
  raises: nothing complained at 11 077 bar, where the lines are enormous. It
  matters because κ_R is a harmonic mean and therefore sensitive to exactly the
  narrow, deep structure that is being under-sampled.

  Scope of the caveat: it touches **only the hypothetical 1 bar reference column**,
  which answers the brief's literal question. It does not touch the `ideal` or
  `aqua` columns (broad, well-resolved lines), and it does not touch the delivered
  tables at all — every one of those is computed at its own cell pressure, so its
  lines are broad. The wstep of 0.01 cm⁻¹ is inherited from
  `generate_linelist_tables.py` deliberately, so the runs stay comparable with what
  is on disk; re-running the 1 bar column at `wstep='auto'` is the follow-up if the
  Extended Data panel is to quote those ratios precisely.

  **(650 K, 300) completes, and the pressure response has an interior minimum:**
  κ_R = 0.7715 (1 bar) → **0.6313** (227.7 bar) → 0.6561 (900 bar). It falls, then
  rises. So the dependence on broadening pressure is not merely non-monotonic in
  sign between points — it is non-monotonic *within* a single point, with a
  turning point somewhere near a few hundred bar. No monotone reasoning about
  broadening pressure survives that, which is worth stating plainly because the
  natural intuition ("more pressure → broader lines → filled windows → higher
  κ_R") is what the brief, and this amendment's first draft, both assumed.

  **Ideal-vs-AQUA splits sharply in two, and the earlier "within ±4%" reading was
  wrong.** At the four moderate points (pressure ratios 0.25×-1.31×) κ_R moves by
  0.962, 0.998, 0.999, 1.011 — negligible. At (2000 K, 1200), pressure ratio
  4.61×, it moves by **×2.218**. An earlier draft of this amendment concluded from
  the first four that "rebuilding the delivered tables on AQUA pressures would not
  move them meaningfully", reasoning from the saturation measured at T = 1000 K
  (d ln κ_R/d ln P = 0.034 between 138 and 599 bar). That reasoning does not
  extrapolate: over the 11 077 → 51 080 bar step the same quantity is **0.521**,
  fifteen times steeper. The conclusion is withdrawn.

  **But the 2.2× is a convention artefact, not a physical result, and it must not
  be quoted as one.** At these pressures the half-width dwarfs the truncation:

  | state | HWHM [cm⁻¹] | × the 50 cm⁻¹ truncation | retained Lorentz area | profile scaled up by | κ_P/κ_R |
  |---|---|---|---|---|---|
  | 1 bar | 0.163 | 0.003 | 0.998 | ×1.002 | 7.18 |
  | ideal, 11 077 bar | 1810 | 36 | 0.0176 | **×56.9** | 4.88 |
  | AQUA, 51 080 bar | 8349 | 167 | 0.0038 | **×262** | **2.11** |

  RADIS is renormalising each profile by up to a factor 262 to put the whole line
  strength inside a ±50 cm⁻¹ box. The spectrum tends to a smooth pseudo-continuum,
  the windows disappear, and κ_R climbs toward κ_P — the contrast falls from 7.2 to
  2.1. That the true κ_R rises steeply here is physically reasonable (a genuine
  γ = 8349 cm⁻¹ really would smear the band into a continuum), but **the magnitude
  is set by an arbitrary 50 cm⁻¹ box, not by the line shape.**

  The defensible statement is therefore: *the broadening convention is worth under
  4% at ρ ≲ 500 and becomes the dominant uncertainty at ρ ≳ 1000, where the
  computed line-list opacity is convention-dominated rather than resolved.* That
  bears directly on the nine `ll_clamped` points of A58's table, which live in
  exactly that regime and carry the largest FP/HITEMP ratios in the manuscript.

  The number the manuscript actually needs from the (650 K, 300) point is the
  **ideal-vs-AQUA systematic: −3.8% on κ_R** (0.6313 / 0.6561) and **+0.1% on κ_P**, despite the
  AQUA pressure being only **0.253×** the ideal-gas one. That is the largest
  ideal-vs-AQUA effect so far (against 0.1-0.2% at the two T = 1000 K points) and
  it is still small. The saturation seen at T = 1000 K — d ln κ_R / d ln P falling
  from 0.156 over the first two decades to 0.034 over the next — is the reason: a
  4× pressure error buys only a few percent in κ_R once past ~100 bar.

  **κ_R is NOT monotonic in broadening pressure, and the sign flips with
  temperature.** Both rows below are direct RADIS runs, like for like:

  | point | 1 bar → ideal P | HWHM at ideal P | HITEMP contrast κ_P/κ_R |
  |---|---|---|---|
  | (1000 K, 30) | 2.051 → 4.432, **×2.16 up** | 32.7 cm⁻¹ | 60 |
  | (1000 K, 129.7) | 2.051 → 4.660, **×2.27 up** | 141 cm⁻¹ | 60 |
  | (650 K, 300) | 0.7715 → **0.6561**, ×0.85 **down** | 267 cm⁻¹ | 464 → 547 |

  At 1000 K, broadening fills the windows and κ_R rises, as expected. At 650 K it
  *falls* by 15%, and the spectral contrast *rises* from 464 to 547 — the
  effective windows get **deeper** as the gas is compressed, which no real line
  shape does. A simple "half-width past the truncation" threshold does not explain
  it on its own: (1000 K, 129.7) also has HWHM ≫ 50 cm⁻¹ and still rises. What
  separates the cases is window depth — at 650 K HITEMP's contrast is ~8× larger,
  so κ_R, a harmonic mean, is far more sensitive to how far-wing absorption is
  treated, and that is exactly what the truncation-plus-renormalisation mangles.

  **No directional claim is made.** The measurement is that the response reverses;
  the mechanism is entangled with the line-shape convention and cannot be
  separated with the runs in hand. This is the concrete reason the ratios are not
  bounds. It also means the *delivered* tables' cold-dense κ_R carries this
  systematic, since they were built with the same truncation at these pressures —
  which lands on the same 12 cold-dense points A58's continuum result lands on.

  **Second, much stronger test of A6's interpolation.** (650 K, 300) is
  interpolated in T across HITEMP's 600-700 K gap, over which κ_R **doubles**
  (0.4599 → 0.9235). Against the direct run:

  | | κ_R | κ_P |
  |---|---|---|
  | interpolated (A6) | 0.65172 | 359.97 |
  | RADIS direct | 0.6561 | 358.8 |
  | difference | **−0.67%** | **+0.33%** |

  Log-linear interpolation in (T, log ρ) holds to under 1% even where the function
  doubles across the interval. Together with the (1000 K, 129.7) ρ-interpolation
  check at 0.044%, A6's scheme is now validated against independently computed
  states in both coordinates.

  A third check at (1500 K, 471.3) gives κ_R **+1.09%**, κ_P −0.03%. It was
  expected to be the hardest ρ test, since 471.3 sits in HITEMP's widest ρ interval
  (300 → 1000); it is not, because κ_R is nearly flat there — only ×1.06 across the
  whole interval — so the 1.09% miss is ~17% of the interval's own spread. Worth
  recording as the least favourable of the three, but the T test above remains the
  stringent one. The genuinely demanding case is still to come: (3000 K, 2229) is
  **above** HITEMP's ρ ceiling and therefore nearest-neighbour clamped.

  **The pressure response tracks window depth, cleanly.** Contrast is quoted at
  1 bar so the four points share a reference:

  | point | contrast κ_P/κ_R at 1 bar | 1 bar → ideal P |
  |---|---|---|
  | (650 K, 300) | **463.6** | ×0.850 **falls** |
  | (1000 K, 30) | 135.3 | ×2.161 rises |
  | (1000 K, 129.7) | 135.3 | ×2.272 rises |
  | (1500 K, 471.3) | 20.1 | ×1.760 rises |

  The reversal happens only at the deep-window end, which is the cold-dense corner
  — the same corner as A58's continuum result and A52's noise-floor question. Three
  independent problems now converge on the same twelve points.

  This one complete point already fixes the shape of the answer, and it is the
  shape the physics predicts: **the broadening convention is a Rosseland-mean
  effect and essentially not a Planck-mean one.** κ_R more than doubles (2.16×)
  between 1 bar and the delivered ideal-gas P; κ_P moves by **0.04%**. Pressure
  broadening fills the transparency windows, the harmonic Rosseland mean lives in
  those windows, the arithmetic Planck mean does not — the same mechanism A53
  identified behind the cold FP/HITEMP ratio. Two readings follow:

  * a 1 bar-reference table, had one been used, would have understated κ_R by ~2×
    already at ρ = 30, and by more at higher density — so the brief's worry was
    well founded even though its premise about the delivered files was not;
  * the systematic *actually* present in the tables is ideal-vs-AQUA, and at this
    point it is **0.1%**, because P_AQUA/P_ideal = 0.98 here. It is the points
    where that pressure ratio departs from 1 — 0.25× at (650 K, 300) and 15.7× at
    (3000 K, 2229) — that will decide whether the delivered tables need rebuilding
    on AQUA pressures, and those are still running.

  **Unplanned by-product: the first direct test of A6's interpolation against
  ground truth.** ρ = 129.7 is not on HITEMP's ρ grid (0.1, 1, 10, 30, 100, 300,
  1000), so the table value there is triangulated in (T, log ρ) — the exact call
  the structure model makes. The RADIS run computes that state *directly*:

  | | κ_R [cm²/g] | κ_P [cm²/g] |
  |---|---|---|
  | interpolated (A6) | 4.6621 | 277.79 |
  | RADIS, computed at (1000 K, 598.6 bar) | 4.660 | 277.8 |
  | difference | **+0.044%** | **−0.005%** |

  Every previous check of the interpolation was either self-consistency or against
  the delivered JSON at points that were *on* the grid; this is the first against
  an independently computed off-grid state. Read it for what it is, though: κ_R
  only moves 2.6% between the bracketing nodes (4.634 at ρ = 100, 4.754 at 300),
  so the scheme is not under strain here. The strong tests are still to come —
  (1500 K, 471.3), which sits in the wide 300-1000 interval, and (3000 K, 2229),
  which is **above** HITEMP's ρ ceiling and therefore nearest-neighbour clamped,
  so it will measure directly the clamp error A53 could only bound at ~30% from
  the delivered analysis.

  Timing note for whoever repeats this: the third spectrum at this point took
  1570 s against ~450 s for the first two at essentially the same pressure, so the
  ~7.3 min/spectrum figure is not reliable and the full run may take considerably
  longer than the 2.2 h it implies. The driver is resumable per point, so this is
  an inconvenience rather than a risk.

  **The answer, in one paragraph.** The broadening-pressure convention is worth
  **under 4% on κ_R at ρ ≲ 500 kg/m³** (ideal-vs-AQUA = 0.962, 0.998, 0.999, 1.011
  at the four points there) and becomes the **dominant uncertainty at ρ ≳ 1000**
  (×2.218 at (2000 K, 1200), ×0.945 at (3000 K, 2229)). The hypothetical 1 bar
  reference the brief asked about would have been wrong by ×0.85-2.27 — but it is
  hypothetical, because the delivered tables do not use it. On the Planck mean the
  convention is negligible everywhere except the most extreme state, where κ_P
  falls 28% (229.4 → 165.2).

  **At ρ ≳ 1000 none of this is resolved physics.** At (3000 K, 2229) the AQUA
  half-width is **64 000 cm⁻¹ — 1281× the truncation, and wider than AER's whole
  0-20 000 cm⁻¹ reference range.** RADIS is packing each line's full strength into
  a ±50 cm⁻¹ box, the spectrum goes grey (κ_P/κ_R falls from 2.75 at 1 bar to 1.68),
  and both means are set by the box rather than by the line shape. The honest
  summary is: **the line-list opacity at ρ ≳ 1000 kg/m³ is convention-dominated**,
  and that is where the nine `ll_clamped` points of A58's table live.

  Regenerate the panel and table with
  `python figures/plot_ed1_broadening.py config/fiducial.yaml`.

## 2026-09-20 (extension opacity campaign: units at source, enlarged domain)

- **A50 — The FP factor-1000 error is fixed IN THE DATA; `sigma`/`sigma_raw`
  were already correct and were NOT scaled.** The extension campaign delivers 21
  new FP spectra in correct units, so a single on-load scale (SPEC §2, A4) could
  no longer serve both halves of the grid. `tools/fix_fp_units.py` therefore
  archives `fp_opacity/` to `fp_opacity_legacy_x1000/`, divides the 27 originals
  at source, merges the 21 new points, stamps every file with a `units_fixed`
  marker, and regenerates `grid.json` over all 48. The registry's FP `scale` is
  now 1.0 and `opacity/loaders.py` **refuses** any FP npz without the marker, so
  a legacy file can never be silently re-ingested.

  **Which arrays were scaled, and the evidence.** The owner's brief listed
  `kappa, alpha, sigma, sigma_raw, planck, rosseland`. Measured on the real
  files, the ×1000 entered at the `alpha` step and not at the cross-sections: in
  all 27 originals `alpha = 1000·n·sigma` (n = ρN_A/M, matched to 10 significant
  figures), and in all 21 extension files `alpha = n·sigma` exactly; `kappa =
  alpha/ρ` holds in both. So `sigma`/`sigma_raw` are already correct
  cm²/molecule. Scaling them would have broken `alpha = n·sigma` in the
  rewritten files and corrupted the cross-sections the dilute noise-floor check
  (A51) reads. **Owner decision (2026-09-20): scale `alpha`, `kappa`, `planck`,
  `rosseland` only; leave `sigma`, `sigma_raw` untouched.** The identity is
  re-asserted on every file before the rewrite and pinned by
  `validation/test_fp_units.py::test_sigma_was_not_scaled`.

  Gates (`validation/test_fp_units.py`, 7 tests): (a) the loader's κ at
  (1000 K, 129.744) is unchanged across the whole operation — κ_R = 35.409291,
  κ_P = 383.437496 cm²/g on [100, 4000] cm⁻¹, captured live immediately before
  the rewrite and now pinned as a constant, since no "before" survives it;
  (b) old-vs-new consistency κ_P(1000, 129.744)/κ_P(1000, 100) = 0.9262, inside
  [0.7, 1.3] and matching the stored full-band scalars' 0.9257; (c) all 48 files
  load with finite positive κ, all marked, all spanning the science band (common
  wn span [1.3359, 6000.11] cm⁻¹ — `T2750_rho1841.02` has the one offset grid).

  The stored `planck`/`rosseland` scalars were scaled for self-consistency only:
  they are FULL-BAND means and no science path reads them (band-limited means
  are recomputed from `kappa(wn)` by the shared table builder, A4).

  **Consequence — the FP grid is now 48 points over T = 650–3000 K,
  ρ = 10–2229 kg/m³**, replacing the 27-point T = 1000–3000 K, ρ = 109–2229 box
  that every "FP floor" number in the repo was written against.

  *SPEC §2 is now stale* (it mandates the ×1e-3 on load). Not edited — it is the
  owner's document; proposed replacement text is in the run summary.

- **A51 — The enlarged grid REVERSES A41: the old fiducial's deep closure now
  rests on interpolated FP data, not on a clamp.** A41's finding was that at
  (M_c = 5 M⊕, f_env = 5%, T_eq = 500 K, T_int = 20 K) every one of the 46 FP
  evaluations closing ρ ∈ [400, 550] was a nearest-neighbour clamp ~300 K below
  the grid's 1000 K floor, which is why the headline was demoted to
  "extrapolated, not evidence". The extension campaign added T = 650 and 800 K
  at ρ = 100-500 — precisely that region — and A48's cold-extension ask is
  therefore satisfied. Re-running the same audit (n_sample = 1500):

  | | A41 (27-point grid) | now (48-point grid) |
  |---|---|---|
  | FP samples in ρ ∈ [400, 550] | 46 | 48 |
  | of those, clamped | **46 / 46**, all below the T floor | **0 / 48**, all interpolated |
  | composite `clamp_fraction_all` | 0.4855 | 0.4562 |
  | envelope mass `clamped_t` | 0.045 (A48) | **0.000209** |
  | envelope mass `in_box_out_of_hull` | 0.954 (A48) | 0.9917 |
  | envelope mass `interpolated` | 0.0009 (A48) | 0.0081 |

  Composite: R_p = 1.71732 R⊕, one RCB at ρ = 215.1, max radiative ρ = 136.8.
  HITEMP is unchanged (its data did not move): R_p = 1.71229 R⊕, three RCBs at
  ρ = 285.3 / 515.3 / 522.4, max radiative ρ = 504.3, `clamp_fraction_all`
  0.3791. **The closure itself still holds** — HITEMP leaves a radiative sliver
  at ρ = 450-504, the composite does not — and it is now made of FP data inside
  the grid.

  *But the old fiducial is still not running on FP data overall.* 99.2% of its
  envelope mass is `in_box_out_of_hull`: inside the grid's (T, ρ) box, outside
  the scattered hull (A6). The campaign closed the "too cold for the grid" ask
  and left the "fill the grid's own interior gaps" ask open — the two A48
  separated. That gap needs no extension of the grid's range at all.

  **Gates updated (5), all consequences of the enlarged grid, none loosened:**
  `test_opacity::test_raw_unscaled_means_trip_the_assertion` (no scale left to
  undo since A50 — the legacy error is now re-introduced explicitly ×1000, and
  the guard is unchanged); `test_s5_clamp_audit::test_fiducial_deep_closure_*`
  (inverted, and renamed to `..._now_rests_on_interpolated_fp_data`);
  `test_s5_clamp_audit::test_clamp_reasons_attributes_to_the_used_subsource`
  (the (700 K, 500) probe is in-grid now; the out-of-grid probe moves to
  500 K); `test_s5_grid_coverage::test_the_old_fiducial_runs_on_extrapolated_
  opacity` (`clamped_t` > 0.01 becomes < 0.01 and > 0, `in_box_out_of_hull`
  > 0.9); `test_s5_grid_coverage::test_extension_nodes_are_reported_separately`
  (bracketing corners may now sit one lattice step above the 650 K floor,
  because the 600-700 K cell straddles it). Suite: 124 → 131 gates.

- **A52 — The dilute FP windows are floor-dominated — but so is nearly the whole
  grid, and the construction contradicts itself.** Task 15 asked whether the nine
  new ρ ≤ 60 cells clear the dipole model's noise floor
  (`figures/plot_s7_dilute_floor.py`, gates in
  `validation/test_s7_dilute_floor.py`). The ARC original,
  `lda_campaign/stepG_dipoles/compare_point_full.py`, is **not on this machine**
  and remains the authority; what is implemented here is:
  σ(ν) ∝ ν(1−e^{−hcν/kT})·D(ν) ⇒ D(ν) = σ/[ν(1−e^{−x})]; the residual is an
  *amplitude* ratio so it enters the power squared, r² = (5.75/12.4)² = 0.215;
  that power is spread flat over the 0.5 fs timestep's Nyquist wavenumber
  (33 356 cm⁻¹) giving D_floor = r²∫D dν/ν_Nyq; mapped back through the same
  prefactor, κ_floor(ν) = D_floor·ν(1−e^{−x}) — which rises across the band, not
  flat in κ. The bandwidth is the dominant assumption and all three choices
  (Nyquist / the computed [1.1, 6000] range / the science band alone) are
  reported.

  **The answer, plainly: yes.** At the nine dilute cells only **0.8-4.8%** of the
  1900-3100 cm⁻¹ transparency window clears 2× the floor (band cores 60-75%).
  And it matters because that window carries the Rosseland mean: the harmonic
  weighting puts **86-93%** of κ_R's integrand inside it at these cells (the
  Planck mean, arithmetic, takes <8% from it). So FP κ_R at dilute densities is
  set by a spectral region the dipole model does not resolve.

  **But the control reverses the framing.** Running all 48 points: the window is
  floor-dominated at *every one of the 34 cells with ρ ≤ 500* (window fraction
  0.000-0.048), including the dense points the existing headline rests on, and
  the nine dilute cells are if anything marginally *better* than the ρ = 100-500
  ones. Only the hot dense corner (ρ ≳ 1300, T ≥ 2500) has windows clearly above
  the floor, and it is also the only place the Rosseland window share drops
  (0.35-0.75 for ρ > 1000). The defensible statement is therefore **not** "the
  dilute points are unreliable" but "the transparency windows are the weak part
  of the FP spectra everywhere" — which does not single out any point for
  exclusion.

  **Open question (blocks no computation; changes the manuscript's wording).** In
  42 of 48 cells the delivered, *smoothed* κ sits BELOW this floor in the window
  (median κ/floor down to 0.06). For an additive white residual that is
  impossible — smoothing removes fluctuation, not bias, so E[κ] ≥ floor always.
  So the data rules out this construction at the quoted power fraction. Two
  limiting models bracket the truth: **(a)** white additive noise in the dipole
  time series — the construction above, which the data contradicts; **(b)** a
  dipole-*surface* fit error correlated with the configuration, whose spectral
  signature tracks the signal, i.e. a ~21.5% **multiplicative** uncertainty at
  every ν, which every cell clears trivially. Under (a) FP κ_R at ρ ≤ 500 is set
  by unresolved windows; under (b) it is not. `compare_point_full.py` knows
  which; this repository cannot. A third possibility is that the delivered
  `sigma` was floor-subtracted before smoothing. **Asked of the owner in the run
  summary; no number in Tasks 16-18 depends on the answer.**

- **A53 — Crossover rebuilt on 48 points; the new cold points expose that the
  FP/line-list Rosseland ratio is a WINDOW-DEPTH difference, not an opacity
  difference.** `opacity/crossover.py` recomputes the M2 ratios through our own
  band-limited table builder over all 48 FP points, writing
  `outputs/crossover_hitemp_v2.json`; `figures/plot_p2_crossover.py` now plots
  that instead of replaying the delivered JSON, which is kept overplotted as the
  independent check. The extension grid is not a subset of HITEMP's (T = 650 and
  ρ = 60/150/500 have no twin), so HITEMP is evaluated through its own
  `GridOpacitySource` — the same triangulated (T, log ρ) call the structure model
  makes — rather than by exact-point matching. Every FP point is interior to
  HITEMP's T range; the nine with ρ > 1000 keep the existing nearest-neighbour
  clamp and are flagged `ll_clamped` everywhere.

  **Agreement with the delivered RADIS analysis, split by how κ_LL was obtained**
  (the split matters: a single "max deviation" would be dominated by the clamps):

  | κ_LL obtained by | n | max \|dev\| Planck | max \|dev\| Rosseland |
  |---|---|---|---|
  | both coords on HITEMP's grid | 3 | 0.38% | **0.24%** |
  | T interpolated (2250, 2750 K) | 15 | 0.63% | 2.98% |
  | ρ clamped (> 1000 kg/m³) | 9 | 2.88% | 29.65% |

  The first row is the real cross-check and it reproduces A25's ~1%. The last is
  expected — RADIS computed those at the true density, we hold κ at ρ = 1000.

  **[Corrected 2026-09-23, A59.]** That last attribution is only partly right. The
  clamp error has since been *measured* against direct RADIS runs at the true
  density: −1.62% at ρ/ρ_ceiling = 1.20 and −5.89% at 2.23, reproducing the
  deviations at those points to <0.3 pp. So the clamp is a ≤6% effect, not a 30%
  one. The two 20-30% outliers are at T = 2250 and 2750 K, the only temperatures
  off HITEMP's T grid, and remain unexplained. See A59.

  **The finding.** At the eight new cold points (T = 650, 800 K) the FP/HITEMP
  Rosseland ratio is **24-213×**, against 1.4-7× everywhere at T ≥ 1000, while
  the Planck ratio stays at 1.3-2.0. The spectral contrast κ_P/κ_R says why:
  HITEMP's runs from 2.4 at 3000 K to **567 at 650 K** — a line list's
  transparency windows deepen enormously as it cools, lines being narrow and
  sparse — while **FP's never leaves 1.0-32 anywhere on the grid**. The FP
  spectrum has no deep windows at all. Since the Rosseland mean is harmonic and
  takes 74-93% of its integrand from the 1900-3100 cm⁻¹ window at ρ ≤ 500 (A52),
  the cold Rosseland ratio is measuring how deep each dataset's windows are, not
  how opaque the water is.

  **Consequence for the headline.** "FP is 2-3× more opaque at ρ ≥ 109, T ≥ 1000"
  survives as a statement about the *legacy* grid. Its extrapolation to the cold,
  dense regime — the regime the whole structural claim lives in, and the one the
  campaign was run to reach — does not: there the disagreement is 24-213× and is
  attributable to FP's filled-in windows under A52's model (a), or is a real and
  extraordinary physical result under model (b). The two cannot be separated in
  this repository. Tasks 17 and 18 report the numbers the model gives either way,
  but the manuscript cannot quote a cold Rosseland ratio as physics until A52's
  open question is answered.

- **A54 — Blend window moved to 10-30 kg/m³ (NOT the proposed 3-10), and the FP
  box corner is a config value.** With FP data now reaching ρ = 10, the owner
  proposed blending 3-10. Measured κ_R seam continuity (FP/HITEMP at the seam,
  T = 1000/1500/2000 K) across the candidates:

  | window | seam ρ | worst seam mismatch | FP at the seam |
  |---|---|---|---|
  | 3-10 | 5.48 | 3.08× | **CLAMPED — below the FP grid's ρ floor** |
  | **10-30** | **17.32** | **3.69×** | interpolated |
  | 30-100 | 54.77 | 7.47× | interpolated |
  | 100-300 (previous default) | 173.21 | 10.50× | interpolated |

  3-10 wins narrowly on continuity but its seam sits *below* ρ = 10, so the
  composite would hand κ to FP exactly where FP has no data — nearest-neighbour
  extrapolation by construction. **Adopted 10-30**: the best continuity among
  windows where FP has data, and it cuts the seam discontinuity from 10.5× to
  3.7×. Pinned by `test_fp_grid_floor::test_the_blend_window_starts_at_or_above_
  the_fp_data_floor`.

  *Caveat, stated rather than buried:* a lower seam hands more of the cold
  envelope to FP, and A52/A53 show that is exactly where FP's κ_R is least
  trustworthy (floor-limited windows; 24-213× disagreement at T < 1000 K). Seam
  continuity alone cannot settle the window; `plot_s1_blend_sensitivity.py`
  remains the scan that bounds it.

  **FP grid floors are now `opacity.fp_grid_t_floor_k` (650) and
  `opacity.fp_grid_rho_floor_kg_m3` (10)**, not module constants.
  `figures/_common.fp_floors(cfg)` is the accessor, `structure.model` keeps
  mirrored constants for the few call sites with no Config, and
  `validation/test_fp_grid_floor.py` asserts the configured values equal the
  built FP table's actual minima — so config cannot drift from data.

  **Gates updated (6), none loosened.** `test_evolution::test_history_comparison_
  diagnostics` — the "FP-relevant epoch" has stopped discriminating: with the ρ
  floor at 10 every RCB in the model qualifies, so the epoch now starts at
  0.16 Gyr instead of >0.5 Gyr; the physics claim it encoded (the *deep* RCB
  regime is late) is re-tested explicitly against the legacy 108.796 threshold.
  `test_history_diagnostics` — the synthetic crossing track starts at ρ = 10, so
  it now passes an explicit threshold rather than the module default.
  `test_s6_candidate` — blend-window pin, and the closure restated (below).
  `test_s4_domain_scan`, `test_s5_grid_coverage` — box-corner constants.
  `test_s5_grid_coverage::test_the_old_fiducial_*` — inverted, see below.

  **Consequence — the old fiducial now runs on FP data.** Envelope mass
  interpolated: **0.0009 (A48) → 0.0081 (A51, after the merge) → 0.9995 (now)**,
  with `in_box_out_of_hull` 0.954 → 0.992 → 0.0004 and `clamped_t` 0.045 →
  0.0002 → 0.0001. The composite's trajectory under the 10-30 window runs
  through the well-covered part of the enlarged grid. A41's central objection to
  the headline case is fully resolved.

  **Consequence — the candidate's closure is now PARTIAL.** A44 could say "zero
  in-box radiative layers" because the box was (1000 K, 108.796), entirely inside
  the FP-dominated region. With the box floor at (650 K, 10) "in box" reaches
  into densities where the composite uses the line list *by design*. At the
  candidate (T_eq = 700 K, T_int = 20 K): HITEMP leaves an in-box layer of
  **3025 ppm** (T 855-1023 K, ρ 10.4-130); the composite leaves **377 ppm**
  (T 855-916 K, ρ 10.1-25.6) — an 8× shrinkage at the box's dilute edge, not an
  elimination. Nothing survives above the blend window's top. R_p 1.75904 →
  1.77269 R⊕.

- **A55 — Scan on the enlarged domain: the FP/line-list dispute AT the boundary
  is 6-54×, against a flip threshold of 1.01-1.21×. The boundaries move; the mass
  does not.** 108 cells rerun with the new composite,
  `outputs/S4_domain_scan_extended.csv`. **All 108 converged inside the AQUA
  table; no skips.**

  *In-box radiative structure.* A47 found ZERO cells with any in-box radiative
  mass against the old (1000 K, 109) box. Against the new (650 K, 10) box most
  cells have some, but it is tiny — the largest is **506 ppm** of envelope mass
  (M_c = 5, wmf = 0.05, T_eq = 500, T_int = 30; layer T 653-731 K, ρ 28.9-60.2).
  **One cell of 108** has the composite close every in-box layer with
  `clamp_fraction_all` < 1% — that same cell, R_p 1.73887 → 1.74763 R⊕ (+0.50%).
  **Zero cells** have a radiative layer lying *entirely* inside the box: every
  layer runs up to ρ ~ 1e-6 at the top of the atmosphere.

  *The Task-17 clamp gate, and why it needs two forms.* As specified — clamped
  fraction where T ≥ 650 and ρ ≥ 10 — it is **unsatisfiable by construction**:
  the corner has no upper bound, so the deep convective interior (past 10 000 K
  against the grid's 3000 K ceiling, past 2229 kg/m³ against its density
  ceiling) counts as in-box and clamps necessarily. All 108 cells violate it.
  `structure.diagnostics.in_box_clamp_fraction` therefore reports both the
  literal corner and the grid's actual **rectangle**, with clamps attributed to
  `too_hot` / `too_dense` / `in_box_out_of_hull`. On the rectangle:

  * **80 of 108 cells pass** at < 1%;
  * **28 fail, and they are a result, not noise**: all 18 cells at T_eq = 300 K
    (clamps **100% `in_box_out_of_hull`** — the cold-dilute corner 650-1000 K,
    10-100 kg/m³, where the campaign has FP at T = 650/800 only for ρ ≥ 100 and
    at ρ = 10-60 only for T ≥ 1000), plus 9 at T_eq = 1300 K and 1 at 500 K over
    small in-box counts;
  * the T_eq = 700-1100 K belt the manuscript lives in is **clean**.

  CLAMP_LIMIT is not relaxed; the violating set is pinned by
  `validation/test_s4_in_box_clamp.py` so filling the gap or regressing both show
  up. **That corner is the campaign's remaining ask.**

  *Task 13 redone at the boundary (`figures/plot_s8_rcb_flip.py`).* A49 measured
  flip thresholds but could not measure the dispute *at* the RCB, because FP was
  clamped from 1000 K everywhere. Now **8 of the 10 strongest cells have FP
  interpolated at the HITEMP RCB** (the other 2 sit at ρ ≈ 54, in the cold-dilute
  hull gap above). At those 8 boundaries (T = 726-1008 K, ρ = 44-124,
  P = 170-283 bar):

  | quantity | range |
  |---|---|
  | HITEMP κ_R at the RCB | 1.06-4.63 cm²/g |
  | FP κ_R at the RCB | 28.7-61.6 cm²/g |
  | **FP / HITEMP** | **6.2-54.3×** |
  | flip threshold ∇_ad/∇_rad | **1.010-1.208×** |
  | boundaries where dispute > threshold | **8 of 8** |

  So the dispute exceeds what the boundary needs by a factor of 5-50, and the
  composite duly moves every boundary — by **0.39-0.85 dex in density** (the
  composite RCB is 2.5-7× less dense). **But the enclosed mass barely moves:
  0.0000-0.0006 of the envelope**, and |ΔR_p|/R_p over the ten strongest cells is
  **0.36-0.89%**. A49's "the boundaries are not robust" is confirmed and now
  measured on real FP data rather than a clamp.

  *The caveat that governs all of it.* Those 6-54× ratios are at T = 726-1008 K,
  which is exactly where A52 finds FP's transparency windows at or below the
  dipole noise floor and A53 finds the FP/HITEMP Rosseland ratio measuring window
  depth rather than opacity. The boundaries move because FP has no transparency
  windows there. Whether that is physics or floor is A52's open question.

- **A56 — Cooling on the new tables: the old fiducial's divergence is enormous
  and radiogenically unreachable; the water worlds' is small and real.**
  `figures/plot_s9_case_evolution.py`, composite vs HITEMP, T_int grid 12-450 K
  at 44 and 88 nodes, `outputs/S9_case_evolution.csv`.

  | case | ΔR@1 | ΔR@5 | ΔR@7 | ΔR@10 | Δs@10 | 44→88 |
  |---|---|---|---|---|---|---|
  | old fiducial (M_c 5, f_env 0.05) | 0.403% | 1.305% | 1.537% | **1.653%** | **35.08%** | converged |
  | M_c 15, wmf 0.50 | −0.000% | −0.0004% | 0.022% | 0.078% | 0.96% | **sign flips — not a result** |
  | M_c 10, wmf 0.50 | 0.000% | 0.010% | 0.045% | 0.121% | 1.35% | converged |
  | M_c 5, wmf 0.50 | 0.000% | 0.030% | 0.103% | 0.237% | 2.15% | converged |

  **Cooling-time ratio** (age to reach the line list's 5 Gyr interior entropy):
  old fiducial **6.92×** — the composite needs 5 Gyr to reach what HITEMP
  reaches at 0.723 Gyr. Water worlds 1.002 / 1.009 / 1.052×. Where the composite
  never reaches the reference entropy inside the T_int grid the ratio is re-posed
  on the composite's own entropy rather than reported as NaN.

  **Convergence.** A44 found a 24→44 node change flipping the sign of a
  sub-percent difference. Three of four cases now converge under 44→88, because
  the differences are far larger than the age integral's noise. The exception is
  M_c = 15, wmf = 0.50, whose ΔR@5 moves −0.0004% → +0.0024%: both are zero to
  within the quadrature, and it is labelled "not a result" rather than quoted.

  **Escape window (A40's null) survives**: max |ΔR/R| 3.1e-5 to 9.0e-5 %, max
  |Δs/s| 1.3e-4 to 1.9e-4 % over the first 100 Myr, all four cases — matching
  A40's 5.7e-5 % / 1.3e-4 % on the old tables. The mass-loss null is untouched by
  the campaign.

  **The decisive cut — radiogenic reachability (A30).** The old fiducial reaches
  its radiogenic T_int floor at **0.386 Gyr**, where ΔR is **0.079%** and Δs
  **0.705%**. Everything that makes it look spectacular — 1.65%, 35%, 6.92× —
  accumulates *below* the floor, at internal temperatures core heating never lets
  it reach. **The water worlds never reach their floor within the track at all**,
  so their whole 10 Gyr history is physically available. This inverts the
  structural ranking: wmf = 0.50 cells show the largest |ΔR_p| at fixed T_int
  (0.67-0.89%, A55) but the smallest *evolutionary* divergence, because a thick
  envelope buries the opacity-sensitive layer in a much larger energy budget —
  and they are the only cases whose divergence is reachable.

- **A57 — A structure build could hang forever; it now fails cleanly.** Task 18's
  first two attempts ran 2.6 and 2.2 CPU-hours without finishing, on a job that
  takes 15 minutes. Two causes, both fixed:

  1. **Duplicated work.** `evolution.history.compare_histories` called
     `_deepest_rcb_rho`, which rebuilt every structure the cooling track had
     just built, purely to read back RCB densities. `CoolingTrack.rcb_rho_si`
     now records them as the track is built; `_deepest_rcb_rho` survives only as
     a fallback for hand-constructed tracks. Halves every history comparison.
  2. **An unbounded ODE solve.** `scipy.solve_ivp` with LSODA can keep halving
     its step on a stiff segment instead of failing, so the call never returns
     and the run hangs with **no diagnostic at all** — no timeout, no iteration
     cap, nothing. `structure.integrator._EvalBudget` now caps right-hand-side
     evaluations per solve (`numerics.max_rhs_evals`, default 200 000 — a healthy
     integration spends under 1000, so 200× headroom) and raises a clear
     `RuntimeError` instead. Every caller already handles a failed build: the
     domain scan skips and logs it, and the shooting solver swallows it on its
     warm-start attempt and falls back to the bracketed solve, so one bad node
     costs seconds rather than a run. With the guard in place the four-case job
     completes in 894 s.

  *Note for anyone who meets it:* `capi_return is NULL / Call-back
  cb_f_in_lsoda__user__routines failed` on stderr is f2py reacting to *any*
  Python exception raised inside the LSODA callback — including a deliberate
  one. It is not itself evidence of a numerical fault.

## 2026-09-14 (relocating the structural demonstration)

- **A41 — At the current fiducial the deep closure is decided by CLAMPED FP
  opacity, not by FP data.** `clamp_fraction_radiative` (A26) judges reliability
  in the radiative zone, which is the right measure for a radiative zone — but
  the headline claim is that the composite *closes* a window the line list opens,
  i.e. it is about a layer the composite has made **convective**, where that
  metric is blind by construction. Added `StructureResult.clamp_fraction_all`
  (clamped fraction over every sampled profile point; the pre-existing
  `clamp_fraction`, which counts the integrator's own steps, is unchanged) and
  `structure/diagnostics.py` (`clamp_audit`, `clamp_reasons`, CLI) which splits
  evaluations by the sub-source actually used and by reason. At the fiducial
  (M_c=5 M⊕, f_env=5%, T_eq=500 K, T_int=20 K, default blend, n_sample=1500):

  | | composite | HITEMP |
  |---|---|---|
  | samples with ρ ∈ [400, 550] | 46 (T = 633–763 K) | 7 (T = 630–637 K) |
  | of those, clamped | **46 / 46, all below FP's 1000 K grid floor** | 0 / 7 |
  | clamp_fraction_all | 0.4855 | 0.3791 |
  | clamp_fraction_radiative | 0.0000 | 0.0000 |
  | whole-profile clamps, by reason | 277 below FP T floor + 1179 in-box-but-outside the scattered hull (A6) | 1137 above HITEMP's ρ ceiling (1000 kg/m³), deep convective interior only |

  **Plainly: the closure of the deep line-list radiative zone at this case rests
  on nearest-neighbour extrapolated FP opacity.** The FP grid starts at 1000 K;
  the envelope there is 630–760 K. The line-list side of the comparison is fully
  in-grid (its 7 samples interpolate; its clamps are all in the deep convective
  interior above ρ = 1000 kg/m³, where opacity is irrelevant). So A28/A35's
  closure is not evidence *from* first-principles data at this case — it is the
  nearest FP grid point (1000 K) applied 300 K too cold. The structural
  demonstration has to move to a case whose line-list radiative zone lies inside
  the FP domain (T ≥ 1000 K, ρ ≥ 108.796), which is the A42 scan.

- **A42 — Scan result: the line-list radiative zone essentially never reaches
  inside the FP data box. One usable cell in 40.** `figures/plot_s4_domain_scan.py`
  builds HITEMP and composite structures over M_c = 5 M⊕ × T_eq ∈ {500, 600,
  700, 800} K × T_int ∈ {20, 30, 40, 50, 75} K × f_env ∈ {5%, 10%} (n_sample =
  1500) and measures each radiative layer's (T, ρ) extent and mass share, the
  radiative mass fraction inside the FP box (T ≥ 1000 K **and** ρ ≥ 108.796
  kg/m³), and the composite's `clamp_fraction_all` (A41). New diagnostics
  `structure.diagnostics.radiative_layers` / `in_box_radiative_mass_fraction`.

  The two requirements pull against each other: a *deep* radiative zone needs a
  cold interior (T_int = 20 K gives max radiative ρ ≈ 100–500) and a *hot* one
  needs strong irradiation, which drives the RCB shallow (T_int ≥ 30 K gives
  ρ ≤ 60; T_int = 75 K gives ρ ≈ 1–2). Results:

  * **38 of 40 cells have exactly zero radiative mass inside the FP box.**
  * The two exceptions are T_eq = 700 K, T_int = 20 K: f_env = 5% (**365 ppm** of
    the envelope mass, layer T 1004–1023 K, ρ 127–130) and f_env = 10% (212 ppm,
    T 1002–1021 K, ρ 119–123). At f_env = 10% the composite still leaves an
    in-box radiative layer and clamps 3.2% of its samples.
  * **No cell has a line-list radiative layer lying entirely inside the box** —
    in every case the zone starts in the dilute, sub-1000 K atmosphere and only
    its base dips in.
  * `clamp_fraction_all` > 1% in 38 of 40 composite runs: unless the deep
    envelope is hot, the composite is extrapolating FP.

  **Candidate fiducial: T_eq = 700 K, T_int = 20 K, f_env = 5%, M_c = 5 M⊕** —
  the only cell where an in-box line-list radiative layer exists *and* the
  composite closes it with `clamp_fraction_all = 0.0000` (every FP evaluation in
  the whole profile interpolated). There: HITEMP radiative T 848–1023 K, ρ up to
  130.5, mass fraction 0.0032, of which 365 ppm in-box; composite radiative
  T 848–1000 K, ρ up to 122.3, zero in-box; R_p 1.75904 → 1.75944 R⊕
  (**+0.023%**). The demonstration is real but thin: it rests on ~4×10⁻⁴ of the
  envelope mass and an 8 kg/m³ slice of density.

- **A43 — G5 was not like-for-like: their T_irr IS our T_eq. Fixing the two
  definition mismatches takes the median residual 16.5% → 8.4%.** Before
  attributing anything to physics, the three matching definitions were checked
  against Aguichine et al. (2021) itself:

  1. **Irradiation (was wrong).** Their Equation (9) defines the irradiance
     temperature T_irr = T_eff·√(R⋆/2a) — the factor 2 is in it, so T_irr *is*
     the zero-albedo, full-redistribution equilibrium temperature our `t_eq_k`
     means. Verified against their own Table 3 (GJ 9827, whose data they take
     from the NASA Exoplanet Archive). With Rodriguez et al. (2018) stellar
     parameters (T_eff = 4340 K, R⋆ = 0.602 R_☉) the formula gives 1181/818/686 K
     for b/c/d against their tabulated **1184/820/686 K** (≤0.3%); with Niraula
     et al. (2017) values (4255 K, 0.651 R_☉) it gives 1204/834/698 K (≤1.8%).
     Dropping the factor 2 gives 1671/1157/969 K — 41% high under either
     parameter set, so the conclusion does not hinge on the stellar radius. A32
     had applied T_eq = T_irr/√2, i.e. modelled every anchor planet one √2 *too
     cold*.
  2. **Radius (was slightly wrong).** Their reported radius is the transiting
     radius: interior boundary plus the atmosphere thickness taken at **0.1 Pa**
     (their atmospheric grid's output). Ours is the top-of-atmosphere radius at
     `boundary.p_top_bar` = 1e-5 bar = 1 Pa — one decade deeper, so slightly
     small. (It is *not* the 10-bar matching radius: matching pressure sets where
     the interior integration starts, not where R_p is measured.)
  3. **Core (already right).** Their x_core = CMF/(1−WMF) is the core fraction of
     the refractory part only and Figure 8 has no metallic core, so x_core = 0 ⇒
     our `cmf = 0` with M_c = (1−WMF)·M_p. No correction.

  Corrections are opt-in `anchor.*` config flags (`t_irr_to_t_eq`,
  `radius_pressure_bar`, `cmf`), **defaults unchanged**, so
  `outputs/g5_residuals.json` still reproduces exactly. Median |ΔR_p/R_p| over
  the eight anchor points (`figures/plot_s5_g5_audit.py`,
  `outputs/g5_residuals_conventions.json`):

  | conventions | median |rel| | vs baseline |
  |---|---|---|
  | baseline (A32) | 0.1653 | — |
  | + irradiation (T_eq = T_irr) | **0.0997** | −0.0656 |
  | + radius level (0.1 Pa) | 0.1534 | −0.0119 |
  | + both | **0.0844** | −0.0809 |

  So **about half of the published G5 residual was a definition mismatch, not
  physics**, and with both corrections the anchor sits at 8.4%, inside SPEC §7's
  5–10% aim. Every corrected residual is still negative (our radii remain small),
  and the WMF trend survives (−3.9 to −9.9% at WMF = 0.2, −6.9 to −17.9% at
  0.5), so the pure-steam-envelope vs layered-water-structure difference (A32)
  remains the residual physics. For the attribution sentence: the 1–100 bar
  matching-pressure scan moves R_p by only **0.349%** at the fiducial
  (1.73943 / 1.73887 / 1.73338 R⊕ at 1 / 10 / 100 bar), so the boundary choice
  cannot account for it either. **Recommendation (owner's call): make
  `anchor.t_irr_to_t_eq = "identity"` the default and re-report G5.**

- **A44 — At the candidate the evolution signal is ~18× smaller than at the old
  fiducial and below the age integral's own grid convergence; the cooling-time
  ratio is 1.000 ± 0.002.** `config/candidate.yaml`
  fixes the A42 case (M_c = 5 M⊕, f_env = 5%, T_eq = 700 K; structural case
  T_int = 20 K) and `figures/plot_s6_candidate_evolution.py` runs the E2
  comparison and the fixed-entropy τ_KH ratio there. The cooling ages are
  integrated on a 24-node T_int grid first and re-run on E2's 44-node grid if the
  A38 age check trips (`evolution.age_grid_problem`, new; `HistoryComparison` now
  carries its two tracks so a caller can check them). The coarse grid passed its
  A38 check, but a `--grid fine` re-run shows the differences at this case are
  **not converged on 24 nodes** — the signal is ~0.1–1%, comparable to the age
  integral's own grid error, and the 5 Gyr difference even changes sign. Both are
  reported; the 44-node numbers are the ones to quote, and the figure and CSV in
  `outputs/` are the fine run:

  | age | Δs [J/kg/K] | Δs/s | ΔR_p [R⊕] | ΔR_p/R_p | (24-node: Δs, Δs/s) |
  |---|---|---|---|---|---|
  | 5 Gyr | +7.91 | +0.082% | +0.000070 | +0.0040% | −13.61, −0.140% |
  | 7 Gyr | +37.68 | +0.398% | +0.000373 | +0.0213% | +22.96, +0.242% |
  | 10 Gyr | +128.32 | +1.388% | +0.001578 | +0.0901% | +109.45, +1.183% |

  FP-relevant epoch 4.86–16.50 Gyr on the fine grid (4.09–16.37 Gyr coarse; the
  deepest RCB first reaches ρ ≥ 108.796 at 3.90 Gyr for the composite, 3.81 Gyr
  for HITEMP — those crossing ages are structural, not integrated, and move by
  <2%); escape window clean as ever (max |ΔR_p/R_p| = 8.61e-5 %, max |Δs/s| =
  1.55e-4 % over the first 100 Myr; 5.45e-5 % / 8.06e-5 % coarse).
  τ_KH(composite)/τ_KH(HITEMP) at fixed interior entropy over T_int,ref ∈ {30, 40}
  K × f_env ∈ {0.02, 0.05, 0.10}: **0.9981–1.0020** — no cooling signal at all.

  Three honest caveats. (i) Compared with the old T_eq = 500 K fiducial (Δs/s =
  +7.2% at 7 Gyr, A37) the candidate's entropy signal at 7 Gyr is **+0.40%**,
  ~18× smaller: moving into the FP domain buys provenance and costs signal. (ii) The candidate is
  clamp-free only at its *structural* case (T_int = 20 K, clamp_fraction_all = 0);
  the fixed-entropy models at T_int,ref = 30–40 K clamp 1.4–35% of their samples,
  so the candidate's cooling numbers do not inherit the structural case's
  provenance. (iii) At this size the differences are at the level of the age
  integral's grid convergence (see the table above), so they should be quoted as
  "no cooling signal resolvable at this case", not as a measured Δs.

- **A45 — The candidate's T_int = 20 K is below the radiogenic floor at every
  age: the in-domain demonstration is not reachable either.** Running the D1/A30
  floor calculation at the candidate (M_c = 5 M⊕, R_p = 1.75944 R⊕, Earth-like
  20 TW per M⊕ at 4.5 Gyr, effective halving every 6 Gyr):

  | age | L_radio | T_int floor |
  |---|---|---|
  | 1 Gyr | 149.8 TW | **35.95 K** |
  | 4.5 Gyr (present day) | 100.0 TW | **32.49 K** |
  | 10 Gyr | 53.0 TW | 27.72 K |

  The candidate's structural case is T_int = 20 K, i.e. **12.5 K below the
  present-day floor and 16 K below the 1 Gyr floor**. v1 omits core luminosity
  (SPEC §4), so the solver reaches 20 K happily; a real 5 M⊕ Earth-like core does
  not let the planet get there. Combined with A42 (the FP box is touched only at
  T_int = 20 K) this closes the loop A30 opened: **the regime where the
  first-principles data actually constrain a line-list radiative zone is, for
  steam sub-Neptunes on a 5 M⊕ core, radiogenically out of reach.** Gate
  `validation/test_s7_radiogenic_floor.py` (the floor stays above 20 K for
  R_p ∈ [1.5, 2.5] R⊕ at 1, 4.5 and 10 Gyr).

- **A46 — The corrected G5 conventions are now the defaults: median residual
  16.5% → 8.4%.** A43 established that Aguichine et al. (2021) Eq. (9)
  T_irr *is* our T_eq and that their reported radius is the transiting radius at
  0.1 Pa, and left both corrections opt-in. They are now the `Anchor` defaults
  (`t_irr_to_t_eq = "identity"`, `radius_pressure_bar = 1e-6` bar); `cmf = 0` was
  already right and is unchanged. `outputs/g5_residuals.json` is regenerated and
  now records the conventions it was computed under.

  | | pre-A46 | A46 default |
  |---|---|---|
  | median \|ΔR_p/R_p\| | 0.1653 | **0.0844** |
  | worst | 0.2279 | 0.1787 |
  | WMF = 0.2 points | −6.2 to −15.0% | **−3.9 to −9.9%** |
  | WMF = 0.5 points | −18.1 to −22.8% | **−6.9 to −17.9%** |

  The G5 gate now *pins* the median at 0.0844 with a 1% tolerance
  (`_MEDIAN_A46`), alongside the existing per-point 30% bound: with the
  definitions corrected the published number is a gate, not a report, so a silent
  change of EOS, opacity or definition cannot pass unnoticed.

  **The residual that remains is attributed to the reference model's layered
  supercritical-water structure.** Every corrected residual is still negative and
  still grows with WMF, which is exactly where our single pure-steam envelope
  departs from their ice/liquid/supercritical water layers. The two competing
  explanations are measured and too small or the wrong sign: the 1–100 bar
  matching-pressure scan moves R_p by **0.349%** (A43), and recomputing the eight
  residuals with full-band HITEMP κ_R instead of the [100, 4000] cm⁻¹ science
  band moves the median **0.0844 → 0.1147**, i.e. band truncation makes the
  agreement *worse*, not better (A40 measured the same effect pre-A46, 0.165 →
  0.172; re-measured here under the corrected conventions).

  Housekeeping: `figures/plot_s5_g5_audit.py` now states both flags explicitly in
  every convention (an empty override dict would silently make its baseline equal
  its endpoint), and `test_repository_defaults_are_unchanged` is renamed
  `test_repository_defaults_are_the_corrected_conventions`. A43's "defaults
  unchanged" is superseded by this entry, not rewritten. **Stale, flagged not
  edited:** `docs/band_limitation.md`'s G5 row still quotes the pre-A46 0.165 →
  0.172.

- **A47 — Extended scan, 108 cells: NO planet in it puts a line-list radiative
  layer inside the FP data box, and no composite run is clamp-free.**
  `figures/plot_s4_domain_scan.py --extended` opens the A42 grid in every
  direction that could plausibly help — T_eq ∈ {300, 500, 700, 900, 1100, 1300} K
  × T_int ∈ {30, 50} K × M_c ∈ {5, 10, 15} M⊕ × wmf ∈ {0.05, 0.20, 0.50}, both
  sources, n_sample = 1500 (5m08s). **All 108 cells converged inside the AQUA
  table; none was skipped.**

  *Mass convention.* `envelope.f_env` is already the envelope's share of the
  *whole* planet (the model sets M_p = M_c/(1−f_env)), so the requested `wmf` and
  the existing `f_env` are the same number. Both are recorded, equal by
  construction, next to the implied M_p — up to 30 M⊕ at M_c = 15, wmf = 0.5.

  | | result |
  |---|---|
  | cells with **any** HITEMP radiative mass inside the FP box | **0 of 108** |
  | same for the composite | 0 of 108 |
  | cells with composite `clamp_fraction_all` < 1% | **0 of 108** (range 0.011–1.000) |
  | HITEMP `clamp_fraction_all` | 0.279–1.000 |

  **Why: the two requirements move in opposite directions across the whole
  space.** Taking the densest radiative zone at each T_eq (M_c = 15 M⊕,
  wmf = 0.05, T_int = 30 K):

  | T_eq | radiative layer T | max radiative ρ | in the FP box? |
  |---|---|---|---|
  | 300 K | 363–364 K | **357** | no — 640 K too cold |
  | 500 K | 606–732 K | **122** | no — 270 K too cold |
  | 700 K | 848–1007 K | 45.8 | no — 2.4× too dilute |
  | 900 K | 1090–1285 K | 32.3 | no — 3.4× too dilute |
  | 1100 K | 1333–1574 K | 27.8 | no — 3.9× too dilute |
  | 1300 K | 1575–1850 K | 26.4 | no — 4.1× too dilute |

  Irradiation heats the radiative zone and thins it at the same rate: by the time
  the layer is hot enough for the FP grid (T ≥ 1000 K, reached at T_eq ≳ 900 K) it
  sits at ρ ≈ 30, a factor 3–4 below the grid's 108.796 floor; where it *is* dense
  enough (T_eq ≤ 500 K) it is 270–640 K too cold. The A42 exceptions were both at
  **T_int = 20 K**, which this grid does not contain: the in-box window is already
  shut at T_int = 30 K everywhere in the extended space, which is consistent with
  A42 and narrows the in-domain demonstration to an even smaller corner than A42
  suggested — one that A45 shows is radiogenically unreachable anyway.

  Because the colour scale (in-box ppm) is zero in all 108 cells, each cell is
  also annotated with how close its line-list radiative zone gets on *both* axes
  (max radiative ρ and that zone's max T, new column
  `{source}_max_radiative_t_k`) — otherwise the figure would be a grid of honest
  zeros saying nothing, and the anti-correlation above is the actual result.

  Mechanics worth keeping: `structure.diagnostics.aqua_coverage` (AQUA *clips*
  out-of-table queries rather than raising, A2, so a scan has to test for it
  explicitly), per-cell skip logging for non-convergence, and the new columns
  (total radiative mass fraction, the envelope's hottest/densest point, and
  `clamp_fraction_all` for **both** sources). Three figures, one per core mass —
  the one-figure-per-script convention is overridden here on request. Gate
  `validation/test_s4_extended_scan.py`.

- **A48 — Coverage, trajectory by trajectory: at the old fiducial 0.09% of the
  envelope mass is FP data; at the candidate 99.9998% is.**
  `figures/plot_s5_grid_coverage.py` draws the 27 FP grid points and their
  interpolation hull on the (T, log ρ) plane with six envelope trajectories over
  it, every sample coloured by what the composite opacity call there actually
  was. `structure.diagnostics.classify_evaluations` assigns one of five
  categories per evaluation and `category_mass_fractions` weights them by
  envelope mass. The fifth category is the one the clamp counters never
  separated: **in-box but outside the hull** — inside the FP grid's (T, ρ) box
  yet outside its scattered convex hull (A6) — which is a gap that can be filled
  without extending the grid's range at all, and it dominates.

  | trajectory (T_eq, T_int, wmf, M_c) | interpolated | in-box, out-of-hull | clamped in T | clamped in ρ |
  |---|---|---|---|---|
  | (a) old fiducial 500/20/0.05/5 | **0.000892** | **0.954237** | 0.044870 (all too cold) | 0 |
  | (b) candidate 700/20/0.05/5 | **0.999998** | 0 | 0 | 0 |
  | (c1) 300/30/0.05/15 | 0.000000 | 0 | 0.734286 (all too cold) | 0.265714 |
  | (c2) 500/30/0.05/15 | 0.193364 | 0.539580 | 0.001346 (too cold) | 0.265709 |
  | (c3) 700/30/0.05/15 | 0.254032 | 0 | 0.745967 (all too **hot**) | 0 |
  | (d) hot 1100/30/0.20/5 | 0.011916 | 0 | 0.988083 (all too **hot**) | 0 |

  The "clamped in T" column is split because it hides two different problems:
  *too cold* for the 1000 K floor (a real gap — it is where the old fiducial's
  closure happens) versus *too hot* for the 3000 K ceiling, which is the deep
  convective interior where opacity does not shape the structure at all. (c)
  could not be "the three cells with the largest in-box radiative fraction"
  because A47 found that fraction is zero in all 108; the tie is broken on max
  radiative ρ, one cell per T_eq, and the script says so when it does it.

  **Campaign specification** (`outputs/S5_grid_coverage.csv`), on the requested
  100 K × 0.25 dex lattice, taking the four bracketing corners of each needed
  sample's cell so enclosure is guaranteed — and then *verified* with the same
  triangulated interpolator the model uses, not assumed:

  * **55 nodes**, T ∈ [1000, 2700] K, ρ ∈ [100, 3162] kg/m³, enclose the in-box
    part of all six trajectories (1898 samples). Enclosure: **verified**.
  * Addendum, reported separately because it is a different ask: **47 nodes**
    over T ∈ [300, 1000] K at the same densities would be needed to cover the
    dense-but-cold segments where the composite currently clamps in T (1337
    samples). This is the extension the old fiducial's headline actually depends
    on.

  Corner nodes may sit up to one lattice step outside the grid's present range;
  that is what enclosing an edge sample requires. Gate
  `validation/test_s5_grid_coverage.py` (categories partition the mass, the clamp
  categories reproduce the source's own clamp mask, and the emitted nodes really
  do enclose the samples they claim to).

- **A49 — Every RCB but two is within 8% of flipping: the boundary positions are
  far more opacity-sensitive than the opacity data are certain.**
  The model's criterion is ∇_rad > ∇_ad with ∇_rad ∝ κ_R, so at fixed structure
  the factor by which κ_R would have to change for a sample to flip is exactly
  ∇_ad/∇_rad there. `structure.diagnostics.rcb_opacity_sensitivity` reports it at
  each radiative layer's deepest sample ("close it": κ up, the layer shrinks by
  one sample) and at the first convective sample below ("extend it": κ down, the
  layer grows by one sample), for the A48 trajectories under HITEMP:

  | trajectory | RCB (T, ρ, P) | κ_R there | close it | extend it |
  |---|---|---|---|---|
  | (a) old fiducial | 629.5 K, 285.3, 180 bar (+2 more at 515/522) | 0.0564 | **×4.098** (layer 1) / ×1.068 (layer 2) | ×0.866 / ×0.824 |
  | (b) candidate | 1025.0 K, 131.7, 615 bar | 0.5103 | **×1.001** | ×0.976 |
  | (c1) 300 K | 364.2 K, 423.6, 0.955 bar | 0.0224 | **×13.1** | ×0.992 |
  | (c2) 500 K | 736.5 K, 124.1, 283 bar | 0.1136 | ×1.030 | ×0.964 |
  | (c3) 700 K | 1014.0 K, 46.7, 216 bar | 0.4737 | ×1.048 | ×0.962 |
  | (d) hot | 1554.8 K, 13.3, 95.4 bar | 2.222 | ×1.075 | ×0.993 |

  Read against the FP/line-list difference of **2–3×** in the dense regime (A28,
  Fx): a boundary needing ×1.001–×1.075 would be moved by an opacity correction
  far smaller than the one actually in dispute, so the *position* of these RCBs
  is not a robust prediction — which is the honest framing for the manuscript's
  RCB numbers. The two exceptions run the other way and are the informative ones:
  the old fiducial's deep radiative window needs **×4.1** to close, i.e. it is
  genuinely a 2–3× -scale question that FP data could settle, and the 300 K case
  needs ×13, i.e. nothing plausible touches it.

  This is a first-order diagnostic: it holds the converged (P, T, ρ) profile
  fixed and asks how close each boundary is to flipping, not where a
  self-consistent model with a different opacity would put it. Printed by
  `figures/plot_s5_grid_coverage.py`; gate
  `validation/test_s6_rcb_sensitivity.py` pins the arithmetic against the
  structure model's own RCBs and layer masses.

## 2026-09-13 (Nature manuscript checks)

- **A35 — Blend prescription is a config knob and a stated systematic; the deep
  radiative zones close under every prescription.** The composite's surviving RCB
  at the cold fiducial sits *inside* the default blend window (ρ ≈ 190–234 vs the
  window 100–300), so its location is set partly by the blend prescription and
  had to be scanned rather than assumed. Added
  `opacity.blend_window_kg_m3` (default `[100, 300]`, unchanged behaviour) as the
  canonical name — `opacity.composite_blend_rho` is kept as a deprecated
  read-only alias so existing figure scripts are untouched — and
  `opacity.blend_mode` (`linear` | `sharp`, default `linear`). **`linear` labels
  the pre-existing C¹ smoothstep log-blend (A12), not a piecewise-linear weight**;
  `sharp` is a hard switch to FP at the window's lower edge with no blending, for
  which a zero-width window (e.g. `[108.796, 108.796]`, the FP grid floor) is
  legal and is refused in `linear` mode. The built-source cache key now carries
  these knobs (`opacity.build_tables.cache_key`), so a scan rebuilds instead of
  silently reusing the first source; the constant-κ injections in the gates and
  `figures/_common` use the same helper. New diagnostic
  `structure.deep_radiative_zones(result, rho_min)` returns the contiguous
  radiative runs denser than a threshold — `[] ` is the closure statement.
  `figures/plot_s1_blend_sensitivity.py` (T_int=20 K, T_eq=500 K, n_sample=1500,
  band [100,4000]) gives, with **HITEMP: 3 RCBs at ρ = 285/515/522, radiative to
  ρ = 504, one deep radiative zone (ρ 450–504)**:

  | case | #RCB | RCB ρ [kg/m³] | max radiative ρ | deep zones (ρ>300) | R_p [R⊕] |
  |---|---|---|---|---|---|
  | [100,200] linear | 1 | 191 | 127 | none | 1.7223 |
  | [100,300] linear (default) | 1 | 218 | 127 | none | 1.7165 |
  | [100,500] linear | 1 | 234 | 172 | none | 1.7139 |
  | sharp @ 108.796 | 3 | 105/116/112 | 109 | none | 1.7257 |

  The RCB density moves 191→234 (±10% about the default) across the blend
  windows, and the sharp switch instead pins the boundary at the switch density
  itself (the discontinuous κ makes the profile flip across it three times) — but
  **no composite prescription leaves a radiative layer at ρ > 300**, so the
  closure claim (A28) is prescription-independent while the RCB *density* is not.
  Gate `validation/test_s1_blend_window.py` asserts closure for every case, with
  HITEMP's deep zone as a positive control so the assertion is not vacuous.

- **A36 — The deep radiative zones are thin slivers at the steam→dense-water
  density jump; report max-radiative-ρ and RCB density separately.** Resolving
  the HITEMP case shows the deep radiative layers sit exactly where AQUA's
  density jumps from ρ = 151 to 538 kg/m³ across P ≈ 176→184 bar at T ≈ 630 K
  (vapour → dense water). ρ is locally non-monotonic there (538, 541, 491, 450,
  504, 541), which is why the deepest *RCB densities* (515, 522) exceed the
  densest *radiative* point (504). Consequence for the manuscript: the numbers
  quoted as "RCB at 453 / 123 kg/m³" are **max-radiative-ρ**, not RCB densities;
  both are now reported side by side in every table. Presence of the zone is
  robust to sampling (3 RCBs and one deep zone at n_sample = 600/1000/1500/2000)
  but its extent is not (max radiative ρ = 439/513/504/477), so gates assert
  presence/absence, never a count or a density.

- **A37 — `opacity.fp_scale`: the conclusions survive a ±20–25% perturbation of
  the FP table.** Added `opacity.fp_scale` (default 1.0), a constant multiplier
  on the FP table's κ_R and κ_P applied at table-build time, i.e. *before* the
  composite blends, and carried in the source cache key. The SPEC §2
  factor-1000 plausibility assertion is evaluated on the **unperturbed** means,
  since it guards a data error and not a deliberate perturbation.
  `figures/plot_s2_opacity_scaling.py` repeats the cold-fiducial structure and
  the E2 cooling comparison (E2's T_int grid; the HITEMP track carries no FP, so
  it is built once and reused) for fp_scale ∈ {0.8, 1.0, 1.25}:

  | fp_scale | #RCB | RCB ρ | max radiative ρ | deep zones | Δs(7 Gyr) [J/kg/K] | Δs/s | ΔR_p(7 Gyr) | ΔR_p/R_p |
  |---|---|---|---|---|---|---|---|---|
  | 0.80 | 1 | 219.0 | 127.9 | none | +533.87 | +7.13% | +0.00467 R⊕ | +0.2725% |
  | 1.00 | 1 | 217.9 | 127.3 | none | +539.56 | +7.21% | +0.00467 R⊕ | +0.2729% |
  | 1.25 | 1 | 216.8 | 126.8 | none | +544.82 | +7.28% | +0.00467 R⊕ | +0.2725% |

  A ±25% change in FP's κ moves the RCB density by ~1% and the 7 Gyr entropy
  bias by ~2% *of itself*; closure holds and the entropy bias keeps its sign
  (positive: the composite stays hotter/higher-entropy at fixed age). The
  composite−HITEMP difference is therefore set by *where* FP takes over, not by
  FP's magnitude. Note the diagnostic asymmetry the manuscript should quote:
  **the interior-entropy bias at 7 Gyr is ~7%, the radius bias only ~0.27%** —
  entropy is the sensitive observable-adjacent quantity, R_p is not. Gate
  `validation/test_s2_fp_scale.py`.

- **A38 — Cooling-track ages need ≳20 T_int nodes; below that the quadrature
  returns negative ages.** Building the fp_scale gate on a 14-node grid gave a
  track spanning [−8.35, 2.79] Gyr: `cumulative_simpson` of 1/L against E
  overshoots when L varies by ~4 decades across too few nodes (A24 already found
  the G6 residual needs n≈60 for 2e-3; the failure mode below ~20 is worse than
  inaccurate, it is unphysical). This is with the real composite/HITEMP opacity,
  where L(T_int) is steepest; the existing constant-κ evolution gates on 12-node
  grids are unaffected, which is why it had not surfaced before. No guard added
  inside `cooling_track` (that
  would change existing gates' behaviour); instead every age-based gate uses
  E2's 44-node grid, and `figures/plot_s2_opacity_scaling._at_age` now **raises
  rather than extrapolating** — `np.interp` would silently clamp "the difference
  at 7 Gyr" to the end of a short track.

- **A39 — Band truncation quantified for HITEMP; it is *not* the cause of the G5
  low bias.** SPEC §3 defers the absolute band systematic to a full-band vs
  common-band comparison for the line list. `opacity/build_tables.py` gained a
  CLI — `python opacity/build_tables.py config/fiducial.yaml --source hitemp
  --band full|common|lo,hi` — writing per-grid-point means (csv + npz);
  `--band full` resolves to the source's **own computed range read from the
  files**, [10, 11000] cm⁻¹ for HITEMP (76 of the 77 (T, ρ) points exist; T=500,
  ρ=0.1 was never delivered), never wider. Results
  (`figures/plot_s3_band_truncation.py`):

  * κ_R(full)/κ_R(common): **min 0.018 (3000 K, 1000 kg/m³), median 0.343, max
    0.994 (500 K, 10 kg/m³)** — i.e. the science band biases κ_R **high**, by
    ~3× at the median and up to ~55× at 3000 K, exactly the direction
    `docs/band_limitation.md` predicted (the harmonic Rosseland mean is
    dominated by the near-IR/visible transparency windows the band excludes).
    The effect is almost entirely a function of T (0.99 at 500 K → 0.02 at
    3000 K), as the Planck peak marches out of the band.
  * κ_P(full)/κ_P(common): min 0.273, median 0.850, max 1.002 — the arithmetic
    mean is far less band-sensitive, which is why Planck is the trustworthy
    headline (A22).
  * `rel_floor` = 1e-6 leaves the **median** ratio unchanged to three digits but
    does bind at individual points (max |Δratio| = 0.137, 30% relative), so the
    floor-free ratio is reported alongside in `outputs/S3_band_truncation.csv`.
  * **G5 median |ΔR_p/R_p| goes 0.1653 (common) → 0.1720 (full)**
    (`outputs/g5_residuals_fullband.json`). The common-band column reproduces
    the committed gate's 0.16530 exactly (same code path). The T_irr = 400 K
    points do not move at all; the T_irr = 1000 K points get worse by 0.7–1.1
    percentage points, since the less opaque full-band tables give a slightly
    smaller radius and our radii were already low. **Conclusion for the
    manuscript: band truncation cannot explain the G5 residual — removing it
    makes agreement slightly worse**, leaving the pure-steam-envelope vs
    full-water-phase model difference (A32) as the dominant cause.

  Diagnostic only: `opacity.band_cm` stays [100, 4000] (A25) and the
  FP-vs-line-list differential stays common-band (A21). Gate
  `validation/test_s3_band_truncation.py`.

- **A40 — Escape-window null quantified; F3 range printed over its actual grid,
  with the (f_env × T_eq) scan opt-in.** Two reporting additions, no default
  changed.

  *(a) Escape window and the FP-floor crossing (E2).* `HistoryComparison` gained
  `escape_window_extrema(window_gyr)` and `first_age_rcb_above(rho_si)` — in the
  module, not the figure, so they are gated (`test_history_diagnostics.py`); the
  window edge is included by interpolation so a sparse age grid cannot miss it,
  and the crossing age is interpolated in log age. At the fiducial (T_eq=500 K,
  M_c=5 M⊕, f_env=5%, E2's 44-node T_int grid), over the first 100 Myr (39 age
  nodes from 9.3e-6 Gyr): **max |ΔR_p/R_p| = 5.74e-5 %, max |Δs/s| = 1.34e-4 %**
  — both *below* the shooting tolerance (`shoot_rtol` = 1e-6), i.e. numerically
  indistinguishable from zero. The escape-window null (A29) is therefore not a
  small difference but no difference. The deepest RCB first reaches ρ ≥ 108.796
  kg/m³ at **1.626 Gyr (composite)** and **1.548 Gyr (HITEMP)**; peak
  |ΔR_p/R_p| = 0.27%, FP-relevant epoch 1.87–12.68 Gyr. (`docs/
  fp_validity_regime.md` still records A29's earlier run at T_eq=600 K — epoch
  2.6–14 Gyr, peak 0.33%; the numbers here are the current E2 default, T_eq=500 K.
  The doc is left as A29's record.)

  *(b) F3 range.* F3's printed summary now reports min/max with the (T_eq, f_env)
  each extreme came from, plus the grid extents. **The plotted figure is
  unchanged**: it is still the f_env scan at the middle T_eq of the config grid
  (700 K), and the default invocation computes exactly what it did before. SPEC
  §8 asks F3 over (f_env, T_eq), so `--t-eq-scan` walks the whole config T_eq
  grid for the *printed* range only — opt-in, because it is ~2.6× the runtime
  (13 min vs 5 min, measured) and the owner's figure should not silently change
  cost. Over
  f_env ∈ [0.01, 0.1] × T_eq ∈ [500, 800] K at T_int,ref = 150 K:

  | ratio to HITEMP | min | max |
  |---|---|---|
  | composite | 0.9979 (T_eq=500, f_env=0.05) | 1.0021 (T_eq=500, f_env=0.10) |
  | POKAZATEL | 1.3670 (T_eq=800, f_env=0.01) | 2.2420 (T_eq=500, f_env=0.01) |
  | FP [extrapolated, unreliable] | 18.07 (T_eq=800, f_env=0.01) | 29.68 (T_eq=500, f_env=0.05) |

  So the manuscript can quote: **composite ≡ HITEMP to ±0.2% across the whole
  grid**, the reliable line-list systematic is **1.37–2.24×** (POKAZATEL vs
  HITEMP), and the pure-FP 18–30× is the opacity-extrapolation artifact already
  flagged in A33/F3's caption. A non-converged fixed-entropy cell is recorded as
  NaN and reported rather than aborting the scan (none occurred).

## 2026-08-04 (deck prep)

- **A32 — G5 anchor is the published Aguichine (2021) fit; coefficients are
  online-only, so the numerical anchor is digitized Figure 8.** Their M–R fit is
  Eq. 29, log₁₀R_p = a·log₁₀M_p + b + exp(−d·(log₁₀M_p + c)), with a,b,c,d
  tabulated per (composition, T_irr) in an online machine-readable table not in
  the PDF. G5 therefore anchors to 8 points read from their Figure 8 (Ma19+, no
  metallic core; read uncertainty ~5%), in
  `validation/data/aguichine2021_fig8_points.csv`. Matching: WMF→f_env
  (pure-steam envelope), T_irr→T_eq via T_irr=√2·T_eq, **cmf=0** (pure-rock core,
  matching their no-metallic-core case — `core.cmf` added to config, A9 default
  0.325), T_int=30 K. Residuals (documented, not tuned, per SPEC §7): our R_p is
  systematically smaller by **6–15% at WMF=20% (meets SPEC's 5–10% aim) and
  18–23% at WMF=50%**. Causes, in order: our pure-steam envelope vs their full
  water-phase (ice/liquid/supercritical) layers — dominant at high WMF; residual
  core-composition difference; cooler mature thermal state + Guillot-radiative vs
  their atmosphere. Gate asserts <30% (structural agreement) and writes
  `outputs/g5_residuals.json`. G5 now runs (was skipping).

- **A33 — POKAZATEL 11×11 regeneration declined (owner).** POKAZATEL stays on the
  27-point FP grid; its role is the overlap-agreement cross-check (banked, A22:
  HITEMP/POKAZATEL agree ~1.5× in κ_R, ~5% in κ_P). **HITEMP is the sole reliable
  structural line list** (full T=500–3000, ρ=0.1–1000 coverage); POKAZATEL and
  pure-FP are flagged unreliable in structure runs (radiative-zone coverage).

- **A34 — Final deck narrative framing (owner-approved, amended).** Headline
  states the failure mechanism and the DFT-MD justification up front, and frames
  the denser-envelope regime as a *mapped applicability boundary of the
  first-principles correction*, not a follow-up commitment: "First-principles
  (DFT-MD) water opacity is 2–3× higher than HITEMP/POKAZATEL at ρ≥109 kg/m³,
  T≥1000 K — the warm-dense regime where the line picture breaks down (discrete
  Voigt lines leave transparent windows a collective continuum fills; far-wing/
  continuum absorption is truncated; the line concept dissolves once collision
  time ≈ vibrational period). For steam sub-Neptunes this regime is convective
  while warm and radiogenically unreachable once old (core ≈ M_c×20 TW floors
  T_int above the ~25 K needed for the deep RCB), so the correction does not
  measurably change their radius (≤1.7%, and unreached) or observable cooling and
  leaves no atmospheric-escape imprint. The mapped applicability boundary of the
  correction is denser envelopes — massive water worlds and magma-ocean-coupled
  layers — whose RCB sits at ρ≥109 while still warm." Reflected in figure
  captions and docs.

## 2026-08-04

- **A30 — Reachability of the FP-relevant regime (checks 1–2). Radiogenic
  heating TRUNCATES it for steam sub-Neptunes.** The FP-relevant regime lives at
  T_int ≲ 25 K (where the RCB reaches ρ≥109). Two reachability checks:
  (1) *Cooling age.* From the tracks (no radiogenics), the fiducial reaches
  T_int=30 at 0.9 Gyr, 20 at 7.4 Gyr, 15 at 9.6 Gyr — old but within Hubble
  time; but the low-mass/high-f_env planets that showed the largest ΔR_p (1.7%,
  D1) reach T_int=20 only at ~61 Gyr, i.e. **never** — their large FP impact is
  unreachable.
  (2) *Radiogenic T_int floor (v1 omission sensitivity).* An Earth-like rocky
  core of mass M_c produces ≈ M_c × 20 TW present-day radiogenic power (Earth ≈
  20 TW, range 15–25; declines to ~0.5× by 10 Gyr via the U/Th/K decay
  ensemble). Setting L_int(T_int) = L_radio gives an equilibrium **T_int floor of
  ≈ 33 K (fiducial, present) → 28 K at 10 Gyr**; ≈ 24–40 K over M_c = 2–20. Since
  the floor (28–40 K for M_c ≳ 3) sits **above** the ~25 K FP threshold,
  radiogenic heating keeps these planets warm enough that the RCB never reaches
  ρ≥109 — **FP never becomes relevant to their structure.** Robust to the
  15–25 TW uncertainty (floor ∝ L^¼ → ±~10%). v1 omits core luminosity (SPEC §4),
  so the code's T_int can go arbitrarily low; this is a stated sensitivity, not a
  code change. Net: for realistic steam sub-Neptunes the FP structural effect is
  **not reached**; it is confined to the denser-envelope follow-up (RCB at ρ≥109
  while warm). See `docs/fp_validity_regime.md`.

- **A31 — Multi-RCB topology is physics, not solver noise (check 3).** For the
  fiducial at T_int=20, pure HITEMP has 3 RCBs (ρ≈[285,509,487]), composite 1
  (ρ≈212). Robust to G7-style refinement (halved tolerances + doubled sampling:
  still 3 RCBs, densities shift <5%). The exact *count* is marginally sensitive
  to the opacity value (a +3% κ_R perturbation flips 3→1), so "3 RCBs" is not a
  robust prediction — but the qualitative FP effect (its 2–3× higher opacity
  closing the deep radiative window the line list opens) far exceeds that ~3%
  flip threshold and is robust. (Realized only below the radiogenic T_int floor,
  A30 — for the fiducial this topology is not reached; marginally reachable for
  the lowest-mass cores.)

- **A29 — FP-relevant domain + history effect (owner scan).** Systematic map of
  where FP propagates to structure/evolution. `D1_fp_domain`: over (M_c 2–20,
  f_env 1–50%) at cold T_int=20, the line-list radiative zone reaches ρ≥109 in
  31/36 cells, but the structural impact |ΔR_p|(composite−line-list) is **≤1.7%**
  everywhere (largest for low-mass/high-f_env; ~0.5% fiducial). `evolution.
  compare_histories` + `E2_history`: matched-IC cooling under composite vs
  HITEMP. Findings — the FP-relevant epoch is **late (~2.6–14 Gyr)**, so a hot
  start has a *shallow* RCB and the **first-100-Myr escape window is unaffected**
  (ΔR_p≈0.00%, RCB ρ≈2–5) — a clean null for mass-loss imprint; the peak radius
  divergence (~0.33% fiducial) is reached late and **persists** to Gyr (no
  reconvergence). The strong-FP regime lies beyond steam sub-Neptunes (denser
  envelopes whose RCB is at ρ≥109 while warm) — the promoted follow-up. See
  `docs/fp_validity_regime.md`. Not asking FP to matter; this is where it does
  (little for radius, qualitatively for RCB topology, not at all for escape).

- **A27 — CORRECTED by A28. (Superseded.)** A27 claimed FP never affects these
  planets because the radiative zone stays at ρ ≲ 7 kg/m³. That was tested only
  at T_int ≥ 50 and was **wrong for cold/old planets**. Retained for the record;
  see A28.

- **A28 — FP matters for OLD/COLD planets; reconciles the ρ≈540 vs ρ≈7 conflict.**
  RCB depth is set by the opacity value and T_int: higher κ or higher T_int →
  convection starts shallower. Sweeping T_int with the *real* opacity: the
  radiative zone is dilute (ρ ≲ 7) for T_int ≳ 50 (young/warm — A27's regime),
  but the RCB **migrates deep into FP territory** as T_int drops, reaching max
  radiative ρ ≈ 60 (T_int=30), 120–450 (T_int=20), 140–680 (T_int=15). This
  reconciles the earlier ρ≈540 (that scan was cold + low constant κ). Since cold
  planets are the observed population, FP is relevant there. **Qualitative
  effect (as the owner predicted):** pure line lists open a *deep radiative
  window* at ρ≥109 that the FP-informed composite **closes** (FP is 2–3× more
  opaque there → those layers convect). At T_int=20 (fiducial) pure HITEMP has 3
  RCBs and a radiative zone to ρ=453; the composite has 1 RCB to ρ=123. The
  instantaneous R_p difference is small (~0.1–0.24%); the structural/RCB
  topology difference is qualitative. Whether it accumulates over cooling
  history (and whether it touches the young escape window) is the A29 evolution
  study. NOTE: the FP-relevant epoch is LATE (cold), so a HOT start has a
  *shallow* RCB — FP relevance grows with age, opposite to a naive expectation.
  See `docs/fp_validity_regime.md`.

- **A25 — Default science band [100, 4000] cm⁻¹ + Rosseland floor; resolves the
  A22 Rosseland discrepancy (it was the BAND, not the floor).** With HITEMP
  delivered, Fx was run FP-vs-HITEMP through our own pipeline. Under the owner's
  convention (band [100, 4000], relative floor 1e-6) our FP/HITEMP κ_R
  reproduces `crossover_hitemp.json` to ~1% (3.43 vs 3.42, 2.11 vs 2.11, 1.52 vs
  1.52); Planck likewise. Decomposition shows the entire ~2× gap was the **band**:
  [50,4500]→[100,4000] moves FP/HT κ_R 6.47→3.43, while the 1e-6 floor is
  **inactive** (all sources' in-band κ_min/κ_max ≥ 1e-6, so it never binds). The
  far-band transparency windows (50–100, 4000–4500 cm⁻¹) dominate the harmonic
  Rosseland mean and inflate the ratio for the fine line lists more than for
  smoothed FP. The default `opacity.band_cm` is therefore [100, 4000] (the
  crossover-consistent, mutually-consistent-across-FP/HITEMP/POKAZATEL core
  band); `rel_floor`=1e-6 is kept for convention parity + robustness (a knob).
  [50,4500] is the reported band systematic (~2× on κ_R). HITEMP and POKAZATEL
  agree on κ_R to ~1.5× (POKAZATEL slightly less opaque). See
  `docs/opacity_crosscheck.md`.

- **A26 — Dilute-limit flat-κ extrapolation + radiative-zone reliability.**
  HITEMP is delivered only for ρ≥0.1 kg/m³ (45 dilute points OOM/rejected in
  RADIS at wstep=0.01 — Doppler widths ~4e-5 cm⁻¹). Below the grid ρ floor, κ is
  held flat in ρ for line-list sources (`SourceSpec.dilute_flat=True`): pressure
  broadening is negligible there, so mass opacity κ = σ·N_A/M becomes
  ρ-independent (cross-section per molecule fixed). These ρ-floored calls are
  counted as `flat_rho` (physically justified) — NOT as unreliable hull clamps.
  Separately, run reliability is now judged on the **radiative-zone** clamp
  fraction, since opacity shapes the structure only there; clamps in the deep
  convective interior (T above the grid ceiling, structure on the adiabat) are
  irrelevant. Effect: composite/HITEMP fiducial runs have radiative clamp = 0.0
  (reliable) with flat_rho ≈ 0.25; pure FP stays 1.0 (unreliable), as it must.
  Flat-extrapolation region and both clamp fractions are in run metadata and
  figure captions. `StructureResult.{clamp_fraction, clamp_fraction_radiative,
  flat_rho_fraction}`.

## 2026-08-01

- **A24 — Phase-2 evolution scheme + G6.** Static-sequence cooling (Lopez &
  Fortney 2014; Ginzburg 2016) parameterized by T_int: a sequence of static
  models gives E(T_int), L(T_int), R_p(T_int); age follows from dt = dE/L,
  integrated as a Simpson cumulative integral of 1/L against E (E monotone in
  T_int). L(T_int) spans ~4 decades, so the age integral is dominated by the
  cool, low-L end; the G6 residual |∫L dt − ΔE|/ΔE converges ~2nd order in grid
  size — Simpson reaches ~2e-3 at n=60 (default), comfortably under the 1% gate
  (trapz needed n≈100). `evolution.cooling_track`; gate G6 in
  `validation/test_evolution.py`. The definitive "stay puffy longer": at fixed
  R_p the more opaque source is older — age_FP/age_POKAZATEL ≈ 15 at R_p≈1.85–1.9
  (T_eq=600), matching the fixed-entropy F3 headline via an independent method.
  Its magnitude inherits the Rosseland ±~2× line-list caveat (A22).

- **A21 — Per-comparison opacity bands (owner decision).** The common
  Rosseland/Planck band is now resolved per *set of sources being compared*, not
  globally: `resolve_band(cfg, sources)`. A pure source uses its own full band
  (`get_source(name)` → that source's band); a figure passes exactly the grid
  sources it plots, so a windowed dataset (POKAZATEL [50,4500]) only constrains
  comparisons it is part of, and the composite (HITEMP+FP) is not capped by it.
  Two consequences the owner wanted: (a) adding a source never silently changes
  previously computed numbers — each result is reproducible on its own band; (b)
  every figure and results file records the band it used (`StructureResult.
  band_cm`, figure captions, `Fx_opacity_compare.csv` band columns). Band-
  agnostic sources (constant-κ) are cached with a band=None key.

- **A22 — Fx opacity comparison: Planck reproduces the crossover, Rosseland runs
  ~2× high.** `Fx_opacity_compare` (FP vs POKAZATEL, exact, unclamped) matches
  `crossover_hitemp.json` (FP vs HITEMP) in trend everywhere and in Planck to
  ~15–20%, but its Rosseland ratios are ~2× the crossover's. Verified not a
  resolution artifact (smoothing POKAZATEL leaves κ_R unchanged) and not an FP
  error — it is a genuine POKAZATEL–HITEMP Rosseland difference and/or a
  band/method difference in the crossover pipeline (Rosseland is window-
  dominated and band-sensitive). Planck is the trustworthy headline; Rosseland
  carries a ±~2× line-list caveat until Fx is run against HITEMP to isolate it.
  See `docs/opacity_crosscheck.md`.

- **A20 — Interim analytic validation G5a (owner-authorized).** The real G5
  (Aguichine 2021 steam-world M–R) needs an external table not on the machine,
  so it skips. Interim, per the owner's request for "a defensible substitute
  from published scaling relations or an analytic limit," two *derivable* limits
  validate the structure integration with no external data (`test_g5a_interim`):
  (1) **bare-core limit** — R_p → R_core (Zeng 2016, already cited) as f_env→0,
  monotone and from above, residual = atmosphere skin; (2) **isothermal scale
  height** — radiative-zone thickness = H·ln(P_base/P_top),
  H = k_B T/(μ_H2O·g), matched to ~2% (asserted at 12%). These are analytic
  limits (like G1), NOT a literature M–R anchor; G5 proper stays open until the
  Aguichine data + a valid opacity (composite/HITEMP) are available.
  Self-consistency bonus: the f_env→0 residual (0.051 R⊕) equals the
  scale-height thickness (0.052 R⊕).

- **A23 — F3 headline is the fixed-INTERIOR-ENTROPY cooling ratio (owner
  decision, resolves A19).** The KH ratio τ_FP/τ_LL is now compared at fixed
  deep-adiabat specific entropy, not fixed T_int. `build_at_entropy` shoots on
  T_int (interior entropy is monotone in T_int) to match a target adiabat. At
  fixed entropy the two sources have essentially the *same R_p* (entropy sets
  the structure), but the more opaque FP throttles the flux to a lower T_int and
  hence lower L, giving a longer cooling time: τ_FP/τ_POKAZATEL ≈ 12 at
  s*=13500 J/kg/K (T_eq=600, f_env=5%) — the "stay puffy longer" headline. The
  fixed-T_int ratio (<1, A19) is retained as a labeled diagnostic panel; the
  distinction is written up to referee standard in `docs/cooling_time_caveat.md`.
  `StructureResult.interior_entropy_j_kg_k` records the deep entropy. The
  definitive version is the Phase-2 radius(age) track.

- **A19 — F3 cooling-time control variable (fork, awaiting decision).** The KH
  ratio τ_FP/τ_LL at *fixed T_int* comes out < 1 (FP cools faster), opposite to
  the "stay puffy longer" headline, because L_int = 4πR²σT_int⁴ is imposed by
  T_int: the puffier FP planet has larger L (∝R²) and is less bound (smaller
  |E|). Fixing T_int hides opacity throttling. The correct comparison is at
  fixed interior entropy (Phase-1-feasible; opaque planet → lower T_int/L →
  longer τ) or the Phase-2 time track. Kept as a labeled fixed-T_int diagnostic
  until the owner picks; not presented as the headline. See
  `docs/cooling_time_caveat.md`.

- **A16 — POKAZATEL tables registered (overlap comparison only).** 27
  ExoMol-POKAZATEL spectra delivered to `~/subneptune_data/pokazatel_opacity/`
  as `T{T}_rho{rho}.npz` (no `_opacity` suffix; keys `wn,kappa,T,rho,rho_kg_m3,
  p_bar,source`; correct units, no ×1000). They sit on the **same (T,ρ) grid as
  FP** (T≥1000 K, ρ≥109 kg/m³) and are **windowed to [50, 4500] cm⁻¹**, so:
  (a) registry globs for hitemp/pokazatel changed to `T*_rho*.npz`; (b) the
  auto common band (A14) correctly re-intersects to **[50, 4499.95] cm⁻¹**
  (FP∩POKAZATEL) whenever POKAZATEL is present — narrower than FP-alone's
  [1.34, 6000], which raises FP κ_R back ~35% (loses the 4500–6000 windows).
  This is SPEC §3-faithful (common band for exact differentials) but means a
  windowed dataset on disk constrains the band globally; when the wide HITEMP
  11×11 tables land we may want per-comparison bands so POKAZATEL's window does
  not cap the composite — FLAGGED for decision then. POKAZATEL, like FP, does
  **not** cover the radiative zone, so POKAZATEL *structure* runs are also
  ~100% clamped; the meaningful POKAZATEL-vs-FP comparison today is at the
  opacity level on the shared grid (figure `Fx_opacity_compare`).

- **A17 — Delivered POKAZATEL confirms ideal-gas broadening pressure.** The
  delivered POKAZATEL `p_bar` equals ρRT/M to <0.1% at every grid point (e.g.
  598.8 bar at T=1000, ρ=129.744), confirming the production used per-cell
  ideal-gas pressure — matching `generate_linelist_tables.py` default
  `PRESSURE_MODE="ideal"`. (The `HitempCalculator` snippet uses a fixed 1-bar
  reference; the actual POKAZATEL production did not, so the generator follows
  the delivered data, not the snippet's p_ref.)

- **A18 — numpy 2.0 compatibility.** The environment numpy was upgraded 1.26→2.4
  between sessions; `np.trapz` was removed (renamed `trapezoid`). Integrations
  now import `trapezoid` with a `trapz` fallback, so the code runs on numpy
  ≥1.22 and ≥2.0 (requirement stays numpy≥1.26).

## 2026-07-31

- **A15 — Default composite blend window moved to ρ = 100–300 kg/m³ (owner
  decision).** The earlier default (30–100, A12) sat below FP's density floor
  (~109 kg/m³), so FP was extrapolated on the dilute side of the window. The
  default is now **100–300**, where FP has real data; the line list carries the
  full dilute side, and the seam (√(100·300) ≈ 173 kg/m³) sits on real FP data.
  ρ = 30–100 is retained as a **documented sensitivity case**, not the default.
  `config.opacity.composite_blend_rho`.

## 2026-07-29

- **A10 — Guillot (2010) boundary conventions.** Outer boundary uses Guillot
  (2010) eq. 49 (collimated two-stream) with planet-average incidence
  μ⋆ = 1/√3 (Eddington isotropic equivalent). The config equilibrium
  temperature T_eq is interpreted as the full-redistribution value, so
  T_irr = 4^(1/4)·T_eq = √2·T_eq (σT_eq⁴ = ¼σT_irr⁴). γ = κ_vis/κ_th is the
  configurable knob (fiducial 1.0 + sensitivity scan, SPEC §10). Behaviour
  verified: γ<1 greenhouse (deep hotter), γ=1 isothermal irradiation term,
  γ>1 inversion; mirror symmetry γ↔1/γ. See `structure/boundary.py`, gate G1.

- **A11 — No FP grid extension; composite opacity instead (owner decision).**
  FP is a condensed-phase method and is invalid in the dilute regime, where
  line lists are the physically correct opacity. Rather than extend FP below
  its density floor (~109 kg/m³), the code uses a *composite* source: a
  line-list source (HITEMP now, POKAZATEL later) for the dilute-to-moderate
  regime (incl. the 500–1000 K radiative zone), the existing FP source for the
  dense regime, blended in log κ over a density window. Diagnosis that forced
  this: for the SPEC §8 fiducials the opacity-sensitive radiative zone sits at
  T≈600–1030 K, ρ≈2×10⁻⁶–540 kg/m³ — almost entirely outside the FP grid
  (T∈[1000,3000], ρ∈[109,2229]); pure-FP structure runs clamp ~100% of steps.

- **A12 — Composite blending.** κ_R and κ_P are blended in log space with a
  smooth weight that runs from the line-list source (dilute) to FP (dense)
  across a configurable density window (default ρ = 30–100 kg/m³; seam =
  geometric-mean centre). Seam location and width are sensitivity-scan
  parameters. The FP/line-list κ ratio at the seam is logged per run as a
  quantified systematic. Pure `fp`/`hitemp`/`pokazatel` sources remain
  selectable for comparison figures. (The default window was updated to
  100–300 in A15 so FP stays above its data floor; 30–100 is now a sensitivity
  case. The per-run seam ratio and clamp fraction expose any FP extrapolation,
  and the window is scannable.)

- **A13 — Adaptive matching just above the RCB (owner decision).** SPEC §5's
  fixed matching pressure (default 10 bar) fails when the radiative–convective
  boundary is shallower than P_match (routine here, and catastrophic for high
  FP opacity: τ~10⁶ at 10 bar → unphysical ~5000 K → collapse). The Guillot
  outer integration now carries a Schwarzschild instability event
  (grad_rad − grad_ad = 0) and matches at min(P_match_config, P_RCB): the
  analytic radiative profile is only used where the atmosphere is actually
  radiative. The configured match pressure remains the *deep* cap and is still
  swept 1–100 bar for sensitivity.

- **A14 — Common band widened + auto-detected.** The Rosseland/Planck common
  band is now the largest range all *available* sources support: [max over
  sources of min(wn), min over sources of max(wn)], verified per file. With FP
  alone this is ≈[1, 6000] cm⁻¹ (FP spectra reach 6000); it is re-intersected
  when line-list spectra land. Mean tables are recomputed on this band. The
  residual truncation of the *true* Rosseland mean (harmonic weighting is
  dominated by transparency windows beyond the band, biasing κ_R high) remains
  a stated systematic — to be quantified against the line-list source once its
  full-band (>6000 cm⁻¹) spectra are available. See `docs/band_limitation.md`.

- **A8 — ρ(P,T) inversion through the vapor–liquid dome.** At fixed T below the
  critical point, the AQUA table has exact P(ρ) plateaus (constant saturation
  pressure across the coexistence density range; up to ~457 flat steps per
  column near 300 K, phase id 1), and small plateaus at high density (ice
  VII/X, superionic; ids −7, −10). P(ρ,T) is thus non-strictly-monotonic, so
  ρ(P,T) is multivalued on the dome. For the inversion the per-column log P
  array is made *strictly* increasing by `maximum.accumulate` plus a cumulative
  ~1e-12 dex nudge applied **only to flat/decreasing entries** — single-phase
  nodes are left unchanged (≤1e-9 relative P shift, far under the G2 0.5%
  tolerance), so single-phase round-trips stay exact. Effect: for P below a
  plateau the vapor branch is returned, for P above it the liquid branch — the
  physical inward-integration path; exactly on P_sat (measure zero) a dome value
  is returned. A `two_phase` diagnostic (phase id) is exposed so structure runs
  report any dome traversal. G2 spot points are taken in single-phase regions.
  The `phase` diagnostic is a nearest-node lookup and is approximate within ~1
  grid cell of a transition; ρ, s, grad_ad remain accurate there. See
  `docs/eos_inversion.md`.

- **A9 — Core mass fraction default 0.325 (Earth-like).** SPEC §4 cites the
  Zeng (2016) "Earth-like rock" M–R relation. Zeng's Earth-like composition is
  CMF = 0.325 (32.5% Fe + 67.5% MgSiO₃), used as the default in the power law
  R/R⊕ = (1.07 − 0.21·CMF)(M/M⊕)^0.27. CMF is a parameter, so pure-rock
  (CMF = 0) or iron-rich cores need no code change. The analytic power law is
  used rather than the tabulated Zeng (2019) grid; valid to a few percent for
  M ≲ 8 M⊕ (fiducial M_c = 5 M⊕ is inside).

- **A1 — Config parsing dependency.** SPEC §1 says "numpy, scipy, matplotlib
  only" but also specifies YAML run configs. PyYAML is added as a *config-only*
  dependency (not a compute dependency); no physics depends on it. Alternative
  considered: stdlib `tomllib`, rejected because SPEC explicitly says YAML.

- **A2 — AQUA table is (ρ,T)-indexed, not (P,T).** SPEC §2/§4 request
  `rho(P,T)`, `s(P,T)`, `grad_ad(P,T)`. The only AQUA table present is
  `aqua_eos_rhot_v1_0.dat`, tabulated on a (ρ,T) grid (1501×301) with P as an
  output column. `rho(P,T)` is therefore obtained by inverting the tabulated
  `P(ρ,T)` along `T = const` (P is monotone increasing in ρ at fixed T for this
  EOS). Inversion scheme and its error will be documented in
  `docs/eos_inversion.md`. The companion AQUA P-T table is not on disk.

- **A3 — `ad_grad` is tabulated directly.** SPEC §2 allows deriving `grad_ad`
  from `s` by finite difference "if not tabulated directly." It *is* tabulated
  (`ad_grad` column, `(dlogT/dlogP)_S`), so it is used directly; no
  finite-difference derivation is needed. `s` and `u` are also tabulated.

- **A4 — FP opacity factor-1000 confirmed.** Loading `T*_opacity.npz` and
  band-integrating `kappa` *raw* gives κ_P ≈ 3×10⁵–4×10⁶ cm²/g; after the
  mandated ×1e-3 scaling κ_P lands at 323–4170 cm²/g, inside SPEC §2's asserted
  1e1–1e4 window. The scaling is applied on load and asserted (hard-fail
  otherwise). FP files also carry `alpha`(=κ·ρ) and precomputed `planck`/
  `rosseland` scalars (in *raw* units); band-limited means are (re)computed
  from spectra on the common band rather than trusting these scalars.

- **A5 — HITEMP/POKAZATEL per-point spectra absent → "wait for real spectra"
  (owner decision, 2026-07-28).** SPEC §2 lists HITEMP and POKAZATEL per-point
  spectra npz as inputs; neither is on disk. Present instead: FP spectra +
  `crossover_hitemp.json` containing *precomputed* HITEMP↔FP Planck/Rosseland
  ratios and M1–M4 crossover metrics per grid point. Owner chose to **wait for
  real spectra** rather than reconstruct HITEMP means as `ratio × FP_mean`,
  because that would break SPEC §3's central guarantee (all sources integrated
  on the *identical* common band "so the differential conclusions are exact") —
  the RADIS ratios were computed over RADIS's own range/method. Consequences:
  the `fp` source is built now; `hitemp`/`pokazatel` are registered but raise a
  clear "spectra not delivered" error until their npz arrive (ingesting them
  then requires only a registry entry, zero code change — SPEC §3);
  `crossover_hitemp.json` is used only as an independent cross-check; and the
  F1/F3 differential figures are deferred until the line-list spectra are
  delivered. (F2, which is per-source, can still run for FP.)

- **A6 — FP opacity grid is scattered, not rectangular; "bilinear" realized as
  linear interpolation on a Delaunay triangulation.** SPEC §3 prescribes
  "bilinear in (T, log rho) on log kappa." The 27 FP points do not form a
  tensor-product grid (e.g. T=1000 K has a single ρ point; ρ sampling differs
  per T), so literal bilinear is undefined. Faithful generalization used:
  `scipy.interpolate.LinearNDInterpolator` on points (T, log₁₀ρ) interpolating
  log₁₀κ (piecewise-linear on triangles = bilinear's scattered-data analogue).
  Extrapolation outside the convex hull is the SPEC §3 gate: clamp to nearest
  grid value (`NearestNDInterpolator`), count every clamped call, flag the run
  unreliable if >20% of steps clamp. See `docs/opacity_interpolation.md`.

- **A7 — Common band fixed at the configured [50, 4500] cm⁻¹, verified
  covered.** SPEC §3 allows "the largest band all sources share — detect and
  record." With only FP present (native [1.1, 6000] cm⁻¹ ⊃ [50, 4500]), the
  configured band is used directly and each source's spectra are asserted to
  span it before its means are computed, guaranteeing an identical band across
  sources when line-list spectra arrive. FP full-band [1.1,6000] vs common-band
  [50,4500] truncation is quantified in `docs/band_limitation.md` as a proxy
  until line-list full-band means are available.
