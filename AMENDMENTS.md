# AMENDMENTS

Every design decision made without a direct citation, or that resolves a
gap/ambiguity in `SPEC.md`, is logged here with a one-line rationale
(SPEC §10). Dated; newest first.

## 2026-09-25 (Extended Data; the Fig. 4 categories; what the ARC pull does and does not contain)

- **A67 — Extended Data: 8 figures and 2 tables, built from what the repository
  and the ARC pull actually hold.** `figures/ed_items.py`, `figures/ed_tables.py`,
  `figures/ed_index.py` → `figs/ed/`. Same `save_pdf` guard as the main figures, so
  the vocabulary rule, the 7 pt floor and the exact width are enforced here too —
  it caught four 6.9 pt labels of mine in ED Fig. 5 before they reached a PDF.

  **What the ARC pull contains is mostly code, not data.** 92 files: scripts,
  READMEs, campaign documentation and 33 PNGs, and exactly **four** data files —
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
  landed, and `tools/grid_extension_effect.py` refuses to run on a partial set.

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
