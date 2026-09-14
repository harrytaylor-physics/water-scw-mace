# AMENDMENTS

Every design decision made without a direct citation, or that resolves a
gap/ambiguity in `SPEC.md`, is logged here with a one-line rationale
(SPEC §10). Dated; newest first.

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
