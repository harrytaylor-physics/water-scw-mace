# AMENDMENTS

Every design decision made without a direct citation, or that resolves a
gap/ambiguity in `SPEC.md`, is logged here with a one-line rationale
(SPEC §10). Dated; newest first.

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
