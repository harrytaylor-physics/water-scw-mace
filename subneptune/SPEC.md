# SPEC.md — `subneptune`: steam-envelope sub-Neptune structure & evolution

Permanent research code (GitHub-grade), not a throwaway. Built to answer,
first, one question for a slide deck: **how does replacing line-list water
opacity (HITEMP / ExoMol-POKAZATEL) with first-principles (FP) opacity
change the interior structure and cooling of steam-envelope sub-Neptunes?**
— and then to remain the group's steam-world model afterwards.

Owner: Harry Taylor (Oxford). Reviewer-of-record for physics choices:
assume Raymond Pierrehumbert will read the outputs; every equation must be
traceable to a citation or an explicit derivation in docs/.

## 0. Scope

v1 = pure-steam (H2O) envelope on an inert rocky core. No H/He, no
mixtures, no tidal/radiogenic heating, no atmospheric escape. These are
documented extensions (section 9), not silent omissions.

Phase 1 (deck-critical, build first): opacity tables -> static structure
-> T(P)/RCB comparison across opacity sources + Kelvin-Helmholtz cooling
timescale ratio.
Phase 2: thermal evolution -> radius(age) tracks per opacity source.
Phase 3+: extensions (section 9).

## 1. Repository layout

    subneptune/
      README.md  SPEC.md  CITATIONS.md  AMENDMENTS.md
      config/           yaml run configs (fiducial planet + grids)
      data/             NOT committed; symlink or path-config to
                        ~/subneptune_data (AQUA tables, fp_opacity/)
      eos/              AQUA loader + thermodynamic interface
      opacity/          table builder + interpolators (FP, HITEMP, POKAZATEL)
      structure/        hydrostatic integrator, boundary conditions
      evolution/        thermal evolution driver (phase 2)
      validation/       gates (pytest) + literature benchmarks
      figures/          plotting scripts (plain style; each fig one script)
      docs/             derivations + assumption ledger

Python >= 3.11; numpy, scipy, matplotlib only. Type hints throughout.
pytest for every gate. No notebooks in the repo. Every physics constant
from scipy.constants or CODATA, cited once in docs/constants.md.

## 2. Data inputs (already on this machine)

- `~/subneptune_data/AQUA/` — AQUA EOS tables (Haldemann, Alibert, Mordasini
  & Benz 2020, A&A 643, A105). Use the tabulated rho(P,T), s(P,T), and
  adiabatic gradient; if grad_ad is not tabulated directly, derive from s
  and document the finite-difference scheme.
- `~/subneptune_data/fp_opacity/T*_opacity.npz` — 27 first-principles
  opacity spectra on the (T, rho) grid, keys: wn [cm^-1], alpha [cm^-1],
  kappa [cm^2/g], T [K], rho_kg_m3.
  **CRITICAL: these carry a legacy factor-1000 error — multiply kappa and
  alpha by 1e-3 on load. Assert plausibility after scaling (band-integrated
  kappa_P must land in ~1e1–1e4 cm^2/g across the grid) and hard-fail
  otherwise.**
- `~/subneptune_data/fp_opacity/crossover_hitemp.json` — per-point HITEMP
  Planck means / ratios (schema: inspect and adapt; keys include T, rho,
  and a planck ratio or means).
- HITEMP and (when delivered) POKAZATEL per-point spectra npz — same key
  convention as FP, correct units (no 1e-3).

## 3. Opacity module

One interface, N sources: `OpacitySource.kappa_R(T, rho)` and
`.kappa_P(T, rho)`.

- `build_tables.py`: from each source's per-point spectra, compute
  band-limited Rosseland and Planck means over the COMMON band
  [50, 4500] cm^-1 (or the largest band all sources share — detect and
  record). Same band for every source so the *differential* conclusions
  are exact; absolute band-truncation is a stated limitation (see docs/
  band_limitation.md: at 1000–3000 K the Planck function peaks at
  2900/T um, i.e. inside the band, but Rosseland weighting has tails —
  quantify the truncation by comparing full-band vs common-band means for
  the line-list source, where full-band is available).
- Interpolation: bilinear in (T, log rho) on log kappa. Extrapolation
  policy is a GATE, not a guess: clamp to nearest grid value, count and
  log every clamped call, and report the fraction of structure-integration
  steps that touched extrapolation. > 20% of steps clamped => the run is
  flagged unreliable in its output metadata.
- Sources registered by name in config: `fp`, `hitemp`, `pokazatel`.
  Adding one must require zero code changes outside a registry entry.

## 4. EOS module

AQUA everywhere in the envelope (it spans vapor/liquid/supercritical/ice
phases with consistent thermodynamics — the reason it is the v1 choice).
Interface: rho(P,T), grad_ad(P,T) [= (dlnT/dlnP)_S], s(P,T), plus
sound-check helpers. Validation against the paper's own published values
at spot points (gate G2). Design the interface so a future custom EOS
(e.g. the group's MLIP-derived surface) drops in beside AQUA.

Core: inert rock, mass M_c, radius from the Zeng, Sasselov & Jacobsen
(2016) / Zeng et al. (2019) Earth-like rock mass-radius relation
(tabulate or use their fitted power law; cite which). No core luminosity
in v1 (documented).

## 5. Structure module (phase 1 core)

Inputs: M_c, M_env (or f_env = M_env/M_p), T_int, T_eq (irradiation),
opacity source name. Output: full profile r(m), P(m), T(m), rho(m),
diagnostics (R_p, RCB location(s), photospheric quantities), serialized
npz + json summary.

Equations (interior): hydrostatic equilibrium dP/dm = -Gm/(4 pi r^4);
mass continuity dr/dm = 1/(4 pi r^2 rho); temperature gradient
grad = min(grad_ad, grad_rad) (Schwarzschild criterion), with
grad_rad = 3 kappa_R L P / (16 pi G m sigma T^4)  [Kippenhahn, Weigert &
Weiss 2012; Pierrehumbert 2010 ch. 4-5]. L = L_int constant through the
envelope in v1 (no internal sources), L_int = 4 pi R_p^2 sigma T_int^4.

Outer boundary: Guillot (2010, A&A 520, A27) analytic irradiated
two-stream T(tau) profile, using kappa_th = kappa_R(T, rho) from the
selected source and a configurable visible/thermal ratio gamma
(default: document a fiducial choice and its sensitivity; gamma is a
knob Ray will ask about). Match interior integration to the Guillot
profile at a configurable pressure (default 10 bar; sensitivity-check
1-100 bar and record in docs/boundary_sensitivity.md).

Numerics: integrate inward from the matching point with scipy solve_ivp
(LSODA or Radau), shooting on R_p until r(M_c) = R_c(M_c) to tolerance;
document the alternative (relaxation) as future work. All tolerances in
config, not hard-coded.

## 6. Evolution module (phase 2)

Standard static-sequence thermal evolution [Lopez & Fortney 2014, ApJ
792, 1; Ginzburg, Schlichting & Sari 2016, ApJ 825, 29 for analytic
cross-checks]: parameterize the envelope by internal entropy (or T_int),
E(sequence) from integrated internal + gravitational energy, advance
dE/dt = -L_int(t) (+ constant irradiation held in T_eq; no escape v1).
Output: R_p(t), L(t), RCB(t) per opacity source, and the headline
comparison: cooling-time ratio between opacity sources at fixed (M_c,
f_env), and inferred-f_env shift at fixed (R_p, age).

## 7. Validation gates (pytest; all must pass before any science figure)

G1  Isothermal-limit recovery: kappa -> 0 => radiative zone approaches
    the Guillot skin solution analytically; assert to tolerance.
G2  EOS spot checks: AQUA loader reproduces >= 10 published table values
    to <= 0.5%.
G3  Hydrostatic residual: post-hoc |dP/dm + Gm/(4 pi r^4)| integrated
    error < 1e-6 relative, every converged model.
G4  Opacity-swap null test: two sources pointing at identical tables =>
    identical structure to machine precision.
G5  Literature anchor: static M-R for steam worlds vs Aguichine, Mousis,
    Deleuil & Marcq (2021, ApJ 914, 84) at matched (M, f_env, T_eq):
    agree to ~5-10% in R_p; document residuals and causes (band-limited
    opacity, boundary choice) rather than tuning to match.
G6  (phase 2) Energy bookkeeping: |Delta E - integral L dt| / E < 1%
    over the full track.
G7  Convergence: R_p change < 0.1% under halved tolerances and doubled
    table resolution, for the fiducial planet.

## 8. Phase-1 deliverables (the deck figures)

Fiducial planet: M_c = 5 M_earth, f_env = 5% (steam), T_eq = 500-800 K
(config grid), T_int grid 30-300 K. For each opacity source:
F1  T(P) profiles overlaid (FP vs HITEMP vs POKAZATEL), RCB marked.
F2  RCB pressure and R_p vs T_int per source.
F3  Kelvin-Helmholtz cooling-time ratio tau_FP / tau_LL vs (f_env, T_eq)
    — the "planets stay puffy longer" quantification.
Figures: plain matplotlib, colorblind-safe if colored, every figure
produced by one script in figures/ taking only a config path.

## 9. Documented extensions (create stubs + docs, do not implement)

H/He-H2O mixed envelopes (ties to the group's miscibility programme);
custom-functional EOS and opacity swap (LDA / hybrid surfaces exist);
atmospheric escape; radiogenic/tidal heating; full non-gray radiative
transfer replacing Guillot.

## 10. Working agreements

- No silent assumptions: anything chosen without a citation goes in
  AMENDMENTS.md with a one-line rationale.
- Every module lands with its gate tests in the same commit.
- Commit granularity: one module or one gate per commit, imperative
  messages.
- When a design decision materially affects Phase-1 figures (band choice,
  gamma, matching pressure), implement the default AND the sensitivity
  scan; never just the default.
