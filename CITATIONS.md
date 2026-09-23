# CITATIONS

Primary references for physics choices. Every equation in the code traces to
one of these or to an explicit derivation in `docs/`.

## Equations of state
- **AQUA** — Haldemann, Alibert, Mordasini & Benz (2020), *A&A* 643, A105.
  H₂O EOS spanning vapor/liquid/supercritical/ice with consistent
  thermodynamics. Used everywhere in the envelope (SPEC §4).

## Core mass–radius
- **Zeng, Sasselov & Jacobsen (2016)**, *ApJ* 819, 127; **Zeng et al. (2019)**,
  *PNAS* 116, 9723. Earth-like rock mass–radius relation for the inert core
  (SPEC §4).

## Opacity
- **HITEMP** — Rothman et al. (2010), *JQSRT* 111, 2139. Line-list water opacity.
- **ExoMol POKAZATEL** — Polyansky et al. (2018), *MNRAS* 480, 2597. Line-list
  water opacity.
- First-principles (FP) opacity — group data (see `docs/` for provenance).
- **MT_CKD water-vapour continuum** — Mlawer, Turner, Masiello, Alvarado et al.
  (2023), *JGR Atmospheres* 128, e2023JD038473 (the v4.0+ formulation: Eqs. 1-4,
  the (T_ref/T)^n self-continuum temperature dependence, the reference-density
  scaling and the radiation term); Mlawer et al. (2012), *Phil. Trans. R. Soc. A*
  370, 2520 (AER's own general reference, and the one their data file cites).
  Implementation is **release 4.3** of AER's distributed model and coefficients
  (<https://github.com/AER-RC/MT_CKD>), used for the `hitemp_mtckd` source
  (A58). **Self-continuum only** — the foreign term carries (1 − x_H2O) and the
  envelope is pure steam — and **extrapolated far outside its fitted range** in
  both T and density; see `opacity/mtckd.py` and AMENDMENTS A58.
- **RADIS** — Pannier & Laux (2019), *JQSRT* 222, 12. Line-by-line code used to
  generate the HITEMP spectra and the A59 broadening-convention study.

## Radiative structure & boundary
- **Guillot (2010)**, *A&A* 520, A27. Analytic irradiated two-stream T(τ)
  profile for the outer boundary (SPEC §5).
- **Kippenhahn, Weigert & Weiss (2012)**, *Stellar Structure and Evolution*.
  Radiative gradient, Schwarzschild criterion.
- **Pierrehumbert (2010)**, *Principles of Planetary Climate*, ch. 4–5.

## Evolution (phase 2)
- **Lopez & Fortney (2014)**, *ApJ* 792, 1. Static-sequence thermal evolution.
- **Ginzburg, Schlichting & Sari (2016)**, *ApJ* 825, 29. Analytic cross-checks.

## Literature anchor (G5)
- **Aguichine, Mousis, Deleuil & Marcq (2021)**, *ApJ* 914, 84. Steam-world
  mass–radius for the static structure benchmark.

## Constants
- CODATA via `scipy.constants` — see `docs/constants.md`.
