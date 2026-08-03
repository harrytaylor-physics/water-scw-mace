# Does FP opacity affect these steam worlds? (Headline caveat)

## Finding

For the SPEC §8 fiducial class (M_c = 5 M⊕, f_env = 1–10% steam, T_eq 500–800 K,
T_int 30–300 K), **FP opacity does not change the interior structure or cooling
at all.** The composite (FP-dense + line-list-dilute) is identical to the pure
line list:

| T_eq | T_int | R_p (composite) | R_p (HITEMP) | max radiative-zone ρ [kg/m³] |
|------|-------|-----------------|--------------|------------------------------|
| 500  | 50    | 1.759           | 1.759        | 6.70 |
| 500  | 150   | 1.849           | 1.849        | 0.04 |
| 500  | 300   | 1.917           | 1.917        | 0.00 |
| 800  | 50    | 1.815           | 1.815        | 2.98 |
| 800  | 150   | 1.888           | 1.888        | 0.09 |
| 800  | 300   | 1.958           | 1.958        | 0.01 |

## Why

Opacity shapes the structure only in the **radiative zone**; the convective
interior follows the adiabat regardless of κ. For these planets the radiative
zone is the outer atmosphere, at **ρ ≲ 7 kg/m³** (down to ~10⁻² at high T_int).
FP is a condensed-phase method with a grid floor at **ρ = 109 kg/m³** — 1.5 to 4
orders of magnitude denser than anywhere opacity matters here. So the composite
uses the line list throughout the radiative zone, and FP enters only the deep
convective interior (ρ > 100), where it does not affect the structure.

## The earlier "~15× stay puffy longer" was an extrapolation artifact

Pure-FP runs give age_FP/age_LL ≈ 15–28 at fixed R_p, but that model
nearest-neighbor-clamps FP from its coldest/least-dense grid point
(T=1000 K, ρ=109) down into the (T≈650 K, ρ≈0.01) radiative zone — physically
unjustified, and now flagged unreliable (radiative clamp fraction = 1.0, A26).
It must **not** be the deck headline.

## What FP vs line-list opacity *does* say (exact, reliable)

At the overlap where FP is valid (dense, hot; ρ ≥ 109, T ≥ 1000), FP is **2–3×
more opaque** than HITEMP/POKAZATEL in the Rosseland mean on the [100, 4000]
band (`Fx_opacity_compare`, reproduces the crossover to ~1%). This is a real,
quantified opacity difference — it simply does not propagate to the structure of
*this* planet class, because that (dense, hot) regime is convective here. It
would matter for denser/hotter envelopes whose RCB sits at ρ ≳ 100 (more massive
planets, magma-ocean-coupled layers) — a documented extension, not v1.

## The real opacity systematic for these planets: line-list choice

The dominant, *reliable* opacity uncertainty is between line lists, not FP vs
line list. HITEMP vs POKAZATEL: κ_R differs by ~1.5× (POKAZATEL less opaque),
giving a **cooling-age ratio ≈ 2×** at fixed R_p. Plus the band systematic
(~2× on κ_R between [100,4000] and [50,4500]).

## Deck-ready statement

- FP water opacity is 2–3× higher than line-list at the overlap (dense, hot) —
  exact, crossover-validated.
- For steam sub-Neptunes, opacity sets the structure in the dilute (ρ ≲ 10),
  cool (600–1000 K) radiative zone — outside FP's condensed-phase validity — so
  **FP does not change their radius or cooling** (composite ≡ line list here).
- The reliable opacity uncertainty is the **line-list choice (~2× in cooling
  age)**, not FP vs line list.
- The "~15× puffier with FP" number is an artifact of extrapolating FP beyond
  its grid and is excluded.
