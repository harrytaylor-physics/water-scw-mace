# Band truncation of opacity means

Rosseland/Planck means are integrated over a **common band shared by every
source**, so the *differential* FP-vs-line-list conclusions are exact (the band
is identical for both). The price is an **absolute** band-truncation error,
quantified here. It is a stated systematic, not a bug.

## Common band (auto-detected, A14)

The band is the largest range all *available* sources support:

    band = [ max over all spectra files of min(wn),
             min over all spectra files of max(wn) ]

verified per file (`opacity.detect_common_band`). With only FP present this is
**≈ [1.34, 6000] cm⁻¹** (FP spectra span [1.11, 6000.34]; the min-of-maxes is
6000.11). It is re-intersected automatically when line-list spectra land — e.g.
if HITEMP starts at 10 cm⁻¹ the low edge moves to 10. Widening from the earlier
[50, 4500] lowered the FP Rosseland means by ~35% (median), as the 4500–6000
cm⁻¹ transparency windows enter the harmonic mean.

## Residual truncation beyond the band (the real systematic)

The band's **upper edge is capped at 6000 cm⁻¹ by FP** (a condensed-phase
dataset that cannot be extended). But the Planck function's peak (per unit
wavenumber) marches to higher wavenumber with temperature — wn_peak ≈ 2.82
k_BT/hc, which is **~5870 cm⁻¹ already at 3000 K** — so substantial radiative
weight sits beyond the band at high T:

| T (K) | Planck energy fraction in [1.34, 6000] (vs a 15000 cm⁻¹ reference) |
|-------|-------------------------------------------------------------------|
| 1000  | 0.975 |
| 1500  | 0.839 |
| 2000  | 0.656 |
| 2500  | 0.502 |
| 3000  | 0.393 |

This matters most for the **Rosseland** mean. Rosseland is a *harmonic* mean
weighted by dB/dT and is dominated by the most transparent frequencies —
the near-IR/optical **transparency windows** between water bands, largely at
wn > 6000 cm⁻¹. Excluding them **biases κ_R high**, which in turn makes
atmospheres artificially optically thick (a direct cause of the very high FP
κ_R and the need for adaptive matching, A13).

## To quantify (needs line-list full-band spectra)

The line-list source can be computed to much higher wavenumber than FP. When
HITEMP spectra are delivered (ideally to ~15000–20000 cm⁻¹), quantify the
systematic by comparing, for the line-list source only:

- κ_R and κ_P on the common band [1.34, 6000] cm⁻¹, vs
- κ_R and κ_P on the line-list's full band,

and record the ratio here per (T, ρ). This isolates the absolute bias without
touching the FP-vs-line-list differential (which stays common-band and exact).

Reproduce band figures: `compute_point_means("fp", cfg, band)` for any band;
`detect_common_band(cfg)` for the auto band.
