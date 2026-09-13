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

## Quantified with the delivered HITEMP spectra (A39)

HITEMP was delivered on [10, 11000] cm⁻¹ — wider than FP but short of the
15000–20000 cm⁻¹ originally hoped for — at 76 of the 77 (T, ρ) grid points
(T=500, ρ=0.1 is missing). Comparing its means on that **full** band with the
**science** band [100, 4000] cm⁻¹ (A25) isolates the absolute bias without
touching the FP-vs-line-list differential, which stays common-band and exact
(A21):

    python opacity/build_tables.py config/fiducial.yaml --source hitemp --band full
    python figures/plot_s3_band_truncation.py config/fiducial.yaml

κ_R(full) / κ_R(common), per (T, ρ) — **< 1 everywhere, i.e. the science band
biases κ_R high**, exactly as predicted above:

| T (K) | ρ=0.1 | 1 | 10 | 30 | 100 | 300 | 1000 |
|-------|------|------|------|------|------|------|------|
| 500  | —    | 0.98 | 0.99 | 0.99 | 0.99 | 0.99 | 0.99 |
| 700  | 0.90 | 0.89 | 0.85 | 0.84 | 0.83 | 0.83 | 0.82 |
| 1000 | 0.59 | 0.43 | 0.36 | 0.35 | 0.34 | 0.33 | 0.32 |
| 1500 | 0.12 | 0.10 | 0.09 | 0.08 | 0.08 | 0.08 | 0.08 |
| 2000 | 0.05 | 0.04 | 0.04 | 0.04 | 0.04 | 0.04 | 0.04 |
| 3000 | 0.02 | 0.02 | 0.02 | 0.02 | 0.02 | 0.02 | 0.02 |

(Full table, including κ_P and the floor-free variant, in
`outputs/S3_band_truncation.csv`.) Median 0.343 over the grid; the dependence is
almost entirely on T, as the Planck peak marches out of the band. The Planck
mean is far less affected (median 0.850, min 0.273) — the harmonic Rosseland
mean is the window-sensitive one, which is why Planck carries the headline
(A22).

**Effect on the G5 anchor.** Rebuilding the Aguichine (2021) comparison with
full-band HITEMP tables moves the median |ΔR_p/R_p| from 0.165 to 0.172
(`outputs/g5_residuals_fullband.json`): the less opaque full-band tables give
slightly *smaller* radii, and ours were already low. **Band truncation therefore
does not explain the G5 residual** — the pure-steam-envelope vs full-water-phase
model difference (A32) remains the dominant cause.

Reproduce band figures: `compute_point_means(name, cfg, band)` for any band;
`detect_common_band(cfg, [name])` for a source's own full band.
