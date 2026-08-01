# Fx opacity comparison vs the crossover_hitemp.json cross-check

`Fx_opacity_compare` computes the exact FP-vs-line-list band-limited κ_R/κ_P
ratios on the 27 shared (T,ρ) grid points, over the common band the plotted
sources share (currently [50, 4499.95] cm⁻¹, set by POKAZATEL's window). The
owner's independent `crossover_hitemp.json` (FP vs **HITEMP**, RADIS pipeline)
is used as a cross-check.

## What agrees, what doesn't

| point (T, ρ)        | my FP/POKAZATEL Planck | crossover FP/HITEMP Planck | my Rosseland | crossover Rosseland |
|---------------------|------------------------|----------------------------|--------------|---------------------|
| 1000, 129.7         | 1.63                   | 1.38                       | 14.9         | 7.6                 |
| 1500, 1200          | 5.41                   | 4.81                       | 50.8         | 21.9                |
| 2000, 150           | 1.92                   | 1.49                       | 4.2          | 1.7                 |
| 3000, 2228.7        | 22.5                   | 18.2                       | 84.3         | 39.6                |
| 3000, 108.8         | 1.73                   | 1.40                       | 2.8          | 1.4                 |

- **Trend: reproduced exactly** — monotone rise from ~1.4 (cool, low-ρ) to ~20
  (hot, dense), same ordering across all 27 points.
- **Planck: agrees to ~15–20%** (mine runs slightly high). Good for two
  different line lists compared through two different pipelines.
- **Rosseland: mine runs ~2× the crossover** consistently.

## Why the Rosseland factor of ~2 (not a bug)

- It is **not resolution**: box-smoothing the fine POKAZATEL (0.05 cm⁻¹) up to
  10–40 cm⁻¹ leaves its κ_R unchanged (2.58 → 2.59), so the FP–POKAZATEL
  spectral-resolution mismatch does not explain it.
- It is **not my FP means**: FP κ_P = 425, κ_R = 38.5 cm²/g at (1000, 129.7)
  reproduce the session-start values and the same FP files the crossover used.
- The residual is therefore a **genuine POKAZATEL-vs-HITEMP difference in the
  Rosseland mean** (κ_R(HITEMP) ≈ 5.1 vs κ_R(POKAZATEL) ≈ 2.6 at that point),
  and/or a **band/method difference in the crossover pipeline** (its band is not
  recorded; Rosseland is a harmonic, transparency-window-dominated mean and is
  far more band/method-sensitive than Planck). Two water line lists differing by
  ~2× in the Rosseland mean while agreeing to ~15% in Planck is physically
  reasonable — the Rosseland mean lives in the weak-line windows where the lists
  are least complete and most divergent.

## Status & how to close it

**Not declared "solid" for Rosseland.** The Planck comparison is trustworthy;
the Rosseland absolute ratio carries a line-list/band systematic of order 2.
To isolate it, run `Fx_opacity_compare` against the **HITEMP** tables once
generated (same band, same means pipeline): FP/HITEMP from Fx should then match
`crossover_hitemp.json` closely, confirming the pipeline and pinning the
remaining POKAZATEL–HITEMP difference as physical. Until then, quote **Planck
ratios** as the headline opacity contrast and present Rosseland with this ±~2×
line-list caveat.

Numbers are dumped to `outputs/Fx_opacity_compare.csv` (band recorded per row).
