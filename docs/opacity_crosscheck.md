# Fx opacity comparison — cross-check vs crossover_hitemp.json (RESOLVED)

`Fx_opacity_compare` computes the exact band-limited κ_R/κ_P ratios FP vs line
list on shared (T,ρ) points. With HITEMP delivered, we ran FP-vs-HITEMP through
our own pipeline against the owner's independent `crossover_hitemp.json`
(FP-vs-HITEMP, RADIS) to close the ~2× Rosseland discrepancy flagged earlier
(A22). **It is closed: the gap was the band convention, not the floor.**

## The resolution

Under the owner's convention — band **[100, 4000] cm⁻¹**, relative floor
1e-6 — our FP/HITEMP means reproduce the crossover to ~1%:

| point (T, ρ)  | FP/HT κ_R (ours) | crossover κ_R | FP/HT κ_P (ours) | crossover κ_P |
|---------------|------------------|---------------|------------------|---------------|
| 1500, 300     | 3.43             | 3.42          | 1.81             | 1.80          |
| 2000, 300     | 2.11             | 2.11          | 1.55             | 1.55          |
| 2500, 300     | 1.52             | 1.52          | 1.55             | 1.54          |

## Decomposition: it was the BAND, not the floor

FP/HITEMP κ_R at (1500, 300), varying band and floor independently:

| band          | floor | FP/HT κ_R |
|---------------|-------|-----------|
| [50, 4500]    | 0     | 6.47      |
| [50, 4500]    | 1e-6  | 6.47      | ← floor inactive
| [100, 4000]   | 0     | 3.43      | ← matches crossover
| [100, 4000]   | 1e-6  | 3.43      |

- The **band** drives it: [50,4500]→[100,4000] moves 6.47→3.43. The far-band
  transparency windows (50–100 and 4000–4500 cm⁻¹) dominate the harmonic
  Rosseland mean and inflate the FP/line-list ratio for the fine line lists more
  than for smoothed FP.
- The **1e-6 floor is inactive**: every source's in-band κ_min/κ_max
  (FP 1.8e-3, POKAZATEL 3e-4, HITEMP 6e-6) is ≥ 1e-6, so the floor never binds.
  (It binds and matters at ≥1e-4; kept as a robustness knob.)

## wstep is not a systematic

HITEMP is at wstep 0.01 cm⁻¹, POKAZATEL at 0.05, FP smoothed to 20. Downsampling
HITEMP 0.01→0.05 changes κ_R by **−0.0%**, so the resolution difference does not
bias the Rosseland comparison — both line lists resolve the in-band windows.

## Line-list agreement and the residual systematic

- HITEMP and POKAZATEL agree on κ_R to ~1.15–1.46× (POKAZATEL slightly less
  opaque) and on κ_P to ~5% — a genuine, modest line-list difference.
- The **band choice is the real ~2× systematic** on the *absolute* Rosseland
  mean (hence on the cooling ratio): [50,4500] gives ~2× larger FP/LL κ_R than
  [100,4000]. We adopt **[100, 4000]** as the default science band (A25): it is
  crossover-consistent, mutually consistent across FP/HITEMP/POKAZATEL, and the
  well-resolved core. The [50,4500] result is reported as the systematic bound.

## Confidence statement

The Rosseland numbers can be stated with confidence **on the [100, 4000] band**:
FP/HITEMP reproduces the independent crossover to ~1%, resolution is irrelevant,
and the two line lists agree to ~1.5×. The remaining uncertainty is the band
choice itself (~2× between [100,4000] and [50,4500]), which propagates to the
headline cooling ratio and is quoted as such.

Numbers: `outputs/Fx_opacity_compare.csv` (band recorded per row). Reproduce the
convention variants with `compute_point_means(..., rel_floor=...)` over any band.
