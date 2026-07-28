# Band truncation of opacity means

SPEC §3 fixes a **common band** [50, 4500] cm⁻¹ for every opacity source so the
*differential* conclusions between sources are exact. The price is an
**absolute** band-truncation error, which this note quantifies. It is a stated
limitation, not a bug: it cancels in FP-vs-line-list comparisons because both
use the identical band.

## Planck energy captured in-band

The Planck function peaks (per unit wavenumber) near wn ≈ 1.62 k_B T / (h c),
which for 1000–3000 K sits inside [50, 4500] cm⁻¹, but the distribution has
tails that leave the band — increasingly so as T rises and the peak marches
toward the upper edge:

| T (K) | Planck energy fraction in [50, 4500] cm⁻¹ |
|-------|-------------------------------------------|
| 1000  | 0.918 |
| 1500  | 0.778 |
| 2000  | 0.681 |
| 2250  | 0.647 |
| 2500  | 0.620 |
| 2750  | 0.598 |
| 3000  | 0.580 |

## Common-band vs full-band means (FP source)

Line-list full-band spectra are not yet available (AMENDMENTS A5), so FP
(native span [1.11, 6000] cm⁻¹) is used as the proxy for full-band. Ratio =
common-band mean / full-band mean, across the 27 FP grid points:

| Mean       | min ratio | median | max ratio |
|------------|-----------|--------|-----------|
| Planck     | 1.017     | 1.269  | 1.494     |
| Rosseland  | 1.018     | 1.527  | 2.144     |

Common-band means run **above** full-band means because the excluded regions
(< 50 and 4500–6000 cm⁻¹) are relatively transparent water windows; dropping
them removes low-opacity weight. Rosseland is affected most (its harmonic
weighting is dominated by transparency windows).

## Consequences

- Absolute R_p / RCB from band-limited means carry this offset; it is reported
  as a limitation on any absolute figure, not tuned away.
- The FP↔line-list *differential* (F1/F3) is unaffected: identical band, so the
  truncation is common-mode and cancels.
- **To revisit once line-list spectra arrive:** recompute this table with the
  line-list source (SPEC §3 asks for full-band vs common-band on the *line-list*
  source specifically) and, if the largest shared band across all delivered
  sources exceeds [50, 4500], widen the common band and re-record (AMENDMENTS A7).

Reproduce: `compute_point_means("fp", cfg, band)` for `band=(50,4500)` and
`band=(2,6000)`.
