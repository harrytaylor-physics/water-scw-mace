# Line-list opacity tables — generation spec (HITEMP now, POKAZATEL later)

What to generate with RADIS so the loader/registry ingest them with zero code
change, the composite source (A11/A12) works, and pure line-list vs FP
comparison figures are exact on a common band.

## Where & filenames

    ~/subneptune_data/hitemp_opacity/T{T}_rho{rho}_opacity.npz
    ~/subneptune_data/pokazatel_opacity/T{T}_rho{rho}_opacity.npz   (later)

Directory names match the registry (`opacity/registry.py`); the loader globs
`T*_opacity.npz` and reads (T, ρ) from the file contents, so the exact filename
text is free as long as it matches the glob.

## npz schema (per file)

Same convention as the FP files, **correct units — no ×1000**:

| key          | units   | shape   | notes                                  |
|--------------|---------|---------|----------------------------------------|
| `wn`         | cm⁻¹    | (N,)    | wavenumber grid (monotone)             |
| `kappa`      | cm²/g   | (N,)    | opacity per unit mass                  |
| `T`          | K       | scalar  | temperature of this spectrum           |
| `rho_kg_m3`  | kg/m³   | scalar  | density of this spectrum               |
| `alpha`      | cm⁻¹    | (N,)    | optional; = κ·ρ(cgs). Recomputed if absent |

Extra keys are fine (ignored). Do **not** apply the legacy 1e-3 factor — that
is tied to the existing FP files only. The loader asserts each spectrum spans
the common band and that κ > 0.

## (T, ρ) grid — 11 × 11 = 121 spectra

Chosen so the line list covers the opacity-sensitive radiative zone
(T≈600–1030 K, ρ≈2×10⁻⁶–540 kg/m³) that FP cannot reach, the composite blend
window (default ρ=100–300; A15), and overlaps FP (T≥1000, ρ≥100) for exact
comparison.

- **T (K):** 500, 600, 700, 800, 900, 1000, 1250, 1500, 2000, 2500, 3000
- **ρ (kg/m³):** 1e-5, 1e-4, 1e-3, 1e-2, 0.1, 1, 10, 30, 100, 300, 1000

Rectangular (every ρ at every T) keeps the interpolation hull convex and
minimizes clamping. If cheap, extend ρ up to 2000 and T to 4000 for fuller deep
coverage (harmless where the interior is convective; reduces reliability flags).

## Wavenumber range & resolution

- **Range:** at least **[≤50, ≥6000] cm⁻¹** so the common band stays wide
  (it is capped at min-of-maxes = FP's ~6000). Ideally extend to **~20000 cm⁻¹**:
  the line list is the only source that can, and it lets us quantify the
  Rosseland band-truncation bias (`docs/band_limitation.md`) — the harmonic
  Rosseland mean is dominated by the near-IR/optical transparency windows above
  6000 cm⁻¹.
- **Resolution:** native/fine (≤1 cm⁻¹). Do not heavily smooth — the Rosseland
  mean is sensitive to unresolved transparency windows. Means are computed at
  the resolution you provide. (The FP files were smoothed to 20 cm⁻¹; the
  FP↔line-list resolution difference is itself physical, so we do not force the
  line list to match — happy to revisit.)

## Pressure for broadening

Use AQUA P(T, ρ) for consistency with the structure EOS (the earlier crossover
set used IAPWS/ideal). Values [bar] on the grid above:

| ρ \ T   |  500  |  600  |  700  |  800  |  900  | 1000  | 1250  | 1500  | 2000  | 2500  | 3000  |
|---------|-------|-------|-------|-------|-------|-------|-------|-------|-------|-------|-------|
| 1e-5    |2.3e-5 |2.8e-5 |3.2e-5 |3.7e-5 |4.2e-5 |4.6e-5 |5.8e-5 |6.9e-5 |1.0e-4 |2.6e-4 |4.1e-4 |
| 1e-4    |2.3e-4 |2.8e-4 |3.2e-4 |3.7e-4 |4.2e-4 |4.6e-4 |5.8e-4 |6.9e-4 |9.6e-4 |1.7e-3 |3.7e-3 |
| 1e-3    |2.3e-3 |2.8e-3 |3.2e-3 |3.7e-3 |4.2e-3 |4.6e-3 |5.8e-3 |6.9e-3 |9.4e-3 |1.3e-2 |2.6e-2 |
| 1e-2    |2.3e-2 |2.8e-2 |3.2e-2 |3.7e-2 |4.2e-2 |4.6e-2 |5.8e-2 |6.9e-2 |9.3e-2 |0.12   |0.18   |
| 0.1     |0.23   |0.28   |0.32   |0.37   |0.42   |0.46   |0.58   |0.69   |0.93   |1.2    |1.5    |
| 1       |2.3    |2.8    |3.2    |3.7    |4.1    |4.6    |5.8    |6.9    |9.2    |12     |14     |
| 10      |21     |26     |31     |36     |41     |46     |58     |69     |92     |1.2e2  |1.4e2  |
| 30      |26*    |69     |87     |1.0e2  |1.2e2  |1.4e2  |1.7e2  |2.1e2  |2.8e2  |3.5e2  |4.2e2  |
| 100     |26*    |1.2e2  |2.2e2  |2.9e2  |3.5e2  |4.4e2  |5.8e2  |6.9e2  |9.2e2  |1.2e3  |1.4e3  |
| 300     |26*    |1.2e2* |3.6e2  |6.2e2  |8.8e2  |1.2e3  |1.7e3  |2.2e3  |3.2e3  |4.1e3  |4.9e3  |
| 1000    |3.6e3  |5.8e3  |8.0e3  |1.0e4  |1.2e4  |1.4e4  |1.9e4  |2.3e4  |3.0e4  |3.7e4  |4.3e4  |

`*` = inside the vapor–liquid coexistence dome (P pinned at P_sat, e.g. ~26 bar
at 500 K). There the gas-phase line list at that ρ/P_sat is an approximation
(the true state is a two-phase mix); these points sit at moderate ρ and are
down-weighted by the composite blend. Regenerate this table for any grid change:
see the snippet in the commit adding this file.

## After delivery

Drop the files in `hitemp_opacity/`; then `get_source("hitemp", cfg)` and
`get_source("composite", cfg)` work with no code change, the common band
re-intersects automatically (A14), and pure `hitemp` / composite runs become
selectable. POKAZATEL: identical spec in `pokazatel_opacity/`.
