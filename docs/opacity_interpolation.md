# Opacity interpolation & the extrapolation gate

## Grid is scattered, not rectangular

SPEC §3 prescribes "bilinear in (T, log ρ) on log κ." The 48 FP grid points do
**not** form a tensor-product grid: T = 1000 K has a single ρ sample, and the ρ
sampling differs at every T (see `~/subneptune_data/fp_opacity/grid.json`).
Literal bilinear interpolation is therefore undefined here.

## Method used (AMENDMENTS A6)

Linear interpolation on a Delaunay triangulation of the scattered points —
the faithful scattered-data generalization of bilinear (piecewise-linear on
each triangle):

- Coordinates: (T [K], log₁₀ρ [ρ in kg m⁻³]).
- Interpolated quantity: log₁₀κ (SI, m² kg⁻¹), for both Rosseland and Planck
  means. Interpolating in log–log keeps the strongly varying opacity smooth.
- `scipy.interpolate.LinearNDInterpolator(..., rescale=True)`. `rescale`
  normalizes the disparate T (~10³) and log₁₀ρ (~2–3) magnitudes to a unit box
  before triangulation so the triangles are geometrically sensible.
- At a grid vertex the interpolant returns that vertex's value exactly (tested
  in `validation/test_opacity.py::test_interpolation_recovers_grid_values`).

## Extrapolation is a gate, not a guess (SPEC §3)

Queries outside the convex hull of the grid are **clamped to the nearest grid
value** (`NearestNDInterpolator`, same rescaled coordinates), and every clamped
call is counted:

- `OpacityTable.clamp` (a `ClampCounter`) tracks total and clamped evaluations;
  `source.clamp_fraction` reports the clamped fraction.
- The structure driver resets the counter before an integration and reads the
  fraction afterward. **> 20% of steps clamped ⇒ the run is flagged unreliable
  in its output metadata** (`opacity.clamp_fraction_threshold`, default 0.20).

### Expected heavy clamping with the current FP grid

The FP grid covers only T ∈ [1000, 3000] K and ρ ∈ [109, 2229] kg m⁻³. A
sub-Neptune envelope extends to much cooler, lower-density outer layers
(T down to T_eq ≈ 500–800 K), which fall outside the hull and clamp. Phase-1
FP structure runs are therefore expected to trip the clamp-fraction flag over
parts of the envelope; this is surfaced honestly rather than hidden, and is a
direct argument for extending the FP grid to lower T/ρ.
