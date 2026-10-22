# AQUA EOS: (ρ,T) → (P,T) inversion

## Why invert

The delivered AQUA table `aqua_eos_rhot_v1_0.dat` is tabulated on a regular
(ρ, T) grid (1501 ρ log-spaced ∈ [1e-10, 1e5] kg m⁻³ × 301 T log-spaced ∈
[100, 1e5] K; T is the fastest axis) with P, grad_ad, s, u, c, phase as output
columns. The structure integrator works in (P, T) (hydrostatic dP/dm), so it
needs ρ(P,T), grad_ad(P,T), s(P,T). The companion AQUA (P,T) table is not on
disk (AMENDMENTS A2), so these are obtained by inverting P(ρ,T) along T = const.

## Interpolation

Bilinear in (log P, log T) on the tabulated quantities:

1. Locate the query in the native log T grid → bracketing columns j0, j1 and a
   linear weight w in log T (clamped to [0,1] outside the grid).
2. Within each column invert P(ρ) at fixed T: `np.interp(logP, logP_col, ·)`
   gives log ρ, s, grad_ad, u, c at the query pressure. logP_col is strictly
   increasing (below).
3. Linearly blend the two columns in log T.

log ρ is interpolated in log–log; s, grad_ad, u, c linearly against log P.

## Coexistence plateau (AMENDMENTS A8)

Below the critical point (647 K) the table has **exact P(ρ) plateaus**: across
the vapor–liquid dome, at fixed T, P holds at the saturation pressure while ρ
spans [ρ_vapor, ρ_liquid] (up to ~457 flat steps per column near 300 K, phase
id 1). Smaller plateaus occur at high density (ice VII/X, superionic; ids −7,
−10). P(ρ,T) is therefore not strictly monotonic and ρ(P,T) is multivalued on
the dome.

For inversion the per-column log P is made **strictly increasing** by
`numpy.maximum.accumulate` followed by a cumulative ~1e-12 dex nudge applied
*only to flat/decreasing entries*:

```
x = maximum.accumulate(logP_col)
dup = [False, *(diff(x) <= 0)]
x = x + cumsum(dup) * 1e-12
```

Single-phase nodes are left unchanged (a downstream single-phase node is
shifted by at most ~n_dome·1e-12 dex ≈ 1e-9 in relative P — negligible vs the
G2 0.5% tolerance), so **single-phase round-trips stay exact** (verified to
~1e-9 in `validation/test_eos.py`). On the dome the nudge makes ρ(P) single
valued and branch-correct:

- P below the plateau → **vapor** branch,
- P above the plateau → **liquid** branch,
- P exactly on P_sat (measure zero) → a value inside the dome.

This matches the physical inward-integration path (vapor until P_sat, then a
density jump to liquid). Example at 500 K: ρ(20 bar) ≈ 9.6, ρ(30 bar) ≈ 817
kg m⁻³ — the dome is crossed as a jump, as it should be.

## Diagnostics & validation

- `EOSState.two_phase` (phase id ≤ 2) flags coexistence; structure runs report
  any dome traversal. The phase lookup is nearest-node and approximate within
  ~1 grid cell of a transition; ρ/s/grad_ad remain accurate there.
- Independent physical anchors (not table round-trips): liquid water ρ(1 bar,
  300 K) ≈ 987 kg m⁻³ (IAPWS ≈ 997, ~1% via interpolation off-node); superheated
  vapor ρ(1 bar, 500 K) ≈ 0.435 kg m⁻³ vs ideal-gas 0.433.
