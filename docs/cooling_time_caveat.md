# F3 Kelvin–Helmholtz cooling time: control-variable caveat

## The issue

SPEC §8 F3 is the "planets stay puffy longer" quantification: τ_FP / τ_LL, with
τ_KH = |E_env| / L_int. Computed at **fixed T_int** (the current F3), the ratio
comes out **< 1** (FP ≈ 0.90 × POKAZATEL) — i.e. FP appears to cool *faster*,
the opposite of the intended headline and opposite to F1/F2/Fx (FP is puffier,
higher opacity).

## Why

In this static model the internal luminosity is imposed by the boundary
condition:

    L_int = 4π R_p² σ T_int⁴.

At fixed T_int, the more opaque (FP) planet is puffier (larger R_p), so **L ∝ R²
is larger**, while its envelope is *less* gravitationally bound (larger r
everywhere ⇒ |E| slightly *smaller*). τ = |E|/L therefore decreases. Measured
at (M_c=5, f_env=5%, T_eq=600, T_int=150):

| source     | R_p [R⊕] | L [W]     | \|E\| [J]  | τ_KH [Gyr] |
|------------|----------|-----------|-----------|------------|
| FP         | 1.946    | 5.56e16   | 2.71e32   | 0.155      |
| POKAZATEL  | 1.876    | 5.16e16   | 2.79e32   | 0.172      |
| constant κ | 1.854    | 5.04e16   | 2.82e32   | 0.177      |

T_int is the wrong control variable for the opacity comparison: it *is* the
surface-flux boundary condition, so fixing it hides the opacity throttling that
makes "stay puffy longer" true. Two planets at the same interior entropy but
different opacity have **different T_int** — the opaque one radiates less (lower
L) and cools slower.

## The physically correct comparison (fork — awaiting owner decision)

1. **Fixed interior entropy / adiabat** (Phase-1-feasible): parameterize the
   sequence by the deep convective specific entropy s (available from AQUA),
   match sources at equal s, and compare τ. The opaque planet then has lower
   T_int and L → longer τ → "stays puffy longer". Requires re-labeling the F3
   x/parameter from T_int to s.
2. **Phase-2 time-integrated track** (definitive): integrate dE/dt = −L(t) along
   the sequence and compare radius(age); "stay puffy longer" is unambiguous
   there. This is the SPEC §6 evolution deliverable.

Until the owner chooses, F3 is kept as a **labeled diagnostic at fixed T_int**
(reproducible, but explicitly *not* the headline), annotated on the figure.
