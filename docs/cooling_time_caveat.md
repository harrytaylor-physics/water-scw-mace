# F3 Kelvin–Helmholtz cooling time: control variable and the "stay puffy longer" result

This note explains why the F3 cooling-time ratio must be taken at fixed interior
entropy, not fixed T_int, and is written to survive an expert referee.

## Definitions

The internal (cooling) luminosity of the static model is

    L_int = 4π R_p² σ T_int⁴,                                            (1)

the Kelvin–Helmholtz time is

    τ_KH = |E_env| / L_int,   E_env = ∫ (u − Gm/r) dm,                   (2)

and the "planets stay puffy longer" claim is that replacing line-list opacity
with the (higher) first-principles opacity **lengthens** the cooling time at
otherwise fixed conditions: τ_FP / τ_LL > 1.

## Why fixed T_int cannot show the effect (it is structural, not numerical)

T_int is not an independent physical property of the planet — via (1) it *is*
the surface-flux boundary condition, i.e. the luminosity itself. Holding it
fixed while swapping the opacity forces both planets to radiate through the same
L at the same T_int, which is precisely the quantity opacity is supposed to
throttle. What then differs is only the structure: the more opaque (FP) planet
has a steeper radiative gradient, a hotter deep adiabat, a larger radius, and —
because it is more extended — a slightly *less* bound envelope (smaller |E|).
With L ∝ R_p² larger and |E| smaller, (2) gives τ_FP < τ_LL. Measured at
(M_c=5 M⊕, f_env=5%, T_eq=600 K, T_int=150 K):

| source     | R_p [R⊕] | L [W]     | \|E\| [J]  | τ_KH [Gyr] | s_int [J/kg/K] |
|------------|----------|-----------|-----------|------------|----------------|
| FP         | 1.946    | 5.56e16   | 2.71e32   | 0.155      | 13991          |
| POKAZATEL  | 1.876    | 5.16e16   | 2.79e32   | 0.172      | 13247          |

So fixing T_int gives τ_FP/τ_LL ≈ 0.90 < 1. This is not a bug and not a small
correction that could be tuned away — it is structural: T_int bakes in L, so the
comparison is ill-posed for the question being asked. Note the interior
entropies differ (13991 vs 13247): the two planets are **not in the same
thermodynamic state**, so comparing their cooling times is not meaningful.

## The correct control variable: interior (deep-adiabat) specific entropy

The cooling history of a convective envelope is a sequence of states labelled by
the specific entropy of the deep adiabat, s_int (Lopez & Fortney 2014; Ginzburg,
Schlichting & Sari 2016). Two planets at the same s_int are in the same interior
thermodynamic state and have (to the accuracy of the outer boundary) the same
radius. Opacity enters by setting **how fast that state radiates away**: a more
opaque envelope transmits less flux, so at fixed s_int it sits at a lower T_int
and lower L, and therefore takes longer to cool. This is the physically correct
"stay puffy longer" statement.

Concretely, `structure.build_at_entropy` shoots on T_int (s_int is monotone in
T_int) to place each source on the same adiabat. At s* = 13500 J/kg/K
(T_eq=600 K, f_env=5%):

| source     | T_int [K] | R_p [R⊕] | L [W]     | τ_KH [Gyr] |
|------------|-----------|----------|-----------|------------|
| FP         | 103.4     | 1.898    | 1.19e16   | 0.735      |
| POKAZATEL  | 191.9     | 1.898    | 1.42e17   | 0.062      |

Now **R_p is identical** (entropy sets the structure), and the opacity shows up
exactly where it should — in the luminosity: FP radiates ~12× less, giving
**τ_FP / τ_POKAZATEL ≈ 12 > 1**. This is the F3 headline (left panel). The
fixed-T_int result is kept as the right-hand diagnostic panel to make the
distinction explicit and pre-empt the referee question.

## Scope and caveats

- With only FP and POKAZATEL (both on the hot/dense grid), the *structure* runs
  used here are opacity-clamped in the radiative zone, so the absolute τ values
  are provisional; the entropy *ratio* is the robust quantity and the mechanism
  is grid-independent. Numbers firm up once the composite/HITEMP tables cover the
  radiative zone.
- The definitive, assumption-free version is the Phase-2 time-integrated track
  (dE/dt = −L(t) along the entropy sequence) giving radius(age) per opacity
  source; that is where "stay puffy longer" is stated without reference to any
  single control variable.
- s_int is taken as the specific entropy at the deepest (base-of-envelope)
  point, which lies on the convective adiabat for these models.
