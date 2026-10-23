# Boundary sensitivity: γ and matching pressure

SPEC §5/§10 require the default AND a sensitivity scan for the outer-boundary
knobs. Fiducial here: M_c=5 M⊕, f_env=5%, T_eq=600 K, T_int=100 K. Values below
use a constant-κ reference (FP is flagged unreliable in the radiative zone;
regenerate against the `composite` source once HITEMP lands — trends are
boundary-driven and largely opacity-independent). Reproduce with
`figures/plot_sensitivity.py`.

## γ = κ_visible / κ_thermal (fiducial γ = 1)

| γ    | R_p [R⊕] | RCB [bar] | T_phot [K] |
|------|----------|-----------|------------|
| 0.25 | 1.8252   | 3.50      | 743        |
| 0.50 | 1.8196   | 2.26      | 740        |
| 1.00 | 1.8162   | 1.59      | 727        |
| 2.00 | 1.8145   | 1.25      | 695        |
| 4.00 | 1.8139   | 1.07      | 658        |

**R_p is robust to γ: ~0.6% across γ ∈ [0.25, 4]** (greenhouse→inversion). The
photosphere temperature and RCB pressure move more (γ<1 warms the deep
radiative layer and pushes the RCB deeper; γ>1 the reverse), but the integrated
radius barely changes. Good news for the "γ is a knob Ray will ask about"
concern — the headline R_p is not γ-sensitive; the RCB/T(P) shape is.

## Matching-pressure cap (fiducial 10 bar, swept 1–100 bar)

With adaptive matching (A13) the model matches at min(cap, P_RCB):

| cap [bar] | R_p [R⊕] | matched at | at RCB? |
|-----------|----------|------------|---------|
| 1         | 1.8214   | 1.00 bar   | no (cap binds) |
| 10        | 1.8162   | 1.69 bar   | yes |
| 100       | 1.8162   | 1.69 bar   | yes |

**Adaptive matching removes the matching-pressure sensitivity.** For any cap at
or below... i.e. any cap ≥ the RCB pressure (~1.69 bar here), the match lands on
the RCB and R_p is *identical* (1.8162). The cap only affects the result if
forced *below* the RCB (cap=1 bar → matches at 1 bar, +0.3% in R_p). This is the
intended fix for the fixed-10-bar failure that collapsed high-opacity (FP)
models: the hand-off now happens where the atmosphere is actually radiative, not
at an arbitrary fixed depth.

## Takeaway

Neither knob materially moves R_p for these steam worlds (each ≲ 1%); they
reshape the T(P)/RCB, which is what F1/F2 display. The absolute-radius
systematics that matter more are the opacity source itself and the band
truncation (`docs/band_limitation.md`).
