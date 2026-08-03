# Where does FP opacity matter for steam sub-Neptunes? (Complete finding)

This note supersedes the earlier "FP never matters" claim (A27, which was tested
only at T_int ≥ 50). The corrected picture (A28/A29) is below.

## RCB depth is set by the opacity value AND T_int

Higher κ or higher T_int → convection starts shallower (grad_rad ∝ κ, L). So:

- **Warm/young planets (T_int ≳ 50 K):** the radiative zone is dilute
  (ρ ≲ 7 kg/m³), far below FP's grid floor (109 kg/m³). FP is irrelevant; the
  composite ≡ pure line list.
- **Cold/old planets (T_int ≲ 30 K):** the RCB migrates deep — max radiative-zone
  ρ ≈ 60 (T_int=30), 120–450 (T_int=20), 140–680 (T_int=15) — into FP territory.

This reconciles the two earlier numbers (ρ≈540 at cold + low κ; ρ≈7 at T_int=50).
Cold planets are the observed population, so FP is relevant for them.

## The qualitative effect: FP closes a deep radiative window

Because FP is 2–3× more opaque than line lists at ρ ≥ 109 (Fx, crossover-
validated), pure line lists **open a deep radiative window** that the FP-informed
composite **closes** (those layers convect under FP). At the fiducial, T_int=20:
pure HITEMP has 3 RCBs and a radiative zone to ρ=453; composite has 1 RCB to
ρ=123. This is a qualitative RCB-topology difference — line lists can spuriously
predict deep radiative zones.

## But the radius impact is small

Across (M_c 2–20 M⊕, f_env 1–50%) at T_eq=500, T_int=20 (`D1_fp_domain`):

- The radiative zone reaches ρ≥109 in **31/36** cells (most of space, at cold T_int).
- The **structural impact |ΔR_p|(composite − line-list) ≤ 1.7%** everywhere,
  largest for low-mass, high-f_env (puffy, low-g) planets; **~0.5% at the
  fiducial**. Closing a near-convective deep window barely moves the radius.

## History effect and the escape window (`E2_history`, A29)

Evolving matched initial conditions (same M_c, f_env, T_eq, T_int sequence)
under composite vs HITEMP, fiducial, T_eq=600:

- **FP-relevant epoch: ≈ 2.6–14 Gyr.** The RCB only reaches ρ≥109 after ~2.6 Gyr;
  FP relevance grows with age (opposite to a naive "hot start = deep RCB").
- **Escape window (first 100 Myr): untouched.** There the RCB is at ρ ≈ 2–5
  kg/m³ (shallow, hot), so composite ≡ line list and **ΔR_p ≈ 0.00%**. FP leaves
  **no imprint on the escape-sculpting epoch** — a clean, reportable null for the
  mass-loss history.
- **Peak radius divergence ≈ 0.33%**, reached late and **persisting** to Gyr
  ages (it does not reconverge — it is driven by the ongoing deep-RCB regime).

## Deck-ready summary

1. FP water opacity is 2–3× higher than line-list at ρ ≥ 109 / T ≥ 1000
   (exact, crossover-validated, Fx).
2. This regime is **convective** in warm/young steam sub-Neptunes, so FP does not
   change their structure; the composite equals the line list there.
3. It becomes the **radiative** deep interior only for **old/cold** planets
   (≳ 2.6 Gyr), where FP **closes deep radiative windows** the line lists open —
   a qualitative RCB-topology change with a **small radius effect (≤ 1.7% across
   planet space, ~0.5% fiducial)** that **persists to Gyr ages**.
4. The **atmospheric-escape window (first 100 Myr) is unaffected** by the FP vs
   line-list opacity choice — no irreversible mass-loss imprint.
5. The **larger structural effect lives beyond steam sub-Neptunes** — denser
   envelopes whose RCB sits at ρ ≥ 109 while still warm (more massive water
   worlds, magma-ocean-coupled layers). That is the natural follow-up where FP's
   opacity difference propagates strongly, promoted from footnote to next study.
