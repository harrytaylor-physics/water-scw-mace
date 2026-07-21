# 05 — post-processing convergence

**What:** stability of the derived quantities against post-processing
choices: Rosseland integration window and floor; Kramers-Kronig band
cutoff (emissivity).
**Criteria:** Planck/Rosseland means stable to <2% under window/floor
variation; KK-derived n stable to <5% at low frequency under band-cutoff
variation (the emissivity module computes the latter automatically per
point as kk_drift).

```
python3 means_stability.py opacity.npz          # wn + kappa fields
```

**Reference results (this project):** Rosseland window 100-5000 cm^-1,
floor 1e-6 relative; Planck over the full band. KK drift <5% at all 27
production points (emissivity_v2 diagnostics). Historical note for the
diary: absolute alpha/kappa in pre-2026-06-29 npz files carry a x1000
unit factor (fixed in code; ratios always unaffected).
