# 04 — spectra convergence (dipole-ACF method)

**What:** the numerical knobs between dipoles.npy and a spectrum:
ACF window (maxlag), total trajectory length, and output smoothing.
**Criteria:** band peak positions stable to <10 cm^-1 and integrated band
intensities to <5% under doubling/halving of each knob. Trajectory length
is the axis that sets production cost — converge it before committing a
grid.

```
python3 spectrum_convergence.py dipoles.npy --dt-fs 5 --temp 3000
```
Produces overlay figures: spectrum vs maxlag {half, ref, double},
vs trajectory fraction {25,50,75,100}%, vs smoothing {10,20,40} cm^-1,
plus a band-metrics table (libration/bend/stretch windows).

**Reference results (this project):** maxlag 8000 samples, Hanning window,
nfft 2^17, 20 cm^-1 output smoothing, 60 ps production @ 5 fs dipole
cadence. Quantum correction w(1-exp(-bhw)). White-noise floor measured in
3900-5500 cm^-1; band peaks quoted vs floor (spectra.py convention).
ACF smoothness gate C(dt)/C(0) >= 0.997.
