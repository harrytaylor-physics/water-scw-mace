# convergence/

Evidence that the numerical parameters are converged. Copy the artefacts
from ARC (~/scw and scw_prod) as listed; each subfolder gets the figure(s),
the generating script, and one paragraph of interpretation.

| Test | Evidence | Source on ARC |
|---|---|---|
| Finite size (32x range) | pressure/RDF/S(k) flat 64 -> 2048 molecules | finite_size.py, plot_finite_size.py |
| Plane-wave cutoff | 600/60 Ry (hybrid campaign benchmarks) | htc_bench.sh outputs |
| Timestep | 0.5 fs production vs coarser (_05fs redo) | run_prod_md_05fs.py, band_compare |
| MD length | 60 ps production; ACF window convergence | spectra.py diagnostics |
| ACF window / nfft | Hanning maxlag 8000, nfft 2^17 | spectra.py |
| Rosseland window/floor | 100-5000 cm^-1, floor 1e-6 | opacity_pipeline.py |
| Low-k S(k) at 673 K | real compressibility physics, not artifact | validation campaign |

Known caveat: the 6400-molecule finite-size runs OOM'd on 44 GB GPUs;
2048 molecules (6144 atoms) is the largest converged size — sufficient
for a 32x volume range over the 64-molecule production box.
