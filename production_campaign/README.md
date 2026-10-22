# production_campaign/

The hybrid (revPBE0-D3) production pipeline that produced the 52-point
first-principles grid, as it ran on ARC: MD, DFT labelling, MACE and dipole-model
training, and the opacity post-processing. **Kept unmodified, as provenance.**
Paths inside the scripts (`/home/shug7909/scw`, `module load ...`) are ARC's and
are left as they ran.

    scw/          72 Python drivers and analysis scripts (from the ARC working
                  directory ~/scw), including opacity_pipeline.py
    arc_scripts/  the Slurm submission scripts (htc_*.sh), the cell lists they
                  iterate over (*.list), and logs/train_v4_7924966.out, the
                  training log of the production potential

| stage | scripts |
|---|---|
| DFT labelling (CP2K) | `arc_scripts/htc_label_*.sh`, `scw/harvest_aimd.py`, `scw/build_trainxyz*.py` |
| MACE training (v4) | `arc_scripts/htc_train_v4.sh`, `htc_train_v4_resume.sh`; log in `arc_scripts/logs/` |
| dipole model | `arc_scripts/htc_dipole_label3*.sh`, `htc_dipole_train_v3.sh`, `scw/build_dipole_xyz*.py` |
| production MD | `arc_scripts/htc_prod_md*.sh`, `htc_md05*.sh`; `scw/run_prod_md*.py` |
| spectra and opacity | `arc_scripts/htc_opac*.sh`, `htc_spectrum05.sh`; `scw/opacity_pipeline.py` |

The production MD loaded `water_scw_v4.model` -- the stage-one model, sha256
`dc4a8977...`, which is [`models/water_scw_v4.model`](../models/). Its error
table is in [`models/mlip_heldout_errors.json`](../models/mlip_heldout_errors.json)
(subneptune A94).

## The delivered spectra are superseded

**`scw/opacity_pipeline.py` line 82 applies the harmonic quantum correction
twice.**

```python
sigma_cm2_raw = (beta * omega**2 / (3.0 * eps0 * c * n_molecules)) * (Smm * qcf) * 1.0e4
```

`beta * omega**2 * Smm` is already the harmonic-corrected absorption, and line 82
multiplies it by `qcf = x / (1 - exp(-x))` again. Every spectrum the pipeline
delivered is too large by that factor -- about 1.2-1.3 in the rotational band and
3.5-7.8 in the vibrational bands at 700 K. The normalisation itself (the
one-sided `/(3 eps0 c N)`) is correct. The pipeline's spectral estimator had two
further defects: a symmetric Hann window on a one-sided autocorrelation, clipped
at zero, and an untapered periodogram whose spectral leakage was mistaken for a
noise floor.

The spectra in [`data/fp_opacity_v5/`](../data/fp_opacity_v5/) replace them.
They were regenerated from the same dipole series by
[`subneptune/opacity/dipole_spectrum.py`](../subneptune/opacity/dipole_spectrum.py)
via [`subneptune/tools/v5_regenerate.py`](../subneptune/tools/v5_regenerate.py):
`sigma = beta w^2 I2 / (6 eps0 c N)` on the two-sided transform, one harmonic
correction, a Hann taper on the dipole **time series** with `sum(w^2)`
normalisation, a continuous `omega^2` prefactor, no zero clip, 10 cm^-1 Gaussian
smoothing. The estimator stays inside `subneptune/` because the planetary model's
gates import it.

| amendment | finding |
|---|---|
| A78 | the double harmonic correction at line 82, proved at the source: `beta <\|Mdot\|^2>` from the dipole series is flat in T while the delivered band integrals fall |
| A79 | the normalisation is correct (factor 1), pinned by a synthetic classical oscillator against the textbook band strength |
| A80 | the estimator replacement (v4) |
| A81 | the apparent noise floor is spectral leakage from the untapered periodogram; a Hann taper on the time series removes it |
| A82 | the smoothing width (10 cm^-1) |
| A83 | v5: all 52 spectra regenerated and propagated downstream |

The amendments are in [`subneptune/AMENDMENTS.md`](../subneptune/AMENDMENTS.md).
Gates that exercise this file as written -- reproducing what it delivered --
import it from here (`subneptune/layout.py: PRODUCTION_SCW`); see
`subneptune/validation/test_a78_dipole_sum_rule.py` and `test_a80_v4_gates.py`.

The raw dipole series these scripts produced are not in the repository; their
SHA-256 fingerprints are in `subneptune/ed_data/arc_dipoles_sha256.json` (the 67
series behind the 52 spectra) and `subneptune/ed_data/arc_a78_dipoles_sha256.json`
(the five cells A78 reproduces).
