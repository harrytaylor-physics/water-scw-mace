# Effect of the XC functional (and dispersion) on opacity

**Question.** How much do the exchange-correlation functional and the
dispersion correction matter for warm-dense-water IR opacity and
emissivity?

**Design.** Fully consistent pipelines compared only at the level of final
spectra on the shared 27-point production grid:
- revPBE0-D3 everywhere (done — the production pipeline)
- LDA everywhere (in progress): LDA force labels -> water_lda_v1 (fresh
  fine-tune from mace-mpa-0-medium, no label mixing), LDA MLWF dipole
  labels -> LDA dipole model, same spectra/opacity post-processing.
- (planned) revPBE0 without D3: isolates dispersion at fixed exchange.

**Pre-registered metrics** (per state point): absorption-edge shift,
Planck-mean ratio, Rosseland-mean ratio, normal-emissivity RMS.

**Status.** Steps A-C of the LDA fast track running (proposal MD done,
selection/labelling in flight). Numerics inherited verbatim from the
hybrid campaign (CUTOFF 600/REL_CUTOFF 60/DZVP-MOLOPT-SR-GTH); only XC
(PADE) and pseudopotentials (GTH-PADE) change. Gate thresholds and
amendments tracked in the campaign's AMENDMENTS.md.
