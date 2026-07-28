# TODO — technical (not physics)

Priority order.

1. **Regenerate stale npz data.** opacity_v3.npz (all 27+ points) and
   scw_port/opacity_surface.npz still carry the x1000 factor; code is fixed.
   Rerun opacity_pipeline over existing dipoles.npy (cheap — no MD needed)
   and retire the --alpha-scale workaround.
2. **Verify eos reference value** in benchmarks/reference/reference.json
   against iapws before first group release (see _note).
3. **Populate reference/holdout_v4.xyz** (held-out labelled split) so the
   parity benchmark actually runs; port eval_split.sh output.
4. **git-lfs setup + push models** (5 files, see models/README.md).
5. **Copy convergence artefacts** from ARC per convergence/README table.
6. **Copy deck figures** into docs/figures/ (deck.html lists exact sources
   per slide; broken images show the source path as fallback text).
7. **Pin versions**: freeze mace-torch/ase/torch versions that v4 was
   trained/evaluated with into requirements.txt (currently unpinned).
8. **CI**: GitHub Action running pytest (code-only tests pass without lfs).
9. **float32 defaults** for all MD drivers (propose/production); float64
   only for parity evaluation.
10. **Fix main_3000K_aimd_db_compare.py line ~198** — same x1000 bug,
    never patched.
11. **Selection/collection skip-markers** consistency pass across pipeline
    scripts (propose has DONE; selection now cached; labels rerun-safe).
12. **n_inf upgrade**: CP2K polarizability single points per state point
    (and per functional) to replace the Lorentz-Lorenz literature alpha_e
    in the emissivity module.

## Post-LDA-campaign (27 Jul 2026)

Done this campaign: LDA label/train/gate chain; 21-point LDA production;
LDA MLWF dipole model; 3-point spectral comparison; LDA emissivity at
compared points; EOS functional comparison + paper figure; dispersion
decomposition; AWE deck.

- [ ] Dense-branch LDA potential (harvest labelled and banked;
      build_v2_dataset.py ready) -> dense spectra + emissivity comparison
- [ ] Larger MLWF dipole campaign to lower the ACF noise floor
      (data lever at ~640 configs documented as spent)
- [ ] H-bonded dimer geometry fix (current scan misses the canonical
      well; D3 difference curve quoted instead)
- [ ] ExoMol/BT2 cross-check: resolve RADIS<->server chunk-naming
      mismatch (first chunk 404s; exomol_selfheal.py has the scaffold)
- [ ] Merge docs/AMENDMENTS_APPEND.md into AMENDMENTS.md
- [ ] EOS figure panel-(a) marker/legend consistency (cosmetic)
