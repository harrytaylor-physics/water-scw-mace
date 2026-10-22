# subneptune

Steam-envelope sub-Neptune interior structure & evolution.

Pure-H₂O (steam) envelope on an inert rocky core, built to answer one
question first: **how does replacing line-list water opacity (HITEMP /
ExoMol-POKAZATEL) with first-principles (FP) opacity change the interior
structure and cooling of steam-envelope sub-Neptunes?** — and to remain the
group's steam-world model afterwards.

Owner: Harry Taylor (Oxford). Physics choices are traceable to a citation or
an explicit derivation in `docs/` (see `CITATIONS.md`). Any choice made
without a citation is logged in `AMENDMENTS.md`.

See `SPEC.md` for the authoritative build specification.

## Status

Phase 1 (deck-critical): opacity tables → static structure → T(P)/RCB
comparison across opacity sources + Kelvin–Helmholtz cooling-timescale ratio.
Phase 2: thermal evolution → radius(age) tracks. Phase 3+: extensions (SPEC §9).

Build order: scaffolding/config → opacity → EOS → structure → literature
anchor (G5) → Phase-1 figures. Each module lands with its gate tests. Phase 2
starts only after every Phase-1 gate passes.

## Layout

    config/       yaml run configs + typed loader (config/schema.py)
    layout.py     where things live in water-scw-mace (spectra, figures, release, models)
    eos/          AQUA loader + thermodynamic interface
    opacity/      table builder + interpolators (FP, HITEMP, POKAZATEL)
    structure/    hydrostatic integrator, boundary conditions
    evolution/    thermal evolution driver (phase 2)
    validation/   gates (pytest) + literature benchmarks
    figures/      plotting scripts (one script per figure, config-path only)
    docs/         derivations + assumption ledger
    outputs/      run artifacts (gitignored)

## Requirements

Python ≥ 3.11; numpy, scipy, matplotlib. PyYAML for config parsing (config
only — see `AMENDMENTS.md`). Install for development:

    pip install -e .

## Data

This package lives inside `water-scw-mace` (see `../MIGRATION.md`). Locations
outside it are set in `layout.py`:

    ../data/fp_opacity_v5/T*_opacity.npz   the 52 FP spectra, v5 estimator (A83), git-LFS,
                                           with grid.json and crossover_hitemp.json
    ../data/release/                       release tables (A73), written by tools/release_tables.py
    ../docs/figures/subneptune/            Figs 1-4 and Extended Data, written by figures/

Raw inputs are not committed (SPEC §1). Configs read them from `paths.data_root`
(default `~/subneptune_data`):

    AQUA/Tables/aqua_eos_rhot_v1_0.dat
    hitemp_opacity/, pokazatel_opacity/     line-list spectra
    mtckd/mtckd_v4.3_self_ref.npz           MT_CKD 4.3 self-continuum reference
                                            coefficients (tools/extract_mtckd_ref.py, A58)

The raw dipole series (`ed_data/arc_dipoles/`, `ed_data/arc_a78/T*/`) are
gitignored and fingerprinted in `ed_data/*_sha256.json`; the gates that need
them skip without them. Overrides: `SUBNEPTUNE_FP_OPACITY_DIR`,
`SUBNEPTUNE_FIGS_DIR`, `SUBNEPTUNE_RELEASE_DIR`.

## Gates

Run all validation gates with `pytest validation/`. No science figure is
produced until the relevant gates pass (SPEC §7).
