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
    data/         NOT committed; path-config or symlink to ~/subneptune_data
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

Not committed. Configs default to `~/subneptune_data`:

    ~/subneptune_data/AQUA/Tables/aqua_eos_rhot_v1_0.dat
    ~/subneptune_data/fp_opacity/T*_opacity.npz  (27 FP spectra)
    ~/subneptune_data/fp_opacity/grid.json
    ~/subneptune_data/fp_opacity/crossover_hitemp.json

Override the root with `paths.data_root` in a config file.

## Gates

Run all validation gates with `pytest validation/`. No science figure is
produced until the relevant gates pass (SPEC §7).
