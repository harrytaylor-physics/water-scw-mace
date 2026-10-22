# MIGRATION — the subneptune planetary model into water-scw-mace

On 2026-10-22 the steam-envelope sub-Neptune interior and evolution model
(`subneptune`, previously a local repository with no remote) and the products of
its v5 opacity work were brought into this repository on the branch
`integrate-subneptune-v5`. This file records what went where, what was left out,
and how the history was handled.

## Layout after the migration

    subneptune/              the planetary model, self-contained, with its history
    production_campaign/     the ARC production pipeline (MD, labels, training, opacity), verbatim
    data/fp_opacity_v5/      the 52 v5 first-principles spectra (git-LFS)
    data/release/            MESA-format kappa_R / kappa_P tables
    docs/figures/subneptune/ Figs 1-4 and the ten Extended Data items (PDF + 200-dpi PNG)
    models/mlip_heldout_errors.json   the v4 potential's error tables

Every location outside `subneptune/` that its code reads or writes is named once,
in [`subneptune/layout.py`](subneptune/layout.py). Provenance strings in the
figure ledgers that name `benchmarks/`, `data/`, `docs/`, `models/` or
`production_campaign/` are relative to the repository root. All other provenance
strings are relative to `subneptune/`.

## Mapping

Source paths are in `~/subneptune` unless marked `~/subneptune_data` (the data
root outside the repository, SPEC §1).

| source | destination | reason |
|---|---|---|
| the repository: `config/`, `constants.py`, `atomicio.py`, `eos/`, `structure/`, `evolution/`, `opacity/`, `tools/`, `figures/` (scripts), `validation/`, `conftest.py`, `pyproject.toml`, `docs/`, `monomer_test/`, `generate_linelist_tables.py`, `AMENDMENTS.md`, `SPEC.md`, `CITATIONS.md`, `README.md` | `subneptune/` (same relative paths) | one self-contained top-level component, named after the package (as `lda_campaign/` and `convergence/` are) |
| — | `subneptune/layout.py` (new) | the host-repository locations, in one place |
| `opacity/dipole_spectrum.py`, `tools/v5_regenerate.py` | `subneptune/opacity/`, `subneptune/tools/` (unchanged) | the v5 estimator stays with the code whose gates import it (8 importers). The A78-A83 fix is documented beside the pipeline it supersedes, in `production_campaign/README.md` |
| `ed_data/arc_a78/scw/*.py` (72 files) | `production_campaign/scw/` | the production MD, labelling, training and opacity scripts that produced the 52-point grid; kept unmodified as provenance |
| `ed_data/arc_scripts/` (Slurm scripts, cell lists, `logs/train_v4_7924966.out`) | `production_campaign/arc_scripts/` | where this repository keeps HPC scripts: beside the Python they run, per campaign |
| `ed_data/arc/water-scw-mace/docs/METHODS.md` (A90-corrected) | `docs/METHODS.md` | the A90 corrections applied to this repository's own file; that copy and this file differed in nothing else |
| `ed_data/mlip_heldout_errors.json` | `models/mlip_heldout_errors.json` | with the model documentation; paths now repository-relative |
| `figs/fig_numbers.json`, `figs/ed/{ed_index,ed_numbers}.json`, `figs/ed/ed_table{1,2}.{csv,tex}` | `docs/figures/subneptune/` | with this repository's figures; the producers now write there (`layout.FIGS_DIR`) |
| `figs/*.pdf` (gitignored there) | `docs/figures/subneptune/` as `.pdf` + `.png` (200 dpi) | regenerated from the committed scripts; `save_pdf` writes both |
| `data/release/` | `data/release/` | release tables with the repository's data (`layout.RELEASE_DIR`) |
| `~/subneptune_data/fp_opacity_v5/` (52 `.npz`, `grid.json`, `crossover_hitemp.json`) | `data/fp_opacity_v5/` (git-LFS) | the v5 spectra, byte-identical copies; read from here (`layout.spectra_dir`) |
| `ed_data/arc_dipoles_sha256.json` | `subneptune/ed_data/` (unchanged) | beside the location the gitignored bundle it fingerprints is read from |
| — | `subneptune/ed_data/arc_a78_dipoles_sha256.json` (new) | fingerprint of the five A78 cells, which are no longer committed |
| `ed_data/arc/scw/` (`lagtest_v2.json`, two CSVs, 33 audit PNGs) | `subneptune/ed_data/arc/scw/` (unchanged) | ARC-pull data the Extended Data reads |
| `monomer_test/water_dipole_v3.model` (LFS) | `subneptune/monomer_test/` (LFS) | `models/` holds `water_dipole_v2` only, so there is nothing to point to instead |

## Left out, and why

| path | reason |
|---|---|
| `ed_data/arc_a78/T*/` — 7 cells: `dipoles.npy`, `opacity_v3.npz` | raw dipole series and superseded v3 spectra; the two `*_BAD64` cells are unused. Stripped from history. The five used cells are fingerprinted in `subneptune/ed_data/arc_a78_dipoles_sha256.json` and gitignored |
| `ed_data/arc_dipoles/` (67 series, 41 MB) | raw production output; never committed (SPEC §1); fingerprinted |
| `ed_data/arc/water-scw-mace/` (56 files) | a stale copy of this repository inside the old one. Stripped from history. Its readers now read this repository's own `benchmarks/reference/reference.json`, `docs/LDA_CAMPAIGN.md` and `docs/METHODS.md` |
| `outputs_v3/` (v3-era outputs, the slide deck `.pptx`), `DECK_OUTLINE.md` | superseded v3 products and presentation material. Stripped from history |
| `figs_v3/`, `data_release_v3/` | superseded v3 products; never committed |
| `outputs/`, including `outputs/manuscript_numbers.json` (the numbers ledger) | run outputs, gitignored; the ledger is rebuilt by `tools/manuscript_numbers.py` |
| `~/subneptune_data/fp_opacity/` (v3), `fp_opacity_v4/`, `fp_opacity_legacy_x1000/` | superseded spectra |
| the rest of `~/subneptune_data/` (AQUA, HITEMP, POKAZATEL, MT_CKD, PALEOS, Aguichine 2025 tracks) | external raw inputs (SPEC §1); still read from `paths.data_root` in the config |
| any manuscript text, bibliography, cover letter or proof | never part of the `subneptune` repository — the history scan below found none |

`figures/build_deck.py` (code that builds a slide deck into the gitignored
`outputs/`) came across as part of the package; the deck itself did not.

## History

The planetary model was imported **with its history**, after filtering:

1. Scan of the whole `~/subneptune` history (115 commits, then 116 with A94).
   - No manuscript files: no `\documentclass`, abstract, `.bib`, proof or cover
     letter. The only `.tex` files are the generated Extended Data tables.
   - No blob over 50 MB (the largest was 2.8 MB).
   - No credential-like path or secret-shaped content.
   - But the history did contain the excluded items above.
2. On a throwaway clone, `git filter-repo --invert-paths` removed from every
   commit: the seven `ed_data/arc_a78/T*/` cells, `ed_data/arc/water-scw-mace/`,
   `outputs_v3/`, `figs_v3/`, `data_release_v3/` and `DECK_OUTLINE.md`. All 116
   commits survive, and a check of every object in the result found none of the
   excluded files.
3. `git subtree add --prefix=subneptune` of the filtered history (commit
   `d8475a7`). Later moves out of `subneptune/` used `git mv`, so `git log
   --follow` traces them.

**Commit hashes.** The first 27 commits kept their hashes. Everything from
`25cdad5` (2026-08-04, "Add deck figures M1-M4 + DECK_OUTLINE.md") onward was
rewritten. `subneptune` A94 is `dc8d47c` in the original repository and
`588318a` here. Hashes that `subneptune/AMENDMENTS.md` cites from after that
point refer to the original repository.

## Running it in place

    cd subneptune
    python -m pytest            # the planetary model's gates (SPEC §8)

- The FP spectra are read from `data/fp_opacity_v5/`.
- Line lists, AQUA, MT_CKD and the other raw inputs are read from `paths.data_root`
  in `config/*.yaml` (default `~/subneptune_data`).
- Gates whose inputs are absent skip cleanly. These are the line lists, the dipole
  bundle (`ed_data/arc_dipoles/`) and the A78 cells (`ed_data/arc_a78/T*/`).
- Environment overrides: `SUBNEPTUNE_FP_OPACITY_DIR`, `SUBNEPTUNE_FIGS_DIR`,
  `SUBNEPTUNE_RELEASE_DIR`.

To rebuild every product: `bash tools/rebuild_producers.sh`. It writes the
figures to `docs/figures/subneptune/`, the release tables to `data/release/` and
the ledger to `subneptune/outputs/`.
