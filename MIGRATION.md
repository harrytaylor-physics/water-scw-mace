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
| `monomer_test/water_dipole_v3.model` (LFS) | `models/water_dipole_v3.model` (LFS, one copy) | the dipole model behind the 52-point spectra, beside the earlier `water_dipole_v2`; the monomer test reads it through `layout.DIPOLE_MODEL` |

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
   `outputs_v3/` and `DECK_OUTLINE.md`. The rules for `figs_v3/` and
   `data_release_v3/` matched nothing, because neither was ever committed. All
   116 commits survive, and a check of every object in the result found none of
   the excluded files.
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

## Commit hash map

Every `~/subneptune` commit (old) and its counterpart on this branch under
`subneptune/` (new), oldest first, from `git filter-repo`'s own
`.git/filter-repo/commit-map`. Checked: all 116 old hashes are commits of
`~/subneptune`, and all 116 new hashes are ancestors of this branch.
27 kept their hash; 89 were rewritten. Use it to follow a hash that
`subneptune/AMENDMENTS.md` cites (`git show <new>`).

| # | author date | old (`~/subneptune`) | new (this repository) | | subject |
|---|---|---|---|---|---|
| 1 | 2026-07-28 16:20 | `a7b9f6d969f1aadfd9837891bb5f9da2843e63a5` | `a7b9f6d969f1aadfd9837891bb5f9da2843e63a5` | same | Scaffold repo, config system, constants, and docs |
| 2 | 2026-07-28 16:35 | `d696d2286a8c1c3254aaa41a5fc03a64abcaf065` | `d696d2286a8c1c3254aaa41a5fc03a64abcaf065` | same | Add opacity module: FP source, band-limited means, interpolation gate (G4) |
| 3 | 2026-07-28 18:12 | `efd2fb9cab6e9142bd73505f8024727c9e245acc` | `efd2fb9cab6e9142bd73505f8024727c9e245acc` | same | Add AQUA EOS module with (P,T) inversion and core M-R (gate G2) |
| 4 | 2026-07-29 13:25 | `4b2a5d7c2ed0071faca3674e4ebb8d36405c222e` | `4b2a5d7c2ed0071faca3674e4ebb8d36405c222e` | same | Add structure module: Guillot boundary, shooting integrator (gates G1/G3/G4/G7) |
| 5 | 2026-07-29 17:06 | `e8c1772bc4e32ed0693da22354186749e0edbb91` | `e8c1772bc4e32ed0693da22354186749e0edbb91` | same | Add adaptive matching at the RCB (AMENDMENTS A13) |
| 6 | 2026-07-29 17:13 | `0a9909d4d4cd9ebf91657de5b6d57808dd7af34b` | `0a9909d4d4cd9ebf91657de5b6d57808dd7af34b` | same | Auto-detect and widen the common opacity band (AMENDMENTS A14) |
| 7 | 2026-07-29 17:52 | `75491227c6e0a644f082114e5d60be2e80f74416` | `75491227c6e0a644f082114e5d60be2e80f74416` | same | Add composite opacity source: line-list + FP, log-blended (A11, A12) |
| 8 | 2026-07-29 18:26 | `38434f64b41c6d81587817445b542012bd9b1820` | `38434f64b41c6d81587817445b542012bd9b1820` | same | Add line-list (HITEMP/POKAZATEL) table generation spec |
| 9 | 2026-07-31 08:27 | `691c2767ecec1eb78b68a134e23f3da6c8173390` | `691c2767ecec1eb78b68a134e23f3da6c8173390` | same | Move default composite blend window to rho=100-300 (AMENDMENTS A15) |
| 10 | 2026-07-31 08:34 | `771a9c37ef6b5b7f68dd9ae478de9a2b893ed731` | `771a9c37ef6b5b7f68dd9ae478de9a2b893ed731` | same | Add standalone RADIS line-list table generator for the cluster |
| 11 | 2026-07-31 08:53 | `3efdadd0b14092cacb3b199ac52a9fe6f55958e5` | `3efdadd0b14092cacb3b199ac52a9fe6f55958e5` | same | Add Phase-1 figures F1/F2/F3 + envelope energy / KH time |
| 12 | 2026-07-31 08:57 | `a9e4156a8bb0211169b39cd80b6ff7595d877fa5` | `a9e4156a8bb0211169b39cd80b6ff7595d877fa5` | same | Add G5 literature-anchor harness (Aguichine 2021), skips pending data |
| 13 | 2026-07-31 09:00 | `47f9b82c63f8ef9f480badcd4fcc07595680d918` | `47f9b82c63f8ef9f480badcd4fcc07595680d918` | same | Add boundary sensitivity scans (gamma, matching pressure) -- SPEC 5/10 |
| 14 | 2026-08-01 13:35 | `9d091b4cff41f9750651952ac1d2c49b4fd94625` | `9d091b4cff41f9750651952ac1d2c49b4fd94625` | same | Register POKAZATEL tables; numpy 2.0 compat (A16-A18) |
| 15 | 2026-08-01 13:43 | `d7ad31a7e4b62f5a65ca60eb747b926e1be5f582` | `d7ad31a7e4b62f5a65ca60eb747b926e1be5f582` | same | Lock line-list generator's RADIS block to HitempCalculator |
| 16 | 2026-08-01 13:59 | `f972e777fd799581aabad2c17c8b6f65dd78b711` | `f972e777fd799581aabad2c17c8b6f65dd78b711` | same | Re-run F1/F2/F3 with POKAZATEL; add exact opacity-comparison figure |
| 17 | 2026-08-01 14:10 | `c41f0cad893d23d03be16eb2b36f29ee7eb6ef95` | `c41f0cad893d23d03be16eb2b36f29ee7eb6ef95` | same | Add interim analytic anchor G5a pending Aguichine data (A20) |
| 18 | 2026-08-01 14:54 | `f1f6f59ffd93508b1b2c8a5017a33aa5cf3d9ff9` | `f1f6f59ffd93508b1b2c8a5017a33aa5cf3d9ff9` | same | Per-comparison opacity bands + Fx ratio CSV & crossover cross-check (A21, A22) |
| 19 | 2026-08-01 15:57 | `9f63e1c2c79799e38ece49c1c8c73c9dd53d8c59` | `9f63e1c2c79799e38ece49c1c8c73c9dd53d8c59` | same | Fixed-entropy F3 headline: stay-puffy-longer via build_at_entropy (A23) |
| 20 | 2026-08-01 21:05 | `0c3628aadf9ea6727d0fa72b5fae3d1a3dfd5141` | `0c3628aadf9ea6727d0fa72b5fae3d1a3dfd5141` | same | Add Phase-2 thermal evolution + gate G6 (A24) |
| 21 | 2026-08-03 21:21 | `e594c67a9efc474aabeedce0e06c9ee0cad84ed7` | `e594c67a9efc474aabeedce0e06c9ee0cad84ed7` | same | Resolve Rosseland discrepancy (band, not floor); dilute-flat + radiative reliability |
| 22 | 2026-08-03 21:52 | `a737e266b84aa3ea4a15ea998d28e6c653c24c8f` | `a737e266b84aa3ea4a15ea998d28e6c653c24c8f` | same | Warm-start the R_p shooting via secant seed (1.5x on sequences) |
| 23 | 2026-08-03 21:57 | `490a0e151e1ecb9cf8fb0bea316373a1d333c9b3` | `490a0e151e1ecb9cf8fb0bea316373a1d333c9b3` | same | Document FP-validity finding: FP doesn't affect these planets (A27) |
| 24 | 2026-08-03 23:28 | `c28af7cf1356979061e544b098a2a5fe27edf8c5` | `c28af7cf1356979061e544b098a2a5fe27edf8c5` | same | FP-relevant domain map + cooling-history study (A28/A29) |
| 25 | 2026-08-04 09:29 | `a38f124b56f032d191f4636ccff0699d307faa05` | `a38f124b56f032d191f4636ccff0699d307faa05` | same | Robustness checks: radiogenic floor truncates FP regime; topology is physics (A30/A31) |
| 26 | 2026-08-04 09:53 | `b2c3037169a70dede4f213d017e08a5656d6911d` | `b2c3037169a70dede4f213d017e08a5656d6911d` | same | Deck figures with corrected framing + radiogenic-floor annotations |
| 27 | 2026-08-04 11:51 | `e9e2ed67d4d485cb9afdb3fc10e8d5637e740f57` | `e9e2ed67d4d485cb9afdb3fc10e8d5637e740f57` | same | Implement gate G5 vs Aguichine 2021 fit; sign-offs (A32-A34) |
| 28 | 2026-08-04 12:14 | `6f430985b91ca8a6a3eddaaececd7aa996693841` | `25cdad5fee34d4a76d3e534facad45ab5ef5d89e` | rewritten | Add deck figures M1-M4 + DECK_OUTLINE.md |
| 29 | 2026-08-04 14:53 | `a8911f8c8ca2fe4efa67eea875cf9352a4a154ca` | `4135b78c3ff31dd9599383c9e0f2000e69aaa9ed` | rewritten | Build reader-mode slide deck (17 slides) from DECK_OUTLINE.md |
| 30 | 2026-08-04 16:14 | `89a3ad5aa72334a852cb5802c1b7042a3e7ad08e` | `5b737fdc743d0635adcf4f770e83c86b3c821280` | rewritten | Regenerate F1/F2/D1/E2/Fx to reader-deck specs; rebuild deck |
| 31 | 2026-09-13 16:40 | `124874f90f8b53bc07191ff331142cc401a69e80` | `c036bdbbc6453ad0290dbc9f5106683d25c4c450` | rewritten | Add blend-window and blend-mode options with the S1 sensitivity scan |
| 32 | 2026-09-13 16:54 | `b3eefa40232c242e8b21664b95bf13be93f67d5f` | `1ef8bc3b08dccda0599003738cc2ba520e9c6013` | rewritten | Add fp_scale opacity perturbation with the S2 robustness scan |
| 33 | 2026-09-13 17:00 | `09987b110c5bb10c3405fb14eb3e1a5beb074871` | `0ba60479b3376512123989fc6d2b0633673f6f8b` | rewritten | Quantify HITEMP band truncation and its effect on the G5 anchor |
| 34 | 2026-09-13 17:02 | `7c413dcf02d09d0a263302b5d60cd072dc6785f1` | `4e8f754ae46d98400e8a041e6dd899e4d17fc81f` | rewritten | Record the measured band truncation in the limitation doc |
| 35 | 2026-09-13 17:13 | `807bbd94fe077f5c250281b13bcba75af475dff5` | `fc651f02f84a112c86fe6cdba59b694c7a9d874e` | rewritten | Report the escape-window null and the F3 cooling-time range |
| 36 | 2026-09-13 17:16 | `1a717c8ea6235fc4498630515221397f7274e607` | `81337fa669e107d14af003c0e1c1b78d0140da2a` | rewritten | Correct two measured numbers in the new amendments |
| 37 | 2026-09-14 12:09 | `b00dcc2ff17ab1155b7d20ed6c9018d2a4b76408` | `1c71a92f77be776377a58a706b2e79282ce88d97` | rewritten | Audit where the deep-closure opacity actually comes from |
| 38 | 2026-09-14 12:17 | `f91c03fecc1158bbc55c57de4bb64760eaab4a78` | `e3585c3fc856ae7fb7c84c5eaf221a208c29feab` | rewritten | Scan for line-list radiative zones inside the FP data box |
| 39 | 2026-09-14 12:23 | `030eb1050600451434cb7badaae15d07b2c14380` | `f5237ba7d66f7f8838a16d956fb7b47674be97f7` | rewritten | Evolve the candidate fiducial and measure its cooling-time ratio |
| 40 | 2026-09-14 12:28 | `4834c2267d69254fe79e8118c9c1ec83c7b8829e` | `9a880ce7624b582daedf4983321e9566c6c0fb39` | rewritten | Report the radiogenic T_int floor at the candidate core |
| 41 | 2026-09-14 12:28 | `a2d463216dad3746a42fdf6936111b6d4f14a129` | `5ef71ec2f1da3d5d292f276f29d41d20a6ed82f6` | rewritten | Audit the G5 anchor definitions before blaming physics |
| 42 | 2026-09-14 12:36 | `4bb9570be2be16ae8c0291636af456d756fa86d4` | `5685c2d0b37157736f12bec3f64508bc63fac7ca` | rewritten | Report the fine-grid S6 numbers and pin the G5 stellar attribution |
| 43 | 2026-09-14 16:02 | `4d6415f4aef1e12a7540dbc372637ee92fc68309` | `1da16ed80bbfbe37878b1709a5ca2df75ef745e7` | rewritten | Adopt the corrected G5 anchor conventions as defaults |
| 44 | 2026-09-14 20:34 | `0cbfc6276ae4886af28896ef87f70f36003bf0f5` | `ba9f366f7d7ad9132f9ab8b181f9aabb508dadc0` | rewritten | Scan the extended planet domain for an in-box radiative layer |
| 45 | 2026-09-14 20:41 | `4c022d757e2f90b9dd0a97236ce56e9699e50b1d` | `60bd3fd0048516b51d40c2bb691d181ea705bae9` | rewritten | Map the FP grid's coverage of envelope trajectories |
| 46 | 2026-09-14 20:42 | `9b0455a405241fd33f419cbfcbcdfc2ad865ca6b` | `bf772158bf642da858f436d85faf807254ec8c5f` | rewritten | Measure how far each RCB is from flipping, in units of opacity |
| 47 | 2026-09-14 20:45 | `60163fb5f7c6239ba02859a3efa606d0bf8b2749` | `dc2051187b7c32acc14080055915a98e4d2ce8fb` | rewritten | Correct the A49 count, the scan timing and stale-CSV handling |
| 48 | 2026-09-14 20:51 | `8134199d09dbc105b17c452ef9c3a9d6d01b1a15` | `e189eded6f4f5ae1cd0100ed775ceaf670e906c7` | rewritten | Annotate the scan cells with their closest approach to the FP box |
| 49 | 2026-09-20 17:15 | `1f95b93703879fbfb37084f42a9a13cb95b2779e` | `b15dbc5f6c9d2c9d0809104990de112d0457a3b1` | rewritten | Fix the FP opacity units at source and merge the extension grid |
| 50 | 2026-09-20 17:17 | `6aa866daa8ccf90c8fb7da984e36027566cfe5b3` | `7444a4a2a81c4591c46005b56d4c9def3ac41ae5` | rewritten | Check the dilute FP spectra against the dipole-model noise floor |
| 51 | 2026-09-20 17:22 | `33a47820c8573e35351d7bedee5c5a86ada0285f` | `0d952860372b4b677dd66391ce909c4a48244640` | rewritten | Rebuild the FP/HITEMP crossover on all 48 grid points |
| 52 | 2026-09-20 18:05 | `0ed64e7f80077132c3c1a8d56651c3741608d6ad` | `75470c36f6a970f453093b4275b49f4f338c60af` | rewritten | Rescan the domain and remeasure the RCB dispute on the enlarged grid |
| 53 | 2026-09-21 08:51 | `aa3c6049336360d46cda18185174114f1cf2250b` | `75a9acd6b039cc361c8a2af15ad1426a1da27e28` | rewritten | Correct the FP grid counts and floor in the docs |
| 54 | 2026-09-21 14:02 | `efa70527ff04d146ee090810b2af974f1dad1776` | `d94046195d7680f0449a0b4099e05e4a52c8dc56` | rewritten | Compare cooling at the strongest cases and stop a structure build hanging |
| 55 | 2026-09-23 17:03 | `95d6b0e71f1e5530174497a604b4f526db8a6417` | `d33c1eb0292fe8c6d1635b1897121902e26a4022` | rewritten | Add the MT_CKD 4.3 self-continuum as a fourth opacity source |
| 56 | 2026-09-23 17:06 | `5aeccf605b5763a1b6bcae01c3a88391ad95633e` | `f74dec26b8d17921aa9de58a502d55bbed75dbb5` | rewritten | Report Fig. 2b against both references and rebuild it on 48 points |
| 57 | 2026-09-23 17:13 | `9fabe9234d414c8fb3eb78bf3aa359e0c8057606` | `75b4441f3c810da4a749aad5e2e975b269455de5` | rewritten | Gate the continuum's density scaling against cntnm at envelope densities |
| 58 | 2026-09-23 17:18 | `01eeb03f510f5b84842f3532d75f3aa27b82f867` | `16b6b494cd1c738afa87ffe43b2e1b6d15a32802` | rewritten | Correct the clamped-point continuum scale and set up Task 21 |
| 59 | 2026-09-23 17:21 | `83c18c82f4dad25b3a1349020fe4c696060073e6` | `0056e529fc781de965f3c64185956d1549788a0e` | rewritten | Keep one failed RADIS point from deleting a whole series in ED1 |
| 60 | 2026-09-23 17:26 | `621c7af4910d190ca2b84942f99a2bdb9476938c` | `97e128eacf53c8375db469bc6c0b477c5db5de95` | rewritten | Record where the line/continuum decomposition stops being well posed |
| 61 | 2026-09-23 17:28 | `a7c1482d31251b12df630df73cc06135cdf50153` | `58a2af9d31f953c0209a23fef830c9c5f0fee715` | rewritten | Say which exponent the half-width estimate uses, and bound the share both ways |
| 62 | 2026-09-23 17:29 | `dee3d3bad24feb67c5253f49bab4718a75a6bdfb` | `dab442e8c0262aab8972ad91bef9a87f60ca0860` | rewritten | Bound the cold-dense continuum share from above as well |
| 63 | 2026-09-23 17:30 | `85a34d0b8d213eb51e931f886c3ddbb07b674ccd` | `573fd7a9bfabff8168cc81993bc7847e41ef43a2` | rewritten | Record the first Task 21 measurement in A59 |
| 64 | 2026-09-23 18:03 | `08ede78287600c0aabe996cfa0da161c80610dae` | `bd0bd5c4372028e53983830b908259c0dc08ddd3` | rewritten | Complete the first Task 21 point in A59 |
| 65 | 2026-09-23 18:16 | `97edc5a8f2011dee0c11f1ae87100a5c80d857cc` | `1b787938a75afa8b6970c4537ed220f16a4d4027` | rewritten | Note the first direct test of the opacity interpolation |
| 66 | 2026-09-23 18:34 | `604bb7e6dc6f564ede4334869e5a395ef44d80de` | `aa53d25bbf9548cb7b09810c2a76df899966a360` | rewritten | Correct the truncation claim: RADIS renormalises, it does not discard |
| 67 | 2026-09-23 18:36 | `302e8dfaade1d197a1edff57ae181ac3bac49833` | `0a93e9bbe97aa430fd61003610f0761201dac298` | rewritten | Record the non-monotonic pressure response and the T-interpolation check |
| 68 | 2026-09-23 18:41 | `783ff3de1e6b43f766c08114233213b9b9c05505` | `580310ba47f0fc354ad8a99c761f240752ecc30f` | rewritten | Complete (650 K, 300): the pressure response has an interior minimum |
| 69 | 2026-09-23 18:58 | `1dddc20fc0c27a525f994db1838e140479e91ba5` | `1a19b4adc74882674a40caf048b5cd70de6f4630` | rewritten | Add the (1500 K, 471.3) rows and the window-depth pattern |
| 70 | 2026-09-23 19:05 | `518f6fcc5584a3a720a0f2570cadd1c470454e81` | `867cd428a94a68266dc6c911060c093bd84982df` | rewritten | Complete (1500 K, 471.3): ideal-vs-AQUA is +1.1% there |
| 71 | 2026-09-23 19:35 | `a326ace202c26c9e0d062bcd0ad06e8cfd66f10c` | `9247c4331f50d07f638a8bdceccfe250973835d1` | rewritten | Record the (2000 K, 1200) rows and the 1 bar resolution caveat |
| 72 | 2026-09-23 19:39 | `e9ff78280ac883319652ec4bcf405f627b334b43` | `c239ccacd83b2362b58b9e1280383b59fca5683a` | rewritten | Withdraw the 'AQUA pressures would not move the tables' conclusion |
| 73 | 2026-09-23 20:06 | `aa20d3185c6848cb9784286ce6b345b8bc5c2c08` | `736b3c6262064f7fafda5083528f8aa298974942` | rewritten | Measure the density clamp; A53's 30% was mis-attributed |
| 74 | 2026-09-23 20:17 | `38ac4c80bb43651850270e96856d650c437a5383` | `67c5fa3130cf033b5b5c8f47d3f52508c969c9b8` | rewritten | Task 21 complete: all 18 spectra, no failures |
| 75 | 2026-09-23 20:18 | `b6861e71ac374fbbd0bb775a563121532811595c` | `553a2f32b5d90b2c2c73059271da5a635ca3b623` | rewritten | Anchor both ED1 panels at unity instead of zero |
| 76 | 2026-09-24 11:16 | `00398fdceda0dd51c31c8b0897be9a7547fba178` | `f37a4d0f959fd8f595a462c07cd602d335d4d702` | rewritten | Split the dilute FP/line-list ratio band by band (Task 22) |
| 77 | 2026-09-24 11:36 | `e5bd46d2bb4dcc2f93139fd1f89ffb93b53ba2f0` | `187b8b2be018b1aedec9111d5d266e484a65bb20` | rewritten | Extend the HITEMP grid, regroup the points physically, add the manuscript figures |
| 78 | 2026-09-24 11:38 | `4fc85e51f9b05f6f789e03c9fd99750b22acbb67` | `9821f6116ecc9d7dbde68c3a15e270aae6e25c7f` | rewritten | Gate the manuscript figures and the one-grid rule |
| 79 | 2026-09-24 11:41 | `8ed0b7c33594c0e50e4862853fb190a25df0d482` | `154599171b628d5800afe640c53a8ff14f81d99d` | rewritten | Update the band gate for the extended grid, and drop a campaign reference |
| 80 | 2026-09-24 11:43 | `73a2b560275d4b2edf96847c4cde18e179681419` | `0f078644e0838588a1d8efc5432840d6224f87f7` | rewritten | Measure, rather than assume, the drift in the stored boundary inputs |
| 81 | 2026-09-24 11:54 | `0df107c5b7541ac4343d85dafd0596b57fe92b51` | `7ab11fccbffe798e7bd3bec3e0cb3e460439953f` | rewritten | Refuse to write a figure with text below 7 pt |
| 82 | 2026-09-24 12:08 | `6f092c0e4cca13b633505ce388135be2ad0f57de` | `16a1300906239a533c4bbd70fb64de54b91987c4` | rewritten | Make the figures actually 183 mm wide, and gate it on the file |
| 83 | 2026-09-24 12:13 | `fdcf0d16b2f192398500d32d67d2348c81c1ab9b` | `d50d52c673c9eb73d4f8e09ea1c6111199716186` | rewritten | Write the numbers ledger in a deterministic order |
| 84 | 2026-09-24 12:18 | `0adaad5eb7c3b86d4900297f35cfeeb4bde71f68` | `1eccb3dae9c8eeb5b761b0a98fb6aae325c8e088` | rewritten | Record that the two line lists differ in isotopologue coverage |
| 85 | 2026-09-24 13:07 | `9076270d6cc5f0697e912a5c2728adc08954dcf0` | `bfc722c38664b2eee305fbcc5353077f9d6b7b5b` | rewritten | Task 24: ARC's MT_CKD reference is bit-for-bit the file A58 built |
| 86 | 2026-09-24 14:11 | `badbdeb7f81c2b573c08415d71a55fb85a6691c9` | `e9818d7a03282781308288ef62ae8665c8938ee9` | rewritten | Task 23 complete: the 2250/2750 K rows close the A53/A59 puzzle |
| 87 | 2026-09-24 14:18 | `a821e220dec9578eb2cb16a7352ca3bebbb8e828` | `6517c08f92a50f7abc0d49264ef931dce8154b8d` | rewritten | Update three gates the extended grid invalidated |
| 88 | 2026-09-24 15:28 | `75e98a91aa01736ddc12ea11af5a301059d986de` | `402f5a91e0749b40019d602a8bbf30665a563c5c` | rewritten | Add the isolated-molecule dipole test bundle (Task 26) |
| 89 | 2026-09-24 16:01 | `984b2f10a23cded6f78f1fbb8aee0f949ff834f2` | `3cd7afa5cc08d967a808e8e1b53481f68d54a2bf` | rewritten | Task 28: confirm the window value, extend the band table to rho = 60 |
| 90 | 2026-09-24 16:08 | `62eb635d47d1cc685ed53da67d9c560095a35937` | `1241207097534410d9db4ddb1733bd98f97d5081` | rewritten | Task 27: the manuscript-numbers ledger, 303 entries, nothing missing |
| 91 | 2026-09-25 19:41 | `0b79b3c2bf2f739b5420e939c773eb4a8ca9a973` | `3a7ccb0906c62ca6f6c4392bd3402f99cc164185` | rewritten | Extended Data: 8 figures, 2 tables, and the Fig. 4 categories |
| 92 | 2026-09-26 12:26 | `162624b73b1c00aeaae5f9555309391ac077ec0f` | `2e4692eb9e53b0be5fecb184099f1d5cd88e28f3` | rewritten | Housekeeping: track the ARC pull and the model binary; rename one ledger key |
| 93 | 2026-09-26 16:11 | `024e7f92acc525b4ec89f748a531543e55989afa` | `56b08a9b04aafe5797e21ddb1a11933b6c3adf8b` | rewritten | Tasks 38-41: the rho^2 scaling test, the wing cut-off, the AQUA entropy, release tables |
| 94 | 2026-09-29 11:08 | `90032e5d9c047e3a866ce615df6c43b009b24d96` | `3a95f489a03f70649ae0bfb9c856b31d0bf28a5e` | rewritten | Tasks 42-44: four corner cells merged, the grid is 52 points |
| 95 | 2026-09-29 13:04 | `1ef8bd99f92b0b6237d930fd593c22058280d736` | `56533ea85591ae9071eec74cb002d78f7e154237` | rewritten | A75: the dilute noise floor, measured from the lag test rather than assumed |
| 96 | 2026-09-29 13:26 | `800ee75b27b539cfd8d2024bfb21071bf681ba42` | `49453b5a59079a07435ce1faf783dee6466eaa7a` | rewritten | A76: the quantum correction is not the free parameter |
| 97 | 2026-09-29 14:55 | `4bf9407a042f562d6b1c82fc6ad658c87e5e2d55` | `14f8f8d7a9884d2b030e1c5b092225c5dce4ad2d` | rewritten | A77: trace the sum-rule violation — grid-wide, and half of A76's index was unverified |
| 98 | 2026-09-29 16:17 | `0522b2bd73a49bdf4d1cb345b076ecaac3e957f0` | `0db8c8b4597d7680a65d6b417f2e3e06a79f0352` | rewritten | A78: the production formula applies the harmonic correction twice |
| 99 | 2026-09-29 16:50 | `240c8583dad79dbe35984f3c11434a003d9b85ee` | `19a6ee0632de0778555765aafd9de325d8c11d4f` | rewritten | A79: the normalisation is correct; the estimator is not; A80 halted |
| 100 | 2026-09-29 19:41 | `c497ffac6e2b51aa0c64d2c8f95775aaf3ffc9f1` | `780cca746e67361ec6b7900a14f73d612672789c` | rewritten | A80: estimator replaced and gated; propagation held at the noise floor |
| 101 | 2026-09-29 19:42 | `66217fae96eb806542d69e8df72746cecc0f36dd` | `9a822e9a14e976a986067c3a63930c268ac48aec` | rewritten | A80 housekeeping: keep the raw dipole bundle out of git, skip its gates without it |
| 102 | 2026-09-29 20:31 | `8151a4cf979a69403218e889e57607f9aa6889d3` | `1f593b46401149d2328dcc69f2fb7538c8073e67` | rewritten | A81: there is no noise floor -- it was spectral leakage; A80's blocker is lifted |
| 103 | 2026-09-29 21:21 | `f013f5bbe0fdc917c3b30daa8139e7a985b8a305` | `18a355eccc67db26de9151ed1517a22de3c7087a` | rewritten | A82: the smoothing width cannot be chosen from a 30 ps record; stopping before propagation |
| 104 | 2026-09-29 23:44 | `f068649d18ba7df958e7accee17ff08cd79384bd` | `db950469552df03cebc205f3187970df9961591c` | rewritten | A83: v5 is live -- estimator propagated; three physics claims moved, none re-pinned |
| 105 | 2026-09-30 08:02 | `9aefb62bddc7d44b280bc9d2a08d38c115f317bf` | `e353eb5161664165eed542fe58984e1854b7372a` | rewritten | A83 outputs: the propagated products, and ED Table 1's uncertainty columns |
| 106 | 2026-09-30 09:42 | `7979f9a4100ed3b747d7ae5262343b9755c4a3e4` | `0d3c70fdaf52b177e4b559e1881ac1c7dab95577` | rewritten | A84: the 6.07 was mine and the v3 baseline was stale; v5's real change is downward |
| 107 | 2026-09-30 11:20 | `c98562ff3cf3972c755a4be5c682a26a435b9432` | `9d60e7a2802eb69a6f856dc6674819ca2b57fa45` | rewritten | A85: the ledger refuses stale inputs; the definitive table is blocked on one case |
| 108 | 2026-09-30 14:11 | `16b1a08844c0c1c30e941aaa7d14dc4b4606332c` | `fbd1180c761cea15ad0403e5d9566d47bc62f447` | rewritten | A86: the stuck case was the critical point; the definitive table at last |
| 109 | 2026-09-30 17:09 | `797ca9565c16152b8dc428a479366fa2f7147f42` | `ab6e953171b60c5c15b5fd0c47508c900fef6223` | rewritten | A87: A86's table was still v5 against v5 in every figure group |
| 110 | 2026-09-30 21:19 | `5470506750b1c9a7ffdf0035109c839721ef8f20` | `8a34048a0309940a69492f83082ec0f7ad2667d1` | rewritten | A88 parts 1-3: window validation reduced, the bound measured, T_int against age |
| 111 | 2026-10-01 14:12 | `d118fad3a1aebba18cfa372da264bb0839aaef9f` | `e85b2c5c34ea968b4ded7bfe242d5c5145d690df` | rewritten | A89: the rho^2 claim survives v5; A88 finished; four groups move |
| 112 | 2026-10-01 18:49 | `ffc1b64723ba95f09f079a9057229a4c94e3c996` | `aeac5d2b64526147f9a4f281ea76422047d7cea2` | rewritten | A90: the pooling assumption holds where it matters; the held-out errors found |
| 113 | 2026-10-01 20:14 | `691fee022d7fba7cd4f8793026f33f0236fc4350` | `348cda7a720898ee0d95e1e0d78c685a74480b1b` | rewritten | A91: the notch is the wing cut-off, and it does not converge |
| 114 | 2026-10-01 20:32 | `1f225e84a7cd614dc4ca4853513b3276dfbeca07` | `2cc6fc94c25882fdea097c2e85fb1e6e918ff6d4` | rewritten | A92: Fig. 1b shows both cut-offs and quotes medians, not one channel |
| 115 | 2026-10-01 21:32 | `344f2a079dcdc5cd894d2c9f043a8c34f98bd3ad` | `0c29371d5628a6c8f5c5ce651d7ea427999f057d` | rewritten | A93: three figure labels, and the colour-bar extension checked not assumed |
| 116 | 2026-10-22 17:54 | `dc8d47cdf6fd8124d5312098aa0f99a1f6c3ebd6` | `588318a16f562e7c934860ed65aeb7f9ff4965e3` | rewritten | A94: the MLIP error table is the production model's (stage one), not SWA |
