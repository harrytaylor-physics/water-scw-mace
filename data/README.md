# data/

Derived data products of the hybrid (revPBE0-D3) campaign. Raw production
output — trajectories and dipole series — is not in the repository; it is
fingerprinted instead (below).

| Where | What | Format |
|---|---|---|
| [fp_opacity_v5/](fp_opacity_v5/) | the 52 first-principles absorption spectra, v5 estimator (650-3000 K, 10-2229 kg m^-3) | `.npz` (git-LFS), one per state point, + `grid.json`, `crossover_hitemp.json` |
| [release/](release/) | Rosseland and Planck mean opacity tables on a regular (log R, log T) grid, from the spectra above | MESA `kap` plain-text format, CSV, fit JSON |

Both are written by code in [`subneptune/`](../subneptune/):
`tools/v5_regenerate.py` writes `fp_opacity_v5/`, and `tools/release_tables.py`
writes `release/`. The locations are set in `subneptune/layout.py`. Each can be
redirected with `SUBNEPTUNE_FP_OPACITY_DIR` / `SUBNEPTUNE_RELEASE_DIR`.

The spectra are git-LFS objects (`.gitattributes`). Run `git lfs install` once
before cloning, as for `models/`.

## Fingerprints of what is not here

| manifest | covers |
|---|---|
| `subneptune/ed_data/arc_dipoles_sha256.json` | the 67 production dipole series behind the 52 spectra (41 MB) |
| `subneptune/ed_data/arc_a78_dipoles_sha256.json` | the five cells subneptune A78 reproduces: their dipole series and the v3 spectra the original pipeline delivered from them |

Gates that need these files skip cleanly when they are absent. Place the bundle
at the manifest's location to run them; the A78 gate checks the hashes.
