# docs/figures/subneptune/

Figures of the planetary model (`subneptune/`). Each is a vector PDF at final
size (183 mm wide), with a 200-dpi PNG of the same canvas for display. All are
generated, never edited by hand. `save_pdf` (`subneptune/figures/_manuscript.py`)
writes both after its checks:
- the 52 points are one grid, so no figure text may call any of them "original",
  "extension" or "campaign";
- nothing is set below the minimum type size.

Regenerate, from `subneptune/`, with `bash tools/rebuild_producers.sh`, or one
item at a time as below. Every number a figure draws is recorded in a ledger:
`fig_numbers.json` here for Figs 1-4, and `ed/ed_numbers.json` for Extended Data.
`tools/manuscript_numbers.py` and `figures/ed_index.py` copy them into the numbers
ledger under the keys listed. That ledger is `subneptune/outputs/`, which is not
committed.

## Main figures

| file | content | producer | ledger keys |
|---|---|---|---|
| `fig1.pdf` / `.png` | the first-principles pipeline, and the spectrum it produces against the line list at two wing cut-offs | `figures/fig1.py` | `figures.fig1.*` (a-c); the stage captions are `pipeline.*` |
| `fig2.pdf` / `.png` | the first-principles opacity surface, and how it differs from line lists | `figures/fig2.py` | `figures.fig2.*` (a-c) |
| `fig3.pdf` / `.png` | what the opacity difference does to a radiative-convective boundary | `figures/fig3.py` | `figures.fig3.*` (a-c) |
| `fig4.pdf` / `.png` | where the opacity disagreement lands relative to real envelopes | `figures/fig4.py` | `figures.fig4.*` |

Run each as `python figures/figN.py config/fiducial.yaml`.

## Extended Data (ten items)

| file | content | producer | ledger keys |
|---|---|---|---|
| `ed/ed_fig1` | validation of the interatomic potential: force parity by regime, NVE energy drift | `figures/ed_items.py` | `ed.ed_fig1.*` |
| `ed/ed_fig2` | Planck-mean opacity over the 52 points; first-principles / HITEMP ratio | `figures/ed_items.py` | `ed.ed_fig2.*` |
| `ed/ed_fig3` | line-list conventions: kappa_R under three broadening pressures, the MT_CKD self-continuum share, band truncation | `figures/ed_items.py` | `ed.ed_fig3.*` |
| `ed/ed_fig4` | density dependence of the absorption: dilute-limit band ratios (a, b); alpha/rho^2, power-law index, mean-free-path map (c-i) | `figures/ed_items.py` | `ed.ed_fig4.*` |
| `ed/ed_fig5` | dipole-model validation: derivatives against revPBE0, model/DFT ratio per mode, residual autocorrelation | `figures/ed_items.py` | `ed.ed_fig5.*` |
| `ed/ed_fig6` | the 108 scanned planets: boundary positions against the computed grid | `figures/ed_items.py` | `ed.ed_fig6.*` |
| `ed/ed_fig7` | external anchor: model radius against the published relation | `figures/ed_items.py` | `ed.ed_fig7.*` |
| `ed/ed_fig8` | evolution of three water worlds: radius, interior entropy, deepest-boundary density | `figures/ed_items.py` | `ed.ed_fig8.*` |
| `ed/ed_table1.csv` / `.tex` | the 52 state points: Planck and Rosseland means, uncertainties, ratios to each line list | `figures/ed_tables.py` | `ed.ed_table1` |
| `ed/ed_table2.csv` / `.tex` | the ten strongest boundaries | `figures/ed_tables.py` | `ed.ed_table2` |

Run `python figures/ed_items.py config/fiducial.yaml`, then `figures/ed_tables.py`
and `figures/ed_index.py` with the same config.
- `ed/ed_index.json` is the index of the ten items.
- The `.tex` tables are generated `tabular` fragments; the `.csv` holds the same
  numbers.
- **Order matters:** run `ed_index.py` *after* `tools/manuscript_numbers.py`,
  because it merges into the ledger (subneptune A85).
