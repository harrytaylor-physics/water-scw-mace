# Literature reference data (for gate G5)

## `aguichine2021.csv` — NOT committed; provide to enable G5

Static mass–radius benchmark for steam worlds from Aguichine, Mousis, Deleuil &
Marcq (2021, ApJ 914, 84). Not on the machine and not fetchable offline, so
gate G5 (`validation/test_g5_literature.py`) **skips** until this file exists.

Expected CSV columns (header row required):

    m_earth,f_env,t_eq_k,r_earth

- `m_earth` — planet mass in Earth masses
- `f_env`   — steam envelope mass fraction (their water mass fraction)
- `t_eq_k`  — equilibrium temperature [K]
- `r_earth` — published planet radius [Earth radii]

Source options: the paper's machine-readable tables, or the authors' grid /
`MSEI`-style output. A handful of points spanning the fiducial box (M ~ 1–10
M⊕, f_env ~ 1–10%, T_eq ~ 500–800 K) is enough for the ≤5–10% check.

G5 also wants a *physically valid* opacity for a fair comparison; with only FP
(flagged unreliable in the radiative zone) available now, run G5 against the
`composite` source once the HITEMP tables land. The residual report records
which opacity was used.
