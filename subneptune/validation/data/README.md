# Literature reference data (for gate G5)

## `Aguichine_2021_ApJ_914_84.pdf` — the source paper

Aguichine, Mousis, Deleuil & Marcq (2021), ApJ 914, 84. Their steam-world
mass–radius fit is Eq. 29:

    log10 R_p = a·log10(M_p) + b + exp(−d·(log10(M_p) + c))

with a,b,c,d tabulated per (composition, T_irr). That coefficient table is an
**online machine-readable table, not contained in the PDF**, so we anchor G5 to
points digitized from their Figure 8 instead (see below and AMENDMENTS A32).

## `aguichine2021_fig8_points.csv` — the G5 anchor

Eight (WMF, T_irr, M, R) points read from Figure 8 (Ma19+ model, no metallic
core; read uncertainty ~5%). Columns: `wmf, t_irr_k, m_earth, r_earth`. Gate
`test_g5_literature.py` matches each point (WMF→f_env, T_irr→T_eq=T_irr/√2,
pure-rock core cmf=0, mature T_int=30 K), computes R_p, and writes
`outputs/g5_residuals.json`. Residuals (documented, not tuned): ~6–15% at
WMF=20%, ~18–23% at WMF=50%, our R_p smaller — driven by our pure-steam envelope
vs their full water-phase layers, plus core-composition and thermal-state
differences.
