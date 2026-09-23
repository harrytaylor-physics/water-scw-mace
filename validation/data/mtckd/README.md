# MT_CKD reference outputs (gate fixtures for A58)

`WATER.COEF.T{296,650,1000,2000}` are outputs of **AER's own compiled MT_CKD 4.3
continuum program**, `cntnm`, not of anything in this repository. They are the
independent reference `validation/test_mtckd_continuum.py` gates `opacity/mtckd.py`
against, so they are committed rather than regenerated.

Each file's long header is AER's provenance/changelog block, verbatim. The data
block is five columns on AER's own 10 cm⁻¹ grid:

    wavenumber [cm^-1]   self_no_rad   foreign_no_rad   self_with_rad   foreign_with_rad

with the coefficients in cm² molecule⁻¹ at the reference density (1013 mb, 296 K).
The gate uses columns 1, 2 and 4.

## How they were produced

Source: <https://github.com/AER-RC/MT_CKD>, release 4.3, built for macOS/gfortran
(single precision) as `cntnm_v4.3_OS_X_gnu_sgl`, with AER's distributed reference
table `data/absco-ref_wv-mt-ckd.nc` (sha256
`69944eb8b045c268e2daeb2cddf536b99f9e9efe7e91c425067b46a5b287ddb3`) present in the
run directory.

`cntnm` reads its state from **stdin**, not from the `INPUT` file its example ships
(the `INPUT` file is ignored — a trap worth recording). Pressure, temperature, path
length and water VMR, one per line:

    for T in 296 650 1000 2000; do
      mkdir -p run_$T && ln -sf .../absco-ref_wv-mt-ckd.nc run_$T/
      printf '1013.0\n%s\n1.0\n1.0\n' "$T" | (cd run_$T && .../cntnm_v4.3_OS_X_gnu_sgl)
      cp run_$T/WATER.COEF WATER.COEF.T$T
    done

VMR = 1 (pure steam) is what the envelope is, and it also zeroes the foreign term
by construction — AER's foreign coefficient carries a factor (1 − x_H2O).

Note that the *printed* self coefficients are at the reference density: they carry
the `(296/T)^n(ν)` temperature dependence and (in column 4) the radiation term, but
no density scaling. That is exactly the quantity `opacity/mtckd.py`'s
`self_continuum_coefficient` returns, which is why the comparison is direct.

The 296 K file alone would be near-tautological — it is the state the reference
table is defined at. The 650/1000/2000 K files are what make the gate bite.

## The density-scaling fixtures

`CNTNM.OPTDPT.T1000_rho129.744` and `CNTNM.OPTDPT.T650_rho300` are the *optical
depth* output of the same binary at two envelope-like states, and they gate the
part the `WATER.COEF` files cannot: the ρ² density scaling, at **175×** and **404×**
MT_CKD's reference density.

`cntnm` takes a pressure, not a density, and converts it with its own ideal-gas
Loschmidt constant, `WTOT = ALOSMT·(P/1013)·(273/T)·L` with
`ALOSMT = 2.68675e19`. So the pressure fed to it is the one that reproduces the
cell's *number density*, `P = n k_B T`:

| T [K] | ρ [kg/m³] | n [cm⁻³] | P fed to cntnm [bar] |
|---|---|---|---|
| 1000 | 129.744 | 4.3371e21 | 598.798 |
|  650 | 300     | 1.0028e22 | 899.970 |

With VMR = 1 and a 1 cm path the printed optical depth is the self-continuum
alone (`W_dry = 0`, so no foreign and no O₂/N₂/CO₂ term), i.e. exactly
`kappa_cont[cm²/g] · ρ[g/cm³] · L[cm]`. Agreement is 8e-4, and the residual is
explained by `ALOSMT` differing from the CODATA value used here by 0.031% — which
the gate asserts rather than absorbs into the tolerance.

    for s in "1000 598.798" "650 899.970"; do
      set -- $s; mkdir -p dens_$1 && ln -sf .../absco-ref_wv-mt-ckd.nc dens_$1/
      printf '%s\n%s\n1.0\n1.0\n' "$(python -c "print($2*1000)")" "$1" \
        | (cd dens_$1 && .../cntnm_v4.3_OS_X_gnu_sgl)
    done
