# Physical constants

All constants are defined once in `constants.py` and sourced from
`scipy.constants` (CODATA 2018) where available. Astronomical/planetary
reference values not in `scipy.constants` use IAU 2015 nominal values and are
cited below. Constants are never re-hardcoded elsewhere in the code.

## From `scipy.constants` (CODATA 2018)

| Symbol      | Name                         | `scipy.constants`      | Value                          | Unit        |
|-------------|------------------------------|------------------------|--------------------------------|-------------|
| `G`         | Newtonian gravitation        | `G`                    | 6.674 30 × 10⁻¹¹               | m³ kg⁻¹ s⁻² |
| `SIGMA_SB`  | Stefan–Boltzmann             | `Stefan_Boltzmann`     | 5.670 374 419… × 10⁻⁸          | W m⁻² K⁻⁴   |
| `K_B`       | Boltzmann                    | `k`                    | 1.380 649 × 10⁻²³              | J K⁻¹       |
| `H_PLANCK`  | Planck                       | `h`                    | 6.626 070 15 × 10⁻³⁴           | J s         |
| `C_LIGHT`   | Speed of light (vacuum)      | `c`                    | 2.997 924 58 × 10⁸             | m s⁻¹       |
| `N_A`       | Avogadro                     | `N_A`                  | 6.022 140 76 × 10²³            | mol⁻¹       |
| `M_U`       | Atomic mass unit             | `m_u`                  | 1.660 539 07 × 10⁻²⁷           | kg          |
| `BAR`       | Bar                          | `bar`                  | 1.000 000 × 10⁵                | Pa          |

## Astronomical reference values (IAU 2015 nominal)

Not provided by `scipy.constants`; taken from IAU 2015 Resolution B3 nominal
values. Earth mass is derived from the nominal terrestrial mass parameter
(GM)⊕ divided by CODATA `G`, so it stays consistent with `G` above.

| Symbol     | Name                              | Value               | Unit  | Source                         |
|------------|-----------------------------------|---------------------|-------|--------------------------------|
| `GM_EARTH` | Nominal terrestrial mass param    | 3.986 004 × 10¹⁴    | m³ s⁻²| IAU 2015 Res. B3               |
| `M_EARTH`  | Earth mass = `GM_EARTH / G`       | 5.972 2 × 10²⁴      | kg    | derived                        |
| `R_EARTH`  | Nominal Earth equatorial radius   | 6.378 1 × 10⁶       | m     | IAU 2015 Res. B3               |

Note: deriving `M_EARTH` from (GM)⊕/`G` (rather than adopting a fixed mass)
means the mass carries `G`'s relative uncertainty (~2×10⁻⁵), which is
negligible for structure work and keeps GM self-consistent.

## H₂O molecular data

| Symbol       | Name                       | Value              | Unit    | Source                    |
|--------------|----------------------------|--------------------|---------|---------------------------|
| `MMW_H2O`    | H₂O mean molecular weight  | 1.801 527 × 10⁻²   | kg mol⁻¹| AQUA table (dilute limit) |

Matches the AQUA table's `mmw` column in the neutral/undissociated limit.
