import os, glob
from ase.io import iread, read, write

HOME   = os.path.expanduser("~")
SNAPS  = os.path.join(HOME, "scw", "snaps")
FFTEST = os.path.join(HOME, "scw", "fftest")
PROD   = os.path.join(os.environ["DATA"], "scw_prod", "T3000_rho108.796_05fs", "prod.xyz")
os.makedirs(SNAPS, exist_ok=True); os.makedirs(FFTEST, exist_ok=True)

# (re)build snapshots if missing
have = sorted(glob.glob(os.path.join(SNAPS, "snap_*.xyz")))
if len(have) < 5:
    K, j = 2000, 0
    for i, at in enumerate(iread(PROD)):
        if i % K == 0:
            write(os.path.join(SNAPS, f"snap_{j:03d}.xyz"), at); j += 1
    print(f"wrote {j} snapshots")
    have = sorted(glob.glob(os.path.join(SNAPS, "snap_*.xyz")))
print(f"{len(have)} snapshots in {SNAPS}")

REVPBE0 = """&GLOBAL
  PROJECT {proj}
  RUN_TYPE ENERGY_FORCE
  PRINT_LEVEL LOW
&END GLOBAL
&FORCE_EVAL
  METHOD QUICKSTEP
  &DFT
    BASIS_SET_FILE_NAME  BASIS_MOLOPT
    BASIS_SET_FILE_NAME  BASIS_ADMM
    POTENTIAL_FILE_NAME  GTH_POTENTIALS
    &MGRID
      CUTOFF 600
      REL_CUTOFF 60
      NGRIDS 5
    &END MGRID
    &QS
      METHOD GPW
      EPS_DEFAULT 1.0E-12
    &END QS
    &AUXILIARY_DENSITY_MATRIX_METHOD
      METHOD BASIS_PROJECTION
      ADMM_PURIFICATION_METHOD MO_DIAG
      EXCH_CORRECTION_FUNC DEFAULT
    &END AUXILIARY_DENSITY_MATRIX_METHOD
    &SCF
      SCF_GUESS ATOMIC
      EPS_SCF 1.0E-6
      MAX_SCF 30
      &OT
        PRECONDITIONER FULL_SINGLE_INVERSE
        MINIMIZER DIIS
      &END OT
      &OUTER_SCF
        EPS_SCF 1.0E-6
        MAX_SCF 20
      &END OUTER_SCF
    &END SCF
    &XC
      &XC_FUNCTIONAL
        &PBE
          PARAMETRIZATION REVPBE
          SCALE_X 0.75
          SCALE_C 1.0
        &END PBE
      &END XC_FUNCTIONAL
      &HF
        FRACTION 0.25
        &SCREENING
          EPS_SCHWARZ 1.0E-6
          SCREEN_ON_INITIAL_P TRUE
        &END SCREENING
        &INTERACTION_POTENTIAL
          POTENTIAL_TYPE TRUNCATED
          CUTOFF_RADIUS 8.0
          T_C_G_DATA t_c_g.dat
        &END INTERACTION_POTENTIAL
        &MEMORY
          MAX_MEMORY 2000
          EPS_STORAGE_SCALING 0.1
        &END MEMORY
      &END HF
      &VDW_POTENTIAL
        DISPERSION_FUNCTIONAL PAIR_POTENTIAL
        &PAIR_POTENTIAL
          TYPE DFTD3
          REFERENCE_FUNCTIONAL revPBE0
          PARAMETER_FILE_NAME dftd3.dat
          R_CUTOFF 15.0
        &END PAIR_POTENTIAL
      &END VDW_POTENTIAL
    &END XC
  &END DFT
  &SUBSYS
    &CELL
      ABC {abc}
      PERIODIC XYZ
    &END CELL
    &TOPOLOGY
      COORD_FILE_NAME {coord}
      COORD_FILE_FORMAT XYZ
    &END TOPOLOGY
    &KIND O
      BASIS_SET DZVP-MOLOPT-SR-GTH
      BASIS_SET AUX_FIT cpFIT3
      POTENTIAL GTH-PBE-q6
    &END KIND
    &KIND H
      BASIS_SET DZVP-MOLOPT-SR-GTH
      BASIS_SET AUX_FIT cpFIT3
      POTENTIAL GTH-PBE-q1
    &END KIND
  &END SUBSYS
  &PRINT
    &FORCES ON
    &END FORCES
  &END PRINT
&END FORCE_EVAL
"""

# GGA = same numerics, ONLY the XC block differs (isolates the functional)
GGA = REVPBE0.replace(
"""      &XC_FUNCTIONAL
        &PBE
          PARAMETRIZATION REVPBE
          SCALE_X 0.75
          SCALE_C 1.0
        &END PBE
      &END XC_FUNCTIONAL
      &HF
        FRACTION 0.25
        &SCREENING
          EPS_SCHWARZ 1.0E-6
          SCREEN_ON_INITIAL_P TRUE
        &END SCREENING
        &INTERACTION_POTENTIAL
          POTENTIAL_TYPE TRUNCATED
          CUTOFF_RADIUS 8.0
          T_C_G_DATA t_c_g.dat
        &END INTERACTION_POTENTIAL
        &MEMORY
          MAX_MEMORY 2000
          EPS_STORAGE_SCALING 0.1
        &END MEMORY
      &END HF
""",
"""      &XC_FUNCTIONAL
        &PBE
          PARAMETRIZATION REVPBE
        &END PBE
      &END XC_FUNCTIONAL
""").replace("REFERENCE_FUNCTIONAL revPBE0", "REFERENCE_FUNCTIONAL revPBE") \
   .replace("BASIS_SET_FILE_NAME  BASIS_ADMM\n    ", "") \
   .replace("""    &AUXILIARY_DENSITY_MATRIX_METHOD
      METHOD BASIS_PROJECTION
      ADMM_PURIFICATION_METHOD MO_DIAG
      EXCH_CORRECTION_FUNC DEFAULT
    &END AUXILIARY_DENSITY_MATRIX_METHOD
""", "").replace("      BASIS_SET AUX_FIT cpFIT3\n", "")

for s in have:
    j   = os.path.basename(s).split("_")[1].split(".")[0]
    L   = read(s).cell.lengths()
    abc = f"{L[0]:.4f} {L[1]:.4f} {L[2]:.4f}"
    with open(os.path.join(FFTEST, f"revpbe0_{j}.inp"), "w") as f:
        f.write(REVPBE0.format(proj=f"revpbe0_{j}", abc=abc, coord=s))
    with open(os.path.join(FFTEST, f"gga_{j}.inp"), "w") as f:
        f.write(GGA.format(proj=f"gga_{j}", abc=abc, coord=s))
print(f"box (first snap): {abc}  ->  L/2 = {L[0]/2:.2f} A (HFX cutoff 8.0 OK)")
print(f"wrote {len(have)} revpbe0_*.inp + {len(have)} gga_*.inp in {FFTEST}")
