import os, re, shutil
SRC = "labels/pool_T2250_rho1.26e+03/f0000"   # insulating, dense, was wrapping
L = float(re.search(r'ABC\s+([\d.]+)', open(f"{SRC}/input.inp").read()).group(1))
os.makedirs("mlwf_test", exist_ok=True)
shutil.copy(f"{SRC}/coord.xyz", "mlwf_test/coord.xyz")
T=f"""&GLOBAL
  PROJECT mlwf
  RUN_TYPE ENERGY
  PRINT_LEVEL LOW
&END GLOBAL
&FORCE_EVAL
  METHOD QUICKSTEP
  &DFT
    BASIS_SET_FILE_NAME BASIS_MOLOPT
    POTENTIAL_FILE_NAME GTH_POTENTIALS
    &MGRID
      CUTOFF 600
      REL_CUTOFF 60
      NGRIDS 5
    &END MGRID
    &QS
      METHOD GPW
      EPS_DEFAULT 1.0E-12
    &END QS
    &SCF
      SCF_GUESS ATOMIC
      EPS_SCF 1.0E-7
      MAX_SCF 60
      &OT
        PRECONDITIONER FULL_SINGLE_INVERSE
        MINIMIZER DIIS
      &END OT
      &OUTER_SCF
        EPS_SCF 1.0E-7
        MAX_SCF 20
      &END OUTER_SCF
    &END SCF
    &LOCALIZE T
      METHOD CRAZY
      MAX_ITER 2000
      &PRINT
        &WANNIER_CENTERS
          IONS+CENTERS
          FILENAME =wannier
        &END WANNIER_CENTERS
      &END PRINT
    &END LOCALIZE
    &XC
      &XC_FUNCTIONAL
        &PBE
          PARAMETRIZATION REVPBE
        &END PBE
      &END XC_FUNCTIONAL
      &VDW_POTENTIAL
        DISPERSION_FUNCTIONAL PAIR_POTENTIAL
        &PAIR_POTENTIAL
          TYPE DFTD3
          REFERENCE_FUNCTIONAL revPBE
          PARAMETER_FILE_NAME dftd3.dat
          R_CUTOFF 15.0
        &END PAIR_POTENTIAL
      &END VDW_POTENTIAL
    &END XC
  &END DFT
  &SUBSYS
    &CELL
      ABC {L:.4f} {L:.4f} {L:.4f}
      PERIODIC XYZ
    &END CELL
    &TOPOLOGY
      COORD_FILE_NAME coord.xyz
      COORD_FILE_FORMAT XYZ
    &END TOPOLOGY
    &KIND O
      BASIS_SET DZVP-MOLOPT-SR-GTH
      POTENTIAL GTH-PBE-q6
    &END KIND
    &KIND H
      BASIS_SET DZVP-MOLOPT-SR-GTH
      POTENTIAL GTH-PBE-q1
    &END KIND
  &END SUBSYS
&END FORCE_EVAL
"""
open("mlwf_test/mlwf.inp","w").write(T)
print(f"src={SRC}  L={L:.4f}A  -> mlwf_test/")
