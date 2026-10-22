#!/usr/bin/env python3
"""Phase 2 prep: turn foundation-MACE geometry pools into revPBE0-D3 single-point jobs.

For each pools/pool_T*_rho*.xyz, pick a sparse, evenly-spaced set of frames and write,
per frame, a CP2K revPBE0-D3+ADMM ENERGY_FORCE input (OT SCF; the insulating bulk regime)
plus a coord.xyz. Each input also prints the Berry-phase cell dipole, so one single point
yields energy + forces + dipole. Paths are collected into frames.list for the job array.

Level of theory is copied verbatim from water64_moderate.inp; only ABC, the HF truncation
radius (must stay < L/2), the MGRID cutoff, the project name and the coord file vary.

Usage (login node, in ~/scw):
    python select_and_stamp.py --per_cell 12
Then it prints the exact `sbatch --array=...` line to run.
"""
import argparse, glob, os, re, numpy as np
from ase.io import read, write

TEMPLATE = """\
&GLOBAL
  PROJECT __PROJECT__
  RUN_TYPE ENERGY_FORCE
  PRINT_LEVEL LOW
&END GLOBAL

&FORCE_EVAL
  METHOD QUICKSTEP
  STRESS_TENSOR ANALYTICAL
  &DFT
    BASIS_SET_FILE_NAME  BASIS_MOLOPT
    BASIS_SET_FILE_NAME  BASIS_ADMM
    POTENTIAL_FILE_NAME  GTH_POTENTIALS
    &MGRID
      CUTOFF      __MGRID__
      REL_CUTOFF  60
      NGRIDS      5
    &END MGRID
    &QS
      METHOD       GPW
      EPS_DEFAULT  1.0E-12
      EPS_PGF_ORB  1.0E-8
    &END QS
    &AUXILIARY_DENSITY_MATRIX_METHOD
      METHOD                    BASIS_PROJECTION
      ADMM_PURIFICATION_METHOD  MO_DIAG
      EXCH_CORRECTION_FUNC      DEFAULT
    &END AUXILIARY_DENSITY_MATRIX_METHOD
    &SCF
      SCF_GUESS  ATOMIC
      EPS_SCF    1.0E-6
      MAX_SCF    30
      &OT
        PRECONDITIONER  FULL_SINGLE_INVERSE
        MINIMIZER       DIIS
      &END OT
      &OUTER_SCF
        EPS_SCF  1.0E-6
        MAX_SCF  20
      &END OUTER_SCF
    &END SCF
    &XC
      &XC_FUNCTIONAL
        &PBE
          PARAMETRIZATION  REVPBE
          SCALE_X          0.75
          SCALE_C          1.0
        &END PBE
      &END XC_FUNCTIONAL
      &HF
        FRACTION  0.25
        &SCREENING
          EPS_SCHWARZ            1.0E-6
          SCREEN_ON_INITIAL_P    TRUE
        &END SCREENING
        &INTERACTION_POTENTIAL
          POTENTIAL_TYPE  TRUNCATED
          CUTOFF_RADIUS   __CUTOFF__
          T_C_G_DATA      t_c_g.dat
        &END INTERACTION_POTENTIAL
        &MEMORY
          MAX_MEMORY             2000
          EPS_STORAGE_SCALING    0.1
        &END MEMORY
      &END HF
      &VDW_POTENTIAL
        DISPERSION_FUNCTIONAL  PAIR_POTENTIAL
        &PAIR_POTENTIAL
          TYPE                  DFTD3
          REFERENCE_FUNCTIONAL  revPBE0
          PARAMETER_FILE_NAME   dftd3.dat
          R_CUTOFF              15.0
        &END PAIR_POTENTIAL
      &END VDW_POTENTIAL
    &END XC
    &PRINT
      &MOMENTS
        PERIODIC .TRUE.
      &END MOMENTS
    &END PRINT
  &END DFT
  &SUBSYS
    &CELL
      ABC       __ABC__
      PERIODIC  XYZ
    &END CELL
    &TOPOLOGY
      COORD_FILE_NAME    coord.xyz
      COORD_FILE_FORMAT  XYZ
    &END TOPOLOGY
    &KIND O
      BASIS_SET          DZVP-MOLOPT-SR-GTH
      BASIS_SET AUX_FIT  cpFIT3
      POTENTIAL          GTH-PBE-q6
    &END KIND
    &KIND H
      BASIS_SET          DZVP-MOLOPT-SR-GTH
      BASIS_SET AUX_FIT  cpFIT3
      POTENTIAL          GTH-PBE-q1
    &END KIND
  &END SUBSYS
  &PRINT
    &FORCES ON
    &END FORCES
  &END PRINT
&END FORCE_EVAL
"""

_OT_BLOCK = """    &SCF
      SCF_GUESS  ATOMIC
      EPS_SCF    1.0E-6
      MAX_SCF    30
      &OT
        PRECONDITIONER  FULL_SINGLE_INVERSE
        MINIMIZER       DIIS
      &END OT
      &OUTER_SCF
        EPS_SCF  1.0E-6
        MAX_SCF  20
      &END OUTER_SCF
    &END SCF"""

_SMEAR_BLOCK = """    &SCF
      SCF_GUESS  ATOMIC
      EPS_SCF    1.0E-6
      MAX_SCF    150
      ADDED_MOS  __ADDED__
      &SMEAR
        METHOD                  FERMI_DIRAC
        ELECTRONIC_TEMPERATURE  __ETEMP__
      &END SMEAR
      &DIAGONALIZATION
        ALGORITHM  STANDARD
      &END DIAGONALIZATION
      &MIXING
        METHOD    BROYDEN_MIXING
        ALPHA     0.4
        NBROYDEN  8
      &END MIXING
    &END SCF"""

# hot/dense cells: same revPBE0-D3+ADMM, but Fermi smearing + diagonalisation
TEMPLATE_HOT = TEMPLATE.replace(_OT_BLOCK, _SMEAR_BLOCK).replace("ADMM_PURIFICATION_METHOD  MO_DIAG", "ADMM_PURIFICATION_METHOD  NONE")
import re as _re
TEMPLATE_HOT = _re.sub(r"[ \t]*&PRINT\s*\n\s*&MOMENTS.*?&END MOMENTS\s*\n\s*&END PRINT\s*\n", "", TEMPLATE_HOT, flags=_re.S)
assert "MOMENTS" not in TEMPLATE_HOT, "hot template still prints moments"

assert TEMPLATE_HOT != TEMPLATE, "hot-template SCF swap failed (text mismatch)"
HOT_T = 2500.0   # cells at/above this temperature use the smeared template

def stamp(atoms, project, outdir, T):
    os.makedirs(outdir, exist_ok=True)
    L = float(atoms.cell.lengths()[0])          # cubic boxes
    cutoff = min(8.0, round(0.49 * L, 1))        # truncated-Coulomb HFX needs < L/2
    rho_mass = sum(atoms.get_masses()) * 1.66054 / atoms.get_volume()  # g/cm^3
    mgrid = 800 if rho_mass > 0.5 else 600
    write(os.path.join(outdir, "coord.xyz"), atoms, format="xyz")
    if T >= HOT_T:
        added = max(40, 3 * (len(atoms) // 3))   # empty states for Fermi smearing
        inp = (TEMPLATE_HOT.replace("__PROJECT__", project)
                           .replace("__ABC__", f"{L:.4f} {L:.4f} {L:.4f}")
                           .replace("__CUTOFF__", f"{cutoff:.1f}")
                           .replace("__MGRID__", str(mgrid))
                           .replace("__ADDED__", str(added))
                           .replace("__ETEMP__", f"{T:.1f}"))
    else:
        inp = (TEMPLATE.replace("__PROJECT__", project)
                       .replace("__ABC__", f"{L:.4f} {L:.4f} {L:.4f}")
                       .replace("__CUTOFF__", f"{cutoff:.1f}")
                       .replace("__MGRID__", str(mgrid)))
    path = os.path.join(outdir, "input.inp")
    with open(path, "w") as f:
        f.write(inp)
    return path

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pools", default="pools/*.xyz")
    ap.add_argument("--per_cell", type=int, default=12)
    ap.add_argument("--max_L", type=float, default=30.0,
                    help="skip cells whose cubic box exceeds this (A); dilute cells "
                         "are handled as gas-phase clusters, not periodic DFT")
    ap.add_argument("--out", default="labels")
    ap.add_argument("--framelist", default="frames.list")
    ap.add_argument("--concurrent", type=int, default=16)
    a = ap.parse_args()

    pool_files = sorted(glob.glob(a.pools))
    if not pool_files:
        raise SystemExit(f"no pools matched {a.pools}")
    entries = []
    skipped = []
    for pf in pool_files:
        stem = re.sub(r"\.xyz$", "", os.path.basename(pf))   # pool_T650_rho1.68
        frames = read(pf, index=":")
        if len(frames) == 0:
            print(f"  WARN empty pool {pf}"); continue
        L = float(frames[0].cell.lengths()[0])
        if L > a.max_L:
            skipped.append((stem, L)); 
            print(f"  SKIP {stem}: box L={L:.1f} A > {a.max_L:.0f} (dilute -> clusters)")
            continue
        Tcell = float(re.search(r"_T(\d+)", stem).group(1))
        idx = np.unique(np.linspace(0, len(frames) - 1, a.per_cell).round().astype(int))
        for i in idx:
            project = f"{stem}_f{i:04d}"
            outdir = os.path.join(a.out, stem, f"f{i:04d}")
            entries.append(stamp(frames[i], project, outdir, Tcell))
        tag = "HOT/smeared" if Tcell >= HOT_T else "OT"
        print(f"  {stem}: L={L:.1f} A, T={Tcell:.0f} [{tag}], {len(frames)} frames -> {len(idx)} labels")

    with open(a.framelist, "w") as f:
        f.write("\n".join(entries) + "\n")
    n = len(entries)
    print(f"\nwrote {n} single-point inputs to {a.framelist}")
    if skipped:
        print(f"skipped {len(skipped)} dilute cells (need cluster treatment): "
              + ", ".join(s for s, _ in skipped))
    print(f"submit with:\n  sbatch --array=0-{n-1}%{a.concurrent} htc_label_array.sh")

if __name__ == "__main__":
    main()
