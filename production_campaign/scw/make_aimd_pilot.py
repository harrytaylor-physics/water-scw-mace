#!/usr/bin/env python3
import os, numpy as np
from ase import Atoms
from ase.io import write
NA=6.02214076e23; M=18.01528e-3
rOH=0.9572; ang=np.deg2rad(104.52)
TEMPL=np.array([[0,0,0],[rOH,0,0],[rOH*np.cos(ang),rOH*np.sin(ang),0]])
def rand_rot():
    q=np.random.normal(size=4); q/=np.linalg.norm(q); w,x,y,z=q
    return np.array([[1-2*(y*y+z*z),2*(x*y-z*w),2*(x*z+y*w)],
                     [2*(x*y+z*w),1-2*(x*x+z*z),2*(y*z-x*w)],
                     [2*(x*z-y*w),2*(y*z+x*w),1-2*(x*x+y*y)]])
def build(nmol, rho):
    L=(nmol*M/NA/rho*1e30)**(1/3)
    ng=int(round(nmol**(1/3)))
    while ng**3 < nmol: ng+=1
    sites=[(i,j,k) for i in range(ng) for j in range(ng) for k in range(ng)][:nmol]
    sym=[]; pos=[]
    for (i,j,k) in sites:
        c=(np.array([i,j,k])+0.5)*L/ng
        m=(rand_rot()@TEMPL.T).T+c; sym+=['O','H','H']; pos+=list(m)
    return Atoms(''.join(sym),positions=np.array(pos),cell=[L,L,L],pbc=True), L
INP="""&GLOBAL
  PROJECT aimd
  RUN_TYPE MD
  PRINT_LEVEL LOW
&END GLOBAL
&MOTION
  &MD
    ENSEMBLE NVT
    STEPS {steps}
    TIMESTEP 0.5
    TEMPERATURE {T}
    &THERMOSTAT
      TYPE NOSE
      &NOSE
        TIMECON 50
      &END NOSE
    &END THERMOSTAT
  &END MD
  &PRINT
    &TRAJECTORY
      &EACH
        MD 20
      &END EACH
    &END TRAJECTORY
    &RESTART
      &EACH
        MD 200
      &END EACH
    &END RESTART
  &END PRINT
&END MOTION
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
      EPS_SCF 1.0E-6
      MAX_SCF 50
      ADDED_MOS {added}
      &SMEAR ON
        METHOD FERMI_DIRAC
        ELECTRONIC_TEMPERATURE {T}
      &END SMEAR
      &DIAGONALIZATION
        ALGORITHM STANDARD
      &END DIAGONALIZATION
      &MIXING
        METHOD BROYDEN_MIXING
        ALPHA 0.4
        NBROYDEN 8
      &END MIXING
    &END SCF
    &XC
      &XC_FUNCTIONAL
        &PBE
          PARAMETRIZATION REVPBE
          SCALE_X 1.0
          SCALE_C 1.0
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
T=4000.0; rho=2940.0; nmol=64; steps=3000
at,L=build(nmol,rho)
os.makedirs("aimd_pilot",exist_ok=True)
write("aimd_pilot/coord.xyz",at)
open("aimd_pilot/aimd.inp","w").write(INP.format(steps=steps,T=int(T),added=2*nmol,L=L))
O=at.positions[at.numbers==8]; d=O[:,None,:]-O[None,:,:]; d-=L*np.round(d/L)
r=np.sqrt((d**2).sum(-1)); np.fill_diagonal(r,9)
print(f"T={T:.0f}K rho={rho:.0f} nmol={nmol} L={L:.3f}A  min O-O={r.min():.2f}A  ADDED_MOS={2*nmol}")
