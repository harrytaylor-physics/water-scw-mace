# run_traj.py -- run MD ONCE, save the trajectory, restart-safe.
#   python run_traj.py <v4|foundation> <T_K> <rho_g_cc> <nmol> <tag> [--prod_ps P] [--equil_ps E] [--dt_fs D]
#
# Design:
#   * All state (positions+velocities+progress) is checkpointed every CHUNK steps to
#     traj/<tag>/<which>_n<nmol>.state.npz  -- a TIMEOUT loses at most one chunk.
#   * Rerunning the same command RESUMES from the checkpoint (equil is skipped if done).
#   * Production frames (positions, every DUMP_FS fs) append to
#     traj/<tag>/<which>_n<nmol>.frames.npz  -- the single source for ALL analyses
#     (RDF, S(k), finite-size, pressure); no analysis ever reruns MD.
#   * Pressure (total = virial + kinetic) is sampled every PSAMPLE steps and stored too.
import sys, os, time, numpy as np
from ase import units
from ase.md.langevin import Langevin
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution, Stationary, ZeroRotation
import bench_common as bc

# ---------------- arguments ----------------
a=sys.argv[1:]
def popopt(name, default, cast=float):
    if name in a:
        i=a.index(name); v=cast(a[i+1]); del a[i:i+2]; return v
    return default
prod_ps = popopt("--prod_ps", None)
equil_ps= popopt("--equil_ps", None)
dt_fs   = popopt("--dt_fs", 0.5)
which=a[0]; T=float(a[1]); rho=float(a[2]); nmol=int(a[3]); tag=a[4]

# defaults scale with size: big boxes self-average -> shorter runs, lighter equil
if equil_ps is None: equil_ps = 40 if nmol<=64 else (10 if nmol<=640 else 5)
if prod_ps  is None: prod_ps  = 60 if nmol<=64 else (12 if nmol<=640 else 6)

DUMP_FS  = 50            # store a frame every 50 fs (plenty for structure)
PSAMPLE  = 10            # sample pressure every 10 steps
CHUNK    = 2000          # max steps between checkpoints (= 1 ps at 0.5 fs)
CKPT_SEC = 300           # also checkpoint at least every 5 min of wall time

out=f"traj/{tag}"; os.makedirs(out, exist_ok=True)
state_f =f"{out}/{which}_n{nmol}.state.npz"
frames_f=f"{out}/{which}_n{nmol}.frames.npz"

dt=dt_fs*units.fs
equil_steps=int(equil_ps*1000/dt_fs); prod_steps=int(prod_ps*1000/dt_fs)
dump_every=max(1,int(DUMP_FS/dt_fs))

at,L=bc.build_box(nmol, rho, seed=0)
isO=np.array([s=="O" for s in at.get_chemical_symbols()])

# ---------------- resume or fresh start ----------------
done_equil=0; done_prod=0
frames=[]; Ps=[]
resumed=False
if os.path.exists(state_f):
    S=np.load(state_f)
    at.set_positions(S["pos"]); at.set_momenta(S["mom"])
    done_equil=int(S["done_equil"]); done_prod=int(S["done_prod"])
    resumed=True
    if os.path.exists(frames_f):
        F=np.load(frames_f)
        frames=[f for f in F["framesO"]]; Ps=list(F["Ps"])
    print(f"RESUME: equil {done_equil}/{equil_steps}  prod {done_prod}/{prod_steps}  "
          f"({len(frames)} frames, {len(Ps)} P samples already stored)")
else:
    MaxwellBoltzmannDistribution(at, temperature_K=T); Stationary(at); ZeroRotation(at)
    print(f"FRESH START: {nmol} mol ({len(at)} atoms) L={L:.2f}A  "
          f"equil {equil_ps} ps + prod {prod_ps} ps at {dt_fs} fs")

at.calc=bc.make_calc(which)

def save_state():
    np.savez(state_f, pos=at.get_positions(), mom=at.get_momenta(),
             done_equil=done_equil, done_prod=done_prod, T=T, rho=rho, L=L, nmol=nmol)
    np.savez(frames_f, framesO=np.array(frames) if frames else np.zeros((0,int(isO.sum()),3)),
             Ps=np.array(Ps), T=T, rho=rho, L=L, nmol=nmol, which=which,
             dump_fs=DUMP_FS, dt_fs=dt_fs)

# ---------------- equilibration (soft start only on genuinely fresh runs) ----------------
if not resumed:
    Langevin(at, 0.1*units.fs, temperature_K=T, friction=0.1).run(2000)   # drain overlaps
t0=time.time()
eq=Langevin(at, dt, temperature_K=T, friction=0.02)
last_ck=time.time()
SUB=100                                        # inner granularity for time-based checkpoints
while done_equil < equil_steps:
    n=min(SUB, equil_steps-done_equil)
    eq.run(n); done_equil+=n
    if done_equil % CHUNK == 0 or time.time()-last_ck > CKPT_SEC or done_equil>=equil_steps:
        save_state(); last_ck=time.time()
        rate=(done_equil*dt_fs/1000)/max(time.time()-t0,1e-9)
        print(f"  equil {done_equil}/{equil_steps} steps  ({rate:.2f} ps/s wall)", flush=True)

# ---------------- production ----------------
dyn=Langevin(at, dt, temperature_K=T, friction=0.01)
step_in_prod=done_prod
def grab():
    global step_in_prod
    step_in_prod+=1
    if step_in_prod % PSAMPLE == 0: Ps.append(bc.pressure_GPa(at))
    if step_in_prod % dump_every == 0: frames.append(at.get_positions()[isO].copy())
dyn.attach(grab, interval=1)
while done_prod < prod_steps:
    n=min(SUB, prod_steps-done_prod)
    dyn.run(n); done_prod+=n
    if done_prod % CHUNK == 0 or time.time()-last_ck > CKPT_SEC or done_prod>=prod_steps:
        save_state(); last_ck=time.time()
        print(f"  prod {done_prod}/{prod_steps} steps  ({len(frames)} frames, {len(Ps)} P samples)", flush=True)

P=np.array(Ps); se=bc.block_se(P,10) if len(P)>=10 else float("nan")
print(f"DONE {tag} {which} N={nmol}: {len(frames)} frames, "
      f"<P>={P.mean():+.3f} +/- {se:.3f} GPa (total)  -> {frames_f}")
