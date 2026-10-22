import numpy as np, sys
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from ase.io import iread
traj=sys.argv[1]; dt=float(sys.argv[2])*1e-15
pos=[]; sym=None
for at in iread(traj):
    pos.append(at.get_positions())
    if sym is None: sym=np.array(at.get_chemical_symbols())
P=np.array(pos); print("frames:",P.shape)
V=(P[2:]-P[:-2])/(2*dt*1e15)   # A/fs, central difference (positions are unwrapped)
def pspec(Vs):
    Vs=Vs-Vs.mean(0)
    n=len(Vs); w=np.hanning(n)[:,None,None]
    F=np.fft.rfft(Vs*w,axis=0)
    return (np.abs(F)**2).sum(axis=(1,2))
freq=np.fft.rfftfreq(len(V),dt); wn=freq/2.99792458e10
SH=pspec(V[:,sym=="H",:]); SO=pspec(V[:,sym=="O",:])
m=(wn>5)&(wn<4500)
plt.figure(figsize=(9,5))
plt.plot(wn[m],SH[m]/SH[m].max(),label="H")
plt.plot(wn[m],SO[m]/SO[m].max(),label="O")
for x in (600,1600,3400): plt.axvline(x,ls="--",c="grey",alpha=.4)
plt.xlabel("cm-1"); plt.ylabel("VDOS (norm.)"); plt.legend()
plt.title("velocity power spectrum (no dipole model)  T=2250 K, 1.26 g/cc")
plt.tight_layout(); plt.savefig(traj.rsplit("/",1)[0]+"/vdos.png",dpi=120)
for lo,hi,nm in [(1200,2000,"bend"),(2800,3900,"stretch")]:
    mm=(wn>lo)&(wn<hi); print(f"{nm}: H-VDOS peak {wn[mm][np.argmax(SH[mm])]:.0f} cm-1, rel {SH[mm].max()/SH[m].max():.3f}")
print("wrote vdos.png")
