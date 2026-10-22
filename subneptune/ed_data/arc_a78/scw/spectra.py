import sys, os, glob, numpy as np, torch
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from ase.io import iread
from mace import data
from mace.tools import torch_geometric, utils, torch_tools

trajdir=sys.argv[1]; dt_fs=float(sys.argv[2]); T=float(sys.argv[3])
mpath=sys.argv[4] if len(sys.argv)>4 else "water_dipole_v2.model"
dt=dt_fs*1e-15
device="cuda"
model=torch.load(f=mpath, map_location=device).to(device); model.eval()
for p in model.parameters(): p.requires_grad=False
z_table=utils.AtomicNumberTable([int(z) for z in model.atomic_numbers])
try: heads=model.heads
except AttributeError: heads=None
rmax=float(model.r_max)
def run(batch):
    cfgs=[data.config_from_atoms(a, head_name="Default") for a in batch]
    ds=[data.AtomicData.from_config(c,z_table=z_table,cutoff=rmax,heads=heads) for c in cfgs]
    dl=torch_geometric.dataloader.DataLoader(dataset=ds,batch_size=4,shuffle=False); o=[]
    for b in dl:
        b=b.to(device)
        with torch.no_grad(): out=model(b.to_dict())
        o.append(torch_tools.to_numpy(out["dipole"])); del b,out; torch.cuda.empty_cache()
    return np.concatenate(o,0)

dipf=os.path.join(trajdir,"dipoles.npy")
if os.path.exists(dipf):
    M=np.load(dipf); print(f"reusing {dipf} {M.shape}")
else:
    dips=[]; buf=[]; n=0
    for at in iread(os.path.join(trajdir,"prod.xyz")):
        buf.append(at)
        if len(buf)>=100: dips.append(run(buf)); n+=len(buf); buf=[]
        if n and n%2000==0 and len(buf)==0: print(f"  {n} frames done", flush=True)
    if buf: dips.append(run(buf)); n+=len(buf)
    M=np.concatenate(dips,0); np.save(dipf,M)
    print(f"dipoles: {M.shape} mean|M|={np.linalg.norm(M,axis=1).mean():.3f} e.A")

M=M-M.mean(0)
def ac(x):
    n=len(x); f=np.fft.rfft(x,n=2*n); a=np.fft.irfft(f*np.conj(f),n=2*n)[:n].real
    a/=(n-np.arange(n)); return a
C=sum(ac(M[:,k]) for k in range(3))
print(f"C(dt)/C(0) = {C[1]/C[0]:.4f}   (>=0.997 smooth; ~0.98 white model-noise floor)")
maxlag=min(len(C)//2,8000); Cw=C[:maxlag]*np.hanning(2*maxlag)[maxlag:]
nfft=2**17; I=np.fft.rfft(Cw,n=nfft).real
freq=np.fft.rfftfreq(nfft,dt); wn=freq/2.99792458e10
hbar=1.054571817e-34; kB=1.380649e-23; w=2*np.pi*freq
alpha=w*(1-np.exp(-hbar*w/(kB*T)))*I
nyq=wn.max()
flm=(wn>min(3900,0.8*nyq))&(wn<min(5500,0.95*nyq))
fl=np.median(np.abs(I[flm]))
print(f"white floor = {fl:.3e} ({100*fl/np.abs(I[(wn>5)]).max():.1f}% of max density)")
for lo,hi,nm in [(50,900,"libration"),(1200,2000,"bend"),(2800,3800,"stretch")]:
    mm=(wn>lo)&(wn<hi)
    if mm.any(): print(f"{nm}: peak {wn[mm][np.argmax(np.abs(I[mm]))]:.0f} cm-1, {np.abs(I[mm]).max()/fl:.2f}x floor")
np.savez(os.path.join(trajdir,"spectrum.npz"), wn=wn, I=I, alpha=alpha, floor=fl, T=T, dt=dt)
fig,ax=plt.subplots(2,1,figsize=(9,7))
m=(wn>5)&(wn<min(4500,0.95*nyq))
ax[0].semilogy(wn[m],np.abs(I[m])); ax[0].axhline(fl,ls=":",c="r"); ax[0].set_ylabel("|I| (log)")
mz=(wn>900)&(wn<min(4200,0.9*nyq)); Iz=np.clip(np.abs(I[mz])-fl,0,None)
ax[1].plot(wn[mz],Iz/(Iz.max() if Iz.max()>0 else 1)); ax[1].set_ylabel("I - floor (vib window)"); ax[1].set_xlabel("cm-1")
for a_ in ax:
    for x in (600,1600,3400): a_.axvline(x,ls="--",c="grey",alpha=.4)
plt.suptitle(f"{os.path.basename(trajdir)}  T={T:.0f} K"); plt.tight_layout()
plt.savefig(os.path.join(trajdir,"spectrum.png"),dpi=120)
print("wrote", os.path.join(trajdir,"spectrum.png"))
