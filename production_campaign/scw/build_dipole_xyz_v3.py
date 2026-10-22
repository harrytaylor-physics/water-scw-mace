import glob, re, os, numpy as np
from ase.io import read, write
def dip(d):
    fn=f"{d}/wannier"
    if not os.path.exists(fn): return None
    L=float(re.search(r'ABC\s+([\d.]+)', open(f"{d}/dip.inp").read()).group(1))
    raw=[l for l in open(fn).read().split('\n') if l.strip()]
    try: n=int(raw[0].split()[0]); body=raw[2:2+n]
    except Exception: return None
    sym=[];pos=[]
    for ln in body:
        p=ln.split()
        if len(p)>=4 and re.match(r'[A-Za-z]',p[0]): sym.append(p[0]); pos.append([float(x) for x in p[1:4]])
    sym=np.array(sym); pos=np.array(pos)
    Oi=np.where(sym=='O')[0]; Hi=np.where(sym=='H')[0]; Xi=np.where(~np.isin(sym,['O','H']))[0]
    if len(Xi)!=4*len(Oi) or len(Hi)!=2*len(Oi): return None
    Op=pos[Oi]
    mic=lambda a,b:(lambda x:x-L*np.round(x/L))(a-b)
    aX=np.array([int(np.argmin((mic(Op,pos[x])**2).sum(1))) for x in Xi])
    Hd=[]
    for hi,h in enumerate(Hi):
        r=np.sqrt((mic(Op,pos[h])**2).sum(1)); Hd.append((float(r.min()),hi,np.argsort(r)))
    Hd.sort(); cap=np.zeros(len(Oi),int); aH=np.full(len(Hi),-1)
    for _,hi,order in Hd:
        for j in order:
            if cap[j]<2: aH[hi]=j; cap[j]+=1; break
    if (cap!=2).any(): return None
    mu=np.zeros(3)
    for o in range(len(Oi)):
        rO=Op[o]
        for hi in np.where(aH==o)[0]: mu+=mic(pos[Hi[hi]],rO)
        for xi in np.where(aX==o)[0]: mu+=-2.0*mic(pos[Xi[xi]],rO)
    return mu, L
out=[]
for d in [l.strip() for l in open("dip3_dirs.list") if l.strip()]:
    if not os.path.exists(f"{d}/dip.out") or "PROGRAM ENDED" not in open(f"{d}/dip.out").read(): continue
    r=dip(d)
    if r is None: continue
    mu,L=r; at=read(f"{d}/coord.xyz"); at.set_cell([L,L,L]); at.set_pbc(True)
    at.info['REF_dipole']=mu.tolist(); at.info['config_type']=d.split('/')[1]; out.append(at)
write("dipole_train_snip.xyz", out)
print(f"dipole_train_snip.xyz: {len(out)} configs with MLWF dipoles")
