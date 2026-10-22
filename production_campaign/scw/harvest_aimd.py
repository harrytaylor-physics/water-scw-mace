import glob, os, re
from ase.io import read, write
os.makedirs("pools_aimd", exist_ok=True)
groups={}
for d in sorted(glob.glob("aimd/cell_T*_rho*")):
    m=re.search(r"cell_T(\d+)_rho(\d+)", d)
    if not m: continue
    groups.setdefault((m.group(1),m.group(2)),[]).append(d)
tot=0
for (T,rho),dirs in sorted(groups.items()):
    out=f"pools_aimd/pool_T{T}_rho{rho}.xyz"
    if os.path.exists(out): os.remove(out)
    n=0
    for d in dirs:
        traj=os.path.join(d,"aimd-pos-1.xyz")
        if not os.path.exists(traj): continue
        L=float(re.search(r"ABC\s+([\d.]+)",open(os.path.join(d,"aimd.inp")).read()).group(1))
        frames=read(traj,index=":")
        keep=frames[50::5]                      # drop ~0.5 ps equil, then every 5th
        for at in keep:
            at.set_cell([L,L,L]); at.set_pbc(True)
            write(out,at,append=True); n+=1
    tot+=n
    print(f"T{T} rho{rho}: {len(dirs)} traj -> {n} frames  (L from inp)")
print(f"total harvested: {tot}")
