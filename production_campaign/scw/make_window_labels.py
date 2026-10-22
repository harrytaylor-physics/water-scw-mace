import os, re
from ase.io import read, write
tmpl=open("mlwf_test/mlwf.inp").read()
frames=read("window_traj_T1000_rho130.xyz",index=":")
assert len(frames)==501 and len(frames[0])==192, (len(frames), len(frames[0]))
L=float(frames[0].cell.lengths()[0])
inp=re.sub(r"ABC\s+[\d. ]+", f"ABC {L:.4f} {L:.4f} {L:.4f}", tmpl)
entries=[]
for i,at in enumerate(frames):
    d=f"diplabels_window/T1000_rho130/f{i:04d}"; os.makedirs(d,exist_ok=True)
    write(f"{d}/coord.xyz",at,format="xyz")
    open(f"{d}/dip.inp","w").write(inp)
    entries.append(f"{d}/dip.inp")
open("dipwin.list","w").write("\n".join(entries)+"\n")
print(f"{len(entries)} MLWF single points -> dipwin.list  L={L:.4f} A")
