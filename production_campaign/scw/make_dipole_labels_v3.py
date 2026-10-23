import glob, os, re
from ase.io import read, write
tmpl=open("mlwf_test/mlwf.inp").read()
entries=[]
for pf in sorted(glob.glob("snippets/snip_T*_rho*.xyz")):
    stem=os.path.basename(pf)[:-4]
    frames=read(pf,index=":")
    L=float(frames[0].cell.lengths()[0])
    inp=re.sub(r"ABC\s+[\d. ]+", f"ABC {L:.4f} {L:.4f} {L:.4f}", tmpl)
    for i,at in enumerate(frames):
        d=f"diplabels_v3/{stem}/f{i:04d}"; os.makedirs(d,exist_ok=True)
        write(f"{d}/coord.xyz",at,format="xyz")
        open(f"{d}/dip.inp","w").write(inp)
        entries.append(f"{d}/dip.inp")
open("dip3.list","w").write("\n".join(entries)+"\n")
print(f"{len(entries)} MLWF single points -> dip3.list")
print(f"submit: sbatch --array=0-{len(entries)-1}%24 htc_dipole_label3.sh")
