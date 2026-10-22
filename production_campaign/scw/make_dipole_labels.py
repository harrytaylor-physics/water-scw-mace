import os, re, glob, shutil
base = open("mlwf_test/mlwf.inp").read()          # the proven &LOCALIZE input
def make_inp(L):
    return re.sub(r'ABC\s+[\d.]+\s+[\d.]+\s+[\d.]+', f'ABC {L:.4f} {L:.4f} {L:.4f}', base)
os.makedirs("diplabels", exist_ok=True); listing=[]
for d in sorted(glob.glob("labels/pool_T*/f*")):
    cell = d.split('/')[1]
    m = re.match(r'pool_T(\d+)_', cell)
    if not m or int(m.group(1)) > 2250: continue   # bulk insulating only, first model
    if not (os.path.exists(f"{d}/input.inp") and os.path.exists(f"{d}/coord.xyz")): continue
    Lm = re.search(r'ABC\s+([\d.]+)', open(f"{d}/input.inp").read())
    if not Lm: continue
    out = f"diplabels/{cell}/{os.path.basename(d)}"; os.makedirs(out, exist_ok=True)
    shutil.copy(f"{d}/coord.xyz", f"{out}/coord.xyz")
    open(f"{out}/dip.inp","w").write(make_inp(float(Lm.group(1))))
    listing.append(out)
open("diplabels/frames.list","w").write("\n".join(listing)+"\n")
print(f"{len(listing)} bulk insulating configs -> diplabels/")
