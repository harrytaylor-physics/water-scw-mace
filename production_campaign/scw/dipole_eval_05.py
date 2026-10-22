import numpy as np, torch
from ase.io import iread
from mace import data
from mace.tools import torch_geometric, utils, torch_tools
device="cuda"
model=torch.load(f="water_dipole_v2.model", map_location=device).to(device); model.eval()
for p in model.parameters(): p.requires_grad=False
z_table=utils.AtomicNumberTable([int(z) for z in model.atomic_numbers])
try: heads=model.heads
except AttributeError: heads=None
rmax=float(model.r_max)
def run(batch):
    cfgs=[data.config_from_atoms(a, head_name="Default") for a in batch]
    ds=[data.AtomicData.from_config(c, z_table=z_table, cutoff=rmax, heads=heads) for c in cfgs]
    dl=torch_geometric.dataloader.DataLoader(dataset=ds, batch_size=4, shuffle=False); o=[]
    for b in dl:
        b=b.to(device)
        with torch.no_grad(): out=model(b.to_dict())
        o.append(torch_tools.to_numpy(out["dipole"])); del b,out; torch.cuda.empty_cache()
    return np.concatenate(o,0)
dips=[]; buf=[]; n=0
for at in iread("md/bulk_05fs/prod.xyz"):
    buf.append(at)
    if len(buf)>=100: dips.append(run(buf)); n+=len(buf); buf=[]
if buf: dips.append(run(buf)); n+=len(buf)
M=np.concatenate(dips,0); np.save("md/bulk_05fs/dipoles.npy", M)
print("dipoles:",M.shape," frames=%d mean|M|=%.3f e.A"%(n,np.linalg.norm(M,axis=1).mean()))
