import os
from ase.io import iread, write
cell = os.path.join(os.environ['DATA'], 'scw_prod', 'T3000_rho108.796_05fs')
os.makedirs("snaps", exist_ok=True)
K, j = 2000, 0          # every 2000 frames (~1 ps apart) -> ~30 decorrelated snapshots
for i, at in enumerate(iread(os.path.join(cell, "prod.xyz"))):
    if i % K == 0:
        write(f"snaps/snap_{j:03d}.xyz", at); j += 1
print(f"wrote {j} snapshots to snaps/")
