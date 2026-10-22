import re, os
from ase.io import read, write
src = "aimd/cell_T6000_rho4894"
frames = read(f"{src}/aimd-pos-1.xyz", index=":")
at = frames[-5]                      # ~100 steps before the crash; settled geometry
L = float(re.search(r'ABC\s+([\d.]+)', open(f"{src}/aimd.inp").read()).group(1))
at.set_cell([L,L,L]); at.set_pbc(True)
out = "aimd/cell_T6000_rho4894_csvr"; os.makedirs(out, exist_ok=True)
write(f"{out}/coord.xyz", at)
inp = open(f"{src}/aimd.inp").read()
inp = re.sub(r'&THERMOSTAT.*?&END THERMOSTAT',
             '&THERMOSTAT\n      TYPE CSVR\n      &CSVR\n        TIMECON 50\n      &END CSVR\n    &END THERMOSTAT',
             inp, flags=re.DOTALL)
inp = re.sub(r'TIMESTEP\s+[\d.]+', 'TIMESTEP 0.25', inp)
open(f"{out}/aimd.inp","w").write(inp)
print(f"restart from step ~{(len(frames)-5)*20}; CSVR + 0.25 fs -> {out}/")
