"""
mace_opacity.py  --  per-cell opacity from MACE-MD dipoles, via the validated backend.

Usage:
    python mace_opacity.py <trajdir> [T_K|auto] [N] [dt_fs] [model] [smooth_cm]

- Reads <trajdir>/dipoles.npy if present  -> pure CPU, instant.
- Else computes it with spectra.py's EXACT dipole eval (GPU), then caches it.
- T and rho are parsed from the dirname (e.g. T3000_rho108.796) unless T given.
- Writes <trajdir>/opacity_v3.npz: wn, sigma, sigma_raw, alpha, kappa, planck, rosseland, T, rho, N, dt.
"""
import sys, os, re, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from opacity_pipeline import spectrum_from_dipoles

trajdir = sys.argv[1].rstrip("/")
base = os.path.basename(trajdir)
mt = re.search(r"T(\d+(?:\.\d+)?)", base)
mr = re.search(r"rho([\d.]+)", base)
T   = float(sys.argv[2]) if len(sys.argv) > 2 and sys.argv[2] not in ("auto", "-") else float(mt.group(1))
N   = int(sys.argv[3])   if len(sys.argv) > 3 else 64
dt_fs = float(sys.argv[4]) if len(sys.argv) > 4 else 1.0
mpath = sys.argv[5] if len(sys.argv) > 5 else "water_dipole_v3.model"
rho_kg_m3 = float(mr.group(1))
density_g_cm3 = rho_kg_m3 / 1000.0
# regime-dependent smoothing, matching the AIMD convention (5-15 / 20 / 40 cm-1)
smooth = float(sys.argv[6]) if len(sys.argv) > 6 else (15.0 if T < 2000 else 20.0 if T < 3500 else 40.0)

dipf = os.path.join(trajdir, "dipoles.npy")
if os.path.exists(dipf):
    M = np.load(dipf)
    print(f"reusing {dipf}  {M.shape}  (CPU-only, no GPU needed)")
else:
    # ---- verbatim dipole eval from spectra.py (only reached if cache is missing) ----
    import torch
    from ase.io import iread
    from mace import data
    from mace.tools import torch_geometric, utils, torch_tools
    device = "cuda"
    model = torch.load(f=mpath, map_location=device).to(device); model.eval()
    for p in model.parameters(): p.requires_grad = False
    z_table = utils.AtomicNumberTable([int(z) for z in model.atomic_numbers])
    try: heads = model.heads
    except AttributeError: heads = None
    rmax = float(model.r_max)
    def run(batch):
        cfgs = [data.config_from_atoms(a, head_name="Default") for a in batch]
        ds = [data.AtomicData.from_config(c, z_table=z_table, cutoff=rmax, heads=heads) for c in cfgs]
        dl = torch_geometric.dataloader.DataLoader(dataset=ds, batch_size=4, shuffle=False); o = []
        for b in dl:
            b = b.to(device)
            with torch.no_grad(): out = model(b.to_dict())
            o.append(torch_tools.to_numpy(out["dipole"])); del b, out; torch.cuda.empty_cache()
        return np.concatenate(o, 0)
    dips = []; buf = []; n = 0
    for at in iread(os.path.join(trajdir, "prod.xyz")):
        buf.append(at)
        if len(buf) >= 100: dips.append(run(buf)); n += len(buf); buf = []
        if n and n % 2000 == 0 and len(buf) == 0: print(f"  {n} frames done", flush=True)
    if buf: dips.append(run(buf)); n += len(buf)
    M = np.concatenate(dips, 0); np.save(dipf, M)
    print(f"dipoles: {M.shape} mean|M|={np.linalg.norm(M, axis=1).mean():.3f} e.A  -> cached {dipf}")

r = spectrum_from_dipoles(M, T_K=T, dt_fs=dt_fs, n_molecules=N,
                          density_g_cm3=density_g_cm3, units="eAng", smooth_sigma_cm=smooth)
out = os.path.join(trajdir, "opacity_v3.npz")
np.savez(out, wn=r["wn_cm-1"], sigma=r["sigma_cm2_per_molecule"], sigma_raw=r["sigma_cm2_per_molecule_raw"],
         alpha=r["alpha_cm-1"], kappa=r["kappa_cm2_g"],
         planck=r["planck_mean_cm2_g"], rosseland=r["rosseland_mean_cm2_g"],
         T=T, rho_kg_m3=rho_kg_m3, N=N, dt_fs=dt_fs, smooth_cm=smooth)
print(f"{base}: T={T:.0f}K rho={rho_kg_m3:.2f}  Planck={r['planck_mean_cm2_g']:.4e}  "
      f"Rosseland={r['rosseland_mean_cm2_g']:.4e} cm2/g  (smooth {smooth:.0f}cm-1)  -> {out}")
