# gen_tasks.py -- write tasks.txt: every isotherm point + finite-size point, one per line.
# Submit as a single SLURM array (see run_array.sh). Edit the lists below to taste.
WHICH      = ["v4","foundation"]
ISOTHERMS  = [(673,"scw673"), (1000,"hot1000")]
ISO_RHO    = [0.05,0.10,0.20,0.30,0.45,0.66,0.80,1.00,1.20]   # g/cc
FS_POINTS  = [(673,0.66,"fs673"), (1000,0.80,"fs1000")]
FS_SIZES   = [64,512,1728]   # 1x, 8x, 27x. Add 4096 (~64x) if the GPU/walltime allow;
                             # 6400 (~100x, ~19k atoms) needs a large-memory GPU + long walltime.

lines=[]
for T,tag in ISOTHERMS:
    for w in WHICH:
        for r in ISO_RHO:
            lines.append(f"python iso_pressure.py {w} {T} {r} {tag}")
for T,r,tag in FS_POINTS:
    for N in FS_SIZES:
        lines.append(f"python finite_size.py v4 {T} {r} {N} {tag}")

open("tasks.txt","w").write("\n".join(lines)+"\n")
print(f"wrote tasks.txt with {len(lines)} tasks")
print(f"submit with:  sbatch --array=1-{len(lines)} run_array.sh")
