# 01 — electronic-structure convergence (CP2K)

**What:** plane-wave CUTOFF (then REL_CUTOFF, then basis) for your system.
**Criterion:** |dE| <= 1 meV/atom and max |dF| <= 10 meV/A vs the highest
cutoff in the scan.
**Protocol:** one representative geometry (your densest/hardest cell is
the binding constraint), single-point ENERGY_FORCE at each cutoff,
identical everything else.

```
python3 make_cutoff_scan.py geom.xyz --cutoffs 400 500 600 700 800 1000 1200 \
    --template my_label_template.inp
# ... run the generated inputs however you like ...
python3 analyze_cutoff_scan.py scan_dir/
```

**Reference results (this project):** hybrid campaign benchmarked 600/60 Ry
with DZVP-MOLOPT-SR-GTH for warm dense water (64 molecules). The LDA
campaign inherited these values verbatim (logged decision, 2026-07) —
rerun this scan at LDA/PADE to verify the inheritance on your hardware.
