# Step E: on-trajectory anchor check (the "skip step D?" arbiter)

**Question.** Parity (12.3/20.2 meV/A) was measured on foundation-model-
proposed geometries. Does v1 stay that accurate on ITS OWN trajectories?
If yes, the re-proposal iteration (step D) buys nothing - skip it and go
straight to production. If on-trajectory errors balloon, the ensembles
differ and D earns its cost.

**Design.** Short v1 MD (5 ps CSVR after 1 ps settle) at 6 spanning
points (dilute-cold, mid, dense branch, all three 3000 K extremes),
3 frames each -> 18 LDA single points -> force RMSE vs v1 predictions.

**Gate GE.1 (pre-registered):** per-regime on-trajectory RMSE <= 1.5x
the parity value (bulk <= 18.5, rcb <= 30.3 meV/A) -> SKIP step D,
proceed to F. Above 1.5x -> run step D as designed. Log either way.

Cost: MD ~30 GPU-min total; labels 18 x ~43 s x 48 cores ~ 10 core-h.

## Run order
    # 1. MD + frame extraction (interactive GPU session, Anaconda3/2025.06-1)
    python3 run_anchor_md.py --model .../models/water_lda_v1_swa.model \
        --workdir /data/phys-plasma-dft/shug7909/lda_pipeline
    # 2. CP2K inputs (login node fine)
    python3 make_anchor_labels.py --workdir ... \
        --template ~/lda_fast_track/stepB_label/label_template.inp
    # 3. submit anchor_labels.slurm (edit account/partition if needed)
    # 4. verdict
    python3 eval_anchor.py --workdir ... \
        --model .../models/water_lda_v1_swa.model
