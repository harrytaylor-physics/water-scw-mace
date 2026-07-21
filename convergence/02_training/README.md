# 02 — training-data convergence (learning curve)

**What:** is the dataset big enough? Train on nested subsets, plot
validation force RMSE vs dataset size.
**Criterion:** plateau — RMSE at 100% within noise of 75%. Still falling
steeply at 100% => data-limited: label more.
**Also:** seed variance (>=2 seeds at one size) separates data-limited
from noise-limited.

```
python3 learning_curve.py train.xyz --fractions 0.25 0.5 0.75 1.0 --seed 7
# trains N models via mace_run_train (GPU strongly advised), then:
python3 learning_curve.py train.xyz --analyze runs/
```

**Reference results (this project):** v4 converged at 1193 configurations
(fine-tuned from mace-mpa-0-medium); lda_v1 used 1080 (40/point x 27).
The fine-tuning regime is data-efficient: hundreds-to-a-thousand labels
suffice because the model learns a correction to the foundation model,
not water from scratch.
