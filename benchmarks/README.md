# benchmarks/

Rerunnable benchmarks with pinned reference values. Anyone in the group
should be able to clone, point at the models, and reproduce these numbers.

```bash
python run_benchmarks.py --model ../models/water_scw_v4.model \
                         --labels reference/holdout_v4.xyz
```

| Benchmark | Reference | Tolerance |
|---|---|---|
| Force parity (bulk regime) | 41 meV/A RMSE | +/- 15% |
| Force parity (RCB regime) | 58 meV/A RMSE | +/- 15% |
| Force parity (hot, pooled) | 119 meV/A RMSE | +/- 15% |
| EOS point 673 K, 1000 kg/m^3 | IAPWS-95 pressure | +/- 10% |
| Timing (192 atoms, GPU fp32) | recorded | — |

`reference/holdout_v4.xyz` is a held-out labelled split from train_v4 —
copy from ARC (or regenerate with the campaign's eval_split.sh). Until it
is present, the parity benchmark reports SKIPPED, not PASS.

The EOS benchmark runs 10 ps NVT at 673 K / 1000 kg/m^3 and compares mean
pressure (virial + ideal-gas kinetic term — the term whose omission was a
real bug in the original workflow) against IAPWS-95.
