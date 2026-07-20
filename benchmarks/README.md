# benchmarks/

Rerunnable benchmarks with pinned reference values. Anyone
should be able to clone, point at the models, and reproduce these numbers.

```bash
python run_benchmarks.py --model ../models/water_scw_v4.model \
                         --labels reference/train_v4.xyz.gz
```

| Benchmark | Reference | Tolerance |
|---|---|---|
| Force parity (bulk regime) | 41 meV/A RMSE | +/- 15% |
| Force parity (RCB regime) | 58 meV/A RMSE | +/- 15% |
| EOS point 673 K, 1000 kg/m^3 | IAPWS-95 pressure | +/- 10% |
| Timing | parity eval ~15 min / 1080 frames (1 GPU, fp64); LDA label SP median 43 s (48 cores, 600 Ry); MACE MD 25k steps ~15-20 min (1 GPU, fp32) | recorded, no gate |

`reference/train_v4.xyz.gz` is the full v4 labelled training set. Reference
RMSEs are training-set parity by regime, matching the original evaluation
methodology (eval_split.sh); optimistic vs true held-out error. T>3000 K frames are excluded (outside the applied envelope).

The EOS benchmark runs 10 ps NVT (CSVR/Bussi thermostat, tau 100 fs) at 673 K / 1000 kg/m^3 and compares mean
pressure (virial + ideal-gas kinetic term) against IAPWS-95.
