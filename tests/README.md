# tests/

Pytest suite. CPU-safe: anything needing a model file skips cleanly if
absent, so code-only CI still passes without git-lfs.

```bash
pytest -m "not slow"        # fast set, < 1 min on CPU with models present
pytest                      # everything (includes 3 ps MD regression)
WATER_SCW_DEVICE=cuda pytest
```

| File | Covers |
|---|---|
| test_potential.py | E-F consistency (autograd identity), dimer binding, close-approach repulsive wall (the deliberately admitted 0.65-0.75 A regime), determinism, slow structure regression |
| test_invariances.py | translation, rotation (energy + force co-rotation), permutation, supercell extensivity with L/2 < r_max (multi-image neighbours, as in the dense production cells) |
| test_dipole.py | water_dipole_v2 loads; dipole finite/sane; rotational equivariance mu(Rx)=R mu(x); translation invariance |
| test_golden.py | pinned reference outputs for both models (self-pinning on first run — commit tests/golden/); catches silent model-file or mace-version drift |
| test_units.py | number-density and opacity unit conversions, including an explicit trap for the historical x1000 bug signature |

Not yet covered (tracked in ../TODO.md): native-vs-compiled model
agreement; fp32-vs-fp64 bounds; direct tests of the opacity/spectra
pipeline once it is ported into the repo.

First run with models present: `pytest tests/test_golden.py` pins golden
values — check the numbers are sane, then commit `tests/golden/*.json`.
