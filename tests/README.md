# tests/

Pytest suite. CPU-safe; every test that needs a model file skips cleanly
if it is absent, so CI without git-lfs still runs the pure-code tests.

```bash
pytest -v                       # from repo root
WATER_SCW_DEVICE=cuda pytest    # to force GPU
```

Covers: model loading, energy-force consistency (finite difference),
translation/rotation/PBC invariance, dimer sanity (binding, no fusion),
determinism, and a fast regression on liquid structure (g_OO first peak).
