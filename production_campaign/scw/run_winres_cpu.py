import torch, runpy, sys
_jl = torch.jit.load
torch.jit.load = lambda f, *a, **k: _jl(f, *a, **{**k, "map_location": "cpu"})
sys.argv = ["window_residual_spectrum.py", "--device", "cpu"]
runpy.run_path("window_residual_spectrum.py", run_name="__main__")
