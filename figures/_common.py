"""Shared plotting style + model helpers for the Phase-1 figures (SPEC §8).

Plain matplotlib, colorblind-safe (Okabe–Ito). Each figure is produced by one
script that takes only a config path. Sources that aren't available yet (e.g.
line lists not delivered) are skipped gracefully so figures auto-fill later.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # headless
import matplotlib.pyplot as plt  # noqa: E402

from config import Config, load_config  # noqa: E402
from eos import AquaEOS  # noqa: E402
from opacity import ConstantOpacitySource, get_source  # noqa: E402
import opacity.build_tables as _bt  # noqa: E402
from structure import StructureResult, build_structure  # noqa: E402

# Okabe–Ito colorblind-safe palette, assigned per source (SPEC §8).
SOURCE_COLOR = {
    "fp": "#0072B2",  # blue
    "hitemp": "#D55E00",  # vermillion
    "pokazatel": "#009E73",  # green
    "composite": "#CC79A7",  # purple
    "constant": "#999999",  # grey (reference)
}
SOURCE_LABEL = {
    "fp": "FP (first-principles)",
    "hitemp": "HITEMP",
    "pokazatel": "POKAZATEL",
    "composite": "Composite (LL+FP)",
    "constant": "constant κ (ref.)",
}


def set_style() -> None:
    plt.rcParams.update({
        "figure.dpi": 130,
        "savefig.dpi": 200,
        "font.size": 11,
        "axes.grid": True,
        "grid.alpha": 0.25,
        "axes.axisbelow": True,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "legend.frameon": False,
        "lines.linewidth": 2.0,
    })


def available_sources(cfg: Config, names) -> list[str]:
    """Names whose opacity source builds now (line lists skipped until delivered)."""
    ok = []
    for name in names:
        try:
            get_source(name, cfg)
            ok.append(name)
        except (FileNotFoundError, AssertionError):
            pass
    return ok


def register_constant_reference(cfg: Config, kappa_si: float = 0.1) -> str:
    """Register a constant-κ reference source under the name 'constant'."""
    key = ("constant", _bt.resolve_band(cfg), str(cfg.data_root_path))
    _bt._cache[key] = ConstantOpacitySource(kappa_si, name="constant")
    return "constant"


def run_t_int_grid(
    cfg: Config, source: str, t_eq_k: float, eos: AquaEOS, n_sample: int = 800
) -> list[StructureResult]:
    """Build converged models across the config T_int grid for one source/T_eq."""
    out = []
    for t_int in cfg.interior.t_int_k:
        out.append(build_structure(cfg, t_eq_k, float(t_int), source, eos=eos, n_sample=n_sample))
    return out


def save(fig, cfg: Config, stem: str) -> Path:
    out_dir = cfg.output_path
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / f"{stem}.png"
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    return path


def load_cfg(argv) -> Config:
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument("config", type=Path, help="run config YAML")
    return load_config(ap.parse_args(argv).config)
