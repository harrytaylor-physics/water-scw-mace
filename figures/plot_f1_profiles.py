"""F1: T(P) profiles overlaid per opacity source, RCB marked (SPEC §8).

    python figures/plot_f1_profiles.py config/fiducial.yaml

Overlays each available opacity source (FP now; HITEMP/POKAZATEL/composite fill
in when their tables land) plus a constant-κ reference, for one (T_eq, T_int).
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # run as script or -m

from eos import AquaEOS
from figures._common import (
    SOURCE_COLOR,
    SOURCE_LABEL,
    available_sources,
    band_caption,
    comparison_band,
    load_cfg,
    register_constant_reference,
    save,
    set_style,
)
from structure import build_structure

import matplotlib.pyplot as plt


def main(argv=None):
    cfg = load_cfg(argv if argv is not None else sys.argv[1:])
    set_style()
    eos = AquaEOS.from_config(cfg)

    t_eq = float(cfg.irradiation.t_eq_k[0])
    t_int_list = cfg.interior.t_int_k
    t_int = float(t_int_list[len(t_int_list) // 2])  # representative

    register_constant_reference(cfg)
    sources = available_sources(cfg, ["fp", "hitemp", "pokazatel", "composite"]) + ["constant"]
    band = comparison_band(cfg, sources)  # one exact band for every curve (A21)

    fig, ax = plt.subplots(figsize=(6.2, 5.2))
    for name in sources:
        res = build_structure(cfg, t_eq, t_int, name, eos=eos, band=band)
        p_bar = res.profile["pressure_pa"] / 1e5
        temp = res.profile["temperature_k"]
        flagged = "" if res.clamp_reliable else "  [clamp>thr]"
        ax.plot(temp, p_bar, color=SOURCE_COLOR[name],
                label=f"{SOURCE_LABEL[name]}  (R_p={res.radius_p_rearth:.2f} R⊕){flagged}")
        # RCB marker(s)
        for i, p_rcb in enumerate(res.rcb_pressures_pa):
            ax.plot(
                _t_at_pressure(res, p_rcb), p_rcb / 1e5, "o",
                color=SOURCE_COLOR[name], ms=6, mec="k", mew=0.5,
                label="RCB" if (name == sources[0] and i == 0) else None,
            )

    ax.set_yscale("log")
    ax.invert_yaxis()  # pressure increases downward (into the planet)
    ax.set_xlabel("Temperature [K]")
    ax.set_ylabel("Pressure [bar]")
    ax.set_title(f"T(P) profiles — $T_{{eq}}$={t_eq:.0f} K, $T_{{int}}$={t_int:.0f} K, "
                 f"$M_c$={cfg.core.mass_mearth:.0f} $M_⊕$, $f_{{env}}$={cfg.envelope.f_env:.0%}\n"
                 f"{band_caption(band)}", fontsize=10)
    ax.legend(fontsize=8.5, loc="upper left")
    path = save(fig, cfg, "F1_tp_profiles")
    print(f"wrote {path}")
    return path


def _t_at_pressure(res, p_pa):
    import numpy as np

    p = res.profile["pressure_pa"]
    t = res.profile["temperature_k"]
    return float(np.interp(np.log(p_pa), np.log(p), t))


if __name__ == "__main__":
    main()
