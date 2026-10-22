"""E1 (Phase 2): radius(age) cooling tracks per opacity source (SPEC §6).

    python figures/plot_e1_radius_age.py config/fiducial.yaml

The definitive "stay puffy longer": R_p(t) for each opacity source at fixed
(M_c, f_env, T_eq). At a fixed age the more opaque (FP) planet is larger; at a
fixed radius it is older. Prints the age ratio at a reference radius and each
track's G6 energy residual. Line lists appear when their tables land.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # run as script or -m

import matplotlib.pyplot as plt
import numpy as np

from eos import AquaEOS
from evolution import cooling_time_ratio, cooling_track
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


def main(argv=None):
    cfg = load_cfg(argv if argv is not None else sys.argv[1:])
    set_style()
    eos = AquaEOS.from_config(cfg)
    t_eq = float(cfg.irradiation.t_eq_k[0])

    register_constant_reference(cfg)
    sources = available_sources(cfg, ["fp", "hitemp", "pokazatel", "composite"]) + ["constant"]
    band = comparison_band(cfg, sources)

    tracks = {s: cooling_track(cfg, s, t_eq, eos=eos, band=band) for s in sources}

    fig, ax = plt.subplots(figsize=(6.6, 5.2))
    for s, tr in tracks.items():
        order = np.argsort(tr.age_gyr)
        reliable = tr.energy_residual < 1e-2
        ax.plot(tr.age_gyr[order], tr.radius_p_rearth[order], color=SOURCE_COLOR[s],
                label=f"{SOURCE_LABEL[s]}  (G6={tr.energy_residual:.1e})")

    ax.set_xscale("log")
    ax.set_xlabel("Age [Gyr]")
    ax.set_ylabel("$R_p$ [$R_⊕$]")
    ax.set_xlim(1e-3, max(tr.age_gyr.max() for tr in tracks.values()))
    ax.set_title(f"Cooling tracks — $M_c$={cfg.core.mass_mearth:.0f} $M_⊕$, "
                 f"$f_{{env}}$={cfg.envelope.f_env:.0%}, $T_{{eq}}$={t_eq:.0f} K\n{band_caption(band)}",
                 fontsize=10)
    ax.legend(fontsize=8.5)
    path = save(fig, cfg, "E1_radius_age")
    print(f"wrote {path}")

    # Headline age ratio at a reference radius both real sources reach.
    if "fp" in tracks:
        lls = [s for s in ("hitemp", "pokazatel") if s in tracks]
        for ll in lls:
            f, g = tracks["fp"], tracks[ll]
            r_ref = 0.5 * (max(f.radius_p_rearth.min(), g.radius_p_rearth.min())
                           + min(f.radius_p_rearth.max(), g.radius_p_rearth.max()))
            print(f"age_FP / age_{ll.upper()} at R_p={r_ref:.3f}: "
                  f"{cooling_time_ratio(f, g, r_ref):.2f}")
    return path


if __name__ == "__main__":
    main()
