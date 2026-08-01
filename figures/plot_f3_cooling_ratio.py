"""F3: Kelvin–Helmholtz cooling-time ratio vs (f_env, T_eq) (SPEC §8).

    python figures/plot_f3_cooling_ratio.py config/fiducial.yaml

Headline (left): τ_FP / τ_LL at **fixed interior entropy** — the "planets stay
puffy longer" quantification. At fixed deep-adiabat entropy the two sources have
the same structure (R_p), but the more opaque FP throttles the flux to a lower
T_int and hence lower L, giving a longer cooling time (ratio > 1).

Diagnostic (right): the same ratio at **fixed T_int**, which comes out < 1
because L = 4πR²σT_int⁴ is imposed by T_int (A19/A23). Kept to make the
control-variable distinction explicit; see docs/cooling_time_caveat.md.

Baseline auto-selects: HITEMP > POKAZATEL > constant-κ placeholder.
"""

from __future__ import annotations

import sys
from dataclasses import replace
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # run as script or -m

import matplotlib.pyplot as plt

from eos import AquaEOS
from figures._common import (
    band_caption,
    comparison_band,
    load_cfg,
    register_constant_reference,
    save,
    set_style,
)
from structure import build_at_entropy, build_structure

F_ENV_GRID = (0.01, 0.02, 0.03, 0.05, 0.10)


def main(argv=None):
    cfg = load_cfg(argv if argv is not None else sys.argv[1:])
    set_style()
    eos = AquaEOS.from_config(cfg)
    t_int_ref = float(cfg.interior.t_int_k[len(cfg.interior.t_int_k) // 2])

    # Baseline: prefer a real line list; else the constant-κ placeholder.
    from figures._common import available_sources

    avail = available_sources(cfg, ["hitemp", "pokazatel"])
    baseline = "hitemp" if "hitemp" in avail else ("pokazatel" if "pokazatel" in avail else "constant")
    if baseline == "constant":
        register_constant_reference(cfg)
    band = comparison_band(cfg, ["fp", baseline])

    fig, (ax_s, ax_t) = plt.subplots(1, 2, figsize=(11.5, 5.0), constrained_layout=True)
    cmap = plt.get_cmap("viridis")
    t_eqs = cfg.irradiation.t_eq_k

    for k, t_eq in enumerate(t_eqs):
        ratio_s, ratio_t = [], []
        for f_env in F_ENV_GRID:
            cfg_f = replace(cfg, envelope=replace(cfg.envelope, f_env=f_env))
            # Baseline defines the reference state (its entropy + its τ at T_int_ref).
            base = build_structure(cfg_f, float(t_eq), t_int_ref, baseline, eos=eos,
                                   n_sample=400, band=band)
            s_star = base.interior_entropy_j_kg_k
            tau_base = base.tau_kh_s
            # Fixed T_int (diagnostic): FP at the same T_int.
            fp_ti = build_structure(cfg_f, float(t_eq), t_int_ref, "fp", eos=eos,
                                    n_sample=400, band=band)
            ratio_t.append(fp_ti.tau_kh_s / tau_base)
            # Fixed entropy (headline): FP matched to the baseline's interior entropy.
            fp_s = build_at_entropy(cfg_f, s_star, float(t_eq), "fp", eos=eos,
                                    band=band, n_sample=400)
            ratio_s.append(fp_s.tau_kh_s / tau_base)
        color = cmap(k / max(1, len(t_eqs) - 1))
        ax_s.plot(F_ENV_GRID, ratio_s, marker="o", ms=5, color=color, label=f"$T_{{eq}}$={t_eq:.0f} K")
        ax_t.plot(F_ENV_GRID, ratio_t, marker="o", ms=5, color=color, label=f"$T_{{eq}}$={t_eq:.0f} K")

    base_lbl = {"hitemp": "HITEMP", "pokazatel": "POKAZATEL",
                "constant": "constant-κ ref (LL placeholder)"}[baseline]
    for ax in (ax_s, ax_t):
        ax.axhline(1.0, color="k", lw=0.8, ls=":")
        ax.set_xscale("log")
        ax.set_xlabel("$f_{env}$")
        ax.set_ylabel(rf"$\tau_{{FP}} / \tau_{{{base_lbl.split()[0]}}}$")
        ax.legend(fontsize=8.5, title="irradiation")

    ax_s.set_title(f"HEADLINE — fixed interior entropy\n(stay puffy longer: ratio > 1)")
    ax_t.set_title("diagnostic — fixed $T_{int}$\n(L imposed by $T_{int}$: ratio < 1, not the story)")
    if baseline in ("pokazatel",):
        ax_s.set_title(ax_s.get_title() + "  [grid overlap; clamped]")

    fig.suptitle(f"Kelvin–Helmholtz cooling-time ratio vs {base_lbl} — "
                 f"$M_c$={cfg.core.mass_mearth:.0f} $M_⊕$, $T_{{int,ref}}$={t_int_ref:.0f} K, "
                 f"{band_caption(band)}", y=1.04)
    path = save(fig, cfg, "F3_cooling_ratio")
    print(f"wrote {path}  (headline=fixed-entropy, baseline={baseline})")
    return path


if __name__ == "__main__":
    main()
