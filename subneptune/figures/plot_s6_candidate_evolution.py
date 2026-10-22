"""S6: evolution and cooling-time ratio at the candidate fiducial (A44).

    python figures/plot_s6_candidate_evolution.py config/candidate.yaml

The candidate (A42: M_c = 5 M⊕, f_env = 5%, T_eq = 700 K) is the one cell of the
S4 scan where a line-list radiative layer lies inside the FP data box and the
composite closes it with every FP evaluation interpolated. This asks what that
buys over a cooling history and over the Kelvin–Helmholtz time:

  (a) the E2 comparison (composite vs HITEMP, matched initial conditions):
      interior-entropy and radius differences at 5, 7 and 10 Gyr, the FP-relevant
      epoch, and the escape-window maxima (A40);
  (b) τ_KH(composite)/τ_KH(HITEMP) at fixed interior entropy, for
      T_int,ref ∈ {30, 40} K over f_env ∈ {0.02, 0.05, 0.10};
  (c) the radiogenic T_int floor for this core (A30/D1), so the manuscript can
      say whether the candidate's T_int is reachable at all.

Cooling ages are integrated on a coarse T_int grid first; if that grid trips the
A38 age check (negative or non-monotone ages) the failure is reported and the run
repeats on the fine grid — never silently skipped.

Outputs ``outputs/S6_candidate_evolution.csv`` and
``outputs/S6_candidate_evolution.png``.
"""

from __future__ import annotations

import csv

import sys
from dataclasses import replace
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # run as script or -m

import atomicio

import matplotlib.pyplot as plt
import numpy as np

from eos import AquaEOS
from evolution import age_grid_problem, compare_histories, t_int_floor_k
from evolution.radiogenic import H_EARTH_NOW_W, radiogenic_power_w
from figures._common import band_caption, comparison_band, load_cfg, save, set_style
from structure import FP_GRID_FLOOR_SI, build_at_entropy, build_structure

AGES_GYR = (5.0, 7.0, 10.0)
ESCAPE_WINDOW_GYR = 0.1
COARSE_GRID_K = np.geomspace(12.0, 450.0, 24)  # tried first
FINE_GRID_K = np.geomspace(12.0, 450.0, 44)  # E2's grid (A29), the fallback
T_INT_REF_K = (30.0, 40.0)
F_ENV_GRID = (0.02, 0.05, 0.10)
LINE_LIST = "hitemp"


def _at_age(age_grid, values, age_gyr):
    """Interpolate on the comparison's age grid; refuse to extrapolate (A38)."""
    if not age_grid[0] <= age_gyr <= age_grid[-1]:
        return float("nan")
    return float(np.interp(age_gyr, age_grid, values))


def history_at_candidate(cfg, eos, t_eq_k, grid_choice: str = "auto"):
    """E2 comparison, coarse grid first, fine grid if the age integral misbehaves.

    ``grid_choice`` "coarse"/"fine" pins one grid instead (a robustness re-run:
    the A38 check only catches an unphysical age integral, not an under-resolved
    one, so the reported differences are worth repeating on the fine grid).
    """
    options = {"auto": (("coarse", COARSE_GRID_K), ("fine", FINE_GRID_K)),
               "coarse": (("coarse", COARSE_GRID_K),),
               "fine": (("fine", FINE_GRID_K),)}[grid_choice]
    for label, grid in options:
        h = compare_histories(cfg, "composite", LINE_LIST, t_eq_k, eos=eos,
                              t_int_grid_k=grid)
        problems = [p for p in (age_grid_problem(h.track_a), age_grid_problem(h.track_b))
                    if p]
        if not problems:
            return h, label, grid.size, []
        print(f"  [A38] the {label} {grid.size}-node T_int grid failed the age check:")
        for p in problems:
            print(f"        {p}")
        if label == options[-1][0]:
            print(f"        the {label} grid is the last one tried — reporting the run "
                  f"as unusable.")
            return h, label, grid.size, problems
        print("        repeating on the fine grid.")
    raise AssertionError("unreachable")


def tau_ratios(cfg, eos, t_eq_k, band):
    """τ_KH(composite)/τ_KH(HITEMP) at fixed interior entropy (the F3 measure)."""
    rows = []
    for t_int_ref in T_INT_REF_K:
        for f_env in F_ENV_GRID:
            c = replace(cfg, envelope=replace(cfg.envelope, f_env=f_env))
            base = build_structure(c, t_eq_k, t_int_ref, LINE_LIST, eos=eos, n_sample=400,
                                   band=band)
            src = build_at_entropy(c, base.interior_entropy_j_kg_k, t_eq_k, "composite",
                                   eos=eos, band=band, n_sample=400)
            rows.append({
                "t_int_ref_k": t_int_ref, "f_env": f_env,
                "tau_kh_hitemp_s": base.tau_kh_s, "tau_kh_composite_s": src.tau_kh_s,
                "ratio": src.tau_kh_s / base.tau_kh_s,
                "t_int_composite_k": src.t_int_k,
                "composite_clamp_fraction_all": src.clamp_fraction_all,
                "composite_clamp_reliable": src.clamp_reliable,
            })
    return rows


def main(argv=None):
    import argparse

    argv = sys.argv[1:] if argv is None else argv
    ap = argparse.ArgumentParser()
    ap.add_argument("config", type=Path, help="run config YAML")
    ap.add_argument("--grid", choices=("auto", "coarse", "fine"), default="auto",
                    help="T_int grid for the cooling ages (default: coarse, falling "
                         "back to fine if the A38 age check trips)")
    args = ap.parse_args(argv)
    cfg = load_cfg([str(args.config)])
    set_style()
    eos = AquaEOS.from_config(cfg)
    t_eq = float(cfg.irradiation.t_eq_k[0])
    band = comparison_band(cfg, ["composite", LINE_LIST])

    print(f"S6 — candidate evolution: {cfg.name}, M_c={cfg.core.mass_mearth:g} M_earth, "
          f"f_env={cfg.envelope.f_env:.0%}, T_eq={t_eq:g} K, {band_caption(band)}")
    h, grid_label, n_nodes, problems = history_at_candidate(cfg, eos, t_eq,
                                                            grid_choice=args.grid)
    ratios = tau_ratios(cfg, eos, t_eq, band)

    # ---------------- figure ----------------
    fig, (ax_h, ax_r) = plt.subplots(1, 2, figsize=(12.6, 5.2),
                                     gridspec_kw={"width_ratios": [1.6, 1.0]},
                                     constrained_layout=True)
    ax_s = ax_h.twinx()
    ax_s.grid(False)
    ax_h.plot(h.age_gyr, h.radius_frac_diff * 100, color="#0072B2", lw=2.0,
              label=r"$\Delta R_p/R_p$  (left)")
    ax_s.plot(h.age_gyr, h.entropy_frac_diff * 100, color="#009E73", lw=1.6, ls="--",
              label=r"interior $\Delta s/s$  (right)")
    for age in AGES_GYR:
        ax_h.axvline(age, color="grey", lw=0.8, ls=":")
    ax_h.axvspan(h.age_gyr.min(), ESCAPE_WINDOW_GYR, color="orange", alpha=0.12,
                 label="first 100 Myr")
    ep = h.fp_relevant_epoch_gyr()
    if ep:
        ax_h.axvspan(ep[0], ep[1], color="green", alpha=0.10,
                     label=f"RCB ρ ≥ {FP_GRID_FLOOR_SI:g} kg/m³")
    ax_h.axhline(0.0, color="k", lw=0.6)
    ax_h.set_xscale("log")
    ax_h.set_xlabel("Age [Gyr]")
    ax_h.set_ylabel(r"$\Delta R_p/R_p$, composite − HITEMP [%]")
    ax_s.set_ylabel(r"interior $\Delta s/s$, composite − HITEMP [%]")
    handles = ax_h.get_legend_handles_labels()
    extra = ax_s.get_legend_handles_labels()
    ax_h.legend(handles[0] + extra[0], handles[1] + extra[1], fontsize=8, loc="upper left")
    ax_h.set_title(f"(a) cooling history, {n_nodes}-node {grid_label} T_int grid", fontsize=9.5)

    x = np.arange(len(F_ENV_GRID))
    for k, t_int_ref in enumerate(T_INT_REF_K):
        # Plotted as a deviation from unity: the whole signal is ±0.2%, which a
        # bar chart anchored at zero would render as four identical bars.
        vals = [(r["ratio"] - 1.0) * 100 for r in ratios if r["t_int_ref_k"] == t_int_ref]
        ax_r.bar(x + (k - 0.5) * 0.35, vals, width=0.35,
                 color=["#0072B2", "#D55E00"][k],
                 label=f"$T_{{int,ref}}$ = {t_int_ref:.0f} K")
    ax_r.axhline(0.0, color="k", lw=0.8)
    ax_r.set_xticks(x, [f"{f:.2f}" for f in F_ENV_GRID])
    ax_r.set_xlabel("$f_{env}$")
    ax_r.set_ylabel(r"$\tau_{KH}$(composite) / $\tau_{KH}$(HITEMP) − 1  [%]")
    ax_r.set_title("(b) cooling-time ratio at fixed interior entropy", fontsize=9.5)
    ax_r.legend(fontsize=8.5)

    fig.suptitle(f"S6 — candidate fiducial: $M_c$={cfg.core.mass_mearth:.0f} $M_⊕$, "
                 f"$f_{{env}}$={cfg.envelope.f_env:.0%}, $T_{{eq}}$={t_eq:.0f} K, "
                 f"{band_caption(band)}", fontsize=10.5)
    path = save(fig, cfg, "S6_candidate_evolution")

    # ---------------- outputs + printed summary ----------------
    out_dir = cfg.output_path
    out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = out_dir / "S6_candidate_evolution.csv"
    with atomicio.open_write(csv_path, newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["quantity", "age_gyr_or_f_env", "t_int_ref_k", "value"])
        for age in AGES_GYR:
            ds = _at_age(h.age_gyr, h.entropy_a - h.entropy_b, age)
            dsf = _at_age(h.age_gyr, h.entropy_frac_diff, age)
            dr = _at_age(h.age_gyr, h.radius_a - h.radius_b, age)
            drf = _at_age(h.age_gyr, h.radius_frac_diff, age)
            w.writerow(["delta_entropy_j_kg_k", age, "", f"{ds:.4f}"])
            w.writerow(["delta_entropy_frac", age, "", f"{dsf:.6e}"])
            w.writerow(["delta_radius_rearth", age, "", f"{dr:.6f}"])
            w.writerow(["delta_radius_frac", age, "", f"{drf:.6e}"])
        for r in ratios:
            w.writerow(["tau_kh_ratio", r["f_env"], r["t_int_ref_k"], f"{r['ratio']:.6f}"])

    print(f"\nwrote {path}")
    print(f"wrote {csv_path}")
    print(f"\n(a) E2 comparison, composite vs {LINE_LIST.upper()} "
          f"({n_nodes}-node {grid_label} T_int grid"
          f"{'; AGE CHECK FAILED' if problems else ''}), "
          f"age range {h.age_gyr.min():.3g}-{h.age_gyr.max():.3g} Gyr")
    head = (f"    {'age [Gyr]':>9s} {'ds [J/kg/K]':>12s} {'ds/s':>10s} "
            f"{'dR_p [R_earth]':>15s} {'dR_p/R_p':>10s}")
    print(head)
    for age in AGES_GYR:
        ds = _at_age(h.age_gyr, h.entropy_a - h.entropy_b, age)
        dsf = _at_age(h.age_gyr, h.entropy_frac_diff, age)
        dr = _at_age(h.age_gyr, h.radius_a - h.radius_b, age)
        drf = _at_age(h.age_gyr, h.radius_frac_diff, age)
        flag = "" if np.isfinite(ds) else "   [outside the track's age range]"
        print(f"    {age:>9.1f} {ds:>12.3f} {dsf * 100:>9.4f}% {dr:>15.6f} "
              f"{drf * 100:>9.4f}%{flag}")
    print(f"    FP-relevant epoch (RCB ρ ≥ {FP_GRID_FLOOR_SI:g}): "
          f"{'none within the track' if ep is None else f'{ep[0]:.3f}-{ep[1]:.3f} Gyr'}")
    esc = h.escape_window_extrema(ESCAPE_WINDOW_GYR)
    cross = h.first_age_rcb_above(FP_GRID_FLOOR_SI)
    print(f"    escape window (first {ESCAPE_WINDOW_GYR * 1000:.0f} Myr, "
          f"{esc['n_grid_points']} nodes from {esc['age_min_gyr']:.2e} Gyr): "
          f"max |ΔR_p/R_p| = {esc['max_abs_radius_frac'] * 100:.3g}%, "
          f"max |Δs/s| = {esc['max_abs_entropy_frac'] * 100:.3g}%")
    print("    deepest RCB first reaches the FP floor: " + ", ".join(
        f"{n} at {('%.3f Gyr' % a) if a is not None else 'never (within the track)'}"
        for n, a in ((h.source_a, cross[h.source_a]), (h.source_b, cross[h.source_b]))))

    print(f"\n(b) τ_KH(composite)/τ_KH({LINE_LIST.upper()}) at fixed interior entropy, "
          f"T_eq={t_eq:g} K")
    head = (f"    {'T_int,ref':>9s} {'f_env':>6s} {'ratio':>9s} {'T_int(comp) [K]':>16s} "
            f"{'comp clamp_all':>14s}")
    print(head)
    for r in ratios:
        print(f"    {r['t_int_ref_k']:>9.0f} {r['f_env']:>6.2f} {r['ratio']:>9.4f} "
              f"{r['t_int_composite_k']:>16.2f} {r['composite_clamp_fraction_all']:>14.4f}")
    lo = min(r["ratio"] for r in ratios)
    hi = max(r["ratio"] for r in ratios)
    print(f"    range over the grid: {lo:.4f} - {hi:.4f}")

    # -- (c) radiogenic T_int floor at the candidate (A30/D1 calculation) --
    t_int_struct = float(min(cfg.interior.t_int_k))
    res = build_structure(cfg, t_eq, t_int_struct, "composite", eos=eos, n_sample=400,
                          band=band)
    m_c = cfg.core.mass_mearth
    print(f"\n(c) radiogenic T_int floor, M_c={m_c:g} M_earth "
          f"(Earth-like: {H_EARTH_NOW_W / 1e12:g} TW per M_earth at 4.5 Gyr, "
          f"effective halving every 6 Gyr), R_p={res.radius_p_rearth:.5f} R_earth")
    for age in (1.0, 4.5, 10.0):
        floor = t_int_floor_k(m_c, res.radius_p_m, age)
        power = radiogenic_power_w(m_c, age)
        tag = {1.0: "  (young)", 4.5: "  (present day)", 10.0: "  (old)"}[age]
        print(f"    age {age:>4.1f} Gyr: L_radio = {power / 1e12:6.2f} TW -> "
              f"T_int floor = {floor:5.2f} K{tag}")
    floor_now = t_int_floor_k(m_c, res.radius_p_m, 4.5)
    floor_1 = t_int_floor_k(m_c, res.radius_p_m, 1.0)
    print(f"    candidate structural case T_int = {t_int_struct:g} K is "
          f"{'ABOVE' if t_int_struct > floor_now else 'BELOW'} the present-day floor "
          f"({floor_now:.2f} K) and "
          f"{'ABOVE' if t_int_struct > floor_1 else 'BELOW'} the 1 Gyr floor "
          f"({floor_1:.2f} K).")
    print(f"    (v1 omits core luminosity, SPEC §4, so the solver can reach "
          f"T_int = {t_int_struct:g} K; the floor is the stated sensitivity, A30.)")
    return path


if __name__ == "__main__":
    main()
