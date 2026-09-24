"""S9: cooling under composite vs line-list opacity at the strongest cases (A56).

    python figures/plot_s9_case_evolution.py config/fiducial.yaml

Task 18. E2's history comparison, run at the three cells the extended domain scan
says the composite changes most, plus the old fiducial (T_eq = 500 K,
T_int = 20 K) rerun on the new 48-point tables so the manuscript can state how
the previously reported numbers moved.

Per case: Δs/s and ΔR_p/R_p at 1, 5, 7 and 10 Gyr, the cooling-time ratio to a
reference entropy (the age each track takes to cool to the same interior
entropy — the "stays hot longer" number, entropy being the state variable the
opacity actually biases, A29), and the escape-window (first 100 Myr) maxima.

The age integral is its own systematic: A44 found a 24- vs 44-node grid flipped
the sign of a sub-percent cooling difference at 5 Gyr. So every case is run at
44 nodes and again at 88, and the 5/7 Gyr numbers are reported side by side. A
difference that does not survive that doubling is not a result and is labelled
so.

Writes ``outputs/S9_case_evolution.csv`` and ``outputs/S9_case_evolution.png``.
"""

from __future__ import annotations

import csv
import sys
from dataclasses import replace
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # run as script or -m

import matplotlib.pyplot as plt
import numpy as np

from eos import AquaEOS
from evolution import compare_histories
from figures._common import load_cfg, save, set_style

N_CASES = 3
N_NODES = 44  # A29/E2's converged grid; doubled to 88 for the convergence check
T_INT_RANGE_K = (12.0, 450.0)
REPORT_AGES_GYR = (1.0, 5.0, 7.0, 10.0)
ESCAPE_WINDOW_GYR = 0.1
SOURCES = ("composite", "hitemp")

# The old fiducial, rerun on the new tables so the change is on the record.
OLD_FIDUCIAL = {"label": "old fiducial", "t_eq_k": 500.0, "t_int_k": 20.0,
                "wmf": 0.05, "mass_core_mearth": 5.0}


def _age_at_entropy(track, entropy) -> float:
    """Age at which a track's interior entropy falls to ``entropy``.

    Entropy decreases with age, so the interpolation is done on the reversed,
    ascending-entropy view. NaN if the track never reaches it.
    """
    s = np.asarray(track.entropy_j_kg_k, float)
    age = np.asarray(track.age_gyr, float)
    idx = np.argsort(s)
    s_sorted, age_sorted = s[idx], age[idx]
    if entropy < s_sorted.min() or entropy > s_sorted.max():
        return float("nan")
    return float(np.interp(entropy, s_sorted, age_sorted))


def cooling_time_ratio_at_entropy(hist, reference_age_gyr: float = 5.0) -> dict:
    """Age ratio to reach the line list's entropy at ``reference_age_gyr``.

    > 1 means the composite (FP-informed) track takes longer to cool to the same
    interior entropy, i.e. the planet stays hot longer under FP opacity.
    """
    s_ref = float(np.interp(reference_age_gyr, hist.age_gyr, hist.entropy_b))
    age_a = _age_at_entropy(hist.track_a, s_ref)
    age_b = _age_at_entropy(hist.track_b, s_ref)
    bound = ""
    if np.isnan(age_a):
        # The composite is the more opaque source, so it cools more slowly and may
        # never reach the line list's reference entropy within the T_int grid.
        # That is a result — "longer than the grid spans" — not a missing number,
        # so it is reported as a LOWER BOUND rather than a NaN. The ratio is then
        # re-posed the other way round, on the composite's own entropy, which the
        # faster-cooling line list does reach (A56).
        s_ref = float(np.interp(reference_age_gyr, hist.age_gyr, hist.entropy_a))
        age_a = reference_age_gyr
        age_b = _age_at_entropy(hist.track_b, s_ref)
        bound = "inverted: age(HITEMP) to the composite's entropy"
    return {
        "reference_age_gyr": reference_age_gyr,
        "reference_entropy_j_kg_k": s_ref,
        "age_composite_gyr": age_a,
        "age_hitemp_gyr": age_b,
        "cooling_time_ratio": (age_a / age_b
                               if age_b and not np.isnan(age_b) else float("nan")),
        "note": bound,
    }


def run_case(cfg, eos, spec, n_nodes: int) -> dict:
    """One case at one age-grid resolution."""
    cfg_c = replace(
        cfg,
        core=replace(cfg.core, mass_mearth=float(spec["mass_core_mearth"])),
        envelope=replace(cfg.envelope, f_env=float(spec["wmf"])),
    )
    grid = np.geomspace(*T_INT_RANGE_K, n_nodes)
    hist = compare_histories(cfg_c, *SOURCES, float(spec["t_eq_k"]), eos=eos,
                             t_int_grid_k=grid)
    out = {"spec": spec, "n_nodes": n_nodes, "hist": hist}
    for age in REPORT_AGES_GYR:
        if age < hist.age_gyr.min() or age > hist.age_gyr.max():
            out[f"dR_{age:g}"] = float("nan")
            out[f"ds_{age:g}"] = float("nan")
            continue
        out[f"dR_{age:g}"] = float(np.interp(age, hist.age_gyr, hist.radius_frac_diff))
        out[f"ds_{age:g}"] = float(np.interp(age, hist.age_gyr, hist.entropy_frac_diff))
    out["escape"] = hist.escape_window_extrema(ESCAPE_WINDOW_GYR)
    out["cooling"] = cooling_time_ratio_at_entropy(hist)
    out["age_range_gyr"] = (float(hist.age_gyr.min()), float(hist.age_gyr.max()))

    # Reachable bounds (A56, A66). A planet stops cooling at its radiogenic floor,
    # so a difference that only accumulates after that age is not one any planet
    # has. The reachable age is therefore min(floor crossing, end of track).
    #
    # Two numbers, because they are not the same and only one of them is the
    # bound A56 quotes. The difference is NOT monotonic in age — the water worlds
    # carry a large early transient — so the maximum over the reachable window
    # (2.1% radius for M_c = 15) is not the difference the planet ends up with
    # (0.08%). The bound is the value AT the reachable age; the maximum is
    # reported beside it as the largest excursion while still cooling.
    # Three things can end the reachable window, and the earliest wins:
    #   - the radiogenic floor, where the planet stops cooling;
    #   - the end of the integrated track;
    #   - the reporting horizon. The tracks run to hundreds of Gyr (the T_int grid
    #     reaches 12 K), so "end of track" is not an age any planet has; the
    #     horizon is the oldest age the manuscript quotes.
    floor_age = hist.radiogenic_floor_age_gyr(float(spec["mass_core_mearth"]))
    track_end = float(hist.age_gyr.max())
    reachable_age = min(track_end, max(REPORT_AGES_GYR))
    if floor_age is not None:
        reachable_age = min(reachable_age, float(floor_age))
    reach = hist.age_gyr <= reachable_age
    out["floor_age_gyr"] = float(floor_age) if floor_age is not None else float("nan")
    out["floor_reached"] = floor_age is not None
    out["reachable_age_gyr"] = reachable_age
    out["track_end_gyr"] = track_end
    out["dR_reachable"] = float(np.interp(reachable_age, hist.age_gyr,
                                          hist.radius_frac_diff))
    out["ds_reachable"] = float(np.interp(reachable_age, hist.age_gyr,
                                          hist.entropy_frac_diff))
    out["dR_max_while_cooling"] = (float(np.nanmax(np.abs(hist.radius_frac_diff[reach])))
                                   if reach.any() else float("nan"))
    out["ds_max_while_cooling"] = (float(np.nanmax(np.abs(hist.entropy_frac_diff[reach])))
                                   if reach.any() else float("nan"))
    return out


def _label(spec) -> str:
    base = spec.get("label", "")
    return (f"{base + ' ' if base else ''}M_c={spec['mass_core_mearth']:.0f} "
            f"wmf={spec['wmf']:.2f} T_eq={spec['t_eq_k']:.0f} "
            f"T_int={spec['t_int_k']:.0f}")


def strongest_specs(cfg, n: int = N_CASES) -> list[dict]:
    """The ``n`` scan cells with the largest |ΔR_p|, as evolution specs."""
    import figures.plot_s8_rcb_flip as s8

    csv_path = cfg.output_path / "S4_domain_scan_extended.csv"
    if not csv_path.exists():
        print(f"  (no {csv_path}; running the old fiducial only)")
        return []
    out = []
    for r in s8.strongest_cases(csv_path, n):
        out.append({"label": "", "t_eq_k": r["t_eq_k"], "t_int_k": r["t_int_k"],
                    "wmf": r["wmf"], "mass_core_mearth": r["mass_core_mearth"]})
    return out


def _write_csv(csv_path, results) -> None:
    """Write every completed case; called after each one so partial runs survive."""
    fields = (["case", "n_nodes", "age_min_gyr", "age_max_gyr"]
              + [f"dR_{a:g}" for a in REPORT_AGES_GYR]
              + [f"ds_{a:g}" for a in REPORT_AGES_GYR]
              + ["cooling_time_ratio", "reference_entropy_j_kg_k",
                 "escape_max_abs_radius_frac", "escape_max_abs_entropy_frac",
                 "floor_age_gyr", "floor_reached", "reachable_age_gyr",
                 "track_end_gyr",
                 "dR_reachable", "ds_reachable",
                 "dR_max_while_cooling", "ds_max_while_cooling"])
    with csv_path.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(fields)
        for r in results:
            w.writerow(
                [_label(r["spec"]), r["n_nodes"], f"{r['age_range_gyr'][0]:.6g}",
                 f"{r['age_range_gyr'][1]:.6g}"]
                + [f"{r[f'dR_{a:g}']:.6e}" for a in REPORT_AGES_GYR]
                + [f"{r[f'ds_{a:g}']:.6e}" for a in REPORT_AGES_GYR]
                + [f"{r['cooling']['cooling_time_ratio']:.6f}",
                   f"{r['cooling']['reference_entropy_j_kg_k']:.6f}",
                   f"{r['escape']['max_abs_radius_frac']:.6e}",
                   f"{r['escape']['max_abs_entropy_frac']:.6e}",
                   f"{r['floor_age_gyr']:.6g}", r["floor_reached"],
                   f"{r['reachable_age_gyr']:.6g}",
                   f"{r['track_end_gyr']:.6g}",
                   f"{r['dR_reachable']:.6e}", f"{r['ds_reachable']:.6e}",
                   f"{r['dR_max_while_cooling']:.6e}",
                   f"{r['ds_max_while_cooling']:.6e}"])


def plot_cases(cfg, results, tag="S9_case_evolution"):
    """Δs/s and ΔR/R vs age, one column per case, both resolutions overlaid."""
    specs = [r["spec"] for r in results if r["n_nodes"] == N_NODES]
    fig, axes = plt.subplots(2, len(specs), figsize=(4.6 * len(specs), 7.0),
                             sharex=True, squeeze=False, constrained_layout=True)
    for col, spec in enumerate(specs):
        for res in [r for r in results if r["spec"] is spec]:
            h = res["hist"]
            ls = "-" if res["n_nodes"] == N_NODES else "--"
            lw = 1.4 if res["n_nodes"] == N_NODES else 1.0
            axes[0][col].plot(h.age_gyr, 100 * h.entropy_frac_diff, ls, lw=lw,
                              color="#CC79A7", label=f"{res['n_nodes']} nodes")
            axes[1][col].plot(h.age_gyr, 100 * h.radius_frac_diff, ls, lw=lw,
                              color="#0072B2", label=f"{res['n_nodes']} nodes")
        for ax in (axes[0][col], axes[1][col]):
            ax.axvspan(1e-4, ESCAPE_WINDOW_GYR, color="orange", alpha=0.12, lw=0)
            ax.axhline(0.0, color="k", lw=0.7, ls=":")
            ax.set_xscale("log")
        axes[0][col].set_title(_label(spec), fontsize=8.5)
        axes[1][col].set_xlabel("age [Gyr]")
    axes[0][0].set_ylabel("Δs/s  [%]   (composite − HITEMP)")
    axes[1][0].set_ylabel("Δ$R_p$/$R_p$  [%]")
    axes[0][0].legend(fontsize=8)
    fig.suptitle("S9 — cooling under composite (FP-informed) vs HITEMP opacity, on the "
                 "48-point FP tables\nshaded: the first 100 Myr escape window; dashed: "
                 "the doubled age grid (convergence check)", fontsize=10)
    return save(fig, cfg, tag)


def main(argv=None):
    import argparse

    argv = sys.argv[1:] if argv is None else argv
    ap = argparse.ArgumentParser()
    ap.add_argument("config", type=Path)
    ap.add_argument("--n-cases", type=int, default=N_CASES)
    args = ap.parse_args(argv)
    cfg = load_cfg([str(args.config)])
    set_style()
    eos = AquaEOS.from_config(cfg)

    # Cheapest case first: the old fiducial is M_c = 5, f_env = 0.05, so it
    # converges in a fraction of the time the wmf = 0.50 water worlds take.
    specs = [OLD_FIDUCIAL] + strongest_specs(cfg, args.n_cases)
    print(f"S9 — cooling comparison (A56) at {len(specs)} cases "
          f"(the {args.n_cases} strongest scan cells + the old fiducial),")
    print(f"  composite vs HITEMP, blend window "
          f"{cfg.opacity.blend_window_kg_m3}, T_int grid "
          f"{T_INT_RANGE_K[0]:g}-{T_INT_RANGE_K[1]:g} K at {N_NODES} and "
          f"{2 * N_NODES} nodes")

    import time

    out_dir = cfg.output_path
    out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = out_dir / "S9_case_evolution.csv"

    results = []
    t_start = time.time()
    total = len(specs) * 2
    for spec in specs:
        for n_nodes in (N_NODES, 2 * N_NODES):
            t0 = time.time()
            print(f"    [{len(results) + 1}/{total}] {_label(spec)} at {n_nodes} "
                  f"nodes ...", flush=True)
            results.append(run_case(cfg, eos, spec, n_nodes))
            # Flush after every case: a stall in a later one must never cost the
            # ones already paid for (this run took hours before, A56).
            _write_csv(csv_path, results)
            print(f"        done in {time.time() - t0:.0f}s "
                  f"(elapsed {time.time() - t_start:.0f}s) -> {csv_path}", flush=True)

    head = (f"\n  {'case':<52s} {'nodes':>5s} | "
            + " ".join(f"{'ΔR@%g' % a:>9s}" for a in REPORT_AGES_GYR) + " | "
            + " ".join(f"{'Δs@%g' % a:>9s}" for a in REPORT_AGES_GYR))
    print(head)
    print("  " + "-" * (len(head) - 3))
    for r in results:
        print(f"  {_label(r['spec']):<52s} {r['n_nodes']:>5d} | "
              + " ".join(f"{100 * r[f'dR_{a:g}']:>8.4f}%" for a in REPORT_AGES_GYR)
              + " | "
              + " ".join(f"{100 * r[f'ds_{a:g}']:>8.4f}%" for a in REPORT_AGES_GYR))

    print(f"\n  CONVERGENCE — does the {N_NODES}->{2 * N_NODES} node doubling preserve "
          f"the 5 and 7 Gyr numbers, in sign and size? (A44 found it did not, on "
          f"the 27-point opacity grid)")
    head = (f"  {'case':<52s} | {'ΔR@5':>10s} {'ΔR@5 (2x)':>10s} {'ΔR@7':>10s} "
            f"{'ΔR@7 (2x)':>10s} | {'verdict':<28s}")
    print(head)
    print("  " + "-" * (len(head) - 2))
    for spec in specs:
        pair = [r for r in results if r["spec"] is spec]
        lo = next(r for r in pair if r["n_nodes"] == N_NODES)
        hi = next(r for r in pair if r["n_nodes"] == 2 * N_NODES)
        signs_hold = all(np.sign(lo[f"dR_{a:g}"]) == np.sign(hi[f"dR_{a:g}"])
                         for a in (5.0, 7.0)
                         if np.isfinite(lo[f"dR_{a:g}"]) and np.isfinite(hi[f"dR_{a:g}"]))
        close = all(abs(lo[f"dR_{a:g}"] - hi[f"dR_{a:g}"])
                    <= 0.25 * max(abs(lo[f"dR_{a:g}"]), 1e-12)
                    for a in (5.0, 7.0)
                    if np.isfinite(lo[f"dR_{a:g}"]) and np.isfinite(hi[f"dR_{a:g}"]))
        # A difference at the 1e-4 % level is a null, not a disagreement: the
        # relative test would otherwise trip on integration noise.
        null = all(abs(lo[f"dR_{a:g}"]) < 1e-6 and abs(hi[f"dR_{a:g}"]) < 1e-6
                   for a in (5.0, 7.0)
                   if np.isfinite(lo[f"dR_{a:g}"]) and np.isfinite(hi[f"dR_{a:g}"]))
        verdict = ("null (|ΔR/R| < 1e-4 %), converged" if null
                   else "converged" if signs_hold and close
                   else "sign holds, size does not" if signs_hold
                   else "SIGN FLIPS — not a result")
        print(f"  {_label(spec):<52s} | {100 * lo['dR_5']:>9.4f}% "
              f"{100 * hi['dR_5']:>9.4f}% {100 * lo['dR_7']:>9.4f}% "
              f"{100 * hi['dR_7']:>9.4f}% | {verdict:<28s}")

    print(f"\n  COOLING-TIME RATIO to the line list's 5 Gyr interior entropy "
          f"(>1 = stays hot longer under FP-informed opacity):")
    for r in results:
        if r["n_nodes"] != N_NODES:
            continue
        c = r["cooling"]
        print(f"    {_label(r['spec']):<52s} s_ref = "
              f"{c['reference_entropy_j_kg_k']:>8.1f} J/kg/K   composite "
              f"{c['age_composite_gyr']:>7.3f} Gyr / HITEMP "
              f"{c['age_hitemp_gyr']:>7.3f} Gyr  =  "
              f"{c['cooling_time_ratio']:.4f}")

    print(f"\n  ESCAPE WINDOW — maxima over the first {ESCAPE_WINDOW_GYR * 1000:.0f} Myr "
          f"(A40's clean null, re-measured on the new tables):")
    for r in results:
        if r["n_nodes"] != N_NODES:
            continue
        e = r["escape"]
        print(f"    {_label(r['spec']):<52s} max |ΔR/R| "
              f"{100 * e['max_abs_radius_frac']:.3e} %   max |Δs/s| "
              f"{100 * e['max_abs_entropy_frac']:.3e} %   "
              f"({e['n_grid_points']} grid points in window)")

    _write_csv(csv_path, results)
    path = plot_cases(cfg, results)
    print(f"\nwrote {path}")
    print(f"wrote {csv_path}")
    return path


if __name__ == "__main__":
    main()
