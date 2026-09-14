"""S4: where does a line-list radiative zone sit INSIDE the FP domain? (A42)

    python figures/plot_s4_domain_scan.py config/fiducial.yaml

A41 showed the fiducial's deep closure is decided by FP opacity clamped below
the FP grid's 1000 K floor, so the demonstration has to move to a case whose
line-list radiative zone lies inside the FP data box (T ≥ 1000 K and
ρ ≥ 108.796 kg/m³). This scans

    M_c = 5 M⊕ × T_eq ∈ {500, 600, 700, 800} K × T_int ∈ {20, 30, 40, 50, 75} K
        × f_env ∈ {0.05, 0.10}

building both the HITEMP and the composite structure at every cell and recording
each radiative layer's (T, ρ) extent and mass share, the in-box radiative mass
fraction, R_p, and the composite's clamp_fraction_all (A41).

The candidate fiducial is the cell with the largest HITEMP in-box radiative mass
fraction that the composite closes while staying honest — composite
clamp_fraction_all < 1%.

Outputs ``outputs/S4_domain_scan.csv`` and ``outputs/S4_domain_scan.png``.
"""

from __future__ import annotations

import csv
import sys
from dataclasses import replace
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # run as script or -m

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle

from eos import AquaEOS
from figures._common import band_caption, comparison_band, load_cfg, save, set_style
from structure import FP_GRID_FLOOR_SI, build_structure
from structure.diagnostics import in_box_radiative_mass_fraction, radiative_layers

T_EQ_GRID = (500.0, 600.0, 700.0, 800.0)
T_INT_GRID = (20.0, 30.0, 40.0, 50.0, 75.0)
F_ENV_GRID = (0.05, 0.10)
N_SAMPLE = 1500
FP_T_FLOOR_K = 1000.0  # FP grid's coldest point; with FP_GRID_FLOOR_SI it is the box
CLAMP_LIMIT = 0.01  # "honest" composite: <1% of profile samples clamped


def _layers_text(layers) -> str:
    return " | ".join(
        f"T {x['t_min_k']:.0f}-{x['t_max_k']:.0f} K, rho {x['rho_min_si']:.3g}-"
        f"{x['rho_max_si']:.3g}, m {x['mass_fraction']:.4f}"
        for x in layers
    )


def run_scan(cfg, eos=None, t_eq_grid=T_EQ_GRID, t_int_grid=T_INT_GRID,
             f_env_grid=F_ENV_GRID, n_sample: int = N_SAMPLE) -> list[dict]:
    """Build HITEMP + composite over the (T_eq, T_int, f_env) grid."""
    eos = eos or AquaEOS.from_config(cfg)
    band = comparison_band(cfg, ["composite", "hitemp"])
    rows = []
    for f_env in f_env_grid:
        cfg_f = replace(cfg, envelope=replace(cfg.envelope, f_env=f_env))
        for t_eq in t_eq_grid:
            guess = {"hitemp": None, "composite": None}
            for t_int in sorted(t_int_grid):  # ascending: R_p rises smoothly
                res = {}
                for name in ("hitemp", "composite"):
                    res[name] = build_structure(
                        cfg_f, t_eq, t_int, name, eos=eos, n_sample=n_sample, band=band,
                        rp_guess_rearth=guess[name],
                    )
                    guess[name] = res[name].radius_p_rearth
                row = {
                    "f_env": f_env, "t_eq_k": t_eq, "t_int_k": t_int,
                    "mass_core_mearth": cfg.core.mass_mearth, "n_sample": n_sample,
                }
                for name in ("hitemp", "composite"):
                    r = res[name]
                    layers = radiative_layers(r)
                    in_box = radiative_layers(r, FP_T_FLOOR_K, FP_GRID_FLOOR_SI)
                    row.update({
                        f"{name}_radius_rearth": f"{r.radius_p_rearth:.5f}",
                        f"{name}_n_radiative_layers": len(layers),
                        f"{name}_radiative_layers": _layers_text(layers),
                        f"{name}_n_rcb": len(r.rcb_densities_si),
                        f"{name}_max_radiative_rho": f"{r.max_radiative_rho_si:.1f}",
                        f"{name}_in_box_radiative_mass_fraction":
                            f"{in_box_radiative_mass_fraction(r, FP_T_FLOOR_K, FP_GRID_FLOOR_SI):.6f}",
                        f"{name}_n_in_box_layers": len(in_box),
                        f"{name}_in_box_layers": _layers_text(in_box),
                        f"{name}_layer_entirely_in_box": any(
                            x["t_min_k"] >= FP_T_FLOOR_K and x["rho_min_si"] >= FP_GRID_FLOOR_SI
                            for x in layers
                        ),
                    })
                row["composite_clamp_fraction_all"] = f"{res['composite'].clamp_fraction_all:.4f}"
                row["composite_clamp_fraction_radiative"] = (
                    f"{res['composite'].clamp_fraction_radiative:.4f}")
                row["hitemp_clamp_fraction_all"] = f"{res['hitemp'].clamp_fraction_all:.4f}"
                rows.append(row)
    return rows


_INT_FIELDS = ("hitemp_n_radiative_layers", "composite_n_radiative_layers",
               "hitemp_n_rcb", "composite_n_rcb", "hitemp_n_in_box_layers",
               "composite_n_in_box_layers", "n_sample")
_FLOAT_FIELDS = ("f_env", "t_eq_k", "t_int_k", "mass_core_mearth")
_BOOL_FIELDS = ("hitemp_layer_entirely_in_box", "composite_layer_entirely_in_box")


def load_rows(csv_path: Path) -> list[dict]:
    """Re-read a previous scan's CSV (the scan itself costs ~12 min)."""
    with Path(csv_path).open() as fh:
        rows = list(csv.DictReader(fh))
    for r in rows:
        for k in _INT_FIELDS:
            r[k] = int(r[k])
        for k in _FLOAT_FIELDS:
            r[k] = float(r[k])
        for k in _BOOL_FIELDS:
            r[k] = r[k] == "True"
    return rows


def _grid_of(rows, f_env, key, cast=float):
    out = np.full((len(T_INT_GRID), len(T_EQ_GRID)), np.nan)
    for r in rows:
        if r["f_env"] != f_env:
            continue
        i = sorted(T_INT_GRID).index(r["t_int_k"])
        j = list(T_EQ_GRID).index(r["t_eq_k"])
        out[i, j] = cast(r[key])
    return out


def main(argv=None):
    import argparse

    argv = sys.argv[1:] if argv is None else argv
    ap = argparse.ArgumentParser()
    ap.add_argument("config", type=Path, help="run config YAML")
    ap.add_argument("--from-csv", action="store_true",
                    help="re-plot from a previous outputs/S4_domain_scan.csv instead of "
                         "rebuilding 80 structures")
    args = ap.parse_args(argv)
    cfg = load_cfg([str(args.config)])
    set_style()
    eos = AquaEOS.from_config(cfg)
    band = comparison_band(cfg, ["composite", "hitemp"])

    out_dir = cfg.output_path
    out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = out_dir / "S4_domain_scan.csv"
    if args.from_csv:
        rows = load_rows(csv_path)
        print(f"re-plotting from {csv_path} ({len(rows)} cells)")
    else:
        rows = run_scan(cfg, eos=eos)
        with csv_path.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)

    # ---------------- figure ----------------
    t_ints = sorted(T_INT_GRID)
    fig, axes = plt.subplots(1, len(F_ENV_GRID), figsize=(6.2 * len(F_ENV_GRID), 5.4),
                             constrained_layout=True)
    axes = np.atleast_1d(axes)
    # The in-box fractions are parts per million of the envelope mass, so annotate
    # and scale in ppm — three decimals of a raw fraction would print 0.000.
    frac_all = np.array([float(r["hitemp_in_box_radiative_mass_fraction"]) * 1e6
                         for r in rows])
    vmax = max(frac_all.max(), 1.0)
    mesh = None
    for ax, f_env in zip(axes, F_ENV_GRID):
        grid = _grid_of(rows, f_env, "hitemp_in_box_radiative_mass_fraction") * 1e6
        comp_box = _grid_of(rows, f_env, "composite_n_in_box_layers", int)
        comp_clamp = _grid_of(rows, f_env, "composite_clamp_fraction_all")
        mesh = ax.pcolormesh(np.arange(len(T_EQ_GRID) + 1), np.arange(len(t_ints) + 1),
                             grid, cmap="viridis", vmin=0.0, vmax=vmax,
                             edgecolors="white", linewidth=0.6)
        for i in range(len(t_ints)):
            for j in range(len(T_EQ_GRID)):
                ax.text(j + 0.5, i + 0.62, f"{grid[i, j]:.0f}", ha="center", va="center",
                        fontsize=8.5,
                        color="white" if grid[i, j] < 0.6 * vmax else "black")
                if comp_box[i, j] > 0:  # composite still radiative inside the FP box
                    ax.add_patch(Rectangle((j, i), 1, 1, fill=False, hatch="//",
                                           edgecolor="#D55E00", lw=0.0))
                if comp_clamp[i, j] > CLAMP_LIMIT:
                    ax.plot(j + 0.5, i + 0.25, marker="x", ms=7, color="#D55E00", mew=1.8)
        ax.set_xticks(np.arange(len(T_EQ_GRID)) + 0.5, [f"{t:.0f}" for t in T_EQ_GRID])
        ax.set_yticks(np.arange(len(t_ints)) + 0.5, [f"{t:.0f}" for t in t_ints])
        ax.set_xlabel("$T_{eq}$ [K]")
        ax.set_ylabel("$T_{int}$ [K]")
        ax.grid(False)
        ax.set_title(f"$f_{{env}}$ = {f_env:.0%}", fontsize=10)
    fig.colorbar(mesh, ax=list(axes),
                 label="HITEMP radiative mass inside the FP box [ppm of envelope]")
    fig.suptitle(
        f"S4 — line-list radiative structure inside the FP data box "
        f"(T ≥ {FP_T_FLOOR_K:.0f} K, ρ ≥ {FP_GRID_FLOOR_SI:g} kg m⁻³): "
        f"$M_c$={cfg.core.mass_mearth:.0f} $M_⊕$, {band_caption(band)}\n"
        f"hatched: the composite still has an in-box radiative layer;  "
        f"×: composite clamp_fraction_all > {CLAMP_LIMIT:.0%} (almost everywhere — "
        f"the composite extrapolates FP unless the deep envelope is hot)", fontsize=10)
    path = save(fig, cfg, "S4_domain_scan")

    # ---------------- printed summary ----------------
    print(f"wrote {path}")
    print(f"wrote {csv_path}")
    print(f"\nS4 — domain scan, M_c={cfg.core.mass_mearth:g} M_earth, n_sample={N_SAMPLE}, "
          f"{band_caption(band)}")
    head = (f"{'f_env':>6s} {'T_eq':>5s} {'T_int':>6s} | {'HT in-box frac':>14s} "
            f"{'HT layers':>9s} {'HT max rad rho':>14s} | {'C in-box layers':>15s} "
            f"{'C clamp_all':>11s} | {'R_p HT':>8s} {'R_p C':>8s}")
    print(head)
    print("-" * len(head))
    for r in rows:
        print(f"{r['f_env']:>6.2f} {r['t_eq_k']:>5.0f} {r['t_int_k']:>6.0f} | "
              f"{float(r['hitemp_in_box_radiative_mass_fraction']):>14.5f} "
              f"{r['hitemp_n_radiative_layers']:>9d} {r['hitemp_max_radiative_rho']:>14s} | "
              f"{r['composite_n_in_box_layers']:>15d} "
              f"{float(r['composite_clamp_fraction_all']):>11.4f} | "
              f"{r['hitemp_radius_rearth']:>8s} {r['composite_radius_rearth']:>8s}")

    closed = [r for r in rows
              if r["composite_n_in_box_layers"] == 0
              and float(r["composite_clamp_fraction_all"]) < CLAMP_LIMIT
              and float(r["hitemp_in_box_radiative_mass_fraction"]) > 0.0]
    print(f"\ncells where the composite closes every in-box radiative layer with "
          f"clamp_fraction_all < {CLAMP_LIMIT:.0%}: {len(closed)} of {len(rows)}")
    if closed:
        cand = max(closed, key=lambda r: float(r["hitemp_in_box_radiative_mass_fraction"]))
        print(f"CANDIDATE FIDUCIAL: T_eq={cand['t_eq_k']:.0f} K, T_int={cand['t_int_k']:.0f} K, "
              f"f_env={cand['f_env']:.2f}, M_c={cand['mass_core_mearth']:.0f} M_earth")
        print(f"  HITEMP in-box radiative mass fraction = "
              f"{float(cand['hitemp_in_box_radiative_mass_fraction']):.5f}; "
              f"in-box layers: {cand['hitemp_in_box_layers'] or 'none'}")
        print(f"  HITEMP radiative layers: {cand['hitemp_radiative_layers']}")
        print(f"  composite: {cand['composite_n_in_box_layers']} in-box radiative layers, "
              f"clamp_fraction_all = {cand['composite_clamp_fraction_all']}, "
              f"clamp_fraction_radiative = {cand['composite_clamp_fraction_radiative']}")
        print(f"  R_p: HITEMP {cand['hitemp_radius_rearth']}, "
              f"composite {cand['composite_radius_rearth']} R_earth")
    else:
        print("CANDIDATE FIDUCIAL: none — no cell has an in-box line-list radiative zone "
              "that the composite closes with clamp_fraction_all < "
              f"{CLAMP_LIMIT:.0%}.")

    entirely = [r for r in rows if r["hitemp_layer_entirely_in_box"]]
    print(f"\ncells with a HITEMP radiative layer lying ENTIRELY inside the FP box "
          f"(no part below {FP_T_FLOOR_K:.0f} K or below ρ = {FP_GRID_FLOOR_SI:g}): "
          f"{len(entirely)} of {len(rows)}")
    for r in entirely:
        print(f"  T_eq={r['t_eq_k']:.0f}, T_int={r['t_int_k']:.0f}, f_env={r['f_env']:.2f}: "
              f"{r['hitemp_radiative_layers']}")
    return path


if __name__ == "__main__":
    main()
