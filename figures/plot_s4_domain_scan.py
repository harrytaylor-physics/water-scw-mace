"""S4: where does a line-list radiative zone sit INSIDE the FP domain? (A42, A47)

    python figures/plot_s4_domain_scan.py config/fiducial.yaml              # A42 grid
    python figures/plot_s4_domain_scan.py config/fiducial.yaml --extended   # A47 grid

A41 showed the fiducial's deep closure was decided by FP opacity clamped below
the then 1000 K grid floor, so the demonstration had to move to a case whose
line-list radiative zone lies inside the FP data box. Since the extension
campaign that box is **T ≥ 650 K, ρ ≥ 10 kg/m³** and is read from config
(``opacity.fp_grid_*``), not hard-coded (A54). The default grid is A42's:

    M_c = 5 M⊕ × T_eq ∈ {500, 600, 700, 800} K × T_int ∈ {20, 30, 40, 50, 75} K
        × wmf ∈ {0.05, 0.10}

and ``--extended`` is A47's, which opens the planet space in every direction
that could plausibly push a radiative zone into the box — hotter and colder
irradiation, heavier cores, and water-world envelope fractions:

    T_eq ∈ {300, 500, 700, 900, 1100, 1300} K × T_int ∈ {30, 50} K
        × M_c ∈ {5, 10, 15} M⊕ × wmf ∈ {0.05, 0.20, 0.50}

``wmf`` is the envelope (water) mass fraction of the *whole* planet. In this
repository ``envelope.f_env`` is already defined that way — the model sets
M_p = M_c/(1−f_env), so f_env = M_env/M_p — hence wmf ≡ f_env and both are
recorded, equal by construction, alongside the implied M_p (A47).

Both sources are built at every cell and we record each radiative layer's (T, ρ)
extent and mass share, the radiative mass fraction inside the FP box, the total
(any-T) radiative mass fraction, the deepest/hottest envelope point, R_p, and
``clamp_fraction_all`` (A41) for *both* sources.

Cells whose structure fails to converge, or whose profile leaves the AQUA table,
are skipped with a logged reason and excluded from the rankings and the figures.

Outputs ``outputs/S4_domain_scan{,_extended}.csv`` and one figure per core mass.
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
from figures._common import (
    band_caption,
    comparison_band,
    fp_floors,
    load_cfg,
    save,
    set_style,
)
from structure import FP_GRID_FLOOR_SI, build_structure
from opacity import get_source
from structure.diagnostics import (
    aqua_coverage,
    in_box_clamp_fraction,
    in_box_radiative_mass_fraction,
    radiative_layers,
)

# --- A42 default grid -------------------------------------------------------
T_EQ_GRID = (500.0, 600.0, 700.0, 800.0)
T_INT_GRID = (20.0, 30.0, 40.0, 50.0, 75.0)
WMF_GRID = (0.05, 0.10)
# --- A47 extended grid ------------------------------------------------------
EXT_T_EQ_GRID = (300.0, 500.0, 700.0, 900.0, 1100.0, 1300.0)
EXT_T_INT_GRID = (30.0, 50.0)
EXT_M_CORE_GRID = (5.0, 10.0, 15.0)
EXT_WMF_GRID = (0.05, 0.20, 0.50)

N_SAMPLE = 1500
# Default FP box corner (schema defaults, A54). main() overrides both from the
# run config, which is authoritative; these keep the module importable alone.
FP_T_FLOOR_K = 650.0
FP_RHO_FLOOR_SI = FP_GRID_FLOOR_SI
CLAMP_LIMIT = 0.01  # "honest" composite: <1% of profile samples clamped
SOURCES = ("hitemp", "composite")


def _layers_text(layers) -> str:
    return " | ".join(
        f"T {x['t_min_k']:.0f}-{x['t_max_k']:.0f} K, rho {x['rho_min_si']:.3g}-"
        f"{x['rho_max_si']:.3g}, m {x['mass_fraction']:.4f}"
        for x in layers
    )


def _source_columns(name: str, res, floors=(FP_T_FLOOR_K, FP_RHO_FLOOR_SI),
                    source=None) -> dict:
    """Everything the scan records about one converged structure."""
    t_floor, rho_floor = floors
    layers = radiative_layers(res)
    in_box = radiative_layers(res, t_floor, rho_floor)
    prof = res.profile
    return {
        f"{name}_radius_rearth": f"{res.radius_p_rearth:.5f}",
        f"{name}_n_radiative_layers": len(layers),
        f"{name}_radiative_layers": _layers_text(layers),
        f"{name}_n_rcb": len(res.rcb_densities_si),
        f"{name}_max_radiative_rho": f"{res.max_radiative_rho_si:.1f}",
        f"{name}_max_radiative_t_k":
            f"{max((x['t_max_k'] for x in layers), default=0.0):.1f}",
        f"{name}_in_box_radiative_mass_fraction":
            f"{in_box_radiative_mass_fraction(res, t_floor, rho_floor):.6f}",
        f"{name}_total_radiative_mass_fraction":
            f"{sum(x['mass_fraction'] for x in layers):.6f}",
        f"{name}_n_in_box_layers": len(in_box),
        f"{name}_in_box_layers": _layers_text(in_box),
        f"{name}_layer_entirely_in_box": any(
            x["t_min_k"] >= t_floor and x["rho_min_si"] >= rho_floor
            for x in layers
        ),
        # deepest/hottest point the envelope reaches (A47)
        f"{name}_t_max_k": f"{float(np.max(prof['temperature_k'])):.1f}",
        f"{name}_rho_max_si": f"{float(np.max(prof['rho_si'])):.1f}",
        f"{name}_clamp_fraction_all": f"{res.clamp_fraction_all:.4f}",
        f"{name}_clamp_fraction_radiative": f"{res.clamp_fraction_radiative:.4f}",
        # Task 17's gate: clamping WHERE THE FP GRID REACHES, not overall (A55).
        # Both boxes plus the reason split, so a violation can be attributed.
        **(_clamp_columns(name, source, res, t_floor, rho_floor)),
    }


def _clamp_columns(name: str, source, res, t_floor: float, rho_floor: float) -> dict:
    """In-box clamp columns for one source, or blanks when no source was passed."""
    keys = ("in_box_clamp_fraction", "n_in_box_samples", "full_box_clamp_fraction",
            "n_full_box_samples", "clamp_too_hot", "clamp_too_dense",
            "clamp_in_box_out_of_hull")
    if source is None:
        return {f"{name}_{k}": "" if "fraction" in k else 0 for k in keys}
    d = in_box_clamp_fraction(source, res, t_floor, rho_floor)
    return {
        f"{name}_in_box_clamp_fraction": f"{d['fraction']:.4f}",
        f"{name}_n_in_box_samples": d["n_in_box"],
        f"{name}_full_box_clamp_fraction": f"{d['full_fraction']:.4f}",
        f"{name}_n_full_box_samples": d["n_full_box"],
        f"{name}_clamp_too_hot": d["too_hot"],
        f"{name}_clamp_too_dense": d["too_dense"],
        f"{name}_clamp_in_box_out_of_hull": d["in_box_out_of_hull"],
    }


def _empty_source_columns(name: str) -> dict:
    """Placeholder columns for a skipped cell, so the CSV stays rectangular."""
    return {
        f"{name}_radius_rearth": "", f"{name}_n_radiative_layers": 0,
        f"{name}_radiative_layers": "", f"{name}_n_rcb": 0,
        f"{name}_max_radiative_rho": "", f"{name}_max_radiative_t_k": "",
        f"{name}_in_box_radiative_mass_fraction": "",
        f"{name}_total_radiative_mass_fraction": "", f"{name}_n_in_box_layers": 0,
        f"{name}_in_box_layers": "", f"{name}_layer_entirely_in_box": False,
        f"{name}_t_max_k": "", f"{name}_rho_max_si": "",
        f"{name}_clamp_fraction_all": "", f"{name}_clamp_fraction_radiative": "",
        **_clamp_columns(name, None, None, 0.0, 0.0),
    }


def run_scan(cfg, eos=None, t_eq_grid=T_EQ_GRID, t_int_grid=T_INT_GRID,
             wmf_grid=WMF_GRID, m_core_grid=None, n_sample: int = N_SAMPLE,
             verbose: bool = True, floors=None) -> list[dict]:
    """Build both sources over the (M_c, wmf, T_eq, T_int) grid.

    A cell whose structure fails to converge, or whose converged profile leaves
    the AQUA table, is recorded with ``skipped = True`` and a ``skip_reason``;
    its source columns are left blank.
    """
    eos = eos or AquaEOS.from_config(cfg)
    band = comparison_band(cfg, list(SOURCES))
    floors = floors or fp_floors(cfg)
    m_core_grid = m_core_grid or (cfg.core.mass_mearth,)
    rows = []
    for m_core in m_core_grid:
        for wmf in wmf_grid:
            cfg_c = replace(
                cfg,
                core=replace(cfg.core, mass_mearth=float(m_core)),
                envelope=replace(cfg.envelope, f_env=float(wmf)),
            )
            for t_eq in t_eq_grid:
                guess = dict.fromkeys(SOURCES)  # warm start resets per (M_c, wmf, T_eq)
                for t_int in sorted(t_int_grid):  # ascending: R_p rises smoothly
                    row = {
                        "mass_core_mearth": float(m_core), "wmf": float(wmf),
                        "f_env": float(wmf),  # identical by construction (A47)
                        "mass_p_mearth": f"{float(m_core) / (1.0 - float(wmf)):.4f}",
                        "t_eq_k": float(t_eq), "t_int_k": float(t_int),
                        "n_sample": n_sample, "skipped": False, "skip_reason": "",
                    }
                    res, reason = {}, ""
                    for name in SOURCES:
                        try:
                            res[name] = build_structure(
                                cfg_c, t_eq, t_int, name, eos=eos, n_sample=n_sample,
                                band=band, rp_guess_rearth=guess[name],
                            )
                            guess[name] = res[name].radius_p_rearth
                        except (RuntimeError, ValueError) as exc:
                            reason = f"{name}: did not converge ({type(exc).__name__}: {exc})"
                            break
                        cov = aqua_coverage(eos, res[name].profile["pressure_pa"],
                                            res[name].profile["temperature_k"])
                        if not cov["inside"]:
                            reason = f"{name}: outside the AQUA table ({cov['summary']})"
                            break
                    if reason:
                        row.update({"skipped": True, "skip_reason": reason})
                        for name in SOURCES:
                            row.update(_empty_source_columns(name))
                        if verbose:
                            print(f"  SKIP M_c={m_core:g} wmf={wmf:.2f} T_eq={t_eq:g} "
                                  f"T_int={t_int:g}: {reason}")
                    else:
                        for name in SOURCES:
                            row.update(_source_columns(
                                name, res[name], floors,
                                source=get_source(name, cfg_c, band=band)))
                    rows.append(row)
    return rows


_INT_FIELDS = ("n_sample",) + tuple(
    f"{s}_{k}" for s in SOURCES
    for k in ("n_radiative_layers", "n_rcb", "n_in_box_layers", "n_in_box_samples",
              "n_full_box_samples", "clamp_too_hot", "clamp_too_dense",
              "clamp_in_box_out_of_hull"))
_FLOAT_FIELDS = ("wmf", "f_env", "t_eq_k", "t_int_k", "mass_core_mearth")
_BOOL_FIELDS = ("skipped",) + tuple(f"{s}_layer_entirely_in_box" for s in SOURCES)


def load_rows(csv_path: Path) -> list[dict]:
    """Re-read a previous scan's CSV (the extended scan itself costs ~5 min)."""
    with Path(csv_path).open() as fh:
        rows = list(csv.DictReader(fh))
    if rows and "wmf" not in rows[0]:
        raise ValueError(
            f"{csv_path} predates A47 (no 'wmf' column): re-run the scan without "
            f"--from-csv to regenerate it"
        )
    for r in rows:
        for k in _INT_FIELDS:
            r[k] = int(r[k] or 0)
        for k in _FLOAT_FIELDS:
            r[k] = float(r[k])
        for k in _BOOL_FIELDS:
            r[k] = r[k] == "True"
    return rows


def _num(row, key, default=np.nan) -> float:
    """A numeric cell, tolerating the blanks a skipped row leaves."""
    v = row[key]
    if v == "" or v is None:
        return default
    return float(v)


def _grid_of(rows, sel, key, t_eq_grid, t_int_grid, cast=float):
    """(T_int × T_eq) array of ``key`` over the rows matching ``sel`` (a dict)."""
    out = np.full((len(t_int_grid), len(t_eq_grid)), np.nan)
    t_ints, t_eqs = sorted(t_int_grid), list(t_eq_grid)
    for r in rows:
        if any(r[k] != v for k, v in sel.items()):
            continue
        i, j = t_ints.index(r["t_int_k"]), t_eqs.index(r["t_eq_k"])
        out[i, j] = np.nan if r["skipped"] else cast(_num(r, key))
    return out


def plot_grid(cfg, rows, m_core, wmf_grid, t_eq_grid, t_int_grid, band, tag,
              floors=None):
    """One figure for a core mass: a (T_eq × T_int) panel per wmf."""
    t_floor, rho_floor = floors or fp_floors(cfg)
    t_ints = sorted(t_int_grid)
    fig, axes = plt.subplots(1, len(wmf_grid),
                             figsize=(5.6 * len(wmf_grid), 5.2),
                             constrained_layout=True)
    axes = np.atleast_1d(axes)
    # In-box fractions are parts per million of the envelope mass: a raw fraction
    # would annotate every cell "0.000".
    ppm = np.array([_num(r, "hitemp_in_box_radiative_mass_fraction", 0.0) * 1e6
                    for r in rows if r["mass_core_mearth"] == m_core and not r["skipped"]])
    vmax = max(ppm.max() if ppm.size else 0.0, 1.0)
    mesh = None
    for ax, wmf in zip(axes, wmf_grid):
        sel = {"mass_core_mearth": m_core, "wmf": wmf}
        grid = _grid_of(rows, sel, "hitemp_in_box_radiative_mass_fraction",
                        t_eq_grid, t_int_grid) * 1e6
        comp_box = _grid_of(rows, sel, "composite_n_in_box_layers", t_eq_grid,
                            t_int_grid, int)
        comp_clamp = _grid_of(rows, sel, "composite_clamp_fraction_all", t_eq_grid,
                              t_int_grid)
        rad_rho = _grid_of(rows, sel, "hitemp_max_radiative_rho", t_eq_grid, t_int_grid)
        rad_t = _grid_of(rows, sel, "hitemp_max_radiative_t_k", t_eq_grid, t_int_grid)
        mesh = ax.pcolormesh(np.arange(len(t_eq_grid) + 1), np.arange(len(t_ints) + 1),
                             grid, cmap="viridis", vmin=0.0, vmax=vmax,
                             edgecolors="white", linewidth=0.6)
        for i in range(len(t_ints)):
            for j in range(len(t_eq_grid)):
                if np.isnan(grid[i, j]):  # skipped cell
                    ax.add_patch(Rectangle((j, i), 1, 1, facecolor="#BBBBBB",
                                           edgecolor="white"))
                    ax.text(j + 0.5, i + 0.5, "skip", ha="center", va="center",
                            fontsize=7.5, color="#333333")
                    continue
                color = "white" if grid[i, j] < 0.6 * vmax else "black"
                ax.text(j + 0.5, i + 0.68, f"{grid[i, j]:.0f}", ha="center", va="center",
                        fontsize=8.5, color=color)
                # How close this cell gets to the box on each axis: without it a
                # grid of honest zeros would say nothing at all (A47).
                ax.text(j + 0.5, i + 0.44,
                        f"ρ≤{rad_rho[i, j]:.0f}\nT≤{rad_t[i, j]:.0f}", ha="center",
                        va="center", fontsize=6.5, color=color, alpha=0.85)
                if comp_box[i, j] > 0:  # composite still radiative inside the FP box
                    ax.add_patch(Rectangle((j, i), 1, 1, fill=False, hatch="//",
                                           edgecolor="#D55E00", lw=0.0))
                if comp_clamp[i, j] > CLAMP_LIMIT:
                    ax.plot(j + 0.5, i + 0.25, marker="x", ms=7, color="#D55E00", mew=1.8)
        ax.set_xticks(np.arange(len(t_eq_grid)) + 0.5, [f"{t:.0f}" for t in t_eq_grid])
        ax.set_yticks(np.arange(len(t_ints)) + 0.5, [f"{t:.0f}" for t in t_ints])
        ax.set_xlabel("$T_{eq}$ [K]")
        ax.set_ylabel("$T_{int}$ [K]")
        ax.grid(False)
        ax.set_title(f"wmf = {wmf:.0%}  ($M_p$ = {m_core / (1 - wmf):.1f} $M_⊕$)",
                     fontsize=10)
    fig.colorbar(mesh, ax=list(axes),
                 label="HITEMP radiative mass inside the FP box [ppm of envelope]")
    fig.suptitle(
        f"S4 — line-list radiative structure inside the FP data box "
        f"(T ≥ {t_floor:.0f} K, ρ ≥ {rho_floor:g} kg m⁻³): "
        f"$M_c$ = {m_core:g} $M_⊕$, {band_caption(band)}\n"
        f"cell: in-box ppm (top), then how close the line-list radiative zone gets — "
        f"max radiative ρ and its max T;  hatched: the composite still has an in-box "
        f"radiative layer;\n×: composite clamp_fraction_all > {CLAMP_LIMIT:.0%};  "
        f"grey: skipped cell", fontsize=10)
    return save(fig, cfg, tag)


def _rank(rows, key, n=10):
    """The ``n`` converged rows with the largest ``key``, descending."""
    live = [r for r in rows if not r["skipped"]]
    return sorted(live, key=lambda r: _num(r, key, 0.0), reverse=True)[:n]


def _cell_label(r) -> str:
    return (f"M_c={r['mass_core_mearth']:>4.1f} wmf={r['wmf']:.2f} "
            f"T_eq={r['t_eq_k']:>6.0f} T_int={r['t_int_k']:>4.0f}")


def print_rankings(rows, floors, n=10) -> None:
    """The two top-ten tables the extended scan is for (A47)."""
    t_floor, rho_floor = floors
    print(f"\ntop {n} cells by HITEMP radiative mass fraction INSIDE the FP box "
          f"(T ≥ {t_floor:.0f} K, ρ ≥ {rho_floor:g} kg/m³):")
    head = (f"  {'cell':<44s} {'HT in-box':>10s} {'C in-box':>9s} {'C clamp_all':>11s} "
            f"{'composite closes it?':<22s}")
    print(head)
    print("  " + "-" * (len(head) - 2))
    for r in _rank(rows, "hitemp_in_box_radiative_mass_fraction", n):
        ht = _num(r, "hitemp_in_box_radiative_mass_fraction", 0.0)
        cb = r["composite_n_in_box_layers"]
        cl = _num(r, "composite_clamp_fraction_all", np.nan)
        if ht <= 0.0:
            verdict = "n/a (nothing in box)"
        elif cb == 0 and cl < CLAMP_LIMIT:
            verdict = "YES, on FP data"
        elif cb == 0:
            verdict = f"yes, but clamped {cl:.0%}"
        else:
            verdict = f"no ({cb} layer(s) left)"
        print(f"  {_cell_label(r):<44s} {ht:>10.6f} {cb:>9d} {cl:>11.4f} {verdict:<22s}")

    print(f"\ntop {n} cells by TOTAL (any-T) HITEMP radiative mass fraction, with each "
          f"layer's extent — where the radiative layers live relative to the FP box:")
    for r in _rank(rows, "hitemp_total_radiative_mass_fraction", n):
        print(f"  {_cell_label(r)}  total radiative m = "
              f"{_num(r, 'hitemp_total_radiative_mass_fraction', 0.0):.4f}, "
              f"max radiative ρ = {r['hitemp_max_radiative_rho']}, "
              f"envelope reaches T = {r['hitemp_t_max_k']} K, ρ = {r['hitemp_rho_max_si']}")
        for layer in (r["hitemp_radiative_layers"] or "none").split(" | "):
            print(f"      {layer}")


def main(argv=None):
    import argparse

    argv = sys.argv[1:] if argv is None else argv
    ap = argparse.ArgumentParser()
    ap.add_argument("config", type=Path, help="run config YAML")
    ap.add_argument("--extended", action="store_true",
                    help="A47 grid: T_eq 300-1300 K x T_int {30,50} x M_c {5,10,15} "
                         "x wmf {0.05,0.20,0.50} (108 cells, ~5 min)")
    ap.add_argument("--from-csv", action="store_true",
                    help="re-plot and re-rank from a previous CSV instead of rebuilding")
    args = ap.parse_args(argv)
    cfg = load_cfg([str(args.config)])
    set_style()
    eos = AquaEOS.from_config(cfg)
    band = comparison_band(cfg, list(SOURCES))
    floors = fp_floors(cfg)
    t_floor, rho_floor = floors

    if args.extended:
        t_eq_grid, t_int_grid = EXT_T_EQ_GRID, EXT_T_INT_GRID
        wmf_grid, m_core_grid = EXT_WMF_GRID, EXT_M_CORE_GRID
        stem = "S4_domain_scan_extended"
    else:
        t_eq_grid, t_int_grid = T_EQ_GRID, T_INT_GRID
        wmf_grid, m_core_grid = WMF_GRID, (cfg.core.mass_mearth,)
        stem = "S4_domain_scan"

    out_dir = cfg.output_path
    out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = out_dir / f"{stem}.csv"
    if args.from_csv:
        rows = load_rows(csv_path)
        print(f"re-plotting from {csv_path} ({len(rows)} cells)")
    else:
        print(f"scanning {len(m_core_grid) * len(wmf_grid) * len(t_eq_grid) * len(t_int_grid)}"
              f" cells x {len(SOURCES)} sources, n_sample={N_SAMPLE} ...")
        rows = run_scan(cfg, eos=eos, t_eq_grid=t_eq_grid, t_int_grid=t_int_grid,
                        wmf_grid=wmf_grid, m_core_grid=m_core_grid, floors=floors)
        with csv_path.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)

    # ---------------- figures: one per core mass (A47) ----------------
    paths = [plot_grid(cfg, rows, m_core, wmf_grid, t_eq_grid, t_int_grid, band,
                       f"{stem}_Mc{m_core:g}" if args.extended else stem,
                       floors=floors)
             for m_core in m_core_grid]

    # ---------------- printed summary ----------------
    for p in paths:
        print(f"wrote {p}")
    print(f"wrote {csv_path}")
    live = [r for r in rows if not r["skipped"]]
    skipped = [r for r in rows if r["skipped"]]
    print(f"\nS4 — domain scan ({'A47 extended' if args.extended else 'A42 default'} grid), "
          f"n_sample={N_SAMPLE}, {band_caption(band)}")
    print(f"  {len(live)} of {len(rows)} cells converged inside the AQUA table; "
          f"{len(skipped)} skipped")
    for r in skipped:
        print(f"    SKIPPED {_cell_label(r)}: {r['skip_reason']}")

    head = (f"  {'M_c':>4s} {'wmf':>5s} {'T_eq':>5s} {'T_int':>5s} | {'HT in-box':>10s} "
            f"{'HT tot rad':>10s} {'HT max rho':>10s} {'HT Tmax':>8s} | "
            f"{'C in-box':>8s} {'C clamp':>8s} {'HT clamp':>8s} | {'R_p HT':>8s} {'R_p C':>8s}")
    print("\n" + head)
    print("  " + "-" * (len(head) - 2))
    for r in rows:
        if r["skipped"]:
            print(f"  {r['mass_core_mearth']:>4.1f} {r['wmf']:>5.2f} {r['t_eq_k']:>5.0f} "
                  f"{r['t_int_k']:>5.0f} |  --- skipped ---")
            continue
        print(f"  {r['mass_core_mearth']:>4.1f} {r['wmf']:>5.2f} {r['t_eq_k']:>5.0f} "
              f"{r['t_int_k']:>5.0f} | "
              f"{_num(r, 'hitemp_in_box_radiative_mass_fraction'):>10.6f} "
              f"{_num(r, 'hitemp_total_radiative_mass_fraction'):>10.4f} "
              f"{r['hitemp_max_radiative_rho']:>10s} {r['hitemp_t_max_k']:>8s} | "
              f"{r['composite_n_in_box_layers']:>8d} "
              f"{_num(r, 'composite_clamp_fraction_all'):>8.4f} "
              f"{_num(r, 'hitemp_clamp_fraction_all'):>8.4f} | "
              f"{r['hitemp_radius_rearth']:>8s} {r['composite_radius_rearth']:>8s}")

    print_rankings(rows, floors)

    closed = [r for r in live
              if r["composite_n_in_box_layers"] == 0
              and _num(r, "composite_clamp_fraction_all", 1.0) < CLAMP_LIMIT
              and _num(r, "hitemp_in_box_radiative_mass_fraction", 0.0) > 0.0]
    print(f"\ncells where the composite closes every in-box radiative layer with "
          f"clamp_fraction_all < {CLAMP_LIMIT:.0%}: {len(closed)} of {len(live)} converged")
    for r in closed:
        print(f"  {_cell_label(r)}: HITEMP in-box "
              f"{_num(r, 'hitemp_in_box_radiative_mass_fraction'):.6f} of envelope mass, "
              f"in-box layers {r['hitemp_in_box_layers'] or 'none'}; "
              f"R_p {r['hitemp_radius_rearth']} -> {r['composite_radius_rearth']} R_earth")
    substantial = [r for r in closed
                   if _num(r, "hitemp_in_box_radiative_mass_fraction", 0.0) > 0.01]
    print(f"  of those, with a SUBSTANTIAL in-box layer (> 1% of envelope mass): "
          f"{len(substantial)}")

    entirely = [r for r in live if r["hitemp_layer_entirely_in_box"]]
    print(f"\ncells with a HITEMP radiative layer lying ENTIRELY inside the FP box "
          f"(no part below {t_floor:.0f} K or below ρ = {rho_floor:g}): "
          f"{len(entirely)} of {len(live)} converged")
    for r in entirely:
        print(f"  {_cell_label(r)}: {r['hitemp_radiative_layers']}")
    return paths[0]


if __name__ == "__main__":
    main()
