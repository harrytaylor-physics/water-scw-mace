#!/usr/bin/env python3
"""Extended Data figures ED1-ED8 (A67).

    python figures/ed_items.py config/fiducial.yaml            # all
    python figures/ed_items.py config/fiducial.yaml --only ed5

Plots and tabulates results this repository already holds; it runs no new
physics. Where a panel's data are not present, the panel is drawn with its axes
and labelled — so the panel lettering the text cites stays stable — and the gap
is recorded in ``figs/ed/ed_numbers.json`` under ``blocked``.

Output goes to ``figs/ed/``: vector PDF at 183 mm (or 89 mm where one panel
suffices), through the same ``save_pdf`` guard as the main figures, so the
vocabulary rule, the 7 pt minimum and the exact width are enforced here too.
"""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
from matplotlib.colors import LogNorm, TwoSlopeNorm
from matplotlib.lines import Line2D

from config import Config, load_config
from figures._manuscript import (
    ED_DIR,
    OKABE_ITO,
    SINGLE_WIDTH_MM,
    FigureNumbers,
    figure,
    mosaic,
    panel_label,
    save_pdf,
    set_manuscript_style,
)

ROOT = Path(__file__).resolve().parents[1]
ARC = ROOT / "ed_data" / "arc"
BAND = (100.0, 4000.0)
LEDGER = "ed_numbers.json"

# Keys are the data's own (opacity/bands.py); DISPLAY is what the manuscript calls
# them. "rotation" and "libration" are the same 300-1000 cm^-1 window — the two
# names are recorded in the ledger so a caption cannot be read as two bands.
BANDS_CM = {"libration": (300.0, 1000.0), "bend": (1400.0, 1800.0),
            "window": (1900.0, 3100.0), "stretch": (3200.0, 4000.0)}
BAND_DISPLAY = {"libration": "rotation", "bend": "bend",
                "window": "window", "stretch": "stretch"}


def _commit(path: str) -> str:
    try:
        out = subprocess.run(["git", "log", "-1", "--format=%h", "--", path],
                             cwd=ROOT, capture_output=True, text=True,
                             check=True).stdout.strip()
        return out or subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=ROOT,
                                     capture_output=True, text=True,
                                     check=True).stdout.strip()
    except Exception:  # noqa: BLE001
        return "unknown"


def _src(nums: FigureNumbers, panel: str, *files: str) -> None:
    nums.panel(panel, source_files=list(files),
               source_commits={f: _commit(f) for f in files})


def _blocked(ax, nums: FigureNumbers, panel: str, what: str, needs: str) -> None:
    """Draw an empty, labelled panel and record why it is empty."""
    ax.set_xticks([])
    ax.set_yticks([])
    for side in ("left", "bottom"):
        ax.spines[side].set_visible(False)
    ax.text(0.5, 0.56, what, ha="center", va="center", transform=ax.transAxes,
            fontsize=8.0, color="#444444")
    ax.text(0.5, 0.40, "data not in the pulled set", ha="center", va="center",
            transform=ax.transAxes, fontsize=7.0, color="#B00020")
    nums.panel(panel, blocked=True, needs=needs, shows=what)


def _json(cfg: Config, name: str):
    p = cfg.output_path / name
    return json.loads(p.read_text()) if p.exists() else None


# ===========================================================================
def ed1(cfg: Config) -> Path:
    """Validation of the interatomic potential."""
    nums = FigureNumbers(figure="ed_fig1", caption_hint=(
        "a, force-parity RMSE by regime against the benchmark tolerance. "
        "b, potential-vs-functional force error (data absent). c, finite-size "
        "check (data absent). d, NVE energy drift against the 2 meV/atom/ps gate."))
    fig, axes = figure(96.0, nrows=2, ncols=2)
    (ax_a, ax_b), (ax_c, ax_d) = axes

    ref_path = ARC / "water-scw-mace" / "benchmarks" / "reference" / "reference.json"
    ref = json.loads(ref_path.read_text())
    parity = ref["parity_meV_A"]
    tol = 0.15
    names = list(parity)
    vals = [float(parity[k]) for k in names]
    x = np.arange(len(names))
    ax_a.bar(x, vals, width=0.5, color=[OKABE_ITO["blue"], OKABE_ITO["vermillion"]],
             edgecolor="black", linewidth=0.4, zorder=3)
    for xi, v in zip(x, vals):
        ax_a.errorbar(xi, v, yerr=v * tol, fmt="none", ecolor="black",
                      elinewidth=0.8, capsize=3, zorder=4)
        ax_a.annotate(f"{v:.0f}", (xi, v * (1 + tol)), textcoords="offset points",
                      xytext=(0, 3), ha="center", fontsize=7.0)
    ax_a.set_xticks(x)
    ax_a.set_xticklabels(["bulk", "boundary"], fontsize=7.0)
    ax_a.set_ylabel("Force RMSE (meV Å$^{-1}$)")
    ax_a.set_ylim(0, max(vals) * 1.45)
    ax_a.annotate("training-set parity,\nnot held out", (0.03, 0.88),
                  xycoords="axes fraction", fontsize=7.0, color="#B00020", va="top")
    nums.axes("a", ax_a, "regime", "Force RMSE (meV/A)")
    nums.panel("a", rmse_meV_per_A=parity, tolerance_frac=tol,
               # reference.json's own `parity_methodology` wording is paraphrased
               # here; only the wording changes, and the file is in source_files.
               methodology=("training-set parity by regime, the same evaluation "
                            "eval_split.sh performs; optimistic vs true held-out. "
                            "Frames with T > 3000 K excluded (outside the applied "
                            "envelope). Paraphrased from reference.json's "
                            "parity_methodology."),
               caveat=("the pulled reference states these are TRAINING-SET parity by "
                       "regime, explicitly 'optimistic vs true held-out error', with "
                       "T > 3000 K frames excluded. The only held-out parity in the "
                       "pull is 12.3 / 20.2 meV/A and belongs to the LDA model, not "
                       "this one."))
    _src(nums, "a", "ed_data/arc/water-scw-mace/benchmarks/reference/reference.json")

    _blocked(ax_b, nums, "b",
             "potential vs its reference functional,\nand revPBE vs revPBE0",
             "per-configuration force differences for the two comparisons; no file "
             "in the pull holds them, and the quoted ~49 and ~376 meV/A appear "
             "nowhere in it")
    _blocked(ax_c, nums, "c", "pressure and g(r) vs molecule count,\n64 - 2048",
             "the finite-size sweep output; convergence/03_md_sampling/finite_size.py "
             "is present but its results are not. The README records only the verdict "
             "('flat over a 32x volume range, 64 -> 2048; 6400 OOM'd')")

    drift_max, gate = 0.017, 2.0
    ax_d.bar([0], [drift_max], width=0.45, color=OKABE_ITO["green"],
             edgecolor="black", linewidth=0.4, zorder=3, label="worst recorded")
    ax_d.axhline(gate, color=OKABE_ITO["vermillion"], lw=1.0, ls="--", zorder=4,
                 label=f"gate {gate:g}")
    ax_d.set_yscale("log")
    ax_d.set_xticks([0])
    ax_d.set_xticklabels(["21 NVE runs"], fontsize=7.0)
    ax_d.set_ylabel("|energy drift| (meV atom$^{-1}$ ps$^{-1}$)")
    ax_d.set_ylim(drift_max / 6, gate * 6)
    ax_d.annotate(f"{drift_max:g}", (0, drift_max), textcoords="offset points",
                  xytext=(0, 4), ha="center", fontsize=7.0)
    ax_d.legend(loc="upper right", fontsize=7.0)
    nums.axes("d", ax_d, "runs", "|drift| (meV/atom/ps)", ylog=True)
    nums.panel("d", worst_drift_meV_per_atom_per_ps=drift_max, gate=gate,
               n_runs=21, n_passed=21,
               note=("a bound, not a distribution: the pull records the worst value "
                     "over 21 LDA NVE runs, not a value per state point"))
    _src(nums, "d", "ed_data/arc/water-scw-mace/docs/LDA_CAMPAIGN.md")

    for ax, letter, dx in ((ax_a, "a", -0.16), (ax_b, "b", -0.16),
                          (ax_c, "c", -0.16),
                          # d's y-label is long and its axis narrow; at -0.16 the
                          # bold letter overprints the "ps^-1" of the label.
                          (ax_d, "d", -0.30)):
        panel_label(ax, letter, dx=dx, dy=1.09)
    return save_pdf(fig, "ed_fig1", nums, out_dir=ED_DIR, ledger_name=LEDGER)


# ===========================================================================
def ed2(cfg: Config) -> Path:
    """Planck-mean surface and the Planck ratio map."""
    from opacity.build_tables import compute_point_means, get_source
    from opacity.crossover import crossover_v3_rows
    from opacity.table import CM2_PER_G_TO_SI

    nums = FigureNumbers(figure="ed_fig2", caption_hint=(
        "a, Planck-mean opacity over the 52 first-principles points. b, the "
        "first-principles / HITEMP Planck ratio on the same axes. c, the "
        "first-principles / POKAZATEL ratio in both means at the 27 states "
        "POKAZATEL covers, which Fig. 2b no longer carries."))
    t_pts, rho_pts, _, kp = compute_point_means("fp", cfg, BAND)
    src = get_source("fp", cfg, band=BAND)
    rows = crossover_v3_rows(cfg, BAND)

    n = 320
    t_lin = np.linspace(t_pts.min(), t_pts.max(), n)
    lr_lin = np.linspace(np.log10(rho_pts.min()), np.log10(rho_pts.max()), n)
    TT, LR = np.meshgrid(t_lin, lr_lin)
    q = np.column_stack([TT.ravel(), LR.ravel()])

    # A91: three panels each carrying a colour bar will not fit across one
    # journal width -- b's colour-bar label ran into c's y-label. Two rows, with
    # the scatter spanning the lower one.
    fig, axd = mosaic(150.0, [["a", "b"], ["c", "c"]])
    ax_a, ax_b, ax_c = axd["a"], axd["b"], axd["c"]
    surf = (10.0 ** src.table._lin(q)[:, 1].reshape(TT.shape)) / CM2_PER_G_TO_SI
    ok = np.isfinite(surf)
    m = ax_a.pcolormesh(10.0**LR, TT, np.ma.masked_invalid(surf), cmap="viridis",
                        norm=LogNorm(np.nanmin(surf[ok]), np.nanmax(surf[ok])),
                        shading="auto", rasterized=True)
    cb = fig.colorbar(m, ax=ax_a, pad=0.012, fraction=0.05, aspect=28)
    cb.set_label("κ$_P$ (cm$^{2}$ g$^{-1}$)", fontsize=7.6)
    cb.ax.tick_params(labelsize=7.0)

    from scipy.interpolate import LinearNDInterpolator

    pts = np.column_stack([np.array([r["T"] for r in rows]),
                           np.log10([r["rho"] for r in rows])])
    ratio = np.array([r["ratio_planck_bare"] for r in rows])
    lin = LinearNDInterpolator(pts, np.log10(ratio), rescale=True)
    rsurf = 10.0 ** lin(q).reshape(TT.shape)
    ok2 = np.isfinite(rsurf)
    m2 = ax_b.pcolormesh(10.0**LR, TT, np.ma.masked_invalid(rsurf), cmap="viridis",
                         norm=LogNorm(np.nanmin(rsurf[ok2]), np.nanmax(rsurf[ok2])),
                         shading="auto", rasterized=True)
    cb2 = fig.colorbar(m2, ax=ax_b, pad=0.012, fraction=0.05, aspect=28)
    cb2.set_label("first principles / HITEMP, κ$_P$", fontsize=7.6)
    cb2.ax.tick_params(labelsize=7.0)

    # -- c: FP / POKAZATEL, both means (A91: moved off Fig. 2b) --------------
    t_pk, rho_pk, kr_pk, kp_pk = compute_point_means("pokazatel", cfg, BAND)
    by_state = {(round(t, 3), round(r, 3)): (kr, kp)
                for t, r, kr, kp in zip(t_pk, rho_pk, kr_pk, kp_pk)}
    pk = []
    for r in rows:
        key = (round(r["T"], 3), round(r["rho"], 3))
        if key in by_state:
            kr_p, kp_p = by_state[key]
            pk.append({"T": r["T"], "rho": r["rho"],
                       "rosseland": r["kappa_R_fp_cm2_g"] / kr_p,
                       "planck": r["kappa_P_fp_cm2_g"] / kp_p})
    if pk:
        cc = [d["T"] for d in pk]
        s_pk = ax_c.scatter([d["rho"] for d in pk], [d["rosseland"] for d in pk],
                            c=cc, cmap="viridis", marker="s", s=19.0,
                            edgecolor="black", linewidth=0.25, zorder=3)
        ax_c.scatter([d["rho"] for d in pk], [d["planck"] for d in pk],
                     facecolors="none",
                     edgecolors=np.array([s_pk.cmap(s_pk.norm(v)) for v in cc]),
                     marker="s", s=30.0, linewidth=0.9, zorder=2)
        ax_c.axhline(1.0, color="black", lw=0.6, ls=":", zorder=1)
        ax_c.set_yscale("log")
        ax_c.set_ylabel("First principles / POKAZATEL")
        ax_c.set_xlim(auto=True)
        cb3 = fig.colorbar(s_pk, ax=ax_c, pad=0.012, fraction=0.05, aspect=28)
        cb3.set_label("Temperature T (K)", fontsize=7.6)
        cb3.ax.tick_params(labelsize=7.0)
        ax_c.legend(handles=[
            Line2D([], [], ls="none", marker="s", mfc="#7f7f7f", mec="black",
                   ms=4.0, label="filled: Rosseland"),
            Line2D([], [], ls="none", marker="s", mfc="none", mec="#7f7f7f",
                   ms=4.6, mew=0.9, label="open: Planck")],
            loc="upper left", fontsize=7.0)
        nums.axes("c", ax_c, "Density (kg/m^3)", "first principles / POKAZATEL",
                  xlog=True, ylog=True)
        nums.panel("c", n_states=len(pk),
                   rosseland_range=[min(d["rosseland"] for d in pk),
                                    max(d["rosseland"] for d in pk)],
                   planck_range=[min(d["planck"] for d in pk),
                                 max(d["planck"] for d in pk)],
                   moved_from="Fig. 2b (A91)",
                   isotopologue_caveat=(
                       "HITEMP spectra carry isotopologues 1, 2 and 3; POKAZATEL "
                       "the main one only, so this ratio is not matched in "
                       "coverage to the HITEMP ratios in a and b"))
        _src(nums, "c", "opacity/crossover.py")
    else:
        _blocked(ax_c, nums, "c", "first principles / POKAZATEL",
                 "no POKAZATEL state matches a first-principles state")

    for ax in (ax_a, ax_b):
        ax.scatter(rho_pts, t_pts, s=8.0, facecolors="none", edgecolors="white",
                   linewidths=0.55, zorder=3)
        ax.scatter(rho_pts, t_pts, s=8.0, facecolors="none", edgecolors="black",
                   linewidths=0.3, zorder=4)
    for ax in (ax_a, ax_b, ax_c):
        ax.set_xscale("log")
        ax.set_xlabel("Density ρ (kg m$^{-3}$)")
    ax_a.set_ylabel("Temperature T (K)")
    ax_b.set_ylabel("Temperature T (K)")

    nums.axes("a", ax_a, "Density (kg/m^3)", "Temperature (K)", xlog=True)
    nums.panel("a", n_computed_points=int(t_pts.size), band_cm=list(BAND),
               kappa_P_range_cm2_g=[float(kp.min()), float(kp.max())],
               colourbar_range_cm2_g=[float(np.nanmin(surf[ok])),
                                      float(np.nanmax(surf[ok]))],
               colour_scale="viridis, log",
               blank_region="outside the convex hull, where the model clamps")
    nums.axes("b", ax_b, "Density (kg/m^3)", "Temperature (K)", xlog=True)
    nums.panel("b", planck_ratio_range=[float(ratio.min()), float(ratio.max())],
               colourbar_range=[float(np.nanmin(rsurf[ok2])),
                                float(np.nanmax(rsurf[ok2]))],
               colour_scale="viridis, log")
    for p in ("a", "b"):
        _src(nums, p, "outputs/crossover_v3.json", "opacity/crossover.py")
    panel_label(ax_a, "a", dx=-0.17, dy=1.06)
    panel_label(ax_b, "b", dx=-0.17, dy=1.06)
    panel_label(ax_c, "c", dx=-0.17, dy=1.06)
    return save_pdf(fig, "ed_fig2", nums, out_dir=ED_DIR, ledger_name=LEDGER)


# ===========================================================================
def ed3(cfg: Config) -> Path:
    """Line-list conventions: broadening, continuum share, band truncation."""
    from opacity.build_tables import compute_point_means

    nums = FigureNumbers(figure="ed_fig3", caption_hint=(
        "a, kappa_R under three broadening conventions at six states. b, MT_CKD "
        "self-continuum share of the line-list kappa_R over the 52 points. "
        "c, band-truncation ratio for HITEMP over its own grid. d, kappa_R "
        "against the line-wing cut-off at the six states of a."))
    # (a) needs the width: six two-line state labels collide in a third of a page.
    # b and c are square (T, rho) maps with colourbars; sharing a row with a
    # line plot leaves dead space on both sides of them, so d gets its own row.
    fig, axd = mosaic(155.0, [["a", "a"], ["b", "c"], ["d", "d"]])
    ax_a, ax_b, ax_c, ax_d = axd["a"], axd["b"], axd["c"], axd["d"]

    doc = _json(cfg, "broadening_sensitivity.json")
    rows = sorted(doc["results"], key=lambda r: (r["T"], r["rho"]))
    convs = [("ref1bar", "1 bar", OKABE_ITO["grey"]),
             ("ideal", "ideal-gas P", OKABE_ITO["blue"]),
             ("aqua", "AQUA P", OKABE_ITO["vermillion"])]
    x = np.arange(len(rows))
    w = 0.26
    for i, (key, lab, col) in enumerate(convs):
        vals = np.array([r.get(f"kappa_R_{key}_cm2_g", np.nan) for r in rows])
        ax_a.bar(x + (i - 1) * w, vals, w, color=col, edgecolor="black",
                 linewidth=0.35, label=lab, zorder=3)
    ax_a.set_yscale("log")
    ax_a.set_xticks(x)
    ax_a.set_xticklabels([f"{r['T']:.0f} K\n{r['rho']:g} kg m$^{{-3}}$" for r in rows],
                         fontsize=7.0)
    ax_a.set_xlabel("T (K)  /  ρ (kg m$^{-3}$)")
    ax_a.set_ylabel("κ$_R$ (cm$^{2}$ g$^{-1}$)")
    ax_a.legend(loc="upper left", fontsize=7.0)
    nums.axes("a", ax_a, "state", "kappa_R (cm^2/g)", ylog=True)
    nums.panel("a", states=[[r["T"], r["rho"]] for r in rows],
               kappa_R_cm2_g={k: [r.get(f"kappa_R_{k}_cm2_g") for r in rows]
                              for k, _l, _c in convs},
               truncation_cm=doc["truncation_cm"])
    _src(nums, "a", "outputs/broadening_sensitivity.json",
         "tools/broadening_sensitivity.py")

    v3 = _json(cfg, "crossover_v3.json")["results"]
    t = np.array([r["T"] for r in v3])
    rho = np.array([r["rho"] for r in v3])
    share = np.array([r["mtckd_share_kappa_R"] for r in v3])
    sc = ax_b.scatter(rho, t, c=share, cmap="viridis", s=26, vmin=0.0, vmax=1.0,
                      edgecolor="black", linewidth=0.3, zorder=3)
    cb = fig.colorbar(sc, ax=ax_b, pad=0.012, fraction=0.05, aspect=28)
    cb.set_label("MT$_-$CKD share of κ$_R$", fontsize=7.6)
    cb.ax.tick_params(labelsize=7.0)
    ax_b.set_xscale("log")
    ax_b.set_xlabel("Density ρ (kg m$^{-3}$)")
    ax_b.set_ylabel("Temperature T (K)")
    nums.axes("b", ax_b, "Density (kg/m^3)", "Temperature (K)", xlog=True)
    nums.panel("b", n_points=int(t.size),
               share_range=[float(share.min()), float(share.max())],
               share_at_650K_500=float(share[(t == 650) & (rho == 500)][0]))
    _src(nums, "b", "outputs/crossover_v3.json", "opacity/mtckd.py")

    # The widest band ALL the line-list files share, detected rather than assumed:
    # the rows added at 2250 and 2750 K start at 50 cm^-1, so it is [50, 11000], not
    # the [10, 11000] the other rows reach down to (A60).
    from opacity.build_tables import detect_common_band

    full = detect_common_band(cfg, ["hitemp"])
    t_h, rho_h, kr_c, _ = compute_point_means("hitemp", cfg, BAND)
    try:
        _, _, kr_f, _ = compute_point_means("hitemp", cfg, full)
        trunc = kr_f / kr_c
        sc2 = ax_c.scatter(rho_h, t_h, c=trunc, cmap="viridis", s=20,
                           norm=LogNorm(trunc.min(), trunc.max()),
                           edgecolor="black", linewidth=0.3, zorder=3)
        cb2 = fig.colorbar(sc2, ax=ax_c, pad=0.012, fraction=0.05, aspect=28)
        cb2.set_label("κ$_R$(full) / κ$_R$(common)", fontsize=7.6)
        cb2.ax.tick_params(labelsize=7.0)
        ax_c.set_xscale("log")
        ax_c.set_xlabel("Density ρ (kg m$^{-3}$)")
        ax_c.set_ylabel("Temperature T (K)")
        nums.axes("c", ax_c, "Density (kg/m^3)", "Temperature (K)", xlog=True)
        nums.panel("c", n_points=int(t_h.size), full_band_cm=list(full),
                   common_band_cm=list(BAND),
                   ratio_range=[float(trunc.min()), float(trunc.max())],
                   note=("full band is DETECTED as the widest all the line-list "
                         "files share. The rows added at 2250 and 2750 K start at "
                         "50 cm^-1, so it is [50, 11000]; the other rows reach "
                         "down to 10. With a 50 cm^-1 line truncation "
                         "nothing below 50 can reach the 100 cm^-1 band edge, so the "
                         "common-band means are unaffected either way"))
        _src(nums, "c", "opacity/build_tables.py")
    except Exception as exc:  # noqa: BLE001
        _blocked(ax_c, nums, "c", "band truncation for HITEMP",
                 f"means on the full band could not be built: {exc}")

    # -- d: the line-wing cut-off (A71) ---------------------------------
    wc = _json(cfg, "wing_cutoff.json")
    if wc and wc.get("results"):
        from matplotlib import colormaps

        cmap = colormaps["viridis"]
        res = sorted(wc["results"], key=lambda r: (r["T"], r["rho"]))
        truncs = sorted(float(t) for t in wc["truncations_cm"])
        rec = {}
        for i, r in enumerate(res):
            col = cmap(i / max(1, len(res) - 1))
            y = [r["by_truncation"][f"{t:g}"]["kappa_R_cm2_g"] for t in truncs]
            ax_d.plot(truncs, y, color=col, lw=1.0, marker="o", ms=3.0, mew=0.3,
                      mec="black", zorder=3,
                      label=f"{r['T']:.0f} K, {r['rho']:g} kg m$^{{-3}}$")
            rec[f"{r['T']:g}_{r['rho']:g}"] = {
                "truncation_cm": truncs, "kappa_R_cm2_g": y,
                "kappa_R_max_over_min": r["kappa_R_max_over_min"],
                "kappa_P_cm2_g": [r["by_truncation"][f"{t:g}"]["kappa_P_cm2_g"]
                                  for t in truncs],
                "reproduces_delivered_to": r.get("reproduces_delivered_to")}
        ax_d.axvline(wc["delivered_truncation_cm"], color="black", lw=0.8, ls="--",
                     zorder=2)
        ax_d.annotate("delivered", (wc["delivered_truncation_cm"], 0.97),
                      xycoords=("data", "axes fraction"), fontsize=7.0, rotation=90,
                      ha="right", va="top")
        ax_d.set_xscale("log")
        ax_d.set_yscale("log")
        ax_d.set_xlabel("Line-wing cut-off (cm$^{-1}$)")
        ax_d.set_ylabel("κ$_R$ (cm$^2$ g$^{-1}$)")
        # lower right: every curve rises to the right, so that corner is clear
        ax_d.legend(loc="lower right", fontsize=7.0, ncol=2, labelspacing=0.25,
                    columnspacing=1.0, borderpad=0.3, handlelength=1.4)
        ax_d.set_xticks(truncs)
        ax_d.set_xticklabels([f"{t:g}" for t in truncs], fontsize=7.0)
        ax_d.minorticks_off()
        nums.axes("d", ax_d, "line-wing cut-off (cm^-1)", "kappa_R (cm^2/g)",
                  xlog=True, ylog=True)
        nums.panel("d", truncations_cm=truncs,
                   delivered_truncation_cm=wc["delivered_truncation_cm"],
                   band_cm=wc["band_cm"], by_point=rec,
                   largest_spread=max(r["kappa_R_max_over_min"] for r in res),
                   renormalisation=wc["renormalisation"],
                   note=("only the truncation varies; the 50 cm^-1 column "
                         "reproduces the delivered kappa_R exactly, which is what "
                         "makes the rest a measurement"))
        _src(nums, "d", "outputs/wing_cutoff.json", "tools/wing_cutoff.py")
    else:
        _blocked(ax_d, nums, "d", "kappa_R vs the line-wing cut-off,\nsix states",
                 "outputs/wing_cutoff.json (python tools/wing_cutoff.py "
                 "config/fiducial.yaml --wavenum-max 4600)")

    panel_label(ax_a, "a", dx=-0.075, dy=1.09)
    panel_label(ax_b, "b", dx=-0.17, dy=1.09)
    panel_label(ax_c, "c", dx=-0.17, dy=1.09)
    panel_label(ax_d, "d", dx=-0.075, dy=1.09)
    return save_pdf(fig, "ed_fig3", nums, out_dir=ED_DIR, ledger_name=LEDGER)


# ===========================================================================
def _scaling_panels(cfg: Config, fig, axd, nums: FigureNumbers,
                    P: dict[str, str]) -> None:
    """Draw the density-scaling panels into ``axd`` under the letters in ``P``.

    A88 folded ED Fig. 9 into ED Fig. 4: the dilute band ratios and the rho^2
    scaling are the same question asked at two densities, so they are one item.
    ``P`` maps this block's own panel letters onto the host figure's, which keeps
    the drawing code readable while the ledger records the letters a reader sees.
    """
    from matplotlib import colormaps

    from opacity.build_tables import compute_point_means

    doc = _json(cfg, "rho_scaling.json")
    if doc is None:
        raise SystemExit("outputs/rho_scaling.json not found — run "
                         "tools/rho_scaling.py first")

    temps = sorted({float(t) for t in doc["sources"]["fp"]})
    cmap = colormaps["viridis"]
    colour = {t: cmap(i / max(1, len(temps) - 1)) for i, t in enumerate(temps)}

    # -- a, b, c: alpha/rho^2 -------------------------------------------
    recorded = {}
    for panel, (key, title) in zip("abc", SCALING_PANELS):
        ax = axd[P[panel]]
        per_panel = {}
        for source, per_src in doc["sources"].items():
            ls, lw, mk = SOURCE_STYLE[source]
            for t_s, iso in per_src.items():
                t_k = float(t_s)
                if source != "fp" and t_k not in colour:
                    continue  # only the isotherms the computed grid also has
                rho = [r["rho_kg_m3"] for r in iso["rows"]]
                q = [r["q"][key] for r in iso["rows"]]
                ax.plot(rho, q, ls=ls, lw=lw, color=colour[t_k],
                        marker=mk if source == "fp" else None, ms=2.6, mew=0.3,
                        mec="black", alpha=1.0 if source == "fp" else 0.65, zorder=3)
                per_panel.setdefault(source, {})[t_s] = {"rho_kg_m3": rho, "q": q}
        # The rho^2 reference: alpha/rho^2 is flat for exactly pair-wise
        # absorption, so a horizontal segment through each isotherm's densest
        # point is what the data would lie on if it were. Drawn per isotherm
        # rather than as one line because the level is not universal.
        refs = {}
        for t_s, iso in doc["sources"]["fp"].items():
            rows = iso["rows"]
            if not rows:
                continue
            anchor = rows[-1]
            ax.plot([rows[0]["rho_kg_m3"], anchor["rho_kg_m3"]],
                    [anchor["q"][key]] * 2, ls=(0, (1, 1)), lw=0.8, alpha=0.6,
                    color=colour[float(t_s)], zorder=2)
            refs[t_s] = {"rho_span_kg_m3": [rows[0]["rho_kg_m3"],
                                            anchor["rho_kg_m3"]],
                         "level": anchor["q"][key]}
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("Density ρ (kg m$^{-3}$)")
        ax.set_ylabel("α / ρ$^2$   (cm$^5$ g$^{-2}$)")
        ax.annotate(title, (0.03, 0.04), xycoords="axes fraction", fontsize=7.0,
                    va="bottom")
        if panel == "a":
            ax.legend(handles=[Line2D([], [], ls=SOURCE_STYLE[k][0],
                                      lw=SOURCE_STYLE[k][1], color="#666666",
                                      marker=SOURCE_STYLE[k][2] if k == "fp" else None,
                                      ms=2.8, mfc="#bbbbbb", mec="black", mew=0.3,
                                      label=SOURCE_LABEL[k])
                               for k in SOURCE_STYLE]
                              + [Line2D([], [], ls=(0, (1, 1)), lw=0.8, alpha=0.6,
                                        color="#666666", label="ρ$^2$ reference")],
                      loc="lower left", fontsize=7.0, handlelength=2.0,
                      labelspacing=0.3, borderpad=0.35,
                      bbox_to_anchor=(0.0, 0.10))
        nums.axes(P[panel], ax, "Density (kg/m^3)", "alpha/rho^2 (cm^5 g^-2)",
                  xlog=True, ylog=True)
        nums.panel(P[panel], diagnostic=key, curves=per_panel,
                   flat_means="exactly pair-wise, alpha proportional to rho^2",
                   rho2_reference_per_isotherm=refs,
                   rho2_reference_note=("dotted: alpha/rho^2 held at each "
                                        "first-principles isotherm's densest "
                                        "point, which is where the data would lie "
                                        "if absorption were exactly pair-wise"))
        _src(nums, P[panel], "outputs/rho_scaling.json", "tools/rho_scaling.py")
        recorded[key] = per_panel

    # Colour is temperature in a-d; a colourbar says so once instead of a
    # seven-entry legend in every panel.
    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import BoundaryNorm

    edges = [temps[0] - 100] + [0.5 * (a + b) for a, b in zip(temps, temps[1:])] \
        + [temps[-1] + 100]
    sm = ScalarMappable(norm=BoundaryNorm(edges, cmap.N), cmap=cmap)
    cb = fig.colorbar(sm, ax=axd[P["c"]], pad=0.015, fraction=0.055, aspect=22,
                      ticks=temps)
    cb.set_label("Temperature T (K)", fontsize=7.6)
    cb.ax.tick_params(labelsize=7.0)

    # -- d: the local index ---------------------------------------------
    ax_d = axd[P["d"]]
    sig_ln = doc["scatter"]["per_diagnostic"]["rosseland"]["sigma_ln_alpha"]
    for t_s, iso in doc["sources"]["fp"].items():
        x = [i["rho_geom_kg_m3"] for i in iso["index"]]
        y = [i["index"]["rosseland"] for i in iso["index"]]
        e = [float(np.sqrt(2) * sig_ln / i["dln_rho"]) for i in iso["index"]]
        ax_d.errorbar(x, y, yerr=e, ls="-", lw=1.0, color=colour[float(t_s)],
                      marker="o", ms=2.8, mew=0.3, mec="black", elinewidth=0.6,
                      capsize=1.6, zorder=3)
    ax_d.axhline(2.0, color="black", lw=0.9, ls="--", zorder=4)
    ax_d.axhline(1.0, color="#777777", lw=0.7, ls=":", zorder=4)
    ax_d.annotate("ρ$^2$  pair-wise", (0.02, 0.97), xycoords="axes fraction",
                  fontsize=7.0, va="top", ha="left")
    ax_d.annotate("ρ$^1$  dilute lines", (0.02, 0.03), xycoords="axes fraction",
                  fontsize=7.0, va="bottom", ha="left")
    ax_d.set_xscale("log")
    ax_d.set_xlabel("Density ρ (kg m$^{-3}$)")
    ax_d.set_ylabel("d ln α / d ln ρ  (κ$_R$)")
    ax_d.set_ylim(-3.5, 7.0)
    nums.axes(P["d"], ax_d, "Density (kg/m^3)", "d ln alpha / d ln rho", xlog=True)
    nums.panel(P["d"], diagnostic="rosseland",
               sigma_ln_alpha_used=sig_ln,
               error_bar=("+-sqrt(2) sigma_ln_alpha / dln_rho, with "
                          "sigma_ln_alpha from the 729.344/810.382 pair at 3000 K "
                          "as the brief specifies"),
               index_by_isotherm={t: {"rho_geom_kg_m3": [i["rho_geom_kg_m3"]
                                                         for i in iso["index"]],
                                      "index": [i["index"]["rosseland"]
                                                for i in iso["index"]]}
                                  for t, iso in doc["sources"]["fp"].items()})
    _src(nums, P["d"], "outputs/rho_scaling.json")

    # -- e, f: the fitted slope, window and Rosseland ---------------------
    fit = doc["fitted_slope_dense"]["per_source"]
    pooled = doc["scatter_pooled"]["per_diagnostic"]
    pair = doc["scatter"]["per_diagnostic"]
    scale = {k: pooled[k]["sigma_ln_alpha"] / pair[k]["sigma_ln_alpha"] for k in pooled}
    src_col = {"fp": OKABE_ITO["blue"], "hitemp": OKABE_ITO["vermillion"],
               "hitemp_mtckd": OKABE_ITO["green"]}
    for panel, key, title in (("e", "window_2550", "2550 cm$^{-1}$ window"),
                              ("f", "rosseland", "κ$_R$, 100-4000 cm$^{-1}$")):
        ax = axd[P[panel]]
        rec = {}
        for source, per_iso in fit.items():
            ts = sorted(float(t) for t in per_iso)
            y = [per_iso[f"{t:g}"][key]["slope"] for t in ts]
            e = [per_iso[f"{t:g}"][key]["standard_error"] * scale[key] for t in ts]
            ax.errorbar(ts, y, yerr=e, ls=SOURCE_STYLE[source][0], lw=1.0,
                        color=src_col[source], marker=SOURCE_STYLE[source][2], ms=3.0,
                        mew=0.3, mec="black", elinewidth=0.7, capsize=1.8, zorder=3,
                        label=SOURCE_LABEL[source])
            rec[source] = {"T_K": ts, "slope": y, "standard_error_pooled": e}
        ax.axhline(2.0, color="black", lw=0.9, ls="--", zorder=4)
        ax.axhline(1.0, color="#777777", lw=0.7, ls=":", zorder=4)
        ax.set_xlabel("Temperature T (K)")
        ax.set_ylabel("fitted d ln α / d ln ρ")
        ax.annotate(f"{title}\nfitted over ρ ≥ 100", (0.97, 0.97),
                    xycoords="axes fraction", fontsize=7.0, ha="right", va="top")
        ax.set_ylim(0.3, 3.5)
        if panel == "e":
            ax.legend(loc="upper left", fontsize=7.0, labelspacing=0.3,
                      borderpad=0.35)
        nums.axes(P[panel], ax, "Temperature (K)", "fitted d ln alpha / d ln rho")
        nums.panel(P[panel], diagnostic=key, rho_min_kg_m3=100.0, by_source=rec,
                   error_bars=("pooled sigma_ln_alpha (23 dof) rather than the "
                               "single-pair estimate, which has one"),
                   sigma_ln_alpha_pair=pair[key]["sigma_ln_alpha"],
                   sigma_ln_alpha_pooled=pooled[key]["sigma_ln_alpha"])
        _src(nums, P[panel], "outputs/rho_scaling.json")

    # -- g: the breakdown map -------------------------------------------
    ax_f = axd[P["g"]]
    bd = doc["breakdown_map"]
    for lam, style in (("1", dict(color="black", lw=1.1, ls="-")),
                       ("0.1", dict(color="black", lw=0.9, ls="--"))):
        rho_l = bd["rho_at_lambda_kg_m3"][lam]
        ax_f.axvline(rho_l, zorder=2, **style)
        ax_f.annotate(f"Λ = {lam}", (rho_l, 0.985),
                      xycoords=("data", "axes fraction"), fontsize=7.0, rotation=90,
                      ha="right", va="top")
    t_fp, rho_fp, _, _ = compute_point_means("fp", cfg, BAND)
    ax_f.scatter(rho_fp, t_fp, s=9.0, facecolor="none", edgecolor="#555555",
                 linewidth=0.45, zorder=3)
    bnds = boundary_table(cfg)
    for cat, (col, filled) in CATEGORY_STYLE.items():
        sel = [b for b in bnds if b["category"] == cat]
        if sel:
            ax_f.scatter([b["rho"] for b in sel], [b["T"] for b in sel], s=7.0,
                         marker="o", facecolor=col if filled else "none",
                         edgecolor=col, linewidth=0.5, zorder=4)
    lim_rec = {}
    for name, lw, col in (("Freedman et al. 2014", 1.0, OKABE_ITO["vermillion"]),
                          ("Kempton et al. 2023", 1.0, OKABE_ITO["orange"])):
        ent = bd["table_limits"].get(name)
        if not ent:
            continue
        ts = sorted(float(t) for t in ent["by_isotherm"])
        rs = [ent["by_isotherm"][f"{t:g}"]["rho_kg_m3"] for t in ts]
        ax_f.plot(rs, ts, color=col, lw=lw, marker="D", ms=3.2, mew=0.3,
                  mec="black", zorder=5,
                  label=f"{name.split(' et')[0]} {ent['pressure_bar']:g} bar")
        # How much of the problem each published table is documented for. Its
        # limit is known at two temperatures only, so it is interpolated linearly
        # in T between them and held flat outside — stated, not hidden.
        def _limit(t_query, ts=ts, rs=rs):
            return float(np.interp(t_query, ts, rs))

        lim_rec[name] = {
            "T_K": ts, "rho_kg_m3": rs, "pressure_bar": ent["pressure_bar"],
            "interpolation": ("linear in T between the two isotherms, flat "
                              "outside them"),
            "n_computed_points_denser": int(sum(1 for tt, rr in zip(t_fp, rho_fp)
                                                if rr > _limit(tt))),
            "n_boundaries_denser": int(sum(1 for b in bnds
                                           if b["rho"] > _limit(b["T"]))),
        }
    ax_f.set_xscale("log")
    ax_f.set_xlabel("Density ρ (kg m$^{-3}$)")
    ax_f.set_ylabel("Temperature T (K)")
    # The boundaries keep their four category colours, but naming all four here
    # would need seven legend rows in a third-width panel; ED Fig. 6a is where the
    # categories are read, and the counts are in this panel's ledger entry.
    handles = [
        Line2D([], [], ls="none", marker="o", mfc="none", mec="#555555", ms=3.4,
               label=f"computed points ({int(np.asarray(t_fp).size)})"),
        Line2D([], [], ls="none", marker="o", mfc=CATEGORY_STYLE["on computed data"][0],
               mec=CATEGORY_STYLE["on computed data"][0], ms=3.4,
               label="boundaries (108), coloured by category as ED Fig. 6a"),
    ] + ax_f.get_legend_handles_labels()[0]
    # upper left: nothing hot sits below 10 kg m^-3
    leg = ax_f.legend(handles=handles, loc="upper left", fontsize=7.0, ncol=2,
                      labelspacing=0.3, handletextpad=0.4, borderpad=0.35,
                      columnspacing=1.0, frameon=True, framealpha=0.93,
                      edgecolor="#bbbbbb")
    leg.get_frame().set_linewidth(0.4)
    ax_f.set_ylim(200.0, 3350.0)
    nums.axes(P["g"], ax_f, "Density (kg/m^3)", "Temperature (K)", xlog=True)
    nums.panel(P["g"], lambda_definition=bd["definition"],
               rho_at_lambda_kg_m3=bd["rho_at_lambda_kg_m3"],
               lambda_depends_on=bd["depends_on"],
               n_computed_points=int(np.asarray(t_fp).size),
               n_boundaries=len(bnds),
               counts_by_category={c: sum(1 for b in bnds if b["category"] == c)
                                   for c in CATEGORY_STYLE},
               table_limits=lim_rec)
    _src(nums, P["g"], "outputs/rho_scaling.json", "eos/aqua.py")


# ===========================================================================
def ed4(cfg: Config) -> Path:
    """Dilute-limit band structure."""
    from tools.text_checks import _load

    nums = FigureNumbers(figure="ed_fig4", caption_hint=(
        "a, first-principles / HITEMP band-integrated absorption ratio at nine "
        "dilute states. b, kappa(nu) for both sources at two of them. c-e, "
        "alpha/rho^2 against density along each isotherm with at least three "
        "densities, at three places in the spectrum; flat is exact pair-wise "
        "absorption. f, the local power-law index. g, h, the index fitted over "
        "rho >= 100 against temperature, by source, at the 2550 cm^-1 window and "
        "in kappa_R. i, the mean-free-path ratio Lambda with the computed points, "
        "the 108 boundaries and the published tables' documented pressure "
        "limits."))
    doc = _json(cfg, "text_checks.json")
    rows = sorted(doc["band_check"], key=lambda r: (r["T"], r["rho"]))
    band_names = list(BANDS_CM)

    # A88: the dilute band ratios (a, b) and the density scaling (c-i) are one
    # item, so they share a figure. Six columns is the common multiple of the top
    # row's two panels and the scaling block's three.
    fig, axd = mosaic(236.0, [["a", "a", "a", "b", "b", "b"],
                              ["c", "c", "d", "d", "e", "e"],
                              ["f", "f", "g", "g", "h", "h"],
                              ["i", "i", "i", "i", "i", "i"]])
    ax_a, ax_b = axd["a"], axd["b"]
    temps = sorted({r["T"] for r in rows})
    rhos = sorted({r["rho"] for r in rows})
    from matplotlib import colormaps

    cmap = colormaps["viridis"]
    colours = {t: cmap(i / max(1, len(temps) - 1)) for i, t in enumerate(temps)}
    markers = {rhos[0]: "o", rhos[1]: "s", rhos[2]: "^"}
    for bi, b in enumerate(band_names):
        for r in rows:
            jitter = (rhos.index(r["rho"]) - 1) * 0.16 + (temps.index(r["T"]) - 1) * 0.05
            ax_a.plot(bi + jitter, r["ratio_bare"][b], ls="none",
                      marker=markers[r["rho"]], ms=4.2, mfc=colours[r["T"]],
                      mec="black", mew=0.35, zorder=3)
    ax_a.axhline(1.0, color="black", lw=0.7, ls=":", zorder=1)
    ax_a.set_yscale("log")
    ax_a.set_xticks(range(len(band_names)))
    ax_a.set_xticklabels([f"{BAND_DISPLAY[b]}\n{BANDS_CM[b][0]:.0f}-{BANDS_CM[b][1]:.0f}"
                          for b in band_names], fontsize=7.0)
    ax_a.set_xlabel("Band (cm$^{-1}$)")
    ax_a.set_ylabel("First principles / HITEMP, band absorption")
    handles = [Line2D([], [], ls="none", marker="o", mfc=colours[t], mec="black",
                      ms=4.2, label=f"{t:.0f} K") for t in temps]
    handles += [Line2D([], [], ls="none", marker=markers[r], mfc="#bbbbbb",
                       mec="black", ms=4.2, label=f"{r:g} kg m$^{{-3}}$") for r in rhos]
    ax_a.legend(handles=handles, loc="upper left", fontsize=7.0, ncol=2)
    nums.axes("a", ax_a, "band", "band-absorption ratio", ylog=True)
    nums.panel("a", bands_cm={b: list(v) for b, v in BANDS_CM.items()},
               band_display_names=BAND_DISPLAY,
               band_naming_note=("'rotation' in the manuscript and 'libration' in "
                                 "the code are the same 300-1000 cm^-1 window"),
               states=[[r["T"], r["rho"]] for r in rows],
               ratio_bare={b: [r["ratio_bare"][b] for r in rows] for b in band_names},
               per_band_range={b: [min(r["ratio_bare"][b] for r in rows),
                                   max(r["ratio_bare"][b] for r in rows)]
                               for b in band_names})
    _src(nums, "a", "outputs/text_checks.json", "tools/text_checks.py",
         "opacity/bands.py")

    shown = []
    for (t_k, rho), col in (((1000.0, 10.0), OKABE_ITO["blue"]),
                            ((1500.0, 30.0), OKABE_ITO["vermillion"])):
        wn_fp, k_fp = _load(cfg, "fp", t_k, rho)
        wn_ll, k_ll = _load(cfg, "hitemp", t_k, rho)
        mf = (wn_fp >= BAND[0]) & (wn_fp <= BAND[1])
        ml = (wn_ll >= BAND[0]) & (wn_ll <= BAND[1])
        ax_b.plot(wn_ll[ml], k_ll[ml], color=col, lw=0.35, alpha=0.55, zorder=2)
        ax_b.plot(wn_fp[mf], k_fp[mf], color=col, lw=1.0, zorder=3)
        shown.append({"T_K": t_k, "rho_kg_m3": rho,
                      "fp_range_cm2_g": [float(k_fp[mf].min()), float(k_fp[mf].max())],
                      "line_list_range_cm2_g": [float(k_ll[ml].min()),
                                                float(k_ll[ml].max())]})
    ax_b.set_yscale("log")
    ax_b.set_xlim(*BAND)
    ax_b.set_xlabel("Wavenumber ν (cm$^{-1}$)")
    ax_b.set_ylabel("Opacity κ (cm$^{2}$ g$^{-1}$)")
    ytop = ax_b.get_ylim()[1]
    for b, (lo, hi) in BANDS_CM.items():
        ax_b.axvspan(lo, hi, color="#000000", alpha=0.055, lw=0, zorder=0)
        ax_b.text(0.5 * (lo + hi), ytop / 2.4, BAND_DISPLAY[b], ha="center",
                  va="bottom", fontsize=7.0, color="#444444", zorder=1)
    ax_b.legend(handles=[
        Line2D([], [], color=OKABE_ITO["blue"], lw=1.0, label="1000 K, 10 kg m$^{-3}$"),
        Line2D([], [], color=OKABE_ITO["vermillion"], lw=1.0,
               label="1500 K, 30 kg m$^{-3}$"),
        Line2D([], [], color="#666666", lw=1.0, label="thick: first principles"),
        Line2D([], [], color="#666666", lw=0.4, label="thin: HITEMP"),
    ], loc="lower left", fontsize=7.0, ncol=2)
    nums.axes("b", ax_b, "Wavenumber (cm^-1)", "Opacity (cm^2/g)", ylog=True)
    nums.panel("b", shown=shown, bands_shaded_cm={b: list(v)
                                                  for b, v in BANDS_CM.items()})
    _src(nums, "b", "opacity/loaders.py")
    _scaling_panels(cfg, fig, axd, nums,
                    {"a": "c", "b": "d", "c": "e", "d": "f", "e": "g", "f": "h",
                     "g": "i"})
    # a's y-label is the longest in the figure; -0.26 puts the bold letter on it.
    for letter, dx in (("a", -0.26), ("b", -0.16),
                       ("c", -0.30), ("d", -0.26), ("e", -0.26),
                       ("f", -0.26), ("g", -0.30), ("h", -0.30), ("i", -0.075)):
        panel_label(axd[letter], letter, dx=dx, dy=1.06)
    return save_pdf(fig, "ed_fig4", nums, out_dir=ED_DIR, ledger_name=LEDGER)




# ===========================================================================
MONOMER_STEP = {"sym": (0.01, "Å"), "asym": (0.01, "Å"), "bend": (1.0, "°")}
MONOMER_COMPONENT = {"sym": 2, "asym": 0, "bend": 2}  # the component that moves
MU_EXPERIMENT_D = 1.855
REF_BASIS = "aug-cc-pvtz"


def _monomer_pairs(mode: str, ref: str):
    """(minus, plus) geometry names for one mode about one reference."""
    order = ["sym", "asym", "bend"]
    base = 1 if ref == "equilibrium" else 7
    i = base + 2 * order.index(mode)
    return f"g{i:02d}", f"g{i + 1:02d}"          # g01 = +, g02 = - ...


# A88: the bands the window-validation reduction uses, as the producing script
# defined them. The "window" band, 1900-3100 cm^-1, is the transparency window the
# conservative bound leans on, and its ratio is the number A88 Task 2 wants.
WINDOW_VALIDATION_BANDS_CM = {
    "libration": (300.0, 1000.0),
    "bend": (1400.0, 1800.0),
    "window": (1900.0, 3100.0),
    "stretch": (3200.0, 4000.0),
    "overtone": (4200.0, 8000.0),
}


def ed5(cfg: Config, residual_npz: Path | None = None) -> Path:
    """Dipole-model validation: derivatives, ratios, and the dipole-model noise."""
    mono = ROOT / "monomer_test" / "pyscf_monomer.json"
    if not mono.exists():
        raise SystemExit(f"{mono} not found")
    doc = json.loads(mono.read_text())
    model_mu, dft_mu = doc["model_mu"], doc["dft_mu"]
    bases = list(dft_mu)

    nums = FigureNumbers(figure="ed_fig5", caption_hint=(
        "a,b dipole component against displacement for the three modes about "
        "r_OH = 0.957 and 1.00 A, model and revPBE0 in three bases. c, model/DFT "
        "derivative ratio. d, dipole-model residual and signal autocorrelation. "
        "e, band-resolved predicted/labelled absorption (pending)."))
    fig, axd = mosaic(150.0, [["a", "b", "c"], ["d", "d", "e"]])

    src_style = {"model": ("-", 2.2, OKABE_ITO["blue"]),
                 "cc-pvdz": (":", 1.4, OKABE_ITO["grey"]),
                 "aug-cc-pvdz": ("--", 1.4, OKABE_ITO["orange"]),
                 REF_BASIS: ("-", 1.4, OKABE_ITO["vermillion"])}
    mode_marker = {"sym": "o", "asym": "s", "bend": "^"}

    span = []
    for panel, ref, label in (("a", "equilibrium", "0.957"), ("b", "stretched", "1.00")):
        ax = axd[panel]
        recorded = {}
        for mode in ("sym", "asym", "bend"):
            gp, gm = _monomer_pairs(mode, ref)
            comp = MONOMER_COMPONENT[mode]
            for sname in ("model", *bases):
                table = model_mu if sname == "model" else dft_mu[sname]
                if gp not in table or gm not in table:
                    continue
                vals = [table[gm][comp], table[gp][comp]]
                if ref == "equilibrium":
                    vals = [table[gm][comp], table["g00"][comp], table[gp][comp]]
                xs = [-1, 0, 1] if len(vals) == 3 else [-1, 1]
                centre = 0.5 * (vals[0] + vals[-1])
                ls, lw, col = src_style[sname]
                ax.plot(xs, [v - centre for v in vals], ls=ls, lw=lw, color=col,
                        marker=mode_marker[mode], ms=3.2, mew=0.3, mec="black",
                        zorder=3)
                recorded.setdefault(mode, {})[sname] = [float(v - centre) for v in vals]
        ax.axhline(0.0, color="black", lw=0.6, ls=":", zorder=1)
        ax.set_xticks([-1, 0, 1])
        ax.set_xlabel("Displacement (steps)")
        ax.set_ylabel("μ component, centred (D)")
        ax.set_title("")  # the reference length goes in the legend, not a title
        if panel == "a":
            # Without this the reader cannot tell the four sources or the three
            # modes apart: colour/style is the source, marker is the mode.
            h = [Line2D([], [], ls=src_style[k][0], lw=src_style[k][1],
                        color=src_style[k][2],
                        label="model" if k == "model" else k)
                 for k in ("model", *bases)]
            h += [Line2D([], [], ls="none", marker=mode_marker[m], ms=3.4,
                         mfc="#bbbbbb", mec="black", mew=0.3, label=m)
                  for m in ("sym", "asym", "bend")]
            ax.legend(handles=h, loc="upper left", fontsize=7.0, ncol=2,
                      handlelength=1.6, columnspacing=0.9, labelspacing=0.3,
                      borderpad=0.35)
        ax.annotate(f"r$_{{OH}}$ = {label} Å", (0.97, 0.03),
                    xycoords="axes fraction", fontsize=7.0, ha="right", va="bottom")
        span.append(max(abs(v) for vs in recorded.values() for w in vs.values()
                        for v in w))
        nums.axes(panel, ax, "displacement (steps)", "centred dipole component (D)")
        nums.panel(panel, reference_r_OH_A=float(label), centred_values_D=recorded,
                   step={m: list(v) for m, v in MONOMER_STEP.items()},
                   component={"sym": "z", "asym": "x", "bend": "z"})
        _src(nums, panel, "monomer_test/pyscf_monomer.json")

    # One scale for both panels, with headroom for a's legend: the pair exists to
    # show that the slopes grow with r_OH, which separate autoscales would hide,
    # and a legend drawn over the data is worse than a wider axis.
    m = max(span)
    for panel in ("a", "b"):
        axd[panel].set_ylim(-1.18 * m, 2.05 * m)

    ax_c = axd["c"]
    derivs = doc["derivatives"]
    labels, ratios, units, refs = [], [], [], []
    for d in derivs:
        ref = "0.957" if abs(d["r_ref"] - 0.9572) < 1e-3 else "1.00"
        r = d["model"] / d["dft"][REF_BASIS]
        labels.append(d["mode"])
        refs.append(ref)
        ratios.append(r)
        units.append(d["unit"])
    x = np.arange(len(labels))
    cols = [OKABE_ITO["blue"] if "sym" == derivs[i]["mode"] else
            OKABE_ITO["vermillion"] if "asym" == derivs[i]["mode"] else
            OKABE_ITO["green"] for i in range(len(derivs))]
    ax_c.bar(x, np.abs(np.array(ratios)) - 1.0, bottom=1.0, width=0.6, color=cols,
             edgecolor="black", linewidth=0.4, zorder=3)
    ax_c.axhline(1.0, color="black", lw=0.8, ls="--", zorder=4)
    for xi, r in zip(x, ratios):
        ax_c.annotate(f"{r:.2f}", (xi, r), textcoords="offset points",
                      xytext=(0, 3 if r >= 1 else -10), ha="center", fontsize=7.0)
    ax_c.set_xticks(x)
    ax_c.set_xticklabels(labels, fontsize=7.0, rotation=35, ha="right",
                         rotation_mode="anchor")
    ax_c.set_ylabel("model / revPBE0 derivative")  # the basis is in the annotation
    ax_c.set_ylim(0.8, max(ratios) * 1.45)
    # the reference length labels the group, not every bar: six two-line ticks
    # across one panel overlap.
    for r in sorted(set(refs)):
        xs = [xi for xi, rr in zip(x, refs) if rr == r]
        ax_c.annotate(f"r$_{{OH}}$ = {r} Å", (float(np.mean(xs)), -0.155),
                      xycoords=("data", "axes fraction"), fontsize=7.0,
                      ha="center", va="top", annotation_clip=False)
    mu_model = float(np.linalg.norm(model_mu["g00"]))
    mu_dft = float(np.linalg.norm(dft_mu[REF_BASIS]["g00"]))
    # three short lines, left-aligned inside the panel: the one-line form ran off
    # the right edge.
    ax_c.annotate(f"{REF_BASIS}\nμ(eq) D: model {mu_model:.3f},\n"
                  f"revPBE0 {mu_dft:.3f}, expt {MU_EXPERIMENT_D:.3f}",
                  (0.03, 0.97), xycoords="axes fraction", fontsize=7.0, va="top",
                  linespacing=1.35)
    nums.axes("c", ax_c, "mode and reference", "model/DFT derivative ratio")
    nums.panel("c", basis=REF_BASIS,
               ratios={f"{d['mode']}@{'0.957' if abs(d['r_ref']-0.9572)<1e-3 else '1.00'}":
                       d["model"] / d["dft"][REF_BASIS] for d in derivs},
               derivatives_model={f"{d['mode']}@{d['r_ref']:.3f}": d["model"]
                                  for d in derivs},
               derivatives_dft={f"{d['mode']}@{d['r_ref']:.3f}": d["dft"]
                                for d in derivs},
               units=sorted(set(units)),
               permanent_dipole_D={"model": mu_model, "revPBE0_" + REF_BASIS: mu_dft,
                                   "experiment": MU_EXPERIMENT_D},
               intensity_proxy_model_over_dft={
                   "0.96": doc.get("proxy_ratio_0.96"),
                   "1.00": doc.get("proxy_ratio_1.00")},
               xc=doc["xc"])
    _src(nums, "c", "monomer_test/pyscf_monomer.json", "monomer_test/pyscf_monomer.py")

    ax_d = axd["d"]
    lag = json.loads((ARC / "scw" / "lagtest_v2.json").read_text())
    cells = list(lag)
    cmapd = [OKABE_ITO["blue"], OKABE_ITO["vermillion"], OKABE_ITO["green"],
             OKABE_ITO["purple"]]
    rec = {}
    for i, cell in enumerate(cells):
        e = lag[cell]
        n = len(e["acf_signal"])
        lags = np.arange(n, dtype=float)  # fs, per the task brief; the JSON has no unit
        ax_d.plot(lags, e["acf_signal"], color=cmapd[i], lw=1.1, marker="o", ms=2.6,
                  zorder=3)
        ax_d.plot(lags, e["acf_resid"], color=cmapd[i], lw=1.1, ls="--", marker="s",
                  ms=2.6, zorder=3)
        rec[cell] = {"acf_signal": e["acf_signal"], "acf_resid": e["acf_resid"],
                     "pairs": e["pairs"], "resid_mD": e["resid_mD"],
                     "signal_mD": e["signal_mD"],
                     "resid_over_signal": e["resid_mD"] / e["signal_mD"]}
    for j, pair in enumerate(lag[cells[0]]["pairs"]):
        ax_d.annotate(f"{pair}", (j, 1.035), fontsize=7.0, ha="center", color="#555555",
                      annotation_clip=False)
    ax_d.annotate("pairs:", (-0.55, 1.035), fontsize=7.0, ha="right", color="#555555",
                  annotation_clip=False)
    ax_d.set_xlabel("Lag (fs)")
    ax_d.set_ylabel("Normalised autocorrelation")
    ax_d.set_ylim(-0.05, 1.10)
    handles = [Line2D([], [], color=cmapd[i], lw=1.1,
                      label=cells[i].replace("snip_", "").replace("_", " "))
               for i in range(len(cells))]
    handles += [Line2D([], [], color="#666666", lw=1.1, marker="o", ms=2.6,
                       label="signal"),
                Line2D([], [], color="#666666", lw=1.1, ls="--", marker="s", ms=2.6,
                       label="residual")]
    ax_d.legend(handles=handles, loc="lower left", fontsize=7.0, ncol=3)
    nums.axes("d", ax_d, "Lag (fs)", "normalised autocorrelation")
    nums.panel("d", cells=rec, n_lags=len(lag[cells[0]]["acf_signal"]),
               lag_unit_note=("the lag axis is in fs as the task brief states; the "
                              "JSON carries no lag metadata"))
    _src(nums, "d", "ed_data/arc/scw/lagtest_v2.json")

    ax_e = axd["e"]
    npz = residual_npz or (ROOT / "outputs" / "window_residual.npz")
    if npz.exists():
        # A88: this panel used to expect pre-reduced ``band_centre_cm`` and
        # ``A_pred_over_A_label`` arrays, which the delivered file does not carry
        # -- it carries the spectra (nu, A_label, A_pred, A_resid) and both dipole
        # series (M_label, M_pred, 501 frames at 2 fs = 1 ps). The panel therefore
        # could not build at all, which is why A87 had to exclude the whole
        # Extended Data group from its comparison. The reduction is done here.
        #
        # The band ratio is a PLAIN CHANNEL SUM, not a trapezoid. On the file's
        # uniform 33.29 cm^-1 grid the two differ only by the endpoint half
        # weights, but that is enough to matter where a band edge lands on a steep
        # part of the spectrum: sums reproduce the producing script's five ratios
        # to <= 0.0004, trapezoids miss the overtone band by 0.086.
        centres, ratio, resid_frac = [], [], []
        with np.load(npz) as z:
            nu = np.asarray(z["nu"], float)
            a_label = np.asarray(z["A_label"], float)
            a_pred = np.asarray(z["A_pred"], float)
            a_resid = np.asarray(z["A_resid"], float)
            for _tag, (_lo, _hi) in WINDOW_VALIDATION_BANDS_CM.items():
                m = (nu >= _lo) & (nu <= _hi)
                centres.append(0.5 * (_lo + _hi))
                ratio.append(float(a_pred[m].sum() / a_label[m].sum()))
                resid_frac.append(float(a_resid[m].sum() / a_label[m].sum()))
        centres = np.asarray(centres, float)
        ratio = np.asarray(ratio, float)
        ax_e.plot(centres, ratio, ls="none", marker="o", ms=4.0,
                  mfc=OKABE_ITO["blue"], mec="black", mew=0.4, zorder=3)
        ax_e.axhline(1.0, color="black", lw=0.7, ls=":", zorder=1)
        nums.panel("e", band_centre_cm=centres.tolist(),
                   A_pred_over_A_label=ratio.tolist(),
                   residual_over_label=resid_frac,
                   bands_cm={k: list(v) for k, v in WINDOW_VALIDATION_BANDS_CM.items()},
                   reduction="plain channel sum over each band (A88)",
                   source_file=str(npz))
    else:
        ax_e.set_xlim(0, 4000)
        ax_e.set_ylim(0, 2)
        ax_e.axhline(1.0, color="black", lw=0.7, ls=":", zorder=1)
        ax_e.text(0.5, 0.62, "band-resolved A$_{pred}$/A$_{label}$\n1 ps labelled "
                             "trajectory,\n1000 K, 130 kg m$^{-3}$",
                  ha="center", va="center", transform=ax_e.transAxes, fontsize=7.4,
                  color="#444444")
        ax_e.text(0.5, 0.30, "pending", ha="center", va="center",
                  transform=ax_e.transAxes, fontsize=8.0, color="#B00020")
        nums.panel("e", blocked=True, needs=f"{npz} (band_centre_cm, "
                                            f"A_pred_over_A_label)",
                   shows="band-resolved predicted/labelled absorption")
    ax_e.set_xlabel("Band centre ν (cm$^{-1}$)")
    ax_e.set_ylabel("A$_{pred}$ / A$_{label}$")
    nums.axes("e", ax_e, "Band centre (cm^-1)", "A_pred / A_label")

    for letter in "abcde":
        panel_label(axd[letter], letter, dx=-0.20, dy=1.08)
    return save_pdf(fig, "ed_fig5", nums, out_dir=ED_DIR, ledger_name=LEDGER)




# ===========================================================================
CONSERVATIVE_DIVISOR = 4.4  # stated in the task brief, not derived here
SCAN_CSV = "S4_domain_scan_extended.csv"


def boundary_table(cfg: Config) -> list[dict]:
    """The 108 line-list boundaries, with their category against the grid.

    Category is decided exactly as the model decides it: a point is "on computed
    data" when the opacity table's own Delaunay triangulation contains it, and
    clamped otherwise. The three clamp causes are then read off the grid box.
    """
    from opacity.build_tables import get_source
    from tools.broadening_sensitivity import aqua_pressure_pa

    src = get_source("fp", cfg, band=BAND)
    lin = src.table._lin
    scale, offset = np.asarray(lin.scale), np.asarray(lin.offset)
    T_g = np.asarray(src.table.temperature_k, float)
    R_g = np.asarray(src.table.rho_si, float)
    t_lo, t_hi, r_lo, r_hi = T_g.min(), T_g.max(), R_g.min(), R_g.max()

    rows = [r for r in csv.DictReader((cfg.output_path / SCAN_CSV).open())
            if r.get("skipped") != "True"]
    out = []
    for r in rows:
        t_b = float(r["hitemp_max_radiative_t_k"])
        rho_b = float(r["hitemp_max_radiative_rho"])
        if t_b <= 0 or rho_b <= 0:
            continue
        q = np.array([(t_b - offset[0]) / scale[0],
                      (np.log10(rho_b) - offset[1]) / scale[1]])
        inside = int(lin.tri.find_simplex(q)) >= 0
        if inside:
            cat = "on computed data"
        elif t_b < t_lo:
            cat = "T below the grid"
        elif rho_b < r_lo:
            cat = "ρ below the grid"
        else:
            cat = "between computed points"
        p_b, phase = aqua_pressure_pa(cfg, t_b, rho_b)
        out.append({"t_eq_k": float(r["t_eq_k"]),
                    "mass_core_mearth": float(r["mass_core_mearth"]),
                    "wmf": float(r["wmf"]), "t_int_k": float(r["t_int_k"]),
                    "T": t_b, "rho": rho_b, "P": p_b, "aqua_phase": phase,
                    "category": inside and "on computed data" or cat,
                    "inside": inside})
    return out


CATEGORY_STYLE = {
    "on computed data": (OKABE_ITO["blue"], True),
    "T below the grid": (OKABE_ITO["vermillion"], False),
    "ρ below the grid": (OKABE_ITO["orange"], False),
    "between computed points": (OKABE_ITO["purple"], False),
}


def _hull_xy(t_pts, rho_pts):
    from scipy.spatial import ConvexHull

    pts = np.column_stack([np.log10(rho_pts), t_pts])
    h = ConvexHull(pts)
    idx = list(h.vertices) + [h.vertices[0]]
    return 10.0 ** pts[idx, 0], pts[idx, 1]


def _boundary_scatter(ax, bnds, t_pts, rho_pts, markers):
    ax.plot(*_hull_xy(t_pts, rho_pts), color="black", lw=0.9, zorder=4,
            label="grid extent")
    ax.scatter(rho_pts, t_pts, s=9.0, facecolors="none", edgecolors="#555555",
               linewidths=0.5, zorder=3)
    for cat, (col, filled) in CATEGORY_STYLE.items():
        for teq, mk in markers.items():
            sel = [b for b in bnds if b["category"] == cat and b["t_eq_k"] == teq]
            if not sel:
                continue
            ax.scatter([b["rho"] for b in sel], [b["T"] for b in sel], marker=mk,
                       s=20.0, facecolor=col if filled else "none", edgecolor=col,
                       linewidth=0.8, zorder=5)


def ed6(cfg: Config) -> Path:
    """The 108 planets against the computed grid."""
    from opacity.build_tables import compute_point_means, get_source
    from opacity.table import CM2_PER_G_TO_SI

    nums = FigureNumbers(figure="ed_fig6", caption_hint=(
        "a, line-list boundary positions over the computed grid, by category and "
        "equilibrium temperature. b, dispute ratio against the factor that would "
        "flip the boundary, with a conservative bound."))
    bnds = boundary_table(cfg)
    t_pts, rho_pts, _, _ = compute_point_means("fp", cfg, BAND)
    teqs = sorted({b["t_eq_k"] for b in bnds})
    marker_cycle = ["o", "s", "^", "D", "v", "P"]
    markers = {t: marker_cycle[i % len(marker_cycle)] for i, t in enumerate(teqs)}

    fig, (ax_a, ax_b) = figure(86.0, ncols=2)
    _boundary_scatter(ax_a, bnds, t_pts, rho_pts, markers)
    ax_a.set_xscale("log")
    ax_a.set_xlabel("Density ρ (kg m$^{-3}$)")
    ax_a.set_ylabel("Temperature T (K)")
    counts = {c: sum(1 for b in bnds if b["category"] == c) for c in CATEGORY_STYLE}
    handles = [Line2D([], [], ls="none", marker="o",
                      mfc=CATEGORY_STYLE[c][0] if CATEGORY_STYLE[c][1] else "none",
                      mec=CATEGORY_STYLE[c][0], ms=4.4,
                      label=f"{c} ({counts[c]})") for c in CATEGORY_STYLE]
    handles += [Line2D([], [], ls="none", marker=markers[t], mfc="none", mec="#555555",
                       ms=4.4, label=f"T$_{{eq}}$ {t:.0f} K") for t in teqs]
    handles += [Line2D([], [], color="black", lw=0.9, label="grid extent")]
    # shifted right of centre so the panel label has the top-left corner to itself
    ax_a.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.56, 1.003),
                fontsize=7.0, ncol=3, borderaxespad=0.0, columnspacing=0.9)
    nums.axes("a", ax_a, "Density (kg/m^3)", "Temperature (K)", xlog=True)
    nums.panel("a", counts_by_category=counts, n_boundaries=len(bnds),
               n_computed_points=int(t_pts.size),
               T_range_K=[min(b["T"] for b in bnds), max(b["T"] for b in bnds)],
               rho_range_kg_m3=[min(b["rho"] for b in bnds),
                                max(b["rho"] for b in bnds)],
               t_eq_levels_K=teqs,
               category_rule=("'on computed data' means the opacity table's own "
                              "Delaunay triangulation contains the point; the other "
                              "three are the clamp causes read off the grid box"))
    _src(nums, "a", "outputs/S4_domain_scan_extended.csv", "figures/ed_items.py")

    flip = list(csv.DictReader((cfg.output_path / "S8_rcb_flip.csv").open()))
    flip.sort(key=lambda r: float(r["fp_over_hitemp"]))
    ratio = np.array([float(r["fp_over_hitemp"]) for r in flip])
    thresh = np.array([float(r["threshold_factor"]) for r in flip])
    cons = ratio / CONSERVATIVE_DIVISOR
    x = np.arange(len(flip))
    ax_b.vlines(x, thresh, ratio, color=OKABE_ITO["grey"], lw=0.8, zorder=2)
    ax_b.plot(x, thresh, ls="none", marker="_", ms=10, mew=1.8,
              color=OKABE_ITO["sky"], zorder=4, label="flip threshold ∇$_{ad}$/∇$_{rad}$")
    ax_b.plot(x, ratio, ls="none", marker="o", ms=4.2, mfc=OKABE_ITO["blue"],
              mec="black", mew=0.4, zorder=5, label="dispute ratio κ$_R$(FP)/κ$_R$(LL)")
    ax_b.plot(x, cons, ls="none", marker="v", ms=4.0, mfc="none",
              mec=OKABE_ITO["vermillion"], mew=0.9, zorder=5,
              label=f"conservative bound (÷{CONSERVATIVE_DIVISOR:g})")
    ax_b.axhline(1.0, color="black", lw=0.7, ls=":", zorder=1)
    ax_b.set_yscale("log")
    ax_b.set_xticks(x)
    ax_b.set_xticklabels([f"{float(r['t_k']):.0f}\n{float(r['rho_si']):.0f}"
                          for r in flip], fontsize=7.0)
    ax_b.set_xlabel("Boundary:  T (K)  /  ρ (kg m$^{-3}$)")
    ax_b.set_ylabel("Factor")
    ax_b.legend(loc="lower center", bbox_to_anchor=(0.5, 1.003), fontsize=7.0,
                ncol=1, borderaxespad=0.0)
    nums.axes("b", ax_b, "boundary", "factor", ylog=True)
    nums.panel("b", n_shown=len(flip),
               dispute_ratio_range=[float(ratio.min()), float(ratio.max())],
               flip_threshold_range=[float(thresh.min()), float(thresh.max())],
               conservative_divisor=CONSERVATIVE_DIVISOR,
               conservative_range=[float(cons.min()), float(cons.max())],
               n_conservative_still_above_threshold=int((cons > thresh).sum()),
               scope_note=("the flip threshold exists only for the boundaries the "
                           "boundary analysis profiled, because it needs grad_ad and "
                           "grad_rad at the point and the scan CSV stores neither. "
                           "The dispute ratio is available for all boundaries on "
                           "computed data; the threshold is not, so panel b shows "
                           f"the {len(flip)} analysed boundaries rather than all "
                           "of them"),
               conservative_divisor_note=("4.4 is taken from the task brief; it is "
                                          "not derived in this repository"))
    _src(nums, "b", "outputs/S8_rcb_flip.csv", "figures/plot_s8_rcb_flip.py")
    # a's legend occupies the strip above the axes, shifted right; the letter goes
    # in the corner it leaves free.
    panel_label(ax_a, "a", dx=-0.17, dy=1.36)
    panel_label(ax_b, "b", dx=-0.17, dy=1.30)
    return save_pdf(fig, "ed_fig6", nums, out_dir=ED_DIR, ledger_name=LEDGER)


# ===========================================================================
def ed7(cfg: Config) -> Path:
    """External anchor: model radius against the published mass-radius relation."""
    nums = FigureNumbers(figure="ed_fig7", caption_hint=(
        "a, model radius against the published anchor at the matched points. "
        "b, residual against water mass fraction, with the matching-pressure "
        "sensitivity annotated."))
    doc = _json(cfg, "g5_residuals.json")
    conv = _json(cfg, "g5_residuals_conventions.json")
    res = doc["residuals"]
    ours = np.array([r["r_ours"] for r in res])
    theirs = np.array([r["r_aguichine"] for r in res])
    rel = np.array([r["rel"] for r in res])
    wmf = np.array([r["wmf"] for r in res])

    fig, (ax_a, ax_b) = figure(80.0, ncols=2)
    lo = min(ours.min(), theirs.min()) * 0.96
    hi = max(ours.max(), theirs.max()) * 1.04
    ax_a.plot([lo, hi], [lo, hi], color="black", lw=0.8, ls="--", zorder=2)
    sc = ax_a.scatter(theirs, ours, c=wmf, cmap="viridis", s=32, edgecolor="black",
                      linewidth=0.35, zorder=3)
    cb = fig.colorbar(sc, ax=ax_a, pad=0.012, fraction=0.05, aspect=28)
    cb.set_label("Water mass fraction", fontsize=7.6)
    cb.ax.tick_params(labelsize=7.0)
    ax_a.set_xlim(lo, hi)
    ax_a.set_ylim(lo, hi)
    ax_a.set_xlabel("Published radius (R$_⊕$)")
    ax_a.set_ylabel("Model radius (R$_⊕$)")
    nums.axes("a", ax_a, "published radius (R_earth)", "model radius (R_earth)")
    nums.panel("a", n_points=len(res), anchor=doc["anchor"], source=doc["source"],
               t_int_k=doc["t_int_k"], cmf=doc["cmf"],
               published_range_rearth=[float(theirs.min()), float(theirs.max())],
               model_range_rearth=[float(ours.min()), float(ours.max())],
               points=[{k: r[k] for k in ("wmf", "t_irr_k", "m_earth",
                                          "r_aguichine", "r_ours", "rel")}
                       for r in res])

    ax_b.axhline(0.0, color="black", lw=0.7, ls=":", zorder=1)
    ax_b.scatter(wmf, 100.0 * rel, s=32, c=OKABE_ITO["blue"], edgecolor="black",
                 linewidth=0.35, zorder=3)
    med = float(doc["median_abs_rel"])
    worst = float(np.abs(rel).max())
    # One line, not a +/- pair: every residual is negative, and a symmetric pair
    # would draw the eye to scatter about the anchor that is not there.
    med_signed = float(np.median(rel))
    ax_b.axhline(100 * med_signed, color=OKABE_ITO["green"], lw=0.8, ls="--",
                 zorder=2)
    # annotated rather than put in a legend: the only free corner is the one the
    # zero line runs through, and a legend box there sits on top of it.
    ax_b.annotate(f"median {100 * med_signed:+.1f}%\n"
                  f"({int((rel < 0).sum())}/{rel.size} below zero)",
                  (0.97, 0.90), xycoords="axes fraction", fontsize=7.0,
                  ha="right", va="top", color=OKABE_ITO["green"])
    ax_b.set_xlabel("Water mass fraction")
    ax_b.set_ylabel("(model − published) / published  (%)")
    spread = float(conv["matching_pressure_rp_spread_frac"])
    ax_b.annotate(f"matching pressure 1–100 bar:\nR$_p$ spread {100 * spread:.2f}%",
                  (0.03, 0.04), xycoords="axes fraction", fontsize=7.0, va="bottom")
    nums.axes("b", ax_b, "water mass fraction", "relative residual (%)")
    n_neg = int((rel < 0).sum())
    nums.panel("b", median_abs_rel=med, worst_abs_rel=worst,
               median_signed_rel=float(np.median(rel)),
               n_residuals=int(rel.size), n_negative=n_neg,
               sign_note=(f"{n_neg} of {rel.size} residuals are negative: the model "
                          f"is systematically SMALLER than the published anchor, so "
                          f"the quoted median |rel| understates a one-sided bias "
                          f"rather than describing scatter about it"),
               residual_range=[float(rel.min()), float(rel.max())],
               matching_pressure_scan_bar=[c["cap_bar"] for c
                                           in conv["matching_pressure_scan"]],
               matching_pressure_rp_spread_frac=spread,
               matching_pressure_radii_rearth=[c["radius_rearth"] for c
                                               in conv["matching_pressure_scan"]])
    for pl in ("a", "b"):
        _src(nums, pl, "outputs/g5_residuals.json",
             "outputs/g5_residuals_conventions.json",
             "validation/test_g5_literature.py")
    panel_label(ax_a, "a", dx=-0.17, dy=1.06)
    panel_label(ax_b, "b", dx=-0.17, dy=1.06)
    return save_pdf(fig, "ed_fig7", nums, out_dir=ED_DIR, ledger_name=LEDGER)




# ===========================================================================
ED8_CASES = [(5.0, OKABE_ITO["sky"]), (10.0, OKABE_ITO["blue"]),
             (15.0, OKABE_ITO["purple"])]
ED8_WMF, ED8_TEQ, ED8_NODES = 0.50, 500.0, 44
ED8_T_INT_RANGE = (12.0, 450.0)
EARLY_GYR = 0.1


def ed8(cfg: Config) -> Path:
    """Evolution of the water worlds under the two opacities.

    Replays the cooling tracks the boundary-evolution analysis already computes,
    with its settings (44 nodes, T_int 12-450 K) — the same calculation, not a new
    one — because only its per-case summary is stored, not the curves.
    """
    from dataclasses import replace

    from eos import AquaEOS
    from evolution.history import compare_histories

    nums = FigureNumbers(figure="ed_fig8", caption_hint=(
        "Radius, interior entropy and deepest-boundary density against age for "
        "three water worlds under the line list and the composite opacity."))
    eos = AquaEOS.from_config(cfg)
    fig, (ax_r, ax_s, ax_b) = figure(76.0, ncols=3)

    recorded = {}
    for m_core, col in ED8_CASES:
        c = replace(cfg, core=replace(cfg.core, mass_mearth=m_core),
                    envelope=replace(cfg.envelope, f_env=ED8_WMF))
        grid = np.geomspace(*ED8_T_INT_RANGE, ED8_NODES)
        h = compare_histories(c, "composite", "hitemp", ED8_TEQ, eos=eos,
                              t_int_grid_k=grid)
        age = h.age_gyr
        for ax, a, b in ((ax_r, h.radius_a, h.radius_b),
                         (ax_s, h.entropy_a, h.entropy_b),
                         (ax_b, h.rcb_rho_a, h.rcb_rho_b)):
            ax.plot(age, b, color=col, lw=1.2, ls="--", zorder=3)
            ax.plot(age, a, color=col, lw=1.2, ls="-", zorder=4)
        floor = h.radiogenic_floor_age_gyr(m_core)
        if floor is not None:
            for ax in (ax_r, ax_s, ax_b):
                ax.axvline(float(floor), color=col, lw=0.8, ls=":", alpha=0.9, zorder=2)
        recorded[f"M_c={m_core:g}"] = {
            "age_range_gyr": [float(age.min()), float(age.max())],
            "radius_composite_rearth": [float(h.radius_a.min()), float(h.radius_a.max())],
            "radius_line_list_rearth": [float(h.radius_b.min()), float(h.radius_b.max())],
            "entropy_composite_J_kg_K": [float(h.entropy_a.min()),
                                         float(h.entropy_a.max())],
            "entropy_line_list_J_kg_K": [float(h.entropy_b.min()),
                                         float(h.entropy_b.max())],
            "rcb_rho_composite_kg_m3": [float(np.nanmin(h.rcb_rho_a)),
                                        float(np.nanmax(h.rcb_rho_a))],
            "rcb_rho_line_list_kg_m3": [float(np.nanmin(h.rcb_rho_b)),
                                        float(np.nanmax(h.rcb_rho_b))],
            "radiogenic_floor_age_gyr": (float(floor) if floor is not None else None),
            "radius_frac_diff_at_10gyr": float(np.interp(
                10.0, age, h.radius_frac_diff)),
            "entropy_frac_diff_at_10gyr": float(np.interp(
                10.0, age, h.entropy_frac_diff)),
        }

    for ax, ylab, ylog in ((ax_r, "Radius R$_p$ (R$_⊕$)", False),
                           (ax_s, "Interior entropy (J kg$^{-1}$ K$^{-1}$)", False),
                           (ax_b, "Deepest boundary ρ (kg m$^{-3}$)", True)):
        ax.axvspan(ax.get_xlim()[0], EARLY_GYR, color="#000000", alpha=0.07, lw=0,
                   zorder=0)
        ax.set_xscale("log")
        if ylog:
            ax.set_yscale("log")
        ax.set_xlabel("Age (Gyr)")
        ax.set_ylabel(ylab)
        ax.set_xlim(max(1e-4, ax.get_xlim()[0]), 12.0)
    handles = [Line2D([], [], color=c, lw=1.2, label=f"M$_c$ = {m:g} M$_⊕$")
               for m, c in ED8_CASES]
    handles += [Line2D([], [], color="#666666", lw=1.2, ls="-", label="composite"),
                Line2D([], [], color="#666666", lw=1.2, ls="--", label="line list")]
    ax_r.legend(handles=handles, loc="lower left", fontsize=7.0)

    any_floor = any(v["radiogenic_floor_age_gyr"] is not None
                    for v in recorded.values())
    for panel, ax, x, y, ylog in (("a", ax_r, "Age (Gyr)", "R_p (R_earth)", False),
                                  ("b", ax_s, "Age (Gyr)",
                                   "interior entropy (J/kg/K)", False),
                                  ("c", ax_b, "Age (Gyr)",
                                   "deepest boundary rho (kg/m^3)", True)):
        nums.axes(panel, ax, x, y, xlog=True, ylog=ylog)
        nums.panel(panel, cases=recorded, water_mass_fraction=ED8_WMF,
                   t_eq_k=ED8_TEQ, n_t_int_nodes=ED8_NODES,
                   shaded_early_epoch_gyr=EARLY_GYR,
                   radiogenic_floor_marked=any_floor,
                   floor_note=("a dotted vertical line marks the radiogenic floor "
                               "crossing where the track reaches it; these water "
                               "worlds " + ("do" if any_floor else "do not") +
                               " reach theirs inside the track"))
        _src(nums, panel, "evolution/history.py", "figures/plot_s9_case_evolution.py")
    for letter, ax in zip("abc", (ax_r, ax_s, ax_b)):
        panel_label(ax, letter, dx=-0.24, dy=1.07)
    return save_pdf(fig, "ed_fig8", nums, out_dir=ED_DIR, ledger_name=LEDGER)


# ===========================================================================
SCALING_PANELS = (("window_2550", "κ at 2550 cm$^{-1}$"),
                  ("window_1000_1200", "κ over 1000-1200 cm$^{-1}$"),
                  ("rosseland", "κ$_R$, 100-4000 cm$^{-1}$"))
SOURCE_STYLE = {"fp": ("-", 1.1, "o"),
                "hitemp": ("--", 0.9, "s"),
                "hitemp_mtckd": (":", 0.9, "^")}
SOURCE_LABEL = {"fp": "first principles", "hitemp": "HITEMP",
                "hitemp_mtckd": "HITEMP + MT_CKD"}




# A88: nine items, not ten -- ED Fig. 9's scaling panels are now ED Fig. 4c-i.
ITEMS = {"ed1": ed1, "ed2": ed2, "ed3": ed3, "ed4": ed4, "ed5": ed5,
         "ed6": ed6, "ed7": ed7, "ed8": ed8}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("config", type=Path)
    ap.add_argument("--only", nargs="*", default=None)
    args = ap.parse_args(argv)
    cfg = load_config(args.config)
    set_manuscript_style()
    ED_DIR.mkdir(parents=True, exist_ok=True)
    names = args.only or list(ITEMS)
    for name in names:
        if name not in ITEMS:
            raise SystemExit(f"unknown item {name!r}; have {sorted(ITEMS)}")
        print(f"wrote {ITEMS[name](cfg)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
