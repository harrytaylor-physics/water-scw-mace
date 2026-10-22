"""Shared style and bookkeeping for the manuscript figures (A62).

One place for everything the journal fixes, so the four figure scripts carry
physics and nothing else:

* **size** — 183 mm double-column width, given in millimetres and converted once;
* **type** — sans-serif, every label at least 7 pt *at final size* (the figures
  are written at 1:1 scale, so the rcParam point sizes are the final ones);
* **colour** — the Okabe-Ito colourblind-safe palette, and nothing outside it;
* **panel labels** — bold lower-case a/b placed outside the axes;
* **no titles inside panels** — axis labels and the caption carry the meaning;
* **vector output** — PDF with text kept as text, into ``figs/``.

It also collects, per panel, every number a reader could take off the figure —
axis ranges, marked values, point counts — into ``figs/fig_numbers.json``, so a
caption can be checked against the figure it describes rather than against the
author's memory. :class:`FigureNumbers` is that ledger.

The 52 first-principles points are one grid. Nothing here, and nothing in the
figures, distinguishes how or when a point was delivered.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
import os
import sys
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")  # headless
import matplotlib.pyplot as plt  # noqa: E402

# A87: figure and Extended Data ledgers live outside cfg.output_path, so a run
# against a different spectra version could not be isolated -- the v3 column of
# A86's table silently carried v5 figure numbers because of it. The root is now
# redirectable. Since the migration into water-scw-mace it defaults to
# docs/figures/subneptune/, beside the repository's other figures (layout.py).
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from layout import FIGS_DIR  # noqa: E402
ED_DIR = FIGS_DIR / "ed"  # Extended Data items

MM = 1.0 / 25.4  # millimetres -> inches
WIDTH_MM = 183.0  # Nature double column
SINGLE_WIDTH_MM = 89.0  # Nature single column, where one panel suffices
WIDTH_IN = WIDTH_MM * MM

# Okabe-Ito, colourblind-safe. Named by role, not by hue, so a figure that
# changes colour does not change meaning.
OKABE_ITO = {
    "black": "#000000",
    "orange": "#E69F00",
    "sky": "#56B4E9",
    "green": "#009E73",
    "yellow": "#F0E442",
    "blue": "#0072B2",
    "vermillion": "#D55E00",
    "purple": "#CC79A7",
    "grey": "#999999",
}
FP_C = OKABE_ITO["blue"]
HITEMP_C = OKABE_ITO["vermillion"]
POKAZATEL_C = OKABE_ITO["green"]
CONTINUUM_C = OKABE_ITO["purple"]
COMPOSITE_C = OKABE_ITO["black"]

# Minimum type size anywhere in a figure, in points at final size.
MIN_PT = 7.0


def set_manuscript_style() -> None:
    """Journal-fixed rcParams. Called once per figure script, before plotting."""
    plt.rcParams.update({
        "figure.dpi": 200,
        "savefig.dpi": 600,
        "pdf.fonttype": 42,  # TrueType: text stays selectable/searchable, not paths
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Helvetica", "Arial"],
        "font.size": 8.0,
        "axes.labelsize": 8.0,
        "axes.titlesize": 8.0,
        "xtick.labelsize": MIN_PT,
        "ytick.labelsize": MIN_PT,
        "legend.fontsize": MIN_PT,
        "figure.titlesize": 8.0,
        "axes.linewidth": 0.6,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.minor.width": 0.45,
        "ytick.minor.width": 0.45,
        "xtick.major.size": 2.6,
        "ytick.major.size": 2.6,
        "lines.linewidth": 1.1,
        "lines.markersize": 3.5,
        "legend.frameon": False,
        "legend.handlelength": 1.4,
        "legend.handletextpad": 0.5,
        "legend.labelspacing": 0.3,
        "axes.grid": False,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.axisbelow": True,
    })


def panel_label(ax, letter: str, dx: float = -0.085, dy: float = 1.045) -> None:
    """Bold lower-case panel label, placed outside the axes."""
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=9.0,
            fontweight="bold", va="top", ha="left")


def figure(height_mm: float, width_mm: float = WIDTH_MM, **kwargs):
    """A figure of the requested height at a journal column width.

    Constrained layout is the default because the width is a deliverable: the
    content must be fitted into the canvas, not the canvas cropped to the content.
    ``width_mm`` is 183 (double column) or 89 (single, where one panel suffices).
    """
    if width_mm not in (WIDTH_MM, SINGLE_WIDTH_MM):
        raise ValueError(
            f"width {width_mm} mm is neither the double-column {WIDTH_MM:g} nor the "
            f"single-column {SINGLE_WIDTH_MM:g}")
    kwargs.setdefault("layout", "constrained")
    return plt.subplots(figsize=(width_mm * MM, height_mm * MM), **kwargs)


def mosaic(height_mm: float, layout, width_mm: float = WIDTH_MM, **kwargs):
    """A journal-width figure laid out by ``subplot_mosaic``.

    For panels that are not all the same size — a wide one over two narrow ones,
    say. Mixing ``add_subplot`` into a grid built by ``subplots`` fights the
    constrained-layout engine and overlaps the axes; a mosaic does not.
    """
    if width_mm not in (WIDTH_MM, SINGLE_WIDTH_MM):
        raise ValueError(f"width {width_mm} mm is not a journal column width")
    kwargs.setdefault("layout", "constrained")
    fig = plt.figure(figsize=(width_mm * MM, height_mm * MM), **kwargs)
    return fig, fig.subplot_mosaic(layout)


@dataclass
class FigureNumbers:
    """Ledger of every number a reader could read off a figure.

    Captions get checked against this, not against recollection. Each panel
    records its axis ranges and whatever it marks, names or counts.
    """

    figure: str
    caption_hint: str = ""
    panels: dict[str, dict[str, Any]] = field(default_factory=dict)

    def panel(self, letter: str, **entries: Any) -> None:
        self.panels.setdefault(letter, {}).update(entries)

    def axes(self, letter: str, ax, x: str, y: str, xlog=False, ylog=False) -> None:
        """Record an axis' label, scale and limits exactly as drawn."""
        self.panel(letter, **{
            "x_label": x, "y_label": y,
            "x_scale": "log" if xlog else "linear",
            "y_scale": "log" if ylog else "linear",
            "x_range": [float(v) for v in ax.get_xlim()],
            "y_range": [float(v) for v in ax.get_ylim()],
        })


_FORBIDDEN = ("original", "extension", "campaign")


def _scan_text(fig) -> list[str]:
    """Every string drawn into the figure, for the vocabulary check."""
    out = []
    for t in fig.findobj(matplotlib.text.Text):
        s = t.get_text()
        if s:
            out.append(s)
    return out


def _too_small(fig) -> list[tuple[str, float]]:
    """Every drawn string below the minimum type size, with its size."""
    out = []
    for t in fig.findobj(matplotlib.text.Text):
        label = t.get_text()
        if not label:
            continue
        size = float(t.get_fontsize())
        if size < MIN_PT - 1e-9:
            out.append((label.replace("\n", " ")[:40], size))
    return out


def save_pdf(fig, stem: str, numbers: FigureNumbers | None = None,
             out_dir: Path | None = None, ledger_name: str = "fig_numbers.json") -> Path:
    """Write ``figs/<stem>.pdf`` as vector, after the vocabulary and type checks.

    Two rules are enforced here rather than in review, because both are invisible
    once the figure is a PDF and expensive to discover late:

    * the manuscript treats the 52 points as one grid, so a figure that says
      "original", "extension" or "campaign" anywhere is a bug;
    * nothing may be set below :data:`MIN_PT` at final size. Setting the rcParams
      is not enough — a per-call ``fontsize=`` overrides them silently, which is
      exactly how sub-minimum type gets into a figure.
    """
    offenders = [
        s for s in _scan_text(fig)
        if any(word in s.lower() for word in _FORBIDDEN)
    ]
    if offenders:
        raise AssertionError(
            f"{stem}: figure text refers to the delivery history, which the "
            f"manuscript does not: {offenders}"
        )
    small = _too_small(fig)
    if small:
        raise AssertionError(
            f"{stem}: {len(small)} text element(s) below {MIN_PT:g} pt at final "
            f"size: " + ", ".join(f"{lab!r} at {sz:g} pt" for lab, sz in small[:6])
        )
    out_dir = FIGS_DIR if out_dir is None else out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / f"{stem}.pdf"
    # NOT bbox_inches="tight": that crops the page to the drawn content, so the
    # PDF comes out 155-182 mm wide instead of the 183 mm the journal fixes. The
    # canvas size is the deliverable, so the layout engine has to fit the content
    # into it — every figure here is built with layout="constrained".
    fig.savefig(path, format="pdf")
    # A 200-dpi PNG of the same canvas beside it, for display in the repository
    # (water-scw-mace keeps PNGs in docs/figures/). The PDF is the deliverable.
    fig.savefig(path.with_suffix(".png"), format="png", dpi=200)
    plt.close(fig)
    if numbers is not None:
        merge_numbers(numbers, out_dir=out_dir, name=ledger_name)
    return path


def merge_numbers(numbers: FigureNumbers, out_dir: Path | None = None,
                  name: str = "fig_numbers.json") -> Path:
    """Merge one figure's ledger into ``<out_dir>/<name>``."""
    out_dir = FIGS_DIR if out_dir is None else out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / name
    doc = json.loads(path.read_text()) if path.exists() else {
        "_meta": (
            "Every number a reader could take off each panel: axis ranges as "
            "drawn, values marked or annotated, and point counts. Captions are "
            "checked against this file. Widths are 183 mm (double column); type "
            f"is >= {MIN_PT:g} pt at final size; the palette is Okabe-Ito."
        ),
        "width_mm": WIDTH_MM,
        "min_font_pt": MIN_PT,
    }
    doc[numbers.figure] = {
        "caption_hint": numbers.caption_hint,
        "panels": numbers.panels,
    }
    # Deterministic order: the header keys first, then the figures by name, so the
    # file diffs cleanly whatever order the scripts were run in.
    header = {k: doc[k] for k in doc if k.startswith("_") or k in ("width_mm",
                                                                  "min_font_pt")}
    figures = {k: doc[k] for k in sorted(doc) if k not in header}
    path.write_text(json.dumps({**header, **figures}, indent=1) + "\n")
    return path
