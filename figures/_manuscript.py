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

The 48 first-principles points are one grid. Nothing here, and nothing in the
figures, distinguishes how or when a point was delivered.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")  # headless
import matplotlib.pyplot as plt  # noqa: E402

FIGS_DIR = Path(__file__).resolve().parents[1] / "figs"

MM = 1.0 / 25.4  # millimetres -> inches
WIDTH_MM = 183.0  # Nature double column
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


def figure(height_mm: float, **kwargs):
    """A 183 mm-wide figure of the requested height, constrained by default.

    Constrained layout is the default because the width is a deliverable: the
    content must be fitted into the canvas, not the canvas cropped to the content.
    """
    kwargs.setdefault("layout", "constrained")
    return plt.subplots(figsize=(WIDTH_IN, height_mm * MM), **kwargs)


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


def save_pdf(fig, stem: str, numbers: FigureNumbers | None = None) -> Path:
    """Write ``figs/<stem>.pdf`` as vector, after the vocabulary and type checks.

    Two rules are enforced here rather than in review, because both are invisible
    once the figure is a PDF and expensive to discover late:

    * the manuscript treats the 48 points as one grid, so a figure that says
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
    FIGS_DIR.mkdir(parents=True, exist_ok=True)
    path = FIGS_DIR / f"{stem}.pdf"
    # NOT bbox_inches="tight": that crops the page to the drawn content, so the
    # PDF comes out 155-182 mm wide instead of the 183 mm the journal fixes. The
    # canvas size is the deliverable, so the layout engine has to fit the content
    # into it — every figure here is built with layout="constrained".
    fig.savefig(path, format="pdf")
    plt.close(fig)
    if numbers is not None:
        merge_numbers(numbers)
    return path


def merge_numbers(numbers: FigureNumbers) -> Path:
    """Merge one figure's ledger into ``figs/fig_numbers.json``."""
    FIGS_DIR.mkdir(parents=True, exist_ok=True)
    path = FIGS_DIR / "fig_numbers.json"
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
