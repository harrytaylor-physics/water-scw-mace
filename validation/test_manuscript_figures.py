"""Gates for the manuscript figures and the one-grid rule (A62).

These check the things a journal and a reader check — vector text, type size,
width, vocabulary, and whether the caption can be verified against the figure —
without re-rendering, because Fig. 3 builds two structure models and Fig. 2
rebuilds the opacity tables. The figures are produced by
``python figures/figN.py config/fiducial.yaml``; if they have not been produced
on this machine the file-level checks skip and the rules themselves are still
gated.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from config import load_config
from figures._manuscript import (
    MIN_PT,
    OKABE_ITO,
    WIDTH_MM,
    FigureNumbers,
    save_pdf,
    set_manuscript_style,
)

ROOT = Path(__file__).resolve().parents[1]
FIGS = ROOT / "figs"
FIGURES = ("fig1", "fig2", "fig3", "fig4")


@pytest.fixture(scope="module")
def cfg():
    return load_config(ROOT / "config" / "fiducial.yaml")


@pytest.fixture(scope="module")
def numbers():
    path = FIGS / "fig_numbers.json"
    if not path.exists():
        pytest.skip("figs/fig_numbers.json not built on this machine")
    return json.loads(path.read_text())


# -- the rules, gated whether or not the figures exist ---------------------
def test_journal_constraints_are_what_we_think_they_are():
    assert WIDTH_MM == 183.0  # Nature double column
    assert MIN_PT == 7.0


def test_style_puts_nothing_below_the_minimum_type_size():
    import matplotlib.pyplot as plt

    set_manuscript_style()
    for key in ("font.size", "axes.labelsize", "xtick.labelsize",
                "ytick.labelsize", "legend.fontsize", "figure.titlesize"):
        assert float(plt.rcParams[key]) >= MIN_PT, key
    assert plt.rcParams["pdf.fonttype"] == 42  # TrueType: text stays text


def test_the_palette_is_okabe_ito():
    """Eight named entries, all valid hex, no duplicates."""
    assert len(OKABE_ITO) == 9  # eight colours plus black
    assert len(set(OKABE_ITO.values())) == len(OKABE_ITO)
    for name, hexcode in OKABE_ITO.items():
        assert hexcode.startswith("#") and len(hexcode) == 7, name
        int(hexcode[1:], 16)


@pytest.mark.parametrize("word", ["original", "Extension", "CAMPAIGN"])
def test_a_figure_naming_the_delivery_history_is_refused(tmp_path, word):
    """The one-grid rule is enforced where the file is written, not in review."""
    import matplotlib.pyplot as plt

    set_manuscript_style()
    fig, ax = plt.subplots()
    ax.set_xlabel(f"a label mentioning the {word} grid")
    with pytest.raises(AssertionError, match="delivery history"):
        save_pdf(fig, "should_not_exist")
    plt.close(fig)
    assert not (FIGS / "should_not_exist.pdf").exists()


def test_a_clean_figure_is_written_as_vector_with_live_text(tmp_path):
    import matplotlib.pyplot as plt

    set_manuscript_style()
    fig, ax = plt.subplots()
    ax.plot([1, 2], [1, 2])
    ax.set_xlabel("Density ρ (kg m$^{-3}$)")
    nums = FigureNumbers(figure="_gate_probe")
    nums.axes("", ax, "Density (kg/m^3)", "value")
    try:
        path = save_pdf(fig, "_gate_probe", nums)
        head = path.read_bytes()[:5]
        assert head == b"%PDF-"
        body = path.read_bytes()
        assert b"/FontFile2" in body or b"/TrueType" in body  # text, not outlines
    finally:
        (FIGS / "_gate_probe.pdf").unlink(missing_ok=True)
        _drop_probe_from_ledger()


def _drop_probe_from_ledger() -> None:
    path = FIGS / "fig_numbers.json"
    if not path.exists():
        return
    doc = json.loads(path.read_text())
    if doc.pop("_gate_probe", None) is not None:
        path.write_text(json.dumps(doc, indent=1, sort_keys=False) + "\n")


# -- the produced figures, when they are present ---------------------------
@pytest.mark.parametrize("stem", FIGURES)
def test_each_figure_is_a_vector_pdf(stem):
    path = FIGS / f"{stem}.pdf"
    if not path.exists():
        pytest.skip(f"{path} not built on this machine")
    body = path.read_bytes()
    assert body[:5] == b"%PDF-"
    assert b"/FontFile2" in body or b"/TrueType" in body
    assert len(body) > 5_000


@pytest.mark.parametrize("stem", FIGURES)
def test_every_panel_records_the_numbers_a_caption_would_quote(numbers, stem):
    assert stem in numbers, f"{stem} missing from figs/fig_numbers.json"
    entry = numbers[stem]
    assert entry["caption_hint"]
    assert entry["panels"], f"{stem} recorded no panels"
    for letter, panel in entry["panels"].items():
        if "x_range" in panel:  # a data panel, not a schematic
            assert len(panel["x_range"]) == 2 and panel["x_range"][0] < panel["x_range"][1]
            assert len(panel["y_range"]) == 2 and panel["y_range"][0] < panel["y_range"][1]
            assert panel["x_label"] and panel["y_label"]


def test_the_ledger_never_names_the_delivery_history(numbers):
    blob = json.dumps(numbers).lower()
    for word in ("original", "extension", "campaign"):
        assert word not in blob, f"figs/fig_numbers.json says {word!r}"


def test_the_ledger_agrees_with_the_grid_size(numbers):
    """Point counts a caption would quote, checked against the grid itself."""
    if "fig2" in numbers:
        assert numbers["fig2"]["panels"]["a"]["n_computed_points"] == 48
    if "fig4" in numbers:
        assert numbers["fig4"]["panels"][""]["n_first_principles_points"] == 48


# -- the one-grid rule in the outputs, not just the figures ----------------
def test_crossover_rows_carry_no_delivery_history(cfg):
    from opacity.crossover import GROUPS, crossover_v3_rows

    rows = crossover_v3_rows(cfg, (100.0, 4000.0))
    assert len(rows) == 48
    assert set(GROUPS) == {"dilute", "cold_dense", "dense"}
    for r in rows:
        assert "is_extension_point" not in r
        assert "in_delivered_set" not in r
        assert r["group"] in GROUPS


def test_the_groups_are_defined_by_where_a_point_sits(cfg):
    from opacity.crossover import (
        COLD_DENSE_RHO_MAX,
        COLD_DENSE_T_MAX,
        DILUTE_RHO_MAX,
        crossover_v3_rows,
    )

    for r in crossover_v3_rows(cfg, (100.0, 4000.0)):
        if r["rho"] <= DILUTE_RHO_MAX:
            assert r["group"] == "dilute"
        elif r["T"] <= COLD_DENSE_T_MAX and r["rho"] <= COLD_DENSE_RHO_MAX:
            assert r["group"] == "cold_dense"
        else:
            assert r["group"] == "dense"
