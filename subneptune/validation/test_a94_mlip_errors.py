"""The MLIP error numbers are the production model's, read from the log (A94).

The training log ``production_campaign/arc_scripts/logs/train_v4_7924966.out`` prints two
"Error-table on TRAIN and VALID" blocks: one for the stage-one model (epoch 109)
and one for the stage-two SWA model (epoch 186). A90 transcribed the stage-two
table (44.4 meV/A, 7.0 meV/atom held out) as the production potential's error.
Production MD ran the stage-one model (sha256 dc4a8977...), whose held-out error
is 57.8 meV/A and 20.8 meV/atom.

These gates re-read both tables from the log, so a hand edit of
``models/mlip_heldout_errors.json`` that swaps them back -- or a transcription
slip -- fails here rather than reaching Fig. 1a and the ledger.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

from layout import HOST, MLIP_ERRORS  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
HELD = json.loads(MLIP_ERRORS.read_text())
LOG = HOST / HELD["source_log"]
PRODUCTION_SHA256 = "dc4a8977d44e58cf959901b39a594abd3daf69fb09e5622bbe7e8d3e74b7f06c"

_STAGE = re.compile(r"Loaded Stage (one|two) model from epoch (\d+)")
_ROW = re.compile(r"\|\s*(train_Default|valid_Default)\s*\|\s*([\d.]+)\s*\|"
                  r"\s*([\d.]+)\s*\|\s*([\d.]+)\s*\|")


def _tables() -> dict[str, dict]:
    """``{"one": {"epoch": 109, "valid_Default": (E, F, relF), ...}, "two": ...}``."""
    out: dict[str, dict] = {}
    stage = None
    for line in LOG.read_text(errors="ignore").splitlines():
        if m := _STAGE.search(line):
            stage = m.group(1)
            out[stage] = {"epoch": int(m.group(2))}
        elif stage and (m := _ROW.search(line)):
            out[stage][m.group(1)] = tuple(float(x) for x in m.group(2, 3, 4))
    return out


def _as_tuple(d: dict) -> tuple[float, float, float]:
    return (d["energy_meV_per_atom"], d["force_meV_per_A"],
            d["relative_force_rmse_percent"])


def test_the_log_is_in_the_repository():
    assert HELD["source_in_repository"] is True
    assert LOG.is_file(), f"{LOG} is cited but missing"


@pytest.mark.parametrize("tag", ["valid_Default", "train_Default"])
def test_quoted_errors_are_the_stage_one_table(tag):
    t = _tables()
    assert t["one"]["epoch"] == 109
    assert _as_tuple(HELD["errors"][tag]) == t["one"][tag]
    assert HELD["production_model"]["training_stage"] == "stage one"
    assert HELD["production_model"]["sha256"] == PRODUCTION_SHA256


@pytest.mark.parametrize("tag", ["valid_Default", "train_Default"])
def test_stage_two_is_recorded_apart_and_not_quoted(tag):
    t = _tables()
    two = HELD["stage_two_not_used_for_production"]
    assert t["two"]["epoch"] == 186
    assert _as_tuple(two["errors"][tag]) == t["two"][tag]
    assert _as_tuple(HELD["errors"][tag]) != t["two"][tag]


def test_fig1a_caption_quotes_the_production_model():
    pytest.importorskip("numpy")
    from config import load_config
    from tools.pipeline_facts import pipeline_facts, stage_captions

    cfg = load_config(ROOT / "config" / "fiducial.yaml")
    facts = pipeline_facts(cfg)
    if not facts:
        pytest.skip("no FP spectra on this machine; captions have no source")
    cap = stage_captions(facts)[2]
    assert "57.8" in cap and "20.8" in cap
    assert "44.4" not in cap
