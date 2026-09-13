"""Band-truncation gate for the line-list means (SPEC §3; A39).

SPEC §3 states the absolute band systematic as a limitation and requires it to be
quantified by comparing full-band with common-band means for the line-list
source. These checks pin the parts that must not drift:

  * ``--band full`` resolves to HITEMP's *own* computed wavenumber range — the
    intersection of what the files contain, never wider (the loader would refuse
    a wider band, and a silently narrower one would understate the systematic);
  * truncating to the science band biases κ_R **high** at every grid point, which
    is the direction the limitation claims;
  * the default band is unchanged by any of this (diagnostic only);
  * the CLI writes the per-point table it advertises.
"""

from __future__ import annotations

import csv

import numpy as np
import pytest

import opacity.build_tables as bt
from opacity import compute_point_means, detect_common_band, resolve_band
from opacity.registry import get_spec

SOURCE = "hitemp"


@pytest.fixture(scope="module")
def bands(cfg):
    return resolve_band(cfg, [SOURCE]), detect_common_band(cfg, [SOURCE])


@pytest.fixture(scope="module")
def means(cfg, bands):
    common, full = bands
    _, _, kr_c, kp_c = compute_point_means(SOURCE, cfg, common)
    t, rho, kr_f, kp_f = compute_point_means(SOURCE, cfg, full)
    return t, rho, kr_c, kp_c, kr_f, kp_f


def test_full_band_is_exactly_what_the_files_contain(cfg, bands):
    """`--band full` = (max of file minima, min of file maxima); never wider."""
    _, full = bands
    spec = get_spec(SOURCE)
    mins, maxs = [], []
    for path in sorted((cfg.data_root_path / spec.subdir).glob(spec.glob)):
        with np.load(path) as d:
            wn = np.asarray(d["wn"], float)
        mins.append(float(wn.min()))
        maxs.append(float(wn.max()))
    assert full == (max(mins), min(maxs))
    assert full[0] >= min(mins) and full[1] <= max(maxs)  # inside the data, always
    assert bt.resolve_band_spec(cfg, SOURCE, "full") == full
    assert bt.resolve_band_spec(cfg, SOURCE, "common") == resolve_band(cfg, [SOURCE])
    assert bt.resolve_band_spec(cfg, SOURCE, "50,4500") == (50.0, 4500.0)


def test_common_band_biases_rosseland_high(means):
    """Band truncation drops transparency windows from a harmonic mean (A39)."""
    _, _, kr_c, _, kr_f, _ = means
    assert np.all(kr_f < kr_c), "κ_R(full) should be below κ_R(common) at every point"
    # The Planck mean is arithmetic and far less window-sensitive.
    assert np.median(kr_f / kr_c) < 0.5


def test_truncation_grows_with_temperature(means):
    """The Planck peak marches out of the band with T, so the ratio falls with T."""
    t, _, kr_c, _, kr_f, _ = means
    ratio = kr_f / kr_c
    cold = ratio[t == t.min()].mean()
    hot = ratio[t == t.max()].mean()
    assert cold > hot
    assert cold > 0.9  # at 500 K the band holds essentially all the weight


def test_default_band_unchanged_by_the_diagnostic(cfg):
    """The band study must not move the science default (A25)."""
    assert cfg.opacity.band_cm == (100.0, 4000.0)
    assert resolve_band(cfg, [SOURCE]) == (100.0, 4000.0)


def test_cli_writes_a_point_means_table(cfg, repo_root, tmp_path, capsys):
    """`python opacity/build_tables.py <cfg> --source hitemp --band full`."""
    bt._cli([str(repo_root / "config" / "fiducial.yaml"), "--source", SOURCE,
             "--band", "full", "--out-dir", str(tmp_path)])
    out = capsys.readouterr().out
    assert "band=full" in out
    csv_path = tmp_path / f"{SOURCE}_means_full.csv"
    npz_path = tmp_path / f"{SOURCE}_means_full.npz"
    assert csv_path.exists() and npz_path.exists()
    rows = list(csv.DictReader(csv_path.open()))
    full = detect_common_band(cfg, [SOURCE])
    assert len(rows) == np.load(npz_path)["temperature_k"].size
    assert float(rows[0]["band_lo_cm"]) == pytest.approx(full[0])
    assert float(rows[0]["band_hi_cm"]) == pytest.approx(full[1])
    assert all(float(r["kappa_R_cm2_g"]) > 0 for r in rows)
