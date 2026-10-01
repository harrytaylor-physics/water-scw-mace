#!/usr/bin/env python3
"""The pipeline facts Fig. 1a states and the ledger records as ``pipeline.*`` (A89).

One function, two callers. Figure 1a annotates each stage of the pipeline with the
number that validates it, and ``tools/manuscript_numbers.py`` records the same
numbers in the ledger; if each dug them out separately they could drift apart, and
a schematic quoting a number the ledger does not carry is a number no reader can
check. So both call this.

Each fact names the artefact that actually carries it rather than prose about it.
That matters here: ``METHODS.md`` disagrees with the delivered data on two counts
-- it says 60 ps per state point and a 27-point grid, where the dipole series give
30 ps and the spectra directory holds 52 cells -- and where they disagree the data
win. Where only prose exists (the training-set composition, the level of theory)
the fact is marked with that provenance.

The two validation numbers are quoted as the pull states them, which is not how a
draft caption stated them:

  * the force parity is **training-set** parity by regime, 41 meV/A bulk and
    58 meV/A at RCB conditions, and ``reference.json`` itself calls this
    "optimistic vs true held-out error". The only held-out parity anywhere in the
    pull is 12.3 / 20.2 meV/A and belongs to the LDA model, not this one.
  * the NVE check is a **drift rate** bound, 0.017 meV/atom/ps worst over 21 runs
    against a 2 meV/atom/ps gate -- not an energy error per atom.
"""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ARC = ROOT / "ed_data" / "arc"
PARITY_REF = "ed_data/arc/water-scw-mace/benchmarks/reference/reference.json"
METHODS = "ed_data/arc/water-scw-mace/docs/METHODS.md"
CAMPAIGN = "ed_data/arc/water-scw-mace/docs/LDA_CAMPAIGN.md"
REGEN = "tools/v5_regenerate.py"

# From LDA_CAMPAIGN.md, which records the worst value over the 21 runs rather than
# one per state point; figures/ed_items.py panel d plots the same pair.
NVE_DRIFT_MEV_ATOM_PS = 0.017
NVE_DRIFT_GATE = 2.0


def pipeline_facts(cfg) -> dict[str, dict]:
    """``{key: {value, unit, source_script, source_file, note}}``, newest data first.

    Returns ``{}`` when the spectra directory is empty, so a caller on a machine
    without the data gets nothing rather than a half-populated schematic.
    """
    import numpy as np

    spectra = sorted(cfg.fp_opacity_path.glob("T*_opacity.npz"))
    if not spectra:
        return {}
    counts, steps, frames = [], set(), set()
    for f in spectra:
        with np.load(f) as z:
            counts.append(int(z["N"]))
            steps.add(float(z["dt_fs"]))
    dip = ROOT / "ed_data" / "arc_dipoles"
    if dip.is_dir():
        for d in sorted(dip.iterdir()):
            m = d / "dipoles.npy"
            if m.is_file():
                frames.add(int(np.load(m, mmap_mode="r").shape[0]))
    dt = sorted(steps)[0]
    here = str(cfg.fp_opacity_path)

    facts: dict[str, dict] = {
        "n_state_points": dict(
            value=len(spectra), unit="-", source_script=REGEN, source_file=here,
            note=None),
        "molecules_per_cell": dict(
            value=sorted(set(counts)), unit="-", source_script=REGEN,
            source_file=here,
            note=f"{counts.count(64)} cells at 64 and {counts.count(216)} at 216"),
        "md_timestep_fs": dict(
            value=dt if len(steps) == 1 else sorted(steps), unit="fs",
            source_script=REGEN, source_file=here, note=None),
        "equilibration_ps": dict(
            value=10.0, unit="ps",
            source_script="ed_data/arc_a78/scw/run_prod_md_05fs.py",
            source_file=None,
            note=("Langevin at 0.5 fs for 20000 steps (line 48), then "
                  "VelocityVerlet at constant energy (line 50)")),
        "n_training_configurations": dict(
            value=1193, unit="-", source_script=METHODS, source_file=None,
            note=("train_v4.xyz: 800 bulk pools 650-2250 K, 344 RCB fill "
                  "2500-3000 K, 49 above 3000 K")),
        "training_temperature_range_k": dict(
            value=[650, 3000], unit="K", source_script=METHODS, source_file=None,
            note=None),
        "label_level_of_theory": dict(
            value="revPBE0-D3 (CP2K), ADMM, DZVP-MOLOPT-SR-GTH, 600 Ry",
            unit="-", source_script=METHODS, source_file=None, note=None),
    }
    if frames:
        facts["production_ps"] = dict(
            value=round(max(frames) * dt / 1000.0), unit="ps",
            source_script=REGEN, source_file="ed_data/arc_dipoles/*/dipoles.npy",
            note=(f"from the dipole series length ({max(frames)} frames); "
                  f"METHODS.md says 60 ps, which the data do not support"))

    ref = json.loads((ARC / "water-scw-mace" / "benchmarks" / "reference"
                      / "reference.json").read_text())
    facts["force_parity_mev_per_a"] = dict(
        value=ref["parity_meV_A"], unit="meV/A", source_script=PARITY_REF,
        source_file=None,
        note=("TRAINING-SET parity by regime, which reference.json itself calls "
              "'optimistic vs true held-out error'; frames above 3000 K excluded. "
              "The only held-out parity in the pull is 12.3 / 20.2 meV/A and "
              "belongs to the LDA model, not this one"))
    facts["nve_drift_mev_per_atom_per_ps"] = dict(
        value=NVE_DRIFT_MEV_ATOM_PS, unit="meV/atom/ps", source_script=CAMPAIGN,
        source_file=None,
        note=(f"worst over 21 NVE runs against a {NVE_DRIFT_GATE:g} gate; a bound, "
              f"not a distribution"))
    return facts


def stage_captions(facts: dict[str, dict]) -> list[str]:
    """The four validation captions Fig. 1a prints, in stage order."""
    if not facts:
        return [""] * 4
    v = {k: f["value"] for k, f in facts.items()}
    t_lo, t_hi = v["training_temperature_range_k"]
    par = v["force_parity_mev_per_a"]
    prod = v.get("production_ps")
    mol = v["molecules_per_cell"]
    return [
        f"{v['n_state_points']} states, {t_lo:g}-{t_hi:g} K",
        f"{v['n_training_configurations']} configurations\n"
        f"{str(v['label_level_of_theory']).split(' (')[0]}",
        f"force parity {par['bulk']:g}/{par['rcb']:g} meV $\\mathrm{{\\AA}}^{{-1}}$\n"
        f"(bulk/RCB, training set)",
        (f"{min(mol):g}-{max(mol):g} molecules, {v['md_timestep_fs']:g} fs\n"
         f"{prod:g} ps after {v['equilibration_ps']:g} ps; drift "
         f"{v['nve_drift_mev_per_atom_per_ps']:g}"
         if prod else f"{min(mol):g}-{max(mol):g} molecules, {v['md_timestep_fs']:g} fs"),
    ]
