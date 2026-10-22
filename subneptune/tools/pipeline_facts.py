#!/usr/bin/env python3
"""The pipeline facts Fig. 1a states and the ledger records as ``pipeline.*`` (A89).

One function, two callers. Figure 1a annotates each stage of the pipeline with the
number that validates it, and ``tools/manuscript_numbers.py`` records the same
numbers in the ledger; if each dug them out separately they could drift apart, and
a schematic quoting a number the ledger does not carry is a number no reader can
check. So both call this.

Each fact names the artefact that actually carries it rather than prose about it.
That matters here: ``METHODS.md`` disagreed with the delivered data on four
counts -- 60 ps per state point, a 27-point grid, NVT production and 64 molecules,
where the data give 30 ps, 52 cells, Langevin-then-NVE and 64 or 216. A90
corrected the file against the artefacts and preserved what it said as delivered.
Where only prose exists (the training-set composition, the level of theory) the
fact is marked with that provenance.

There are **three** different error numbers here and they are not interchangeable
(A89 found two of them; A90 added the third):

  * the **held-out** split error of the production potential, 57.8 meV/A and
    20.8 meV/atom on ``valid_Default`` (56.9 / 23.5 on ``train_Default``). This is
    the honest generalisation error and it comes from the training job's own log,
    ``production_campaign/arc_scripts/logs/train_v4_7924966.out``, via
    ``models/mlip_heldout_errors.json``. A94: that log holds TWO error tables.
    A90 transcribed the stage-two (SWA, epoch 186) table, 44.4 / 7.0, but
    production MD ran the stage-one model (epoch 109, sha256 dc4a8977...), so
    the stage-one table is the one quoted. The stage-two numbers stay in the JSON,
    labelled as not used for production, and are not read here.
  * the **training-set parity by regime**, 41 meV/A bulk and 58 meV/A at RCB
    conditions, which ``reference.json`` itself calls "optimistic vs true held-out
    error". A different quantity, not a competing estimate; the 58 at RCB remains
    the number for the regime the structural claim lives in. The only held-out
    parity inside the pull is 12.3 / 20.2 meV/A and belongs to the LDA model.
  * the NVE check, which is a **drift rate** bound -- 0.017 meV/atom/ps worst over
    21 runs against a 2 meV/atom/ps gate -- not an energy error per atom.
"""

from __future__ import annotations

import json
from pathlib import Path

import layout

ROOT = Path(__file__).resolve().parents[1]
# Paths naming benchmarks/, docs/, models/ and production_campaign/ are relative to
# the water-scw-mace root (layout.py); the rest are relative to subneptune/.
PARITY_REF = "benchmarks/reference/reference.json"
METHODS = "docs/METHODS.md"
CAMPAIGN = "docs/LDA_CAMPAIGN.md"
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
            source_script="production_campaign/scw/run_prod_md_05fs.py",
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
                  f"METHODS.md said 60 ps until A90 corrected it against this"))

    ref = json.loads(layout.PARITY_REF.read_text())
    facts["force_parity_mev_per_a"] = dict(
        value=ref["parity_meV_A"], unit="meV/A", source_script=PARITY_REF,
        source_file=None,
        note=("TRAINING-SET parity by regime, which reference.json itself calls "
              "'optimistic vs true held-out error'; frames above 3000 K excluded. "
              "Computed with the production (stage-one) model: A94 reproduced "
              "41.0 / 58.3 from it. The only held-out parity in the pull is "
              "12.3 / 20.2 meV/A and belongs to the LDA model, not this one"))
    held = json.loads(layout.MLIP_ERRORS.read_text())
    prod = held["production_model"]
    _hnote = (f"production potential ({prod['training_stage']}, "
              f"{prod['checkpoint'].rsplit('/', 1)[-1]}, sha256 "
              f"{prod['sha256'][:8]}...); from {held['source_log']} lines "
              f"{prod['log_lines']}, recorded in models/mlip_heldout_errors.json. "
              f"Split error -- a different quantity from the training-set parity "
              f"by regime. The log's stage-two (SWA) table is not this model's (A94)")
    for tag, field, key, unit in (
            ("valid_Default", "force_meV_per_A", "mlip_heldout_force_mev_per_a", "meV/A"),
            ("valid_Default", "energy_meV_per_atom", "mlip_heldout_energy_mev_per_atom", "meV/atom"),
            ("train_Default", "force_meV_per_A", "mlip_train_force_mev_per_a", "meV/A"),
            ("train_Default", "energy_meV_per_atom", "mlip_train_energy_mev_per_atom", "meV/atom")):
        facts[key] = dict(
            value=held["errors"][tag][field], unit=unit,
            source_script="models/mlip_heldout_errors.json", source_file=None,
            note=f"{tag}; {_hnote}")
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
        f"held out {v['mlip_heldout_force_mev_per_a']:g} meV "
        f"$\\mathrm{{\\AA}}^{{-1}}$, {v['mlip_heldout_energy_mev_per_atom']:g} meV/atom\n"
        f"parity {par['bulk']:g}/{par['rcb']:g} (bulk/RCB, train)",
        # A91: the NVE drift bound stays in the ledger as
        # pipeline.nve_drift_mev_per_atom_per_ps, but off the figure -- it is a
        # stability check on the integrator, not a validation of the opacity, and
        # it was the one caption number a reader could not place.
        (f"{min(mol):g}-{max(mol):g} molecules, {v['md_timestep_fs']:g} fs\n"
         f"{prod:g} ps after {v['equilibration_ps']:g} ps equilibration"
         if prod else f"{min(mol):g}-{max(mol):g} molecules, {v['md_timestep_fs']:g} fs"),
    ]
