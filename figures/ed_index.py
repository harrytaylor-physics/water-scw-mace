#!/usr/bin/env python3
"""Write figs/ed/ed_index.json and fold the ED numbers into the manuscript ledger.

    python figures/ed_index.py config/fiducial.yaml

``ed_index.json`` maps each Extended Data label to a one-line description, so the
text can cite items by number without opening the figures. It also records, per
item, whether it is built, partly built (some panel has no data), or blocked —
because a citation to a panel that is a placeholder needs to say so.

The same numbers are then merged into ``outputs/manuscript_numbers.json`` under
``ed.*`` keys, so the single source of truth for the read-through covers the
Extended Data as well as the main figures.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from config import load_config
from figures._manuscript import ED_DIR

ROOT = Path(__file__).resolve().parents[1]

DESCRIPTIONS = {
    "ed_fig1": "Validation of the interatomic potential: force-parity RMSE by "
               "regime, and the NVE energy-drift bound against its gate.",
    "ed_fig2": "Planck-mean opacity over the 48 computed points, and the "
               "first-principles / HITEMP Planck ratio on the same axes.",
    "ed_fig3": "Line-list conventions: kappa_R under three broadening pressures, "
               "the MT_CKD self-continuum share of kappa_R, and the band-truncation "
               "ratio over the line-list grid.",
    "ed_fig4": "Dilute-limit band structure: the band-integrated absorption ratio "
               "at nine states, and kappa(nu) for both sources at two of them.",
    "ed_fig5": "Dipole-model validation: dipole derivatives against revPBE0 at two "
               "bond lengths, the model/DFT ratio per mode, and the dipole-model "
               "residual autocorrelation.",
    "ed_fig6": "The 108 scanned planets: boundary positions against the computed "
               "grid by category, and the dispute ratio against the factor that "
               "would flip each analysed boundary.",
    "ed_fig7": "External anchor: model radius against the published relation at the "
               "matched points, and the residual against water mass fraction.",
    "ed_fig8": "Evolution: radius, interior entropy and deepest-boundary density "
               "against age for three water worlds under both opacities.",
    "ed_fig9_scaling": "Density scaling: alpha/rho^2 against density at three "
                       "places in the spectrum, the local and fitted power-law "
                       "index by source, and the mean-free-path breakdown map "
                       "with the published tables' documented limits.",
    "ed_table1": "The 48 state points: first-principles Planck and Rosseland means "
                 "and the ratios to each line-list reference.",
    "ed_table2": "The ten strongest boundaries: state, both opacities, dispute "
                 "ratio, flip threshold, shift and the provenance of the "
                 "first-principles value.",
}


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


def _status(entry: dict) -> tuple[str, list[str]]:
    panels = entry.get("panels")
    if panels is None:
        return "built", []
    blocked = [k or "single" for k, v in panels.items() if v.get("blocked")]
    if not blocked:
        return "built", []
    if len(blocked) == len(panels):
        return "blocked", blocked
    return "partly built", blocked


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("config", type=Path)
    args = ap.parse_args(argv)
    cfg = load_config(args.config)

    ledger_path = ED_DIR / "ed_numbers.json"
    if not ledger_path.exists():
        raise SystemExit(f"{ledger_path} not found — build the ED items first")
    ledger = json.loads(ledger_path.read_text())

    index = {"_meta": ("Extended Data items: label -> one-line content, status, and "
                       "the file that holds it. 'partly built' means a panel is a "
                       "labelled placeholder because its data are not in the pulled "
                       "set; a citation to that panel has to say so.")}
    for label, desc in DESCRIPTIONS.items():
        entry = ledger.get(label)
        is_table = label.startswith("ed_table")
        if entry is None:
            status, blocked, files = "not built", [], []
        elif is_table:
            status, blocked, files = "built", [], entry.get("files", [])
        else:
            status, blocked = _status(entry)
            files = [f"figs/ed/{label}.pdf"]
        index[label] = {"description": desc, "status": status, "files": files}
        if blocked:
            index[label]["placeholder_panels"] = blocked
            index[label]["placeholder_needs"] = {
                (k or "single"): entry["panels"][k].get("needs")
                for k in entry["panels"] if entry["panels"][k].get("blocked")}
    (ED_DIR / "ed_index.json").write_text(json.dumps(index, indent=1) + "\n")

    # fold into the manuscript ledger
    mpath = cfg.output_path / "manuscript_numbers.json"
    doc = json.loads(mpath.read_text()) if mpath.exists() else {}
    added = 0
    for label, entry in ledger.items():
        # the ledger also carries scalar header keys (width_mm, min_font_pt)
        if label.startswith("_") or not isinstance(entry, dict):
            continue
        script = ("figures/ed_tables.py" if label.startswith("ed_table")
                  else "figures/ed_items.py")
        commit = _commit(script)
        if label.startswith("ed_table") or label == "functional_sensitivity":
            doc[f"ed.{label}"] = {
                "value": entry, "unit": "-", "source_script": script,
                "source_file": f"figs/ed/{'ed_numbers.json'}", "commit": commit,
                "available": True}
            added += 1
            continue
        for panel, values in entry.get("panels", {}).items():
            for key, val in values.items():
                if key in ("source_files", "source_commits"):
                    continue
                doc[f"ed.{label}.{panel or 'single'}.{key}"] = {
                    "value": val, "unit": "-", "source_script": script,
                    "source_file": "figs/ed/ed_numbers.json", "commit": commit,
                    "available": not values.get("blocked", False)}
                added += 1
    doc["_n_entries"] = sum(1 for k in doc if not k.startswith("_"))
    doc["_n_unavailable"] = sum(1 for k, v in doc.items()
                                if not k.startswith("_") and not v.get("available"))
    mpath.write_text(json.dumps(doc, indent=1) + "\n")

    print(f"\n  ED INDEX — {len([k for k in index if not k.startswith('_')])} items")
    for label, e in index.items():
        if label.startswith("_"):
            continue
        extra = (f"  [placeholder: {', '.join(e['placeholder_panels'])}]"
                 if "placeholder_panels" in e else "")
        print(f"    {label:<12s} {e['status']:<13s} {e['files'][0] if e['files'] else '-'}"
              f"{extra}")
    print(f"\n  wrote {ED_DIR / 'ed_index.json'}")
    print(f"  merged {added} ED entries into {mpath} "
          f"({doc['_n_entries']} total, {doc['_n_unavailable']} unavailable)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
