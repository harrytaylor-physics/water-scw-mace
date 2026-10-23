#!/usr/bin/env python3
"""Model vs DFT dipole and dipole derivatives for the isolated molecule (A64).

    python compare.py

Reads ``model_dipoles.json`` and ``dft_dipoles.json`` and answers one question:
is ``water_dipole_v3``'s stretch dipole derivative right for a *gas-phase*
molecule? A61 found the first-principles spectrum carrying 2.4-3.9x the line
list's stretch-band intensity and 0.35-0.50x its rotational-band intensity, with
a stretch/bend balance 3.3-3.6x the line list's, and the same at 1000, 1500 and
2000 K. Band intensity goes as |dmu/dQ|^2, so a stretch derivative too large by
~1.8-2x would produce the observed excess on its own.

What is printed:

* **mu at equilibrium** — model, DFT, and experiment (1.855 D). The permanent
  dipole is the easy check; a model can get it right and still have the wrong
  derivative, so agreement here is necessary, not sufficient.
* **dmu/dQ by central differences** for symmetric stretch, asymmetric stretch and
  bend, at both references, as vectors and magnitudes, with |model|/|DFT|.
* **the squared-magnitude ratio (stretch/bend)** for each source — the
  intensity-ratio proxy the dilute-limit result is about. Bend derivatives are
  converted from D/deg to D/A using the bond length (r*dtheta is the displacement
  a stretch of dr is comparable to), so the ratio is dimensionless and means what
  it looks like.

Reading the answer: model/DFT stretch-derivative ratio near 1 means the dipole
model is not the cause of the stretch excess; 1.5-2 or more means it is.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
MU_EXPERIMENT_D = 1.855
MODES = ("sym_stretch", "asym_stretch", "bend")
REFERENCES = ("equilibrium", "stretched")
VERDICT_NOT_THE_CAUSE = 1.25
VERDICT_IS_THE_CAUSE = 1.5


def _index(results: list[dict]) -> dict:
    return {(r["reference"], r["mode"], r["sign"]): r for r in results}


def _norm(v) -> float:
    return math.sqrt(sum(x * x for x in v))


def derivative(idx: dict, reference: str, mode: str) -> dict | None:
    """Central difference (mu(+) - mu(-)) / (2*delta), in D/A or D/deg."""
    plus, minus = idx.get((reference, mode, +1)), idx.get((reference, mode, -1))
    if plus is None or minus is None:
        return None
    delta = abs(plus["displacement"])
    vec = [(a - b) / (2.0 * delta) for a, b in zip(plus["mu_debye"], minus["mu_debye"])]
    return {"vector": vec, "magnitude": _norm(vec),
            "unit": f"D/{plus['displacement_unit']}", "delta": delta}


def intensity_proxy(derivs: dict, r0_a: float) -> dict:
    """(|dmu/dQ_stretch| / |dmu/dQ_bend|)^2, both in D/A, so dimensionless.

    The bend derivative arrives in D/deg. The displacement comparable to stretching
    a bond by dr is the arc r*dtheta, so D/deg is converted with
    (180/pi)/r_0 to D/A-equivalent before the ratio is formed.
    """
    bend = derivs.get("bend")
    if bend is None:
        return {}
    bend_per_a = bend["magnitude"] * (180.0 / math.pi) / r0_a
    out = {"bend_magnitude_d_per_a_equivalent": bend_per_a}
    total_sq = 0.0
    for mode in ("sym_stretch", "asym_stretch"):
        d = derivs.get(mode)
        if d is None:
            continue
        out[f"{mode}_over_bend_squared"] = (d["magnitude"] / bend_per_a) ** 2
        total_sq += d["magnitude"] ** 2
    if total_sq:
        out["total_stretch_over_bend_squared"] = total_sq / bend_per_a**2
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--model-json", type=Path, default=HERE / "model_dipoles.json")
    ap.add_argument("--dft-json", type=Path, default=HERE / "dft_dipoles.json")
    ap.add_argument("--geoms", type=Path, default=HERE / "geoms.json")
    ap.add_argument("--out", type=Path, default=HERE / "compare.json")
    args = ap.parse_args(argv)

    missing = [p for p in (args.model_json, args.dft_json) if not p.exists()]
    if missing:
        raise SystemExit(
            "missing " + ", ".join(str(p) for p in missing) +
            "\nRun eval_model.py (compute node) and harvest_dft.py (after the CP2K "
            "array) first; see README.md.")
    model = json.loads(args.model_json.read_text())["results"]
    dft = json.loads(args.dft_json.read_text())["results"]
    geoms = json.loads(args.geoms.read_text())
    r0 = {"equilibrium": geoms["r_eq_a"], "stretched": geoms["r_stretched_a"]}
    mi, di = _index(model), _index(dft)

    print("\n  TASK 26 — isolated water monomer: model vs revPBE0/MLWF\n")
    print("  PERMANENT DIPOLE at equilibrium (r_OH = "
          f"{geoms['r_eq_a']} A, HOH = {geoms['theta_eq_deg']} deg)")
    eq = {}
    for label, idx in (("model", mi), ("DFT", di)):
        r = idx.get(("equilibrium", "reference", 0))
        if r is None:
            continue
        eq[label] = r["mu_norm_debye"]
        print(f"    {label:<11s} |mu| = {r['mu_norm_debye']:.4f} D   "
              f"({r['mu_debye'][0]:+.4f}, {r['mu_debye'][1]:+.4f}, "
              f"{r['mu_debye'][2]:+.4f})")
    print(f"    {'experiment':<11s} |mu| = {MU_EXPERIMENT_D:.4f} D")
    for label, val in eq.items():
        print(f"    {label} / experiment = {val / MU_EXPERIMENT_D:.4f}")

    derivs = {"model": {}, "DFT": {}}
    print(f"\n  DIPOLE DERIVATIVES, central differences")
    head = (f"    {'reference':<12s} {'mode':<13s} {'source':<6s} | "
            f"{'d(mu)/dQ vector':>34s} | {'|d(mu)/dQ|':>11s} {'unit':<7s} | "
            f"{'model/DFT':>9s}")
    print(head)
    print("    " + "-" * (len(head) - 4))
    ratios = {}
    for reference in REFERENCES:
        for mode in MODES:
            dm = derivative(mi, reference, mode)
            dd = derivative(di, reference, mode)
            if dm:
                derivs["model"].setdefault(reference, {})[mode] = dm
            if dd:
                derivs["DFT"].setdefault(reference, {})[mode] = dd
            ratio = (dm["magnitude"] / dd["magnitude"]) if dm and dd and dd["magnitude"] else None
            if ratio is not None:
                ratios[(reference, mode)] = ratio
            for label, d in (("model", dm), ("DFT", dd)):
                if d is None:
                    continue
                v = d["vector"]
                shown = (f"{ratio:>9.3f}" if label == "DFT" and ratio is not None
                         else " " * 9)
                print(f"    {reference:<12s} {mode:<13s} {label:<6s} | "
                      f"({v[0]:+10.5f},{v[1]:+10.5f},{v[2]:+10.5f}) | "
                      f"{d['magnitude']:>11.5f} {d['unit']:<7s} | {shown}")

    print(f"\n  INTENSITY PROXY  (|d(mu)/dQ_stretch| / |d(mu)/dQ_bend|)^2, "
          f"both as D/A")
    proxies = {}
    head = (f"    {'reference':<12s} {'source':<6s} | {'sym/bend':>10s} "
            f"{'asym/bend':>10s} {'total/bend':>11s}")
    print(head)
    print("    " + "-" * (len(head) - 4))
    for reference in REFERENCES:
        for label in ("model", "DFT"):
            d = derivs[label].get(reference)
            if not d:
                continue
            p = intensity_proxy(d, r0[reference])
            proxies.setdefault(reference, {})[label] = p
            if not p:
                continue
            print(f"    {reference:<12s} {label:<6s} | "
                  f"{p.get('sym_stretch_over_bend_squared', float('nan')):>10.3f} "
                  f"{p.get('asym_stretch_over_bend_squared', float('nan')):>10.3f} "
                  f"{p.get('total_stretch_over_bend_squared', float('nan')):>11.3f}")
        both = proxies.get(reference, {})
        if len(both) == 2 and both["DFT"].get("total_stretch_over_bend_squared"):
            r = (both["model"]["total_stretch_over_bend_squared"]
                 / both["DFT"]["total_stretch_over_bend_squared"])
            print(f"    {reference:<12s} {'ratio':<6s} | "
                  f"{'':>10s} {'':>10s} {r:>11.3f}   <- model/DFT")

    stretch_ratios = [v for (ref, mode), v in ratios.items() if "stretch" in mode]
    verdict = "no stretch derivatives available on both sides"
    if stretch_ratios:
        worst = max(stretch_ratios)
        if worst < VERDICT_NOT_THE_CAUSE:
            verdict = (
                f"model/DFT stretch-derivative ratio is {min(stretch_ratios):.2f}-"
                f"{worst:.2f}, near 1 -> THE DIPOLE MODEL IS NOT THE CAUSE of the "
                f"stretch excess; look elsewhere (the dynamics, or the line list).")
        elif worst >= VERDICT_IS_THE_CAUSE:
            verdict = (
                f"model/DFT stretch-derivative ratio is {min(stretch_ratios):.2f}-"
                f"{worst:.2f} -> THE DIPOLE MODEL IS THE CAUSE: it returns a "
                f"liquid-like dmu/dr for a gas-phase stretch, and band intensity "
                f"goes as the square of it.")
        else:
            verdict = (
                f"model/DFT stretch-derivative ratio is {min(stretch_ratios):.2f}-"
                f"{worst:.2f}, between the two thresholds -> partly implicated; "
                f"the squared effect is "
                f"{min(stretch_ratios) ** 2:.2f}-{worst ** 2:.2f}x in intensity, "
                f"against the 2.4-3.9x the spectrum shows, so it cannot be the "
                f"whole story either way.")
    print(f"\n  READING: {verdict}")

    args.out.write_text(json.dumps({
        "meta": ("Isolated-monomer dipole test: water_dipole_v3 against "
                 "revPBE0/MLWF at identical geometries (A64). Derivatives are "
                 "central differences; the intensity proxy converts the bend "
                 "derivative to D/A with the arc length r*dtheta."),
        "mu_experiment_debye": MU_EXPERIMENT_D,
        "equilibrium_mu_norm_debye": eq,
        "derivatives": derivs,
        "model_over_dft_magnitude": {f"{k[0]}/{k[1]}": v for k, v in ratios.items()},
        "intensity_proxy": proxies,
        "verdict": verdict,
    }, indent=1) + "\n")
    print(f"\n  wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
