#!/usr/bin/env python3
"""Compare two manuscript ledgers, group by group (A89).

    python tools/column_diff.py outputs/manuscript_numbers_v3_rebuilt.json \\
                                outputs/manuscript_numbers.json

Both files must be built by running every producer against the spectra of their
own version -- that is what ``tools/rebuild_v3_column.sh`` and
``tools/rebuild_producers.sh`` are for. Comparing a rebuilt column against a
stored one is how A83 reported a stale CSV as a baseline and how A86 reported
"unchanged" for figure groups that were v5 on both sides.

A group is the part of a key before the first dot. For each group this prints how
many keys the two columns share, how many of those are numeric, and how many of
the numeric ones moved by more than ``--tol`` in relative terms. Keys present in
only one column are listed separately rather than counted as unchanged: an empty
group and an unchanged group are different findings and must not look alike.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


def _value(entry):
    return entry.get("value") if isinstance(entry, dict) else entry


def _num(v):
    return isinstance(v, (int, float)) and not isinstance(v, bool) and math.isfinite(v)


def compare(a: dict, b: dict, tol: float) -> dict:
    keys_a, keys_b = set(a), set(b)
    common = sorted(keys_a & keys_b)
    groups: dict[str, dict] = {}
    for k in common:
        g = k.split(".")[0]
        rec = groups.setdefault(g, {"n_common": 0, "n_numeric": 0, "changed": [],
                                    "non_numeric_changed": []})
        rec["n_common"] += 1
        va, vb = _value(a[k]), _value(b[k])
        if _num(va) and _num(vb):
            rec["n_numeric"] += 1
            denom = abs(va) if va else 1.0
            rel = (vb - va) / denom
            if abs(rel) > tol:
                rec["changed"].append({"key": k, "v3": va, "v5": vb, "rel": rel})
        elif va != vb:
            rec["non_numeric_changed"].append({"key": k, "v3": va, "v5": vb})
    only_a, only_b = {}, {}
    for k in sorted(keys_a - keys_b):
        only_a.setdefault(k.split(".")[0], []).append(k)
    for k in sorted(keys_b - keys_a):
        only_b.setdefault(k.split(".")[0], []).append(k)
    return {"groups": groups, "only_in_first": only_a, "only_in_second": only_b,
            "n_first": len(keys_a), "n_second": len(keys_b), "n_common": len(common)}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("first", type=Path, help="the v3 column")
    ap.add_argument("second", type=Path, help="the v5 column")
    ap.add_argument("--tol", type=float, default=1e-6,
                    help="relative change counted as a difference (default 1e-6)")
    ap.add_argument("--show", type=int, default=6,
                    help="changed keys to print per group (default 6)")
    ap.add_argument("--json", type=Path, default=None)
    args = ap.parse_args(argv)

    a = json.loads(args.first.read_text())
    b = json.loads(args.second.read_text())
    res = compare(a, b, args.tol)

    print(f"\n  {args.first.name}: {res['n_first']} keys")
    print(f"  {args.second.name}: {res['n_second']} keys")
    print(f"  shared: {res['n_common']}\n")
    tot_num = tot_chg = 0
    print(f"  {'group':<24} {'shared':>7} {'numeric':>8} {'changed':>8}   state")
    for g in sorted(res["groups"]):
        r = res["groups"][g]
        tot_num += r["n_numeric"]; tot_chg += len(r["changed"])
        if r["n_numeric"] == 0:
            state = "no numeric keys to compare"
        elif r["changed"]:
            state = f"CHANGED ({len(r['changed'])}/{r['n_numeric']})"
        else:
            state = "unchanged"
        if r["non_numeric_changed"]:
            state += f", {len(r['non_numeric_changed'])} non-numeric differ"
        print(f"  {g:<24} {r['n_common']:>7} {r['n_numeric']:>8} "
              f"{len(r['changed']):>8}   {state}")
    print(f"\n  numeric keys compared: {tot_num}; changed: {tot_chg} "
          f"({100.0*tot_chg/tot_num:.1f}%)" if tot_num else "\n  nothing numeric to compare")

    for g in sorted(res["groups"]):
        r = res["groups"][g]
        if not r["changed"]:
            continue
        print(f"\n  --- {g} ---")
        for c in sorted(r["changed"], key=lambda c: -abs(c["rel"]))[:args.show]:
            print(f"    {c['key']:<52} {c['v3']:>14.6g} -> {c['v5']:>14.6g}"
                  f"  {100*c['rel']:+8.2f}%")
        if len(r["changed"]) > args.show:
            print(f"    ... and {len(r['changed']) - args.show} more")

    for name, d in (("only in " + args.first.name, res["only_in_first"]),
                    ("only in " + args.second.name, res["only_in_second"])):
        if d:
            print(f"\n  {name}:")
            for g, ks in sorted(d.items()):
                print(f"    {g:<24} {len(ks):>4} keys   e.g. {ks[0]}")
    if args.json:
        args.json.write_text(json.dumps(res, indent=1) + "\n")
        print(f"\n  wrote {args.json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
