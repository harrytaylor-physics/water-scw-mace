#!/usr/bin/env python3
"""Self-healing ExoMol runner. Wraps radis_kappa's per-point computation in
a retry loop that repairs the two failure modes we hit:
  1. corrupted .trans.bz2 chunk -> delete it (+ cached .h5 siblings), retry
  2. RADIS refuses to re-download a missing chunk -> fetch it directly from
     exomol.com and retry
Run on the login node in tmux (downloads), then the batch job for compute.

    python3 exomol_selfheal.py --grid grid.json --database BT2
"""
import argparse
import json
import pathlib
import re
import sys
import time
import urllib.request

EXOMOL_BASE = "https://exomol.com/db/H2O/1H2-16O/{db}/{fname}"
MAX_RETRIES_PER_FILE = 3
MAX_TOTAL_REPAIRS = 40

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from radis_kappa import steam_pressure_bar, M_H2O, NA, KB  # noqa: E402
import numpy as np  # noqa: E402


def repair(fname_in_err, db, seen):
    p = pathlib.Path(fname_in_err)
    seen[p.name] = seen.get(p.name, 0) + 1
    if seen[p.name] > MAX_RETRIES_PER_FILE:
        raise SystemExit(f"{p.name}: failed {MAX_RETRIES_PER_FILE} repairs - "
                         f"chunk likely absent on exomol.com; escalate "
                         f"(manual listing check / xsec route).")
    # delete the file and any cached derivatives sharing its stem
    stem = p.name.split(".")[0]
    for sib in p.parent.glob(stem + "*"):
        print(f"  repair: deleting {sib.name}")
        sib.unlink()
    # direct fetch (covers the radis-won't-redownload case)
    url = EXOMOL_BASE.format(db=db, fname=p.name)
    print(f"  repair: fetching {url}")
    try:
        urllib.request.urlretrieve(url, p)
        print(f"  repair: fetched {p.stat().st_size/1e6:.1f} MB")
    except Exception as e:
        print(f"  repair: direct fetch failed ({e}) - leaving absent; "
              f"radis will be given the chance to fetch it")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--grid", required=True)
    ap.add_argument("--database", default="BT2")
    ap.add_argument("--numin", type=float, default=50.0)
    ap.add_argument("--numax", type=float, default=4200.0)
    ap.add_argument("--wstep", type=float, default=0.05)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    from radis import calc_spectrum
    out = pathlib.Path(args.out or "ll_exomol")
    out.mkdir(exist_ok=True)
    seen, repairs = {}, 0
    for pt in json.load(open(args.grid)):
        T, rho = pt["T"], pt["rho"]
        f = out / f"T{T}_rho{rho}.npz"
        if f.exists():
            print(f"skip (done): {f.name}")
            continue
        p_bar, psrc = steam_pressure_bar(T, rho)
        while True:
            try:
                s = calc_spectrum(args.numin, args.numax,
                                  molecule="H2O", isotope="1",
                                  pressure=p_bar, Tgas=T, mole_fraction=1.0,
                                  path_length=1.0,
                                  databank=("exomol", args.database),
                                  wstep=args.wstep, verbose=False)
                break
            except Exception as e:
                m = re.search(r"(/\S+\.trans\.bz2)", str(e))
                if m and repairs < MAX_TOTAL_REPAIRS:
                    repairs += 1
                    print(f"[{repairs}] chunk problem: {m.group(1)}")
                    repair(m.group(1), args.database, seen)
                    time.sleep(2)
                    continue
                raise
        wn, absc = s.get("abscoeff", wunit="cm-1")
        n_ideal_cm3 = (p_bar * 1e5) / (KB * T) * 1e-6
        kappa = absc / n_ideal_cm3 * NA / M_H2O
        np.savez(f, wn=wn, kappa=kappa, T=T, rho=rho,
                 p_bar=p_bar, pressure_source=psrc,
                 source=f"RADIS/exomol/{args.database} (self-healed run)")
        print(f"T{T} rho{rho}: p={p_bar:.1f} bar ({psrc}) -> {f.name}")


if __name__ == "__main__":
    main()
