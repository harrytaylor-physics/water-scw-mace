#!/usr/bin/env python3
"""Stream-subsample frames from production trajectories into snippet files
(memory-safe: never loads a 400 MB trajectory at once).

  python3 make_snippets.py --production /data/.../lda_pipeline/production \
      --points T1000_rho129.744 T2000_rho500 T3000_rho810.382 \
      --n-frames 80 --out snippets
"""
import argparse
import pathlib


def stream_sample(traj, n_frames, discard_frac=0.05, offset_frac=0.0):
    with open(traj) as f:
        first = f.readline()
        nat = int(first.split()[0])
        block = nat + 2
        f.seek(0)
        nlines = sum(1 for _ in f)
    total = nlines // block
    start = int(total * discard_frac)
    step = (total - start - 1) / max(n_frames - 1, 1)
    idx = set(start + int((i + offset_frac) * step) for i in range(n_frames))
    out, keep, buf, fi = [], False, [], -1
    with open(traj) as f:
        for ln, line in enumerate(f):
            if ln % block == 0:
                fi += 1
                keep = fi in idx
                buf = []
            if keep:
                buf.append(line)
                if len(buf) == block:
                    out.extend(buf)
    return "".join(out), total


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--production", required=True)
    ap.add_argument("--points", nargs="+", required=True)
    ap.add_argument("--n-frames", type=int, default=80)
    ap.add_argument("--out", default="snippets")
    ap.add_argument("--offset-frac", type=float, default=0.0,
                    help="0.5 = sample midway between an existing 0.0 grid")
    ap.add_argument("--tag", default="",
                    help="suffix for snippet stem, e.g. b -> snip_..._b.xyz")
    args = ap.parse_args()
    outdir = pathlib.Path(args.out)
    outdir.mkdir(exist_ok=True)
    for pt in args.points:
        traj = pathlib.Path(args.production) / pt / "prod.xyz"
        if not traj.exists():
            print(f"{pt}: prod.xyz NOT FOUND - skipping (production not done?)")
            continue
        text, total = stream_sample(traj, args.n_frames,
                                    offset_frac=args.offset_frac)
        suffix = f"_{args.tag}" if args.tag else ""
        f = outdir / f"snip_{pt}{suffix}.xyz"
        f.write_text(text)
        print(f"{pt}: {args.n_frames} frames sampled from {total} -> {f}")


if __name__ == "__main__":
    main()
