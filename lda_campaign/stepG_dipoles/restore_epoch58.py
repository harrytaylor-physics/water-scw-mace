#!/usr/bin/env python3
"""Rebuild a loadable .model from the surviving epoch-58 checkpoint by
loading the current (same-architecture) .model and swapping in the old
weights. Run on a GPU node.

  python3 restore_epoch58.py \
      --current /data/.../lda_pipeline/checkpoints/water_dipole_lda_run-1.model \
      --ckpt    /data/.../lda_pipeline/checkpoints/water_dipole_lda_run-1_epoch-58.pt \
      --out     /data/.../lda_pipeline/water_dipole_lda.model
"""
import argparse
import torch


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--current", required=True)
    ap.add_argument("--ckpt", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    model = torch.load(args.current, map_location="cpu")
    ckpt = torch.load(args.ckpt, map_location="cpu")
    state = None
    if isinstance(ckpt, dict):
        for key in ("model", "state_dict", "model_state_dict"):
            if key in ckpt:
                state = ckpt[key]
                break
        if state is None and all(hasattr(v, "shape") or isinstance(v, dict)
                                 for v in list(ckpt.values())[:3]):
            state = ckpt
    if state is None:
        raise SystemExit(f"unrecognised checkpoint layout: {list(ckpt)[:8]}")
    missing, unexpected = model.load_state_dict(state, strict=False)
    print(f"loaded epoch-58 weights: missing={len(missing)} "
          f"unexpected={len(unexpected)}")
    if len(missing) > 5:
        raise SystemExit(f"too many missing keys - layout mismatch: "
                         f"{missing[:5]}")
    torch.save(model, args.out)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
