#!/usr/bin/env python3
"""Build the LDA MLWF template: keep the hybrid mlwf.inp's localization/print
machinery, transplant the &XC from a converged step-B LDA label.inp, strip
ADMM/HF wiring, swap pseudos to PADE.

  python3 make_lda_mlwf.py --mlwf ~/scw/mlwf_test/mlwf.inp \
      --lda-label /data/.../lda_pipeline/labels/T1000_rho129.744/sel_0000/label.inp \
      --out mlwf_lda.inp
"""
import argparse
import pathlib
import re


def block(text, start, end=None):
    end = end or start
    m = re.search(rf"([ \t]*&{start}\b.*?^[ \t]*&END {end}\b[^\n]*)",
                  text, re.S | re.M)
    return m.group(1) if m else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mlwf", required=True)
    ap.add_argument("--lda-label", required=True)
    ap.add_argument("--out", default="mlwf_lda.inp")
    args = ap.parse_args()
    mlwf = pathlib.Path(args.mlwf).read_text()
    lda = pathlib.Path(args.lda_label).read_text()

    for tag in ("LOCALIZE", "WANNIER"):
        if tag not in mlwf.upper():
            print(f"WARNING: '{tag}' not found in mlwf.inp - check this is "
                  f"the MLWF template")
    lda_xc = block(lda, "XC")
    old_xc = block(mlwf, "XC")
    if lda_xc is None or old_xc is None:
        raise SystemExit("could not extract &XC from one of the inputs")
    out = mlwf.replace(old_xc, lda_xc)

    admm = block(out, "AUXILIARY_DENSITY_MATRIX_METHOD")
    if admm:
        out = out.replace(admm, "")
    out = re.sub(r"^\s*BASIS_SET_FILE_NAME\s+\S*ADMM\S*\n", "", out, flags=re.M)
    out = re.sub(r"^\s*BASIS_SET\s+AUX_FIT\s+\S+\n", "", out, flags=re.M)
    out = out.replace("GTH-PBE", "GTH-PADE")

    for bad in ("&HF", "DFTD3", "VDW_POTENTIAL"):
        if bad in out:
            raise SystemExit(f"'{bad}' survived the transplant - inspect inputs")
    pathlib.Path(args.out).write_text(out)
    print(f"wrote {args.out}")
    print("REVIEW: LDA/PADE XC only, no HF/ADMM/D3, GTH-PADE pseudos, "
          "LOCALIZE + Wannier print sections intact, ABC line present.")


if __name__ == "__main__":
    main()
