# Step G (scoped): LDA dipoles for the 3 meeting points

Chain (mirrors the proven hybrid v3/v3.5 workflow):
  1. make_lda_mlwf.py      - LDA MLWF template from hybrid mlwf.inp + step-B XC
  2. make_snippets.py      - subsample frames from the 3 production trajectories
  3. gen_dipole_inputs.py  - per-frame CP2K MLWF inputs (clone of make_dipole_labels_v3)
  4. dip_lda.slurm         - label array (settings cloned from htc_dipole_label3.sh)
  5. [needs build_dipole_xyz_v3.py + htc_dipole_train_v35.sh from ~/scw]
     collect -> train water_dipole_lda with the identical v3.5 recipe

Default frames/point: 80 (240 MLWF SPs; LDA MLWF is cheaper than hybrid -
budget well under the hybrid campaign's).
