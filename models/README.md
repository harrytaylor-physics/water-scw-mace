# models/

Trained model files. Use git-lfs (`git lfs track "*.model"`) — files are
tens to hundreds of MB.

| File | What | Labels | Notes |
|---|---|---|---|
| water_scw_v4.model | production potential | revPBE0-D3, 1193 cfgs | apply <= 3000 K |
| water_scw_v4_compiled.model | TorchScript export | — | for LAMMPS/inference |
| water_dipole_v2.model | dipole model | MLWF revPBE0-D3 | feeds spectra.py |
| mace-mpa-0-medium.model | foundation model | MPA-0 | fine-tune parent |
| water_lda_v1.model | LDA potential | LDA/PADE | in training (campaign) |

Copy from ARC: `scp shug7909@htc-login:~/scw/{water_scw_v4,water_dipole_v2,mace-mpa-0-medium}.model models/`

Provenance rule: every model added here gets a row with its label level,
training-set file, and the commit of the training script that produced it.
Never mix label levels within one training set.
