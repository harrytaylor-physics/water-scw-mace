# Entries to append to AMENDMENTS.md

- LDA training recipe deltas vs hybrid v4: batch 4 vs 2, 400 vs 200
  epochs, SWA start 300 vs 120. Rationale: throughput; parity gates
  passed.
- GE.1 anchor FAIL at dense-hot; diagnosed as proposal-ensemble coverage
  (fused-O frame at 0.822 A labels at 161 eV/A under LDA). v2 retrain
  shelved for meeting scope; harvest data banked.
- Production dump_every corrected 10 -> 4 (5 fs -> 2 fs) after Nyquist
  truncation of the stretch band; three comparison points re-run; archive
  heterogeneous (18 @ 5 fs, 3 @ 2 fs).
- GF.2 waived at T3000/rho108.8 after three draws (9.9 / 6.4% with
  10 / 20 ps settle); achieved 2807 K used in QCF and Planck weighting
  and stated on figures.
- LDA dipole model: batch 4 vs v3.5's 10 (memory); epoch-58
  best-validation checkpoint selected; 640-config retrain increased ACF
  noise — data lever documented as spent.
- LDA kappa_R and absolute band intensities floor-limited; not quoted.
- EOS pressures: MLIP stress + analytic ideal-gas kinetic term;
  per-functional own-trajectory convention; hybrid anchored to IAPWS at
  the benchmark point.
- H-bonded dimer scan geometry does not reach the canonical well;
  binding curve withheld, D3 difference curve quoted instead.
