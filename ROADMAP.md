# Roadmap — what we want to understand

Not a task list; a list of questions, roughly ordered by leverage.

1. **How much does the functional matter for opacity?** The LDA-everywhere
   vs revPBE0-D3-everywhere comparison (running). Headline numbers:
   Planck/Rosseland-mean ratios and absorption-edge shifts per state point.
2. **How much of that is dispersion specifically?** revPBE0 with D3 toggled
   off — same exchange, dispersion isolated. Cheap once the LDA machinery
   exists; turns "functionals differ" into "here is why".
3. **Where does structure stop mattering and electronic response take
   over?** Cross pipelines (hybrid optics on LDA trajectories) only if
   referees or the numbers demand attribution.
4. **Is emissivity really density-only?** Current result: e_P collapses
   onto n_inf(rho). Partly built in by the single-alpha_e baseline. The
   per-state-point polarizability upgrade answers whether temperature
   genuinely drops out.
5. **Nuclear quantum effects.** The Berens-style VDOS machinery exists in
   the demixing project; applying it to opacity-relevant states would
   bound the classical-nuclei assumption.
6. **Dissociation onset and ionic conductivity.** sigma_DC from the
   emissivity diagnostics vs T,rho — where does molecular water end?
   Connects to the superionic literature and to the >3000 K descope.
7. **Line-list vs MD crossover.** At what density do HITEMP/MT-CKD line
   lists stop being adequate and the MD-based continuum take over? A
   quantitative crossover map would serve the whole community.
8. **Mixtures.** H-He-H2O opacity on top of the demixing project's
   trajectories — the long-range convergence of both DPhil threads.
