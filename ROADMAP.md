# Roadmap — what we want to understand

Not a task list; a list of questions.

1. **How much does the functional matter for opacity?** The LDA
   vs revPBE0-D3 comparison (running). Headline numbers:
   Planck/Rosseland-mean ratios and absorption-edge shifts per state point.
2. **How much of that is dispersion specifically?** revPBE0 with D3 toggled
   off — same exchange, dispersion isolated. Cheap once the LDA machinery
   exists.
3. **Is emissivity really density-only?** Current result: e_P collapses
   onto n_inf(rho). Partly built in by the single-alpha_e baseline. The
   per-state-point polarizability upgrade answers whether temperature
   genuinely drops out.
4. **Nuclear quantum effects.** Do the warm dense nuceli need a quantum treatment? How relevent, to the opacity and EOS, is the classical-nuclei assumption at supercritical conditions?
5. **Dissociation onset and ionic conductivity.** sigma_DC from the
   emissivity diagnostics vs T,rho — where does molecular water end?
   Connects to the superionic literature and to the >3000 K descope.
6. **Line-list vs MD crossover.** At what density do HITEMP/MT-CKD line
   lists stop being adequate and the MD-based continuum take over? A
   quantitative crossover map would serve the whole community.
