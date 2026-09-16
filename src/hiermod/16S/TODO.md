# TODO list

Running list for the 16S lognormal hierarchical modeling work.

## Model 5 (MDLS2)

- [ ] Manuscript methods: short paragraph on the prior sensitivity check
  (5.5_MDLS2_16S_sensitivity.R) -- tight (dexp(2.5)) vs loose (dexp(1.5))
  sigma priors barely moved sigma[cell]/sigma_loc/sigma_tr or the May/July
  contrasts, so the calibrated prior isn't driving the result.
- [ ] Manuscript methods: short paragraph justifying sigma[cell] (residual
  SD split by Season, not just Management) -- see model comparison against
  the sigma[Mg]-only reduced fit in 5.6_MDLS2_16S_comparison.R.

## Control variables (Models 5-6, MDSYC/MDSYCV -- from-scratch rebuild)

- [ ] Manuscript result, likely main text not just supp: the Management gap
  reported by every pre-covariate model (MD/MDv/MDS/MDST/MDSY(z), May-mean
  contrast approx -237) roughly HALVES once weather/sequencing covariates
  are controlled for (MDSYC/MDSYCV, approx -107). This should be framed as
  revising those earlier models' own effect size (they were confounded),
  not as reporting two equally-valid alternative estimates.
- [ ] Mechanism, now confirmed not just suspected: Organic samples were
  sequenced significantly deeper than Conventional (Welch t=-4.95,
  p=1.6e-6, ~0.63 SD gap in seq_depth_z), b_seq's posterior is credibly
  negative (-0.204, 89% CI excludes 0), and raw Hill_1 vs Seq_depth
  correlation is r=-0.22 (p<0.001) -- the covariate both differs by group
  AND predicts the outcome, i.e. real confounding, not a spurious control.
  See src/check_seq_depth_confound.R for the supporting figure.
- [ ] Methods: state explicitly that seq_depth_z is derived from the DADA2
  pipeline's post-filtering (pre-rarefaction) read count -- a proxy for
  true sequencing depth, not a direct measurement.
- [ ] Reporting scale for b_deg/b_precip/b_seq: report primarily as
  multiplicative/percent effects (exp(b) - 1) on the log scale, not a
  single Hill-scale slope -- the log-linear likelihood means the absolute
  Hill-scale effect isn't baseline-independent (a fixed % change is a
  different number of Hill_1 units at different Mg x Mo baselines).
  Illustrative Hill-scale predicted curves at reference cells (via
  means_MDSYC()'s existing back-transform) as a secondary/supplementary
  visual, not the primary reported number.
- [ ] Open, real (not just hypothetical) collinearity concern: deg_h_z and
  precip_72h_z differ hugely by Month/Season (deg_h_z +/-0.73 SD,
  precip_72h_z +/-0.55 SD between July/May) -- expected, since season is
  largely weather, but this is exactly what MDLSYC_model.R's own
  rho_deg_season/rho_seq_mu collinearity-aware SBC stress test (ITS
  lineage, never run for 16S) was built to check. Worth actually running
  once Model 7 (Tree merge) settles.
- [ ] Cultivar (Model 6) does not explain the sigma[Mg] Organic/
  Conventional residual-variance asymmetry (median 0.285 vs MDSYC's own
  0.31 -- essentially unchanged). That asymmetry remains open/unexplained
  by anything modeled so far.
- [ ] Decide: does the seq-depth-confound figure go in the main text or
  supplement? Undecided.
