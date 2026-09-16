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
  contrast approx -237) SHRINKS once weather/sequencing covariates are
  controlled for. Under within-Season covariate centering (current,
  MDSYC/MDSYCV/MDSTYCV real fits rerun 2026-09-16): May-mean approx -184 to
  -190, July-mean approx -22 to -24, seasonal-change-in-gap approx +161 to
  +166 (89% PI excludes 0 in all three models). This is a smaller
  correction than the earlier globally-centered result (approx -107) --
  expected and correct: global centering let b_deg/b_precip absorb part of
  the genuine between-season temperature/precipitation difference, which
  now correctly flows into the season term instead (see the within-Season
  centering entry below). Frame as revising the pre-covariate models' own
  effect size (real sequencing-depth confound, see mechanism below), not
  as two equally-valid alternative estimates.
- [ ] Mechanism, now confirmed not just suspected: Organic samples were
  sequenced significantly deeper than Conventional (Welch t=-4.95,
  p=1.6e-6, ~0.63 SD gap in seq_depth_z), b_seq's posterior is credibly
  negative (approx -0.20 to -0.22 across MDSYC/MDSYCV/MDSTYCV, 89% CI
  excludes 0), and raw Hill_1 vs Seq_depth
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
- [x] Open, real (not just hypothetical) collinearity concern: deg_h_z and
  precip_72h_z differ hugely by Month/Season (deg_h_z +/-0.73 SD,
  precip_72h_z +/-0.55 SD between July/May) -- expected, since season is
  largely weather, but this is exactly what MDLSYC_model.R's own
  rho_deg_season/rho_seq_mu collinearity-aware SBC stress test (ITS
  lineage, never run for 16S) was built to check. Run for MDSYCV
  (6.5_MDSYCV_16S_collinearity_check.R, real measured
  rho_deg_season=0.735/rho_precip_season=0.547/rho_seq_mg=0.310): clean at
  n_sbc=100 (0 divergences, max Rhat 1.005), b_deg~s_conv correlates in the
  posterior (r=-0.61) but doesn't translate into miscalibration. NOTE:
  superseded by the within-Season covariate centering decision below --
  once deg_h_z/precip_72h_z are centered within Season, this specific
  collinearity mostly disappears from the real data by construction; the
  rho_seq_mg=0.310 part of this check remains relevant (seq_depth_z stays
  globally centered).
- [ ] Cultivar (Model 6) does not explain the sigma[Mg] Organic/
  Conventional residual-variance asymmetry (median 0.285 vs MDSYC's own
  0.31 -- essentially unchanged). That asymmetry remains open/unexplained
  by anything modeled so far.
- [ ] Decide: does the seq-depth-confound figure go in the main text or
  supplement? Undecided.
- [x] Design decision, implemented 2026-09-16: deg_h_z/precip_72h_z are
  now centered WITHIN each Season (mean-subtract per Time level in
  `0_SETUP.R`, then scale by the within-season residual SD), not globally
  -- otherwise the Mo/gamma/s_conv/gap_shift terms report the seasonal
  shift net of the portion explained by b_deg/b_precip, i.e. "the seasonal
  change if temperature had been constant," which isn't the estimand of
  interest. Within-Season centering makes the covariates orthogonal to Mo
  by construction (confirmed: cor(deg_h_z, Mo) and cor(precip_72h_z, Mo)
  both ~1e-16, down from 0.735/0.547) while day-to-day sampling jitter
  within a season is still controlled for, and the real between-season
  temperature/precipitation difference now flows into the season term,
  where it belongs. seq_depth_z stays globally centered -- its Management
  confound (cor=0.310, unchanged) is a technical artifact we want removed,
  not a substantive effect to preserve. 5.3/5.4, 6.3/6.4, 7.3/7.4 (real
  fits) rerun clean with the new covariates; SBC calibration was
  unaffected as expected (simulators already drew these covariates
  independently of Mg/Mo). See the Management-gap entry above for the
  resulting effect-size shift.
- [ ] Backlog idea (tabled, not built): sigma[cell] -- residual SD split
  by Management x Season (4 cells), not just Management (2 levels, current
  sigma[Mg]). Would let season-specific variance asymmetry between
  Management groups show up explicitly, if it ever looks substantively
  important. Old MDLS2 lineage had a sigma[cell] precedent (see Model 5
  section above) but for Location x Season, not Management x Season --
  worth revisiting that comparison's own methodology if this gets built.
