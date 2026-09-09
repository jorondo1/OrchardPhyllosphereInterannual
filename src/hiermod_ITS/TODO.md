# TODO list

Running list for the ITS lognormal hierarchical modeling work. Model
reasoning lives in `MODEL_HISTORY.md`.

- [ ] Recalibrate the latest model for bacteria

## Model 7 (MDLSYC) - Control variables

- [x] Remove weather covariates from Model 6, move to Model 7 (`6.1`-`6.3`
  stripped, `7.1`-`7.3` written)
- [x] Add sequencing depth (`Seq_depth` was already in `div`; standardized
  as `seq_depth_z`, logged first, in `0_SETUP.R`)
- [ ] Keep deg-hours? See Doorman et al. 2013 to decide:
  https://doi.org/10.1111/j.1600-0587.2012.07348.x
- [ ] Confirm `log(Seq_depth)` empirically against real richness/discovery
  curves, not just the theoretical rarefaction-curve argument in
  `MODEL_HISTORY.md`
- [ ] Collinearity-aware SBC follow-up: `simulate_from_priors()` draws
  `deg_h_z`/`precip_72h_z`/`seq_depth_z` independently, so SBC tests
  recoverability, not calibration under the real Season/Year correlation
- [ ] Confirm reference-covariate/Cultivar convention reads well in
  writeup (`means_MDLSYC()` defaults all three covariates to 0, omits
  `cv[Cv]` from the headline estimand)
- [ ] Run real fit/SBC/PPC (model code written and smoke-tested, not yet fit)

## Model 6 (MDLSY)

- [ ] Rerun real fit/SBC/PPC: response is now shifted (`Hill_1 - 1`) and
  weather covariates moved out to Model 7, so the numbers in
  `MODEL_HISTORY.md` predate both changes
- [ ] Organic-May under-dispersion (PPC finding): `yrep` too peaked/narrow
  vs broader observed, check `sigma[3]` directly
- [x] `scale_dexp_rate()` variance budget calibration, used from Model 6
  onward (K genuinely grows 3→4 here; no-op for Model 4→5)
- [x] Fixed `idx$Cv` (was 4 blank placeholder levels vs Cultivar's real 5)
- [x] Repo reorg: `ITS/` holds every script as `X.1_model`/`X.2_validation`/
  `X.3_analysis`, retrofitted onto Models 1-5
- [x] `make_index()` gained `palette`/`palette_n()`
- [x] `variance_component_panels()`: one patchwork panel per statistic
- [x] `save_report()` gained `model_name=`
- [x] `variance_partition_MDLSY()`: Bayesian R2 + VPC, verified to sum to
  1 per draw
- [x] Merged `summarize_sbc()`+`save_sbc_report()` into one ggplot-based fn
- [x] `variance_component_panels()` gained `sd_stats=`
- [x] `GLOSSARY.md`: Bayesian/stats term reference for this project

## Model 5 (MDLS2)

- [ ] Compare `postpred_density_MDLS2` vs Model 4: does Conv-July `yrep`
  tighten without hurting the other 3 cells?
- [ ] Tree PPC/SBC budget + Organic-July bimodality (same open items as
  Model 4)
- [ ] Run 5-shifted's real validation (recovery smoke-tested only); decide
  reported model (shifted vs plain vs Model 6)
- [x] Combined-index workaround for `sigma[cell]` (`ulam()` can't
  double-index `dexp`)
- [x] Built `5.1_MDLS2_model.R` + validation/analysis scripts
- [x] Fixed `save_report()` call-site bug (`recovery=` arg)
- [x] Shifted-lognormal variant, bug found and fixed (see
  `MODEL_HISTORY.md`)

## Model 4 (MDLS)

- [ ] PPC/SBC budget for Tree params (~129 `tr[Tree]` vs MDLv's 4-5
  Location; current `n_sbc`/`iter` likely too slow, pick a workable combo)
- [ ] Organic-July bimodality: real, in postpred overlay; needs mixture
  likelihood if pursued, check Location/Year mapping vs small-n (51)
  artifact first
- [ ] Model 3 (MDLv) contrast-stat PPC tail issue never resolved
- [x] Fixed `run_sbc()`'s hardcoded 2-col `contrast_fn` vs MDLS's 4-col
  `mean`, added `contrast_fn=` param
- [x] Rebuilt `sim_div_MDLS()`: Tree-owned Lo/Mg + paired May/July rows
- [x] Shared `estimand_rows()`/`check_recovery()`/`estimand_panels()` for
  Model 5+
- [x] PPC: 4-cell density overlay + 3 contrast-stat checks

## Discussion / writeup notes

- [ ] SBC underconfident but results clean: Model 4's SBC showed a
  hump-shaped rank histogram (posteriors too wide), yet the real fit's
  May/July/Change contrasts all landed with 89% PI clean of zero
- [ ] Raw medians hide the interaction, model recovers it: raw `Hill_1`
  shows a flat Organic season (May 14.9 → July 15.0), but the model
  recovers a real gap-shift. Hypothesis: unbalanced Location/Year sampling
  across cells (Simpson's-paradox-flavoured), confirm before writeup
- [ ] `deg_h`/`precip_72h` collinearity with Season/Year (Model 7):
  confirmed real, `deg_h` r=-0.73 with Season, `precip_72h` r=-0.55; real
  posterior `cor(b_deg,s_conv)=-0.71` (from the earlier weather-in-Model-6
  fit, needs reconfirming in Model 7). Doesn't bias estimates, widens and
  correlates posteriors, can slow sampling

## Potential robustness checks

- [ ] Hill order 2: focus on dominant species (harmonic mean)

## Paper outline (population-parameter results)

Reference bullets for writing this up, not the writeup itself. Pull
numbers from `save_report()`'s `.txt` outputs and figures from
`save_gg()`'s PDFs, don't re-derive.

**Methods**
- Hierarchical (shifted-)lognormal model of Hill_1 (order-1 diversity);
  Management x Season fixed interaction of interest
- Shift: Hill numbers floor at 1, plain lognormal support is (0, Inf);
  `Dv = Hill_1 - 1` fit, floor added back on backtransform
- Estimand definitions: May gap, July gap, seasonal change in gap; each a
  contrast of population-level (posterior mean/median) diversity, not a
  raw group average
- Validation note: prior predictive checks, SBC, posterior predictive
  checks (details in Supplementary)

**Results**
- 3 headline estimands: posterior median + 89% credible interval, in
  text/table
- Main figure: `fit_contrast_mean`/`fit_contrast_median` panels, observed
  reference lines built in
- Optional secondary figure: two-line interaction plot for
  shape-at-a-glance

**Discussion**
- Direction/magnitude of estimated gaps
- Naive-vs-modeled divergence: raw pooled medians suggest a flat Organic
  season, model recovers a real shift; explain why, not just note it
- Robustness: qualitative finding (direction, growth across season) holds
  across Models 3-7, not just the final spec
- Caveats: Saint-Benoit/Windsor Location-Management confound (Model 3),
  Year formally pooled only as of Model 6, `sigma_tr`/`sigma[cell]` only
  weakly separable (~2 obs/tree), Organic-July bimodality if unresolved,
  weather/depth collinearity with Season (Model 7)

**Supplementary**
- PPC diagnostic plots, SBC rank histogram + KS test, parameter recovery
  tables, full `precis()` output
- Brief model-evolution narrative, condense `MODEL_HISTORY.md`, don't redo it
- [x] Build a model-to-LaTeX function: resolved on notation (standard
  hierarchical-model form, `mu_i`/per-observation subscripts, not a
  literal `alist()` transcription); see `model6_MDLSY_latex.txt` for the
  worked Model 6 example
