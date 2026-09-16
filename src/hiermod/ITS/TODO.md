# TODO list

Running list for the ITS lognormal hierarchical modeling work. Model
reasoning lives in `MODEL_HISTORY.md`.

## Conventions for the ITS rebuild

The 16S lognormal model family was rebuilt from scratch, one addition at a
time, each fully SBC-validated (n_sbc=100 then n_sbc=400) before moving on
-- the same discipline applies here once the ITS rebuild starts (mirroring
`src/hiermod/16S/1_MD` through `7_MDSTYCV`, not retrofitting the scripts
below, which predate these conventions). Established along the way, apply
from Model 1 onward:
- Each model's `alist()` object gets a display name as an ATTRIBUTE, not a
  list element (`ulam()` iterates every list element expecting a formula
  quote() and errors on anything else -- same reason a fit's cstanfit
  lives at `attr(fit, "cstanfit")` rather than a slot):
  `attr(model_XXX_ITS, "name") <- "Some Tolkien Name"`, right next to the
  model's own `model_id_XXX <- "XXX"` line. `save_report()` picks this up
  automatically (`attr(model, "name")`); only pass `model_name=` when a
  script genuinely reports two models at once.
- Calibration-stage outputs (every `N.2*` script, plus any later
  collinearity-style stress test) write to `<hiermod_out_dir>/Calibration`;
  fit/analysis scripts (`N.3`/`N.4`/comparison scripts) keep writing
  directly to `<hiermod_out_dir>`.
- Once a model has multiple same-family effect panels (Year, covariates,
  Cultivar, Tree, ...), combine them into one `patchwork`-stacked figure
  (`p_year / p_covariates / p_cultivar`) instead of saving each separately.
- Build contrast/effect tibbles via
  `bind_rows(tibble(group = ..., value = ...), ...) %>% mutate(statistic = "...")`
  instead of repeating `statistic = "..."` on every row.
- `n_keep = 1000` as the standard `plot_mcmc_pairs()` thinning target, for
  both calibration-stage and fit-stage pairs checks.
- Break long calls (`save_report(...)`, `plot_mcmc_pairs(...)`) across
  multiple lines for legibility rather than one long line.
- No leftover interactive-only diagnostic calls (`hist(div$Hill_1, ...)`,
  `hist(dat_sim$Dv_shifted, ...)`) left in calibration scripts.
- Weather covariates (`deg_h_z`/`precip_72h_z`) centered WITHIN each Season
  (mean-subtract per Time level, scale by the within-season residual SD),
  not globally -- see `src/hiermod/16S/TODO.md`'s "Control variables"
  section for the full reasoning. `seq_depth_z` stays globally centered
  (its Management confound is a technical artifact to remove, not a
  substantive effect to preserve).

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
- [ ] Organic-May PPC mismatch: not a dispersion/`sigma[cell]` issue (the
  full `total_var`, residual+Location+Tree+Year, already covers the raw
  empirical spread comfortably) -- traced instead to Windsor supplying
  54% of Organic-May/July data with the strongest Location effect
  (`b[Lo]` median +0.5, a ~1.65x multiplier); population-average estimand
  is correctly discounting this, raw pooled median isn't a fair
  benchmark. See "Location x Season interaction" under Potential
  robustness checks for the fuller writeup and what (if anything) to do
  about it
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
- [ ] Location x Season interaction -- shelved for now, revisit if a
  reviewer questions whether the season-flip finding is really just
  Windsor. Full reasoning below so this doesn't need re-deriving.

  **The finding that raised the question:** Model 6's Organic-May PPC
  mismatch traced to Windsor supplying 54% of Organic-May/July data, with
  a raw May->July change (log-scale) of -0.34 -- the opposite sign from
  Compton (+0.61) and Milton (+1.11), the other two Organic-sampled
  orchards.

  **Why it doesn't actually contradict the headline finding:** benchmarked
  against what Conventional does at each site (Compton -1.92, Saint-Benoit
  -1.60, Milton -1.11 -- steep declines everywhere), even Windsor's -0.34
  is far shallower than Conventional's mildest decline. Rough per-location
  `gap_shift`-equivalents (Organic change minus ~-1.5 average Conventional
  decline): Milton ~+2.65, Compton ~+2.15, Windsor ~+1.20. Same sign at
  all 3 Organic-sampled locations -- the "Organic buffers the May->July
  drop" direction holds -- but a real 2x+ spread in magnitude that the
  single pooled `gap_shift` currently can't show.

  **What `h[Lo]` (a Location-varying Season interaction) would buy:** add
  `gamma <- s_conv + gap_shift*(Mg-1) + h[Lo]*sigma_h` (same non-centered
  pattern as `b[Lo]`/`tr[Tr]`/`yr[Yr]`).
  - May gap: untouched -- `gamma` only multiplies `(Mo-1)`, zero in May.
  - July gap / seasonal change in gap: point estimate likely similar
    (`h[Lo]` is mean-zero), but the credible interval should honestly
    widen to reflect real between-location heterogeneity that's currently
    either silently absorbed into `sigma[cell]` or just missing from
    `gap_shift`'s own posterior spread.
  - Buys a real robustness statement instead of one pooled number:
    `sigma_h`'s posterior (how consistent is the interaction across
    orchards) plus per-location `h[Lo]` estimates.

  **Real limits, not fixable by more model structure:** Saint-Benoit has
  zero Organic data (`h[Saint-Benoit]` fully unidentified, pure
  prior/shrinkage); Windsor has zero Conventional data (`h[Windsor]`
  can't separate an Organic-specific season effect from "something about
  Windsor as a place" -- a data gap, not a modeling one); only 3
  informative locations feed `sigma_h`, expect the same weak-identification
  flavour already seen with `sigma_yr` (3 years).

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
