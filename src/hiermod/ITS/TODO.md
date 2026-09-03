# hiermod TODO — deferred items

Running list of things pinned for later in the ITS lognormal hierarchical
modeling work (models 1-4). Add to it as new deferrals come up; check items
off (or delete them) once addressed. Done items are kept as a one-line
changelog, not a full writeup -- the reasoning lives in the code/commits.

## Model 4 (MDLS) -- open

- [ ] **PPC/SBC strategy for Tree-level parameters.** ~129 `tr[Tree]`
      parameters vs. MDLv's 4-5 Location ones. Each SBC replicate refits the
      whole model, so MDLb's settings (n_sbc=100, iter=10000) are probably
      too slow here. Decide on a workable n_sbc/iter combo before running
      the SBC section for real.

- [x] **Pool `yr[Yr]` eventually.** Done in Model 6 (MDLSY):
      `yr[Yr]*sigma_yr`, non-centered, `means_MDLSY()`'s `total_var` gains
      `sigma_yr^2`. Still applies to Models 1-5, which keep `yr[Yr]` fixed
      -- not retrofitted into their own alists, only Model 6 changes this.

- [ ] **Organic-July bimodality -- watch, don't chase yet.** The postpred
      density overlay shows real bimodality in observed Organic-July data
      that a single lognormal component can't reproduce (unlike the
      Conventional-July case above, this isn't a sigma-sharing fix -- would
      need a mixture likelihood). Could be genuine ecological
      heterogeneity or a small-n (51 obs) artifact; check whether it maps
      onto specific Locations/Years before deciding it's worth modeling.

- [ ] **Double-check Model 3's (MDLv) contrast-statistic PPC tail issue.**
      Observed median/MAD contrast landed in the tails there, never
      resolved. Model 4's own contrast-stat PPC (`postpred_stat_*`) came
      back clean, but confirm Model 3's cause before fully trusting that.

## Model 4 (MDLS) -- done

- [x] Fixed `run_sbc()`'s hardcoded 2-column contrast assumption
      (`contrast_from_means()`) silently mismatching `means_MDLS()`'s
      4-column `mean` -- root cause of the first SBC run's catastrophic
      D=0.44 miscalibration. `run_sbc()` now takes `contrast_fn=`;
      `contrast_may_gap_MDLS()` supplies Model 4's.
- [x] Rebuilt `sim_div_MDLS()`: Tree owns a fixed Mg/Lo + paired May/July
      row; `N_samples` + true params only, `p_dropout` for design
      imbalance; Parameter recovery section checks fixed effects and all
      3 estimands.
- [x] Shared `estimand_rows()`/`check_recovery()`/`estimand_panels()`
      (postcontrast_helpers.R, hiermod_core.R) for Model 5+ to reuse.
- [x] Posterior contrast: mean- and median-scale panels via
      `estimand_panels()`, plus an interaction plot and a forest plot of
      the 3 contrasts as a more efficient headline display (density
      panels were redundant -- all 3 are algebraic functions of the same
      few fixed effects).
- [x] Posterior predictive check: overall density overlay (4-cell
      Mg×Mo via `interaction()`, no need to split calls) + 3 hand-written
      contrast-stat PPCs (May gap/July gap/change).

## Model 5 (MDLS2) -- open

- [ ] **Run the expensive stuff for real.** Model spec + Parameter
      recovery section are built and smoke-tested at small scale; Prior
      predictive check / SBC (n_sbc=100, iter=5000, ~129 Tree params, same
      cost profile as Model 4) / real Model fit / PPC still need a real run.

- [ ] **Compare `postpred_density_MDLS2.pdf` against Model 4's.** The whole
      point of the 4-cell sigma: does Conventional-July's `yrep` now hug
      the observed spike, without degrading the other 3 cells?

- [ ] Same Tree-level PPC/SBC-strategy and Organic-July-bimodality items as
      Model 4 above -- unresolved there, still open here too, same model
      structure otherwise.

- [ ] **Run `5_ITS_lognormal_MDLS2_shifted.R`'s real validation.** Parameter
      recovery smoke-tested at reduced scale only; still need the real
      Prior predictive check / SBC / Model fit / PPC run, and a check of
      whether divergences are meaningfully lower than plain Model 5's once
      the likelihood-misspecification bug is out of the picture. Decide
      whether the shifted variant becomes the reported model or stays a
      robustness check against the plain one.

## Model 5 (MDLS2) -- done

- [x] Verified empirically that `ulam()` can't double-index a
      `dexp`-distributed scale parameter (`sigma[Mg,Mo]`) -- works for
      `dnorm`-type varying effects (`matrix[i,j]:a ~ dnorm(...)`) but not
      here; cmdstan errors either way naive double-indexing or the matrix
      idiom is tried. Used a combined index instead: `cell = (Mg-1)*2 + Mo`,
      `sigma[cell] ~ dexp(3)`.
- [x] Built `models/model_5_MDLS2.R` (sources `model_4_MDLS.R`, patches
      just the likelihood/sigma-prior lines) and `5_ITS_lognormal_MDLS2.R`
      -- see "Code organization" below.
- [x] Found and fixed a real bug in `save_report()`/its Model 4-5 call
      sites: `fit_summary` calls were passing `post_full()`'s raw list
      output where `report_contrasts_full()`/`compute_contrasts()` needs
      the long `statistic`/`group`/`value` tibble -- broke
      `report_contrasts_full(pf)` outright. Fixed both scripts' call
      ordering (build `pc_estimands_means` first) and added a `recovery=`
      argument to `save_report()` so `check_recovery()` tables get logged
      into the `.txt` report, not just printed to console.
- [x] **Shifted-lognormal variant, `5_ITS_lognormal_MDLS2_shifted.R`.**
      Hill numbers can't go below 1, so `Dv` is fit as `Dv_shifted - 1`
      (raw floor-0 lognormal draw) rather than `Dv_shifted` itself
      directly. Two things verified the hard way, not assumed:
      - `Dv - 1 ~ dlnorm(...)` (expression on the likelihood's LHS)
        doesn't compile in ulam -- needs a precomputed plain-symbol
        column instead.
      - An early version fit the likelihood to `Dv_shifted` (the
        floor-1 quantity) directly -- confirmed via a side-by-side sim
        test to systematically bias `sigma[cell]` low (e.g. recovered
        0.43 vs true 0.56) and was very likely the main driver of an 8%
        divergence / 33% max-treedepth / all-chains-bad-E-BFMI blowup in
        Parameter recovery. Fixed by leaving `model$likelihood`
        completely unmodified from Model 5 (`Dv ~ dlnorm(...)`, always
        the unshifted/floor-0 quantity) -- the shift only ever belongs in
        `means_MDLS2()`'s backtransform (`shift=` argument) and in how
        real data is prepared (`Dv = Hill_1 - 1`, subtracting the floor
        rather than adding it).
      - No separate "Model 6" needed *for the shifted variant* in the end
        -- `sim_div_MDLS2()`/`means_MDLS2()`/`simulate_from_priors()` (all
        in `5.1_MDLS2_model.R`) now take an optional `shift=`, so one
        model file serves both the plain and shifted numbered scripts.
        `post_full()` (postcontrast_helpers.R) gained a `...` passthrough
        so `shift=1` can reach `means_fn` through it. (Model 6/MDLSY did
        end up getting built afterward, for unrelated reasons -- pooling
        Year, adding weather covariates and Cultivar -- see below.)

## Model 6 (MDLSY) -- open

- [ ] **Run the expensive stuff for real.** Model spec, `ulam()`
      feasibility (the new `b_deg*deg_h_z` idiom specifically), and the
      full 6.1/6.2/6.3 pipeline are smoke-tested at reduced scale only
      (real data, low iteration); Parameter recovery / Prior predictive
      check / SBC / real Model fit / PPC still need a real run.
- [ ] **Collinearity-aware calibration check, if `b_deg`/`b_precip`'s
      posteriors end up doing real interpretive work.** `simulate_from_priors()`
      draws `deg_h_z`/`precip_72h_z` as independent `rnorm(0,1)` -- this
      validates that the model *can* recover `b_deg`/`b_precip` in
      principle, not how well-calibrated it is under the real design's
      `deg_h`/`precip_72h` ~ Season/Year correlation specifically (see the
      Discussion note below). A follow-up SBC variant that resamples
      `deg_h_z`/`precip_72h_z` from real `div` (possibly jointly with
      `Mo`/`Yr`, to preserve the correlation) would close that gap -- not
      required to trust the current build, but worth doing before leaning
      heavily on the weather coefficients' precision in the writeup.
- [ ] Decide the reference weather/Cultivar convention is actually what
      gets reported: `means_MDLSY()` defaults to `deg_h_z=precip_72h_z=0`
      (this sample's own average weather) and omits `cv[Cv]` from the
      headline estimand entirely (reported separately, same precedent
      `yr[Yr]` set pre-Model-6) -- confirm this reads naturally once real
      numbers are in hand, not just in the abstract.

## Model 6 (MDLSY) -- done

- [x] Fixed `idx$Cv` (`0_SETUP.R`): was `make_index(div$Cultivar, levels =
      c("", "", "", ""))` -- 4 blank placeholders vs. Cultivar's real 5
      levels, silently making every `to_index()` call return `NA`.
- [x] Repo reorg: `hiermod/ITS/` now holds every ITS-specific script
      (numbered scripts + what used to be `models/*.R`, now flat
      `X.1_..._model.R` files, + `archive/`/`TODO.md`/`MODEL_HISTORY.md`/
      `CODE_REVIEW.md`); `hiermod/utils/` stays put (generic, model-
      agnostic). Every model now three files: `X.1_..._model.R`
      (definition only) -> `X.2_..._validation.R` (recovery, prior-PC,
      SBC, the real fit, PPC) -> `X.3_..._analysis.R` (posterior contrast,
      variance components, no refit). Models 1-5 (+shifted) retrofitted
      into this shape in the same pass, not just Model 6 going forward.
      Forest/interaction plots dropped from every `_analysis.R` file
      (density-only views preferred). Old top-level files were git's
      problem to remember, not duplicated into a manual backup.
- [x] `make_index()` (hiermod_core.R) gained `palette`/`palette_n()` --
      Location/Cultivar/Management colour schemes now live on `idx`
      itself instead of being hand-built per plot script.
- [x] New `variance_component_panels()` (postcontrast_helpers.R): one
      patchwork panel per statistic, each with its own right-side legend,
      for panels with many unrelated colours (Location/Year/Cultivar/
      sigma's) where one shared bottom legend got unreadably long.
- [x] `save_report()` gained `model_name=` -- the `.txt` header used to
      always print the literal string `"model"` (every call site's
      variable is named that); every model now gets a real (fittingly
      irreverent) nickname instead.
- [x] `variance_partition_MDLSY()` (6.1_MDLSY_model.R): Bayesian R2 and
      the variance-partition coefficients (VPC) are the same decomposition
      -- "Explained" is `var(fixed-effect part of mu across observations,
      per draw) / total`, the rest is each random effect's population
      variance (+ cell-membership-weighted residual) as a share of that
      same total. All 5 shares sum to 1 per draw (verified in-session, not
      just assumed) -- plotted as a stacked bar (posterior medians) with
      an aligned ridge-density panel above it for the full posterior shape.

## Discussion / writeup notes

- [ ] **SBC said underconfident, results still came out clean.** Model 4's
      SBC (n_sbc=100) showed a hump-shaped rank histogram (posteriors a
      bit too wide -- the conservative failure mode), yet the real fit's
      May/July/Change contrasts each landed with an 89% PI clean of zero.
      A clean result surviving a model that (if anything) overstates
      uncertainty is a point in its favour, not a contradiction -- worth
      including in the writeup.

- [ ] **`deg_h`/`precip_72h` collinearity with Season/Year (Model 6).**
      Confirmed in the real data, not hypothetical: `deg_h` correlates
      -0.73 with Season and varies sharply by Year (940->1060->1559 across
      2022-2024); `precip_72h` correlates -0.55 with Season, also varies
      by Year. What this does and doesn't imply:
      - Doesn't bias `b_deg`/`b_precip`/`s_conv`/`gap_shift`/`yr[Yr]` --
        Bayesian regression with correlated predictors stays valid, just
        less precise (wider, correlated posteriors), not wrong. Not the
        near-rank-deficient case either (r=-0.73/-0.55, not ~1) -- there's
        real separable information, just less of it than if the
        covariates were independent.
      - Can slow/complicate sampling -- correlated predictors create
        harder posterior geometry for HMC (same flavour as `sigma_tr`/
        `sigma[Mg]`'s ~-0.39 correlation already seen in Model 4, just a
        fixed-effect version of it) -- watch divergences/ESS on `b_deg`/
        `b_precip`/`s_conv`/`gap_shift`/`yr[Yr]*sigma_yr` specifically once
        the real fit runs, not just overall.
      - Does leave a real gap between "verified the model can recover
        these parameters" and "verified it's well-calibrated under the
        real design's specific collinearity" -- SBC as built tests the
        former only (independent synthetic `deg_h_z`/`precip_72h_z`, see
        Model 6 -- open above). Worth being explicit about this
        distinction in the writeup rather than treating SBC as having
        fully validated the weather coefficients.
      - Practical read: `deg_h`/`precip_72h` are here to *control for*
        weather so the Management/Season estimand isn't confounded by "May
        just happens to be cooler/wetter than July" -- that's a legitimate
        reason to include a covariate even when it's correlated with
        Season already in the model, standard practice for isolating a
        phenological effect from a weather-driven one. The cost is
        wider/more entangled uncertainty on both sides, not invalidity.

- [ ] **Raw medians hide the interaction; the model recovers it.** Raw
      `Hill_1` medians: Conventional May 17.8 (n=72) -> July 3.41 (n=75);
      Organic May 14.9 (n=50) -> July 15.0 (n=51) -- Organic looks flat
      raw, but the fitted model recovers a real May->July gap-shift.
      Leading hypothesis: unbalanced Location/Year sampling across the
      Mg×Mo cells distorts the raw pooled medians (Simpson's-paradox-
      flavoured), which `b[Lo]`/`yr[Yr]` correct for -- confirm/explain
      properly before writing this up.

## Paper outline (population-parameter results)

Short reference bullets for writing this up -- not the writeup itself.

- [ ] **Write the manuscript Methods/Results/Discussion using the bullets
      below**, pulling numbers from `save_report()`'s `.txt` outputs and
      figures from `save_gg()`'s PDFs rather than re-deriving anything.

**Methods**
- Hierarchical (shifted-)lognormal model of Hill_1 (order-1 diversity);
  Management x Season as the fixed-effect interaction of interest, Year as
  an unpooled fixed control, Location and Tree as partially-pooled random
  intercepts (repeated May/July measures per tree).
- One line on the shift: Hill numbers floor at 1, a plain lognormal has
  support (0, Inf) -- `Dv = Hill_1 - 1` fit, floor added back on
  backtransform.
- Explicit estimand definitions: May gap, July gap, seasonal change in gap
  -- each a contrast of population-level (posterior mean/median) diversity,
  not a raw group average.
- One line noting the model was validated via prior predictive checks,
  simulation-based calibration, and posterior predictive checks (details
  in Supplementary, not Methods itself).

**Results**
- The 3 headline estimands: posterior median + 89% credible interval,
  reported in text/table.
- Main figure: `fit_contrast_mean`/`fit_contrast_median` panels (raw group
  posteriors + Contrast, per estimand), with the observed (raw-data)
  reference lines now built in.
- Optional secondary figure: the two-line interaction plot, for readers
  who want the shape at a glance before the density panels.

**Discussion**
- Interpretation of direction/magnitude of the estimated gaps.
- The naive-vs-modeled divergence as its own point: raw pooled medians
  suggested a flat Organic season, the model recovers a real shift --
  worth explaining *why* (unbalanced Location/Year sampling across cells),
  not just noting that it happened (see Discussion notes below).
- Robustness: the qualitative finding (direction of the gaps, growth
  across season) holds up across Models 3-5, not just the final
  specification -- worth stating explicitly rather than only reporting
  the last model's numbers.
- Caveats: Saint-Benoit/Windsor Location-Management confound (Model 3),
  Year not formally pooled (can't report a year-to-year variance, only 3
  point estimates), `sigma_tr`/`sigma[cell]` only weakly separable given
  ~2 observations/tree, Organic-July bimodality if still unresolved.

**Supplementary**
- PPC diagnostic plots (density overlay, contrast-stat checks), SBC rank
  histogram + KS test, Parameter recovery tables, full `precis()` output.
- Brief model-evolution narrative (`MODEL_HISTORY.md` has the full trace
  already written -- condense, don't redo).

## Cross-cutting -- done

- [x] Removed the `mean_`/`cv_` round-trip in Prior predictive
      check/SBC setup across all scripts (confirmed exact identity);
      simulators take `loga`/`sigma` directly. Generalized `draw_true()`
      into one shared version (`sbc_helpers.R`); `cv_to_sigma()`
      (hiermod_core.R) kept for picking a human-readable CV by hand.
- [x] **New `models/*.R` files, one per numbered script** (`model_1_MD.R`
      .. `model_5_MDLS2.R`): each holds that model's `alist()`, `means_fn()`,
      `sim_div_fn()`, and (from script 2 onward) `simulate_from_priors()`/
      any custom `contrast_fn()` -- narrative/rationale and the actual
      fits/plots stay in the numbered script. Model 5 sources Model 4's
      file and patches only the 2 changed named elements directly
      (`model$likelihood <- quote(...)`; no override-checking helper, not
      worth it for 1-2 lines). Every retrofitted script re-verified
      (sourced + a fast fit) against its pre-extraction behavior.
- [x] **Found and fixed a real, pre-existing bug in `check_recovery()`**
      (hiermod_core.R): it unconditionally `exp()`'d `true` before
      comparing against `post_draws`' PI, but `post_draws` is always on
      whatever native scale the caller already used (log scale for
      loga/s/gap_shift, natural scale for sigma-type parameters) -- a
      scale mismatch that made `covered` meaningless (and the `true`
      column not comparable to `post_median`) for every call site,
      including Model 4's already-run `fixed_recovery` table. Fixed to
      compare like-for-like, no backtransform; re-verified with a real
      fit showing sane `covered=TRUE` results. **Model 4's own
      `fixed_recovery` output should be re-run/re-read** -- prior readings
      of it aren't trustworthy.
