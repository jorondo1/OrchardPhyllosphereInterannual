# Glossary

Terms used across this project's model-building conversations, grouped by
topic. Written to disambiguate pairs that are easy to swap by accident,
not as a general stats reference; entries assume the reader already knows
roughly what a posterior is. Cross-references point to the concrete R
object/function where a term is realized in this codebase.

## Distributions in the pipeline

- **Prior**: the distribution over a parameter *before* seeing data, e.g.
  `sigma_yr ~ dexp(2)`. Encodes what's plausible a priori.
- **Likelihood**: the distribution of the data given parameters, e.g.
  `Dv ~ dlnorm(mu, sigma[cell])`. The piece that actually reads the real
  measurements.
- **Posterior**: the distribution over parameters after combining prior
  and likelihood via Bayes' rule (`extract.samples(fit)`). What all the
  contrast/variance plots in this project ultimately describe.
- **Derived posterior / derived quantity**: something computed from raw
  posterior draws of one or more parameters, not itself a named parameter
  in the model. Example: `b[Lo]*sigma_loc` (a location's log-scale offset)
  is derived from two separate raw parameters; a Management contrast
  (`mean_org - mean_conv`) is derived from `means_fn()`'s output. Most of
  what this project reports (contrasts, medians, variance shares) is a
  derived quantity, not a raw model parameter, and its shape can differ a
  lot from any of its ingredients' own shapes (see "skew transfer" below).
- **Prior predictive distribution**: what the *data* would look like if
  parameters were drawn from the prior and pushed through the likelihood,
  with no real data involved. Checked via `prior_predictive_spaghetti()`
  (`predictive_checks.R`); catches priors that imply absurd outcomes (e.g.
  Hill-diversity values in the thousands) before touching real data.
- **Posterior predictive distribution**: the same idea, but parameters
  drawn from the *posterior* (after fitting real data). Used for posterior
  predictive checks (PPC): does simulated data from the fitted model
  actually resemble the real observed data? (`plot_ppc_overlay()`,
  `plot_ppc_season_contrast_stats()`.)
- **yrep**: shorthand (bayesplot/Stan convention) for one posterior
  predictive draw, a full simulated dataset generated from one posterior
  parameter draw. What PPC density-overlay plots compare against the real
  `y` (observed data).

**Easy to confuse**: *prior predictive* vs *posterior predictive* use the
same mechanism (simulate data through the likelihood) with a different
source of parameters (prior vs. posterior). *Posterior* (a parameter's
distribution) vs *posterior predictive* (a simulated dataset's
distribution, one level downstream): the predictive version always
answers "what would the data look like," never "what is the parameter."

## Estimands and point summaries

- **Estimand**: the specific quantity of scientific interest being
  estimated, e.g. "the May diversity gap between Organic and
  Conventional." Distinct from any one parameter; an estimand is usually a
  *function* of several parameters (see "derived quantity" above).
- **Contrast**: a difference between two group-level estimands, e.g.
  `mean_org - mean_conv`. `compute_contrasts()`, `contrast_from_means()`,
  `contrast_may_gap_MDLS2()` all produce these.
- **Mean vs. median vs. mode (peak)**: three "typical value" summaries
  that coincide only for symmetric distributions.
  - *mode/peak*: the single most probable value (where a density plot is
    tallest).
  - *median*: the 50th percentile; on this project's working scale
    (log-median parameterization, `loga`), this is what
    `lognormal_mean(mu, total_var=0, shift)` returns.
  - *mean*: the probability-weighted average. For a lognormal,
    `lognormal_mean(mu, total_var, shift) = shift + exp(mu + total_var/2)`,
    always >= the median (Jensen's inequality), and it specifically
    requires `total_var` (marginalizing over every random effect) to
    compute correctly. The median needs no `total_var` at all, though
    random effects still govern posterior *width*, which matters for both
    mean and median even when only the mean's point-formula needs it.
  Under a **skewed** posterior, mean and mode pull apart in the direction
  of the heavy tail: the mean is dragged toward it, the mode isn't. This
  is why the variance-component panels (`variance_component_panels()`,
  `postcontrast_helpers.R`) draw a dashed per-group **mean** line: for a
  skewed derived effect (e.g. a non-centered `z*sigma`-style term), reading
  central tendency off the visual peak alone understates the skew's own
  direction.
- **Marginalizing (over a variable)**: integrating/averaging a joint
  distribution over one of its variables to get the distribution of the
  rest, without conditioning on any particular value of the marginalized-
  out variable. Concretely: `total_var = sigma[cell]^2 + sigma_loc^2 +
  sigma_tr^2 (+ sigma_yr^2 in Model 6)` is the marginal variance once
  Location/Tree/(Year) are averaged over, not fixed at any one level: what
  the mean formula needs, and what the median formula ignores.
- **Skew transfer** (informal term used in this project, not standard
  jargon): in a non-centered random effect (`z[i]*sigma`, `z ~ dnorm(0,1)`,
  `sigma` positive and itself right-skewed, e.g. from a `dexp()` prior),
  the product's skew direction is set by the *sign* of that group's
  typical `z`: negative `z` gives a left-skewed product (heavier left
  tail, since large `sigma` draws push it further negative), positive `z`
  gives a right-skewed product. A structural byproduct of non-centering
  crossed with a positivity-constrained scale, visible in Model 6's
  per-year effect panels (2022/2023/2024 show different skew shapes even
  though all are centered near 0).

## Model structure

- **Fixed effect** (this project's usage): an effect estimated
  independently per level with no shared/pooling prior tying levels
  together beyond a common, un-scaled prior shape, e.g. `loga[Mg] ~
  dnorm(...)` for Management. Also called **unpooled**.
- **Random effect / partially pooled effect**: levels share a common
  population-level scale (`sigma_x`) that is itself estimated from data,
  letting levels "borrow strength" from each other (shrinkage toward the
  population mean, more so for levels with less data). `b[Lo]*sigma_loc`,
  `tr[Tr]*sigma_tr`, and (as of Model 6) `yr[Yr]*sigma_yr` are all
  partially pooled this way.
- **Non-centered parameterization**: the coding trick used throughout this
  project for every partially pooled term: instead of sampling the effect
  directly (`x_raw[i] ~ dnorm(0, sigma_x)`), sample a standardized
  `z`-score (`x[i] ~ dnorm(0,1)`) and multiply by the scale separately
  (`x[i]*sigma_x` used inside `mu`). Mathematically identical, but avoids a
  "funnel" geometry in the sampler that causes divergences when `sigma_x`
  is small.
- **Shrinkage / partial pooling**: the tendency of a partially pooled
  level's estimate to be pulled toward the shared population mean, more
  strongly the less data/more noise that level has on its own. This is
  what random effects actually buy you: honest uncertainty and better
  estimates for sparse groups, not just a change to the point estimate.
- **Hierarchical / multilevel model**: a model with at least one
  partially-pooled (random-effect) grouping structure. Every model from
  Model 4 onward in this project (Location/Tree/Year as random effects on
  top of fixed Management/Season effects).
- **Population** (as in "population SD", "population of years"): the
  higher-level group that a partially pooled effect's levels are treated
  as draws from, e.g. "the population of years" is the hypothetical
  broader set of years that 2022/2023/2024 are a sample of 3 draws from.
  `sigma_yr` is the SD of *that* population, distinct from any one year's
  own estimate.
- **Interaction**: when one predictor's effect depends on the level of
  another, e.g. `gamma <- s + gap_shift*(Mg-1)` (Season's effect on
  diversity differs by Management), the McElreath-style two-equation form
  used since Model 4, where `gap_shift` *is* the interaction term
  directly.
- **Collinearity** (fixed-effect side) / **aliasing** (variance-component
  side, informal usage here): two predictors so correlated in the actual
  observed data that the model can't cleanly separate their individual
  contributions, inflating and correlating their posteriors (e.g.
  `cor(b_deg, s_conv) = -0.71`, `deg_h` vs. Season `r = -0.73`). Doesn't
  bias the *joint* fit or any derived estimand that doesn't isolate one
  coefficient alone; dropping one collinear predictor doesn't "clean up"
  the other's estimate, it lets the dropped variable's explanatory power
  silently transfer into whichever predictor remains.

## Variance decomposition / R2D2M2

- **R2D2M2 / R2-D2 priors** (Aguilar & Bürkner 2023 EJS; generalizes
  Zhang/Naughton/Bondell/Reich 2020 JASA and Yanchenko/Bondell/Reich's
  GLMM extension; background: Gelman 2006): the named literature for
  holding a *total* variance budget roughly fixed as more variance
  components are added to a model, by scaling each component's prior rate
  so their combined expected variance doesn't balloon just from adding
  terms.
- **Variance budget calibration (VBC)** (this project's shorthand for
  applying the R2D2M2 idea): `scale_dexp_rate(rate_ref, k_ref, k_new)` in
  `hiermod_core.R`. For `X ~ Exponential(rate)`, `E[X^2] = 2/rate^2`, so
  scaling `rate` by `sqrt(k_new/k_ref)` holds `E[sum(sigma_i^2)]` fixed
  across a change in how many *summed* variance terms feed `total_var`.
  Does **not** address (a) more values within one already-summed slot
  (e.g. `sigma[cell]`'s 4 cells vs. a 2-cell version, a max-of-n-draws
  effect, not a sum-of-K effect), or (b) new additive fixed-effect terms in
  `mu` itself (`cv[Cv]`, `b_deg`, `b_precip`, `b_seq`): these affect mean
  and median directly via `exp(mu)` and are outside this mechanism
  entirely. Used from Model 6 onward; tabled for Model 5 since K didn't
  actually grow there. Model 7 doesn't run a new calibration pass either
  (same reason, its three added covariates are (b) not (a)) -- it just
  hardcodes Model 6's already-calibrated rates into its own priors.
- **Bayesian R2 / variance-partition coefficient (VPC)**: the share of
  total outcome variance attributable to each variance source (fixed
  effects combined = "explained"/R2, plus each random effect's own
  `sigma_x^2`, plus residual), computed *per posterior draw* so it has its
  own full distribution, not just a point estimate. All shares sum to
  exactly 1 per draw (a true composition); `variance_partition_MDLSY()`
  (`6.1_MDLSY_model.R`), visualized as the stacked-bar (posterior medians)
  + ridge-density (full posterior shape) combo in `6.3_MDLSY_analysis.R`.

## Diagnostics (sampling / calibration)

- **Divergence**: a Hamiltonian-Monte-Carlo transition that failed the
  sampler's numerical-accuracy check, usually a sign the posterior
  geometry is hard to explore accurately near that region (classic cause:
  a "funnel" where `sigma_x` is small and its associated `x` term is
  unconstrained). Report a divergence *rate*, not a raw count, relative to
  total sampling transitions (`(iter/2)*chains` per fit here, given this
  codebase's default 50/50 warmup/sampling split), since a raw count alone
  is meaningless without knowing the denominator.
- **R-hat (Rhat)**: a convergence diagnostic comparing between-chain vs.
  within-chain variance for a parameter; should be ~1.00. Values well
  above 1 mean chains haven't mixed/agree.
- **ESS (effective sample size)**, bulk/tail: how many *independent*-ish
  samples the (autocorrelated) MCMC draws are worth, for estimating the
  distribution's center (bulk) or its extremes (tail). Low ESS relative to
  total draws means high autocorrelation, noisier estimates than the raw
  draw count suggests.
- **E-BFMI (Energy-Bayesian Fraction of Missing Information)**: flags when
  the sampler's momentum resampling step isn't efficiently exploring the
  posterior's energy levels. Low E-BFMI often co-occurs with divergences/
  max-treedepth hits and, like divergences, points to a real
  posterior-geometry problem rather than "just run it longer."
- **Max treedepth**: NUTS hit its step-count ceiling before naturally
  U-turning; frequent hits suggest either a genuinely hard-to-explore
  posterior or a treedepth limit set too low for this model's shape.
- **`pairs()` / `mcmc_pairs()` diagnostic plot**: pairwise scatter of
  parameters (often with divergent draws marked) used to *localize* which
  specific parameter pair is driving divergences/E-BFMI/treedepth issues,
  e.g. a tight diagonal ridge between two parameters signals they're only
  weakly separately identified by the data.
- **SBC (Simulation-Based Calibration)**: fit the model to many datasets
  simulated from its own priors, and check whether the true (simulating)
  parameter value's *rank* among the resulting posterior draws is
  uniformly distributed across replicates. A properly calibrated model
  should show a flat rank histogram; a systematic U-shape/skew (not just
  one noisy bin) signals a real calibration problem. `run_sbc()` /
  `save_sbc_report()` (`sbc_helpers.R`). Distinct from divergence checking:
  SBC is about calibration *across replicates in aggregate*, not a
  per-replicate sampler-health check.
- **KS test (Kolmogorov-Smirnov)**: the statistical test SBC's rank
  histogram is checked against (`ks.test(ranks, "punif")`): are the rank
  statistics consistent with a Uniform(0,1) distribution? A low p-value
  means a real deviation from calibration, not just visual noise in one
  bin.
- **Parameter recovery (check)**: fitting the model to one dataset
  simulated from known, hand-picked "true" parameter values and checking
  that the posterior actually recovers something close to those true
  values. Cheaper/faster than full SBC (one replicate, not many), good for
  a first sanity pass, but doesn't test calibration *across* the space of
  plausible parameter values the way SBC does.
- **PPC (Posterior Predictive Check)**: compare data simulated from the
  fitted posterior (`yrep`) against the real observed data, visually
  (`ppc_dens_overlay_grouped()`) or via specific test statistics
  (`ppc_stat()`, `plot_ppc_season_contrast_stats()`). Tests whether the
  *fitted* model, not just its priors, produces realistic data, distinct
  from prior-predictive checking, which never touches real data at all.

## Frequentist/Bayesian framing (background, used once for contrast)

- **P(data | H0)**: the frequentist null-hypothesis-testing quantity (a
  p-value): how surprising the observed data would be *if* a null
  hypothesis were exactly true. Answers a different question than
  estimation.
- **P(parameter | data)**: the Bayesian posterior: given the data actually
  observed, what's the distribution over the parameter's possible values.
  This is what every fit in this project reports directly, via Bayes' rule
  and an explicit prior; no separate "hypothesis test" step is needed to
  get an estimate with uncertainty attached.
