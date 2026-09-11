# Model history — ITS diversity hierarchical models

Retrospective trace of *why* each model exists, in file order. Not a to-do
list (see `TODO.md`) and not a methods section: a record of the reasoning
trail, so a later reader can see that each layer of complexity answered a
specific question rather than chasing a better-looking fit.

## Why each variable is modeled the way it is

- **Management (`loga[Mg]`).** Fixed, unpooled (2 levels, no wider
  population to pool across). It's the comparison of interest, not
  something to shrink toward a common mean.
- **Location (`b[Lo]*sigma_loc`).** Partially pooled. Locations are a
  sample from a wider population of orchards with uneven sample sizes, and
  reporting a "mean diversity" at all requires marginalizing over
  Location-level variation (`total_var` in every `means_*()`).
- **Tree (`tr[Tr]*sigma_tr`).** Partially pooled, added at Model 4. The
  same trees are sampled in both May and July, so without a Tree term
  those rows get wrongly treated as independent, understating uncertainty.
- **Year (`yr[Yr]`).** Fixed/unpooled through Model 5, despite only 3
  years and uneven sampling. Switched to partially pooled
  (`yr[Yr]*sigma_yr`, non-centered) at Model 6, structurally identical to
  `b[Lo]*sigma_loc`.
- **Cultivar (`cv[Cv]`).** Fixed, unpooled, added at Model 6 (5 levels).
  Left out of the headline Management x Season estimand, reported via its
  own posterior panel instead. Expected to show ~no effect, but kept in
  the posterior record rather than assumed away.
- **Weather and sequencing depth (`b_deg`, `b_precip`, `b_seq`).**
  Continuous control covariates, added at Model 7. Standardized once,
  centrally, in `0_SETUP.R`. Weather is known to correlate with Season and
  Year; this is an identifiability limit of the observed design, not
  something any model spec fixes.
- **Season (`gamma <- s_conv + gap_shift*(Mg-1)`).** Not a random effect:
  only 2 levels, each individually meaningful, and the May->July *change*
  is itself the estimand. McElreath's two-equation form makes `gap_shift`
  the interaction directly.
- **`sigma[Mg]` / `sigma[cell]`.** Residual spread, allowed to differ by
  group. A single shared residual SD visibly distorted the fit once groups
  actually differed in spread (seen in Model 2's simulation demo and
  Model 4's postpred mismatch).
- **`g[Lo]*sigma_g`** (Model 3 only). Location-varying Management gap,
  confirmed real but dropped from Model 4 onward to isolate Season first.
- **Non-centered form everywhere.** Standard HMC efficiency fix: the
  centered form creates a sampling funnel when the hyperparameter is
  small; sampling a `dnorm(0,1)` draw and multiplying by the scale avoids
  that.
- **The floor-1 shift (`Hill_1 - 1`).** A support-matching fix, not an
  effect-structure choice: Hill numbers can't go below 1, plain lognormal
  support is `(0, Inf)`.

## Model 1 — `1_ITS_lognormal_MD.R` (MD / MDv)

**Structure:** `Dv ~ dlnorm(mu, sigma)`, `mu <- loga[Mg]`. Two variants:
MD (one shared `sigma`), then MDv (`sigma[Mg]`, one per Management group).

**Motivation:** Does Hill diversity differ between Conventional and
Organic, and is a single shared residual spread even safe? MDv exists
because the constant-variance model visibly distorted the recovered mean
contrast when simulated data had unequal group variances.

**Status:** No SBC section. Both models are simple enough that a full
calibration check wasn't necessary.

## Model 2 — `2_ITS_lognormal_MDL.R` (MDL / MDLb)

**Structure:** + `b[Lo]*sigma_loc` (non-centered, partially-pooled
Location random intercept).

**Motivation:** Observations come from multiple orchards; treating them as
exchangeable ignores real sampling structure. First model where the
fixed-effect estimand needs to marginalize over a population of locations.

**What it revealed:** The exploratory priors implied an implausibly heavy
tail in the prior predictive check. MDLb tightened `sigma[Mg]`
(dexp(1)→dexp(3)) and `sigma_loc` (dexp(1)→dexp(2)) in response, and is
the version actually fit to real data. First model with a real SBC
section, establishing the pattern every later model reuses.

## Model 3 — `3_ITS_lognormal_MDLv.R` (MDLv)

**Structure:** + `g[Lo]*sigma_g*(Mg-1)`: the Organic-Conventional gap
itself varies by Location, on top of Model 2's shared Location intercept.

**Motivation:** Management maps to one specific orchard pairing per
Location, so "does the gap differ by Location" was worth testing directly
rather than assuming Model 2's shared-gap simplification was fine.

**What it revealed:** `sigma_g` sits comfortably above 0: the gap really
does vary by Location. Saint-Benoit and Windsor (Conventional-only /
Organic-only) are a real limitation: `g[Lo]` there is entangled with
`b[Lo]`, not informative on its own.

## Model 4 — `4_ITS_lognormal_MDLS.R` (MDLS)

**Structure:** Refocused on Season instead of Location-varying gap:
`loga[Mg] + gamma*(Mo-1) + b[Lo]*sigma_loc + yr[Yr] + tr[Tr]*sigma_tr`,
`gamma <- s_conv + gap_shift*(Mg-1)`. Adds `yr[Yr]` (fixed, unpooled) and
`tr[Tr]*sigma_tr` (Tree-level random intercept).

**Motivation:** Same trees were sampled in both May and July, so treating
those rows as independent understates uncertainty, a materially worse
problem than the nuisance-variance questions Models 2-3 addressed. Year
entered as fixed because sampling was unbalanced across years/locations.
Model 3's Location-varying gap was dropped here to isolate Season, a
scope decision, not an oversight.

**What it revealed:**
- A real bug in `run_sbc()`'s default `contrast_fn` (assumed a 2-column
  `mean`, broke on this model's 4-column one), caught and fixed.
- SBC showed mild underconfidence (posteriors a bit too wide), the
  conservative failure mode.
- The real-data PPC showed a shape mismatch: Conventional-July's observed
  density is tightly clustered near a low value, but the posterior
  predictive spread was visibly wider. This became the seed for Model 5.
- `sigma_tr` and `sigma[Mg]` showed a moderate (~-0.39) negative
  correlation, an expected residual-vs-random-effect trade-off given only
  ~2 observations/tree.

## Model 5 — `5_ITS_lognormal_MDLS2.R` (MDLS2)

**Structure:** `sigma[Mg]` (2 cells) → `sigma[cell]` (4 cells, Mg x Mo),
via a combined index `cell = (Mg-1)*2 + Mo`. Everything else identical to
Model 4.

**Motivation:** Tests directly whether Model 4's Conventional-July
mismatch was a sigma-sharing artifact: one residual SD per Management
group forces a single spread to fit two seasons that plausibly differ.

**What it took to get there:** `ulam()` can't double-index a
`dexp`-distributed scale parameter (`dnorm` varying effects work this way,
`dexp` doesn't), so a combined index (`cell`) is used instead. First model
built under the current per-model-code `hiermod/Models/CODE_model.R` +
`X.2_calibration`/`X.3_validation`/`X.4_analysis` file split.

**What it revealed:** Real-fit PPC shows the Conventional-July `yrep`
hugging the observed spike much more closely than Model 4's did. SBC's KS
test looks good (D=0.13, p=0.67 at n_sbc=30), but divergences across
replicates are non-trivial, plausibly because each `sigma[cell]` is backed
by roughly half the data `sigma[Mg]` had.

## Model 5, shifted — `5b.2_MDLS2_shifted_calibration.R`

**Structure:** Same `model` as plain Model 5, unmodified. Only data
preparation changes: real `Dv` is fit as `Hill_1 - 1`, and
`sim_div_MDLS2()`/`means_MDLS2()`/`simulate_from_priors()` gained an
optional `shift=` argument (default `NULL`/`0`, so plain Model 5 is
untouched).

**Motivation:** Hill numbers can't go below 1, but a plain lognormal has
support `(0, Inf)`. A small, honest correction for the response's true
floor.

**What it took to get there:** `Dv - 1 ~ dlnorm(...)` doesn't compile in
`ulam()` (needs a precomputed column). An early version fit the
likelihood directly to the shifted quantity instead of the underlying
floor-0 draw, which biased `sigma[cell]` low (0.43 vs true 0.56) and
caused an 8% divergence / 33% max-treedepth blowup. Fixed by never
changing `model$likelihood`: the shift belongs only in the backtransform
and in how real data is prepared, never in the likelihood's own declared
outcome. This convention carries forward to every model after Model 5
(see Model 6's note below, where it was initially missed).

**Status:** Parameter recovery smoke-tested. Real prior predictive check,
SBC, model fit and PPC not yet run (see `TODO.md`).

## Model 6 — `hiermod/Models/MDLSY_model.R` (MDLSY, "The Varietal")

**Structure:** Patches Model 5. `mu <- loga[Mg] + gamma*(Mo-1) +
b[Lo]*sigma_loc + yr[Yr]*sigma_yr + tr[Tr]*sigma_tr + cv[Cv]`. Two
additions over Model 5:
- `yr[Yr]` pooled into `yr[Yr]*sigma_yr` (non-centered): Year moves from
  fixed/unpooled to partially pooled, `sigma_yr ~ dexp(2)`.
- `cv[Cv] ~ dnorm(0,1)`: Cultivar (5 levels), fixed/unpooled.

Weather and sequencing-depth covariates were part of this model
originally, then split out into Model 7 once it became clear they answer
a different question (controlling for nuisance covariates) than Year/
Cultivar do (structural additions). See Model 7 below.

**Motivation:** Partial pooling gives Year proper shrinkage and honest
between-year uncertainty instead of an unregularized fixed estimate per
level. Cultivar was worth having in the posterior record even under a
null-effect expectation, rather than assumed away.

**What it took to get there:**
- Cultivar uses `make_index()`'s default levels, same as `Lo`/`Tr`, not an
  explicit override (an earlier 4-blank-placeholder override silently
  broke every `to_index()` call).
- `yr[Yr]*sigma_yr` needed no new `ulam()` idiom: structurally identical
  to `b[Lo]*sigma_loc`.
- Initially fit on raw `Hill_1` (unshifted), missing that Model 5-shifted's
  floor-1 convention was meant to carry forward. Fixed the same way:
  `dat$Dv` built as `div$Hill_1 - 1`, `shift = 1` threaded through
  simulation and backtransform calls, likelihood untouched.

**Variance budget calibration (VBC):** first model where the number of
summed `total_var` terms genuinely grows (K=3 → K=4, adding `sigma_yr^2`),
unlike Model 5 where K didn't change. See `GLOSSARY.md` for the R2D2M2
background and what `scale_dexp_rate()` does and doesn't cover.

**Status:** Code fit and reviewed once with weather covariates and the
unshifted response both still present; both have since changed (weather
moved to Model 7, response now shifted), so the real fit/SBC/PPC numbers
need regenerating (see `TODO.md`). Durable findings likely to hold across
the refit: Year's partially-pooled posterior with only 3 groups came back
appropriately wide (`sigma_yr` weakly identified by so few units, while
each year's own offset is well pinned down by its own data); per-year
effect panels show different skew shapes as a structural byproduct of
non-centering (see "skew transfer" in `GLOSSARY.md`), which prompted
adding a dashed posterior-mean line to `variance_component_panels()`.
Also confirmed `bayesplot::mcmc_pairs()` (1.15.0 / ggplot2 4.0.3) breaks
with >2 chains, worked around via `thin_for_pairs()` (`hiermod_core.R`).

## Model 7 — `hiermod/Models/MDLSYC_model.R` (MDLSYC, "The Weatherman")

**Structure:** Patches Model 6. Adds three standardized control
covariates as additive fixed slopes: `b_deg*deg_h_z + b_precip*precip_72h_z
+ b_seq*seq_depth_z`, each `~ dnorm(0,1)`.

**Motivation:** Separates Model 6's structural additions (Year pooling,
Cultivar) from covariates whose only job is controlling for nuisance
variation, so each question gets its own model. `deg_h`/`precip_72h`
control for day-to-day weather instead of letting it hide inside Season/
Year's coarser categorical structure. `Seq_depth` (raw pre-rarefaction
read count per sample) is new here: it's a partial proxy for detection
effort, since a sample with more reads is more likely to have discovered
its true ASV diversity before rarefaction. It doesn't correct for the
opposite confound (a low-true-diversity sample is easier to fully
discover with fewer reads), but is judged better than omitting depth
entirely, given DADA2-inferred ASVs can't support coverage-based
rarefaction the way OTUs with intact singleton counts can.

**`Seq_depth` treatment:** logged, then standardized (`log_seq_depth_z`
in `0_SETUP.R`), not a raw z-score. Raw `Seq_depth` spans ~3.1K-61K
reads, right-skewed, and the thing it's meant to proxy (how completely a
sample's diversity was discovered) is a saturating relationship, the
shape of a rarefaction curve, which a linear effect on the log scale
matches far better than a linear effect on raw counts.

**Variance budget calibration:** K stays at 4 here. The three control
covariates are additive `mu`-level fixed effects, the same category
`cv[Cv]` already was in Model 6, not new summed variance terms. Rather
than recomputing a calibration pass in the validation script (Model 6's
own pattern), Model 6's calibrated rates are hardcoded directly into
`hiermod/Models/MDLSYC_model.R`'s priors (`sigma ~ dexp(3.46)`, `sigma_loc`/
`sigma_tr`/`sigma_yr ~ dexp(2.31)`), inherited directly from `model_MDLSY_ITS`
rather than resourced from `6.2_MDLSY_calibration.R`, where they only
existed as a local variable before being promoted. This also let
`7.2_MDLSYC_calibration.R` drop the separate uncalibrated-vs-
calibrated fit comparison Model 6's script has: with the prior already
calibrated from the start, there's nothing to compare against.

**Tooling:** `contrast_recovery()` (`postcontrast_helpers.R`) was
extracted here and retrofitted onto Models 4, 5, 5-shifted, and 6's
validation scripts: the May-gap/July-gap/seasonal-change posterior +
true-value computation was identical boilerplate in every one of them.

**Known limitation:** all three covariates are likely correlated with
Season and/or Year in the real data (already confirmed for weather in an
earlier Model 6 fit: `deg_h` r=-0.73 with Season). SBC validates
recoverability against synthetic data where the three are simulated
independently, not identifiability under the real correlation structure.
Widens and correlates posteriors; doesn't bias the joint fit or any
estimand that isn't isolating one coefficient alone.

**Status:** Model code written and smoke-tested (recovery/means/variance-
partition functions all verified to run). Real fit, SBC and PPC not yet
run (see `TODO.md`).
