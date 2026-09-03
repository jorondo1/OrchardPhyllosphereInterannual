# Model history — ITS diversity hierarchical models

Retrospective trace of *why* each model exists, in file order. This is not
a to-do list (see `TODO.md` for open items) and not a methods section --
it's a record of the reasoning trail, so a later reader (including us) can
see that each layer of complexity answered a specific question rather than
being added to chase a better-looking fit. See also the "Discussion /
writeup notes" section of `TODO.md` for the standing question of how to
present this sequence honestly (robustness-across-specifications, not just
the final model's numbers).

## Why each variable is modeled the way it is

Quick reference for *why* each effect got the treatment it did (fixed vs.
partially-pooled vs. hand-built interaction vs. not modeled at all) --
collected here since each script only justifies its own choice, not the
full set at once.

- **Management (`loga[Mg]`)** -- fixed, unpooled (2 levels, no wider
  population to pool across): it's the actual comparison of interest, not
  something to shrink toward a common mean.
- **Location (`b[Lo]*sigma_loc`)** -- partially pooled: Locations are a
  sample from a wider population of orchards with uneven sample sizes, and
  reporting a "mean diversity" at all requires marginalizing over
  Location-level variation (see `total_var` in every `means_*()`).
- **Tree (`tr[Tr]*sigma_tr`)** -- partially pooled, added at Model 4: the
  same trees are sampled in both May and July, so without a Tree term those
  two rows get (wrongly) treated as independent, understating uncertainty
  (pseudoreplication).
- **Year (`yr[Yr]`)** -- fixed, unpooled, unlike Location/Tree: only 3
  years, sampled very unevenly (many Locations missing from 2022) -- too
  few and too unbalanced to estimate a population-level `sigma_yr`
  usefully. Treated as a control variable, not a random sample from "the
  population of years."
- **Season (`gamma <- s_conv + gap_shift*(Mg-1)`)** -- not a random effect:
  only 2 levels (May/July), each individually meaningful rather than
  exchangeable draws from a season population, and the May->July *change*
  is itself the estimand. McElreath's two-equation form makes `gap_shift`
  the interaction directly, instead of deriving it after fitting separate
  main effects.
- **`sigma[Mg]` / `sigma[cell]`** -- residual (within-cell) spread, kept
  separate from Location/Tree variance and allowed to differ by group:
  Model 2's simulation demo and Model 4->5's postpred mismatch both showed
  a single shared residual SD visibly distorts the fit once groups actually
  differ in spread.
- **`g[Lo]*sigma_g`** (Model 3 only) -- Location-varying Management gap,
  confirmed real (`sigma_g` > 0) but dropped from Model 4 onward as a
  deliberate scope decision (isolate Season first), not because it stopped
  mattering.
- **Non-centered form everywhere** (`b[Lo]*sigma_loc`, `tr[Tr]*sigma_tr`,
  `g[Lo]*sigma_g`) -- standard HMC efficiency fix: the centered form
  (`b[Lo] ~ dnorm(0, sigma_loc)`) creates a sampling funnel when the
  hyperparameter is small; multiplying a `dnorm(0,1)` draw by it avoids
  that.
- **The Model-5-shifted floor (`Hill_1 - 1`)** -- not an effect-structure
  choice, a support-matching one: Hill numbers can't go below 1, plain
  lognormal support is `(0, Inf)`.

## Model 1 -- `1_ITS_lognormal_MD.R` (MD / MDv)

**Structure:** `Dv ~ dlnorm(mu, sigma)`, `mu <- loga[Mg]`. Two variants:
MD (one shared `sigma`), then MDv (`sigma[Mg]`, one per Management group).

**Motivation:** The baseline question -- does Hill diversity differ at all
between Conventional and Organic, and is a single shared residual spread
even a safe assumption? MDv exists because the constant-variance model
visibly distorted the recovered mean contrast when simulated data had
unequal group variances (demonstrated in-script before fitting the real
data) -- i.e. it was motivated by a *known failure mode*, not a fishing
expedition.

**Status:** No SBC section -- both models are simple enough (no pooling,
no hierarchical structure) that a full calibration check wasn't judged
necessary.

## Model 2 -- `2_ITS_lognormal_MDL.R` (MDL / MDLb)

**Structure:** + `b[Lo]*sigma_loc` (non-centered, partially-pooled
Location random intercept).

**Motivation:** Observations come from multiple orchards/locations;
treating them as exchangeable ignores real structure in the sampling
design. This is the first model where the *fixed-effect estimand* really
requires marginalizing over a population of locations to report a
sensible "mean diversity" -- see `lognormal_mean()`'s `total_var`.

**What it revealed:** The exploratory (MDL) priors implied an implausibly
heavy tail in the prior predictive check; MDLb tightened `sigma[Mg]`
(dexp(1)->dexp(3)) and `sigma_loc` (dexp(1)->dexp(2)) in response, and is
the version actually fit to real data. First model with a real SBC section
(`run_sbc()`/`summarize_sbc()`), establishing the pattern every later model
reuses.

## Model 3 -- `3_ITS_lognormal_MDLv.R` (MDLv)

**Structure:** + `g[Lo]*sigma_g*(Mg-1)` -- the Organic-Conventional gap
itself varies by Location (non-centered), on top of Model 2's shared
Location intercept.

**Motivation:** At each Location, Management maps to one specific orchard
pairing (Compton = PMB/Conventional + COM/Organic, etc.), so "does the
orchard matter beyond its Location" and "does the gap differ by Location"
are the same question here -- worth testing directly rather than assuming
Model 2's shared-gap simplification was fine.

**What it revealed:** `sigma_g` sits comfortably above 0 -- the gap really
does vary by Location, not just by chance. Saint-Benoit and Windsor
(Conventional-only / Organic-only) flagged as a real limitation: `g[Lo]`
there is entangled with `b[Lo]`, not informative on its own.

## Model 4 -- `4_ITS_lognormal_MDLS.R` (MDLS)

**Structure:** Refocused on Season instead of Location-varying gap:
`loga[Mg] + gamma*(Mo-1) + b[Lo]*sigma_loc + yr[Yr] + tr[Tr]*sigma_tr`,
`gamma <- s_conv + gap_shift*(Mg-1)` (McElreath-style two-equation
interaction). Adds `yr[Yr]` (fixed, unpooled) and `tr[Tr]*sigma_tr`
(Tree-level random intercept).

**Motivation:** Same trees were sampled in both May and July within a
year -- treating those rows as independent understates uncertainty
(pseudoreplication), a materially worse problem than the nuisance-variance
questions Models 2-3 were addressing. Year entered as a fixed effect
because sampling was unbalanced across years/locations (many locations not
sampled in 2022). Model 3's Location-varying gap was deliberately dropped
here ("add back later") to isolate the Season question -- documented as a
scope decision, not an oversight.

**What it revealed:**
- A real bug in `run_sbc()`'s default `contrast_fn` (assumed a 2-column
  `mean`, silently broke on this model's 4-column one) -- caught and fixed
  before trusting any SBC result on this model.
- SBC (n_sbc=100) showed mild *underconfidence* (posteriors a bit too wide)
  -- the conservative failure mode, not overconfidence.
- The real-data PPC showed a real shape mismatch: Conventional-July's
  observed density is tightly clustered near a low value, but the
  posterior predictive spread for that cell was visibly wider -- the
  seed for Model 5.
- `sigma_tr` and `sigma[Mg]` showed a moderate (~-0.39) negative
  correlation in the posterior -- an expected residual-vs-random-effect
  trade-off given only ~2 observations/tree, not a bug.

## Model 5 -- `5_ITS_lognormal_MDLS2.R` (MDLS2)

**Structure:** `sigma[Mg]` (2 cells) -> `sigma[cell]` (4 cells, Mg x Mo),
via a combined index `cell = (Mg-1)*2 + Mo`. Everything else identical to
Model 4.

**Motivation:** Directly tests whether Model 4's Conventional-July
mismatch was a sigma-sharing artifact: one residual SD per Management
group forces a single spread to fit two seasons that plausibly have
different true variance.

**What it took to get there:** `ulam()` can't double-index a
`dexp`-distributed scale parameter -- verified empirically (works for
`dnorm`-type varying effects via `matrix[i,j]:a`, not for this); used the
combined-index workaround instead. This was also the first model built
under the new `models/*.R` file-per-model organization (retrofitted onto
Models 1-4 too), since the delta from Model 4 was genuinely small (2
changed lines in the `alist`, 3 functions needing a 4-cell-aware rewrite).

**What it's revealed so far:** Real-fit PPC suggests the Conventional-July
`yrep` now hugs the observed spike much more closely than Model 4's did
(the thing this model set out to test). SBC's KS test looks good (D=0.13,
p=0.67 at n_sbc=30), but total divergences across replicates (~450/replicate)
are non-trivial -- plausibly because each `sigma[cell]` is now backed by
roughly half the data `sigma[Mg]` had, making the variance components
harder for the sampler to cleanly separate. Still open: the real SBC/PPC
run at full replicate count, and the Tree-level PPC/SBC budget question
inherited unresolved from Model 4 (see `TODO.md`).

## Model 5, shifted -- `5_ITS_lognormal_MDLS2_shifted.R`

**Structure:** Not a new Stan model -- `model` is used completely
unmodified from `models/model_5_MDLS2.R`. Only how the data is prepared
changes: real `Dv` is fit as `Hill_1 - 1`, and `sim_div_MDLS2()`/
`means_MDLS2()`/`simulate_from_priors()` all gained an optional `shift=`
argument (default `NULL`/`0`, so plain Model 5 is untouched) so one model
file serves both variants.

**Motivation:** Hill numbers can't go below 1, but a plain lognormal has
support `(0, Inf)` -- a small, honest correction for the response's true
floor, worth checking as a robustness variant against plain Model 5's
numbers rather than assuming it doesn't matter.

**What it took to get there, the hard way:**
- `Dv - 1 ~ dlnorm(...)` (an expression on the likelihood's left-hand
  side) doesn't compile in `ulam()` -- confirmed by trying it. Needs a
  precomputed, plain-symbol column instead.
- An early version fit the likelihood directly to the *shifted* quantity
  (`Dv_shifted ~ dlnorm(...)`, i.e. the floor-1 value) instead of the
  underlying floor-0 draw. This produced an 8% divergence / 33%
  max-treedepth / all-chains-bad-E-BFMI blowup in Parameter recovery.
  Side-by-side simulation test (same seed, same true params) confirmed
  this was a real, biased misspecification -- recovered `sigma[cell]`
  values were systematically low (e.g. 0.43 vs true 0.56) -- not just a
  tuning problem. Fixed by never changing `model$likelihood` at all: the
  shift belongs only in `means_MDLS2()`'s backtransform and in how real
  data is subtracted down to floor-0 before fitting, never in the
  likelihood's own declared outcome.
- `post_full()` (postcontrast_helpers.R) gained a `...` passthrough so
  `shift=1` can reach `means_fn` through it, for the raw mean/median
  panels (contrasts never needed it -- a shared additive constant cancels
  in any group1-group2 difference).

**Status:** Parameter recovery smoke-tested at reduced scale, divergences
back down to a normal range for this model family. Real Prior predictive
check / SBC / Model fit / PPC not yet run (see `TODO.md`).
