# hiermod_core.R -- shared low-level building blocks used by every hiermod
#
# - script and every other utils/*.R file: category<->index codebooks and the
# - lognormal mu/sigma <-> raw mean/CV conversions.


# Category <-> index codebook. make_index(x) derives levels from the column
# (or pass levels= to force an order, e.g. a reference level first).
#   $to_index(x)    category values -> integer index (for ulam() data lists)
#   $to_label(i)    integer index   -> factor with the original labels
#   $to_label_n(i)  integer index   -> factor with "<label> (n=<count>)"
#   $palette        named vector, colour per level (auto hue_pal() unless
#                   palette= is passed -- named or plain, either works)
#   $palette_n()    same colours, re-keyed to the to_label_n() labels (what
#                   a `group` column built from to_label_n() actually holds)
make_index <- function(x, levels = sort(unique(as.character(x))), palette = NULL){
  counts   <- as.integer(table(factor(x, levels = levels)))
  labels_n <- paste0(levels, " (n=", counts, ")")
  if (is.null(palette)) palette <- scales::hue_pal()(length(levels))
  if (is.null(names(palette))) names(palette) <- levels
  list(
    levels     = levels,
    counts     = counts,
    palette    = palette,
    to_index   = function(x) match(as.character(x), levels),
    to_label   = function(i) factor(levels[i], levels = levels),
    to_label_n = function(i) factor(labels_n[i], levels = labels_n),
    palette_n  = function() setNames(palette, labels_n)
  )
}

# raw-scale CV implied by a lognormal sigma
sigma_to_cv <- function(sigma) {sqrt(exp(sigma^2) - 1)}

# inverse of sigma_to_cv(): the lognormal sigma implied by a target raw-scale CV
cv_to_sigma <- function(cv) {sqrt(log(1 + cv^2))}

# E[exp(X)] = exp(mu + var/2) for X ~ Normal(mu, var): the one calculation
# every raw-scale mean in this workflow reduces to. total_var is whatever sum
# of sigma^2/sigma_loc^2/... applies for that model/group, built inline in
# each model's own means_*(post) function.

# shift (default 0) is for a shifted-lognorma response (model 5 onwards)
lognormal_mean <- function(mu, total_var, shift = 0) {shift + exp(mu + total_var / 2)}

# Fixed-effect recovery table for a Parameter recovery section (model 4
# onward, alongside the mean-scale estimand check via post_full()/
# estimand_rows()). true: named vector of true values (e.g. c(loga1 = ..,
# s_conv = .., gap_shift = ..)), on whatever scale post_draws is already on
# (e.g. log scale for loga/s_conv/gap_shift, natural scale for sigma-type
# parameters) -- true and post_draws must match scale, this does no
# backtransform of its own. post_draws: named list of posterior draw
# vectors, same names as true. One row per parameter, with the true value
# flagged in/out of its own 89% PI.
check_recovery <- function(true, post_draws){
  purrr::imap_dfr(post_draws, function(draws, name) tibble(
    param = name,
    true = true[[name]],
    post_median = median(draws),
    lo89 = PI(draws)[1],
    hi89 = PI(draws)[2]
  )) %>% mutate(covered = true >= lo89 & true <= hi89)
}

# Rate for sigma ~ dexp(rate) that holds E[total_var] = E[sum(sigma_i^2)]
# roughly fixed as more independent variance components get summed into it
# (e.g. means_MDLS2()'s sigma[cell]^2+sigma_loc^2+sigma_tr^2 growing an extra
# +sigma_yr^2 term in Model 6) -- without this, each added component just
# stacks its own independent chance of a large draw on top of the existing
# ones, and the lognormal mean (exp(mu+total_var/2)) is exponentially
# sensitive to that sum, so the prior predictive mean/SD creep up with every
# addition even though no single component's own prior changed.
#
# For X ~ Exponential(rate), E[X^2] = 2/rate^2. Holding sum(E[sigma_i^2])
# fixed across a k_ref -> k_new change in how many terms are summed means
# scaling every component's rate by sqrt(k_new/k_ref) -- this is the
# two-line version of the same idea behind R2D2M2-style variance-
# decomposition priors for multilevel models (Aguilar & Bürkner 2023,
# "Intuitive joint priors for Bayesian linear multilevel models: the
# R2D2M2 prior", Electronic Journal of Statistics -- implemented as
# brms::R2D2(); generalizes Zhang/Naughton/Bondell/Reich 2020's R2-D2
# shrinkage prior and Yanchenko/Bondell/Reich's GLMM extension): put a
# prior on the TOTAL variance you find plausible, then split it across
# components, instead of letting each new component add its own
# independent chunk on top. See also Gelman (2006, Bayesian Analysis)
# on hierarchical variance-parameter priors more generally.
#
# Doesn't address a separate, related effect: more VALUES drawn from one
# shared-form slot (e.g. sigma[cell]'s 4 cells vs sigma[Mg]'s 2) raises the
# chance that *some* cell draws large across a full prior-predictive
# ensemble, even though any single observation's own total_var still only
# sums the same number of terms -- that one is still best caught by the
# ordinary prior predictive check, not this formula.
scale_dexp_rate <- function(rate_ref, k_ref, k_new) rate_ref * sqrt(k_new / k_ref)

# Caps an mcmc_pairs() plot's draw count instead of asking ggplot to render
# tens of thousands of point objects (a vector PDF of all of them can hit
# 40MB+; points also overlap heavily well before that many, so nothing is
# visually lost by capping). Plain uniform random subsample, not a
# divergent-preserving one -- forcing every divergent draw in while capping
# the rest inflates their visual share relative to the bulk, distorting the
# exact ratio a pairs plot is usually read for. 10,000 draws (the default)
# comfortably still shows divergences at the rates seen so far in this
# project (1-9% -> 100-900 points).
#
# Also caps chains at 2 (max_chains): bayesplot::mcmc_pairs() (bayesplot
# 1.15.0 / ggplot2 4.0.3, confirmed) has a real bug in its default
# multi-chain conditioning -- with >2 chains, the divergence-colour overlay
# throws a ggplot2 aesthetics-length error every time, reproduced even on
# full, unthinned draws with no subsetting involved at all. Exactly 2
# chains is the one setting that reliably works, so this keeps only the
# first 2 -- still shows the same within-chain funnel/ridge geometry (a
# property of the posterior, not particular to which chain), just from
# fewer of the real fit's 6 chains at once. Revisit once bayesplot ships a
# fix; not something to work around further than this at the call site.
#
# Uses cs$draws()'s plain draws_array (iteration x chain x variable) and
# base `[` array indexing throughout, not posterior::subset_draws() on a
# draws_df -- the latter silently renumbers .iteration on the kept object,
# which breaks matching it back against nuts_params()'s original Chain/
# Iteration numbering (found the hard way: every np row got dropped).
# Plain array indexing has no such bookkeeping to go stale.
#
# cs: a fit's compiled cmdstanr object (attr(fit, "cstanfit")). Returns
# list(draws, np), matched to the same kept (chain, iteration) pairs --
# pass straight to bayesplot::mcmc_pairs(thinned$draws, np = thinned$np).
thin_for_pairs <- function(cs, variables, n_keep = 10000, max_chains = 2){
  arr <- cs$draws(variables = variables)
  np  <- bayesplot::nuts_params(cs)

  chain_idx <- seq_len(min(dim(arr)[2], max_chains))
  arr <- arr[, chain_idx, , drop = FALSE]

  n_iter    <- dim(arr)[1]
  per_chain <- min(n_iter, n_keep %/% length(chain_idx))
  keep_idx  <- sample(seq_len(n_iter), per_chain)
  arr       <- arr[keep_idx, , , drop = FALSE]

  keep_iter_labels <- as.integer(dimnames(arr)$iteration)
  list(
    draws = arr,
    np    = np[np$Chain %in% chain_idx & np$Iteration %in% keep_iter_labels, ]
  )
}
