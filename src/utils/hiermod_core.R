# hiermod_core.R -- category<->index codebooks and lognormal mean/CV
# conversions shared by every hiermod script.

# Builds a category<->index codebook for one grouping variable (Management,
# Location, ...): to_index() for ulam() data lists, to_label()/to_label_n()
# for readable labels in plots, and a colour palette. Needed so every model
# script converts the same string columns the same way, without repeating
# the match()/factor() logic (or drifting: this once caused a flipped
# Conventional/Organic index across scripts).
#   $to_index(x)    category values -> integer index
#   $to_label(i)    integer index   -> factor with the original labels
#   $to_label_n(i)  integer index   -> factor with "<label> (n=<count>)"
#   $palette        named vector, colour per level
#   $palette_n()    same colours, re-keyed to the to_label_n() labels
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

# Raw-scale CV implied by a lognormal sigma.
sigma_to_cv <- function(sigma) {sqrt(exp(sigma^2) - 1)}

# Inverse of sigma_to_cv(): the lognormal sigma implied by a target raw-scale CV.
cv_to_sigma <- function(cv) {sqrt(log(1 + cv^2))}

# Mean of a lognormal(mu, total_var), optionally shifted. Every model's
# means_*() backtransform reduces to this one formula; total_var is
# whichever sum of sigma^2/sigma_loc^2/... applies, built inline per model.
lognormal_mean <- function(mu, total_var, shift = 0) {shift + exp(mu + total_var / 2)}

# Compares posterior draws against known true values (parameter recovery):
# one row per parameter with its 89% PI and whether the true value falls
# inside it. true and post_draws must already be on the same scale.
check_recovery <- function(true, post_draws){
  purrr::imap_dfr(post_draws, function(draws, name) tibble(
    param = name,
    true = true[[name]],
    post_median = median(draws),
    lo89 = PI(draws)[1],
    hi89 = PI(draws)[2]
  )) %>% mutate(covered = true >= lo89 & true <= hi89)
}

# Rescales an Exponential prior's rate so E[sigma^2] stays fixed when more
# independent variance components get summed together (R2D2M2-style budget
# calibration, see GLOSSARY.md). Needed because summing more components
# would otherwise inflate the total implied variance even though no single
# component's own prior changed. Doesn't address more *values* within one
# already-summed slot (e.g. sigma[cell]'s 4 cells vs. 2) -- that's a
# separate effect, best caught by the ordinary prior predictive check.
scale_dexp_rate <- function(rate_ref, k_ref, k_new) round(rate_ref * sqrt(k_new / k_ref), 2)

# Subsamples an mcmc_pairs() plot's draws (capped at 2 chains, n_keep total)
# so the output stays a reasonable size and sidesteps a bayesplot bug where
# mcmc_pairs() errors past 2 chains (bayesplot 1.15.0 / ggplot2 4.0.3,
# confirmed). Plain random subsample, not divergence-preserving.
# At project typical divergence rates (1-9%) that's still hundreds of
# divergent points shown. cs: a fit's compiled cmdstanr object
# (attr(fit, "cstanfit")). Returns list(draws, np) for
# bayesplot::mcmc_pairs(thinned$draws, np = thinned$np).
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
