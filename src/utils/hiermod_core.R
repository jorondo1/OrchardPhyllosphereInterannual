# hiermod_core.R -- shared hiermod helpers: category codebooks, lognormal
# conversions, recovery checks, pairs plots.

# Category <-> index codebook for one grouping variable
# - one place for the string -> integer mapping (avoids flipped indices across scripts)
#   $to_index(x)    category values -> integer index (for ulam data)
#   $to_label(i)    integer index   -> factor with original labels
#   $to_label_n(i)  integer index   -> factor "<label> (n=<count>)"
#   $palette        named colour vector
#   $palette_n()    same colours, keyed to to_label_n() labels
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

# Lognormal sigma -> raw-scale CV
sigma_to_cv <- function(sigma) {sqrt(exp(sigma^2) - 1)}

# Raw-scale CV -> lognormal sigma (inverse of the above)
cv_to_sigma <- function(cv) {sqrt(log(1 + cv^2))}

# Mean of lognormal(mu, total_var), optionally shifted
# - total_var = whichever sum of sigma^2 terms the model has
lognormal_mean <- function(mu, total_var, shift = 0) {shift + exp(mu + total_var / 2)}

# Parameter recovery: posterior median + 89% PI per parameter, true value covered?
# - true and post_draws on the same scale
check_recovery <- function(true, post_draws){
  purrr::imap_dfr(post_draws, function(draws, name) tibble(
    param = name,
    true = true[[name]],
    post_median = median(draws),
    lo89 = PI(draws)[1],
    hi89 = PI(draws)[2]
  )) %>% mutate(covered = true >= lo89 & true <= hi89)
}

# Rescale an Exponential prior rate so E[sum of sigma^2] stays fixed when
# more variance components are summed (R2D2M2-style budget, see GLOSSARY.md)
scale_dexp_rate <- function(rate_ref, k_ref, k_new) round(rate_ref * sqrt(k_new / k_ref), 2)

# Subsample draws for mcmc_pairs(): <= 2 chains, n_keep draws total
# - keeps file size down; bayesplot errors beyond 2 chains
# - random subsample (divergences not preferentially kept)
# - cs = attr(fit, "cstanfit"); returns list(draws, np)
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

# Full pairwise posterior grid, divergences highlighted
# - variables: exact Stan names, e.g. "loga[1]", "sigma[2]", "cv_1"
# - save as PNG (save_gg(..., type = "png")): PDFs of many points get huge
#   p <- plot_mcmc_pairs(fit_sim, variables = c("loga[1]", "loga[2]", ...))
#   save_gg("mcmc_pairs", model_id, p, width = 14, height = 14, type = "png")
plot_mcmc_pairs <- function(fit, variables, n_keep = 500, max_chains = 2,
                             point_size = 0.3, point_alpha = 0.3){
  cs <- attr(fit, "cstanfit")
  thinned <- thin_for_pairs(cs, variables = variables, n_keep = n_keep, max_chains = max_chains)
  bayesplot::mcmc_pairs(
    thinned$draws, np = thinned$np,
    off_diag_args = list(size = point_size, alpha = point_alpha))
}

# "Management Season" group label for PPC panels; needs `idx` in scope (0_INDEX.R)
ppc_group <- function(dat){
  interaction(idx$Mg$to_label(dat$Mg), idx$Mo$to_label(dat$Mo), sep = " ")
}
