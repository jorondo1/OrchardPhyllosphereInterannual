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
