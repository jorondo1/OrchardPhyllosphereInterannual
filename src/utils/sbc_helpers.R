# sbc_helpers.R -- simulation-based calibration (SBC): repeatedly simulate
# from the priors, refit, and check whether the true value's rank in the
# resulting posterior is uniformly distributed across replicates.

# Pulls the i-th prior draw out of extract.prior()'s output as a plain
# named list, for feeding to a model's sim_div_*() as one "true" parameter
# set (an indexed parameter comes back as a matrix row, a scalar as one
# vector element).
draw_true <- function(priors, i){
  purrr::imap(priors, function(x, name){
    x <- as.matrix(x)
    if (ncol(x) == 1) x[i, 1] else x[i, ]
  })
}

# Runs the full SBC loop: simulate data from prior draw i, refit the
# model, and record where the true (simulating) value ranks in the
# resulting posterior, for n_sbc replicates.
#   model_fit   : a compiled ulam fit -- supplies both the model spec and
#                 the prior draws (extract.prior()), so they can't drift apart.
#   simulate_fn : function(true_params) -> data.frame matching the model's data
#   means_fn    : function(post) -> list(mean = cbind(v1, v2), ...)
#   contrast_fn : function(post, true_params, means_fn) -> list(post_contrast,
#                 true_contrast); defaults to a 2-column mean, must be
#                 supplied when a model's estimand doesn't reduce to that
#   extra_estimands : optional named list of function(post, true_params) ->
#                 list(post_draws, true_value) -- one entry per additional
#                 parameter to rank-check beyond the main contrast (e.g.
#                 list(sigma_loc = function(post, tp) list(post_draws = post$sigma_loc,
#                 true_value = tp$sigma_loc))). A single derived contrast can
#                 be well-calibrated while another parameter isn't -- this
#                 checks as many as you like inside the same (expensive)
#                 refit loop instead of needing a separate SBC run per
#                 parameter. Defaults to none, so every existing caller's
#                 behavior is unchanged.
#   n_parallel  : replicates run concurrently via parallel::mclapply()
#                 (safe with cmdstanr's subprocess-based chains)
# Returns a list of per-replicate results, for save_sbc_report().
run_sbc <- function(
    model_fit, simulate_fn, means_fn, contrast_fn = contrast_from_means,
    extra_estimands = list(),
    n_sbc = 8, iter = 1000, chains = 2, cores = chains, refresh = 0,
    control = NULL, n_parallel = 1) {

  model  <- model_fit@formula
  message('Extracting parameters from priors...')
  priors <- extract.prior(model_fit, n = n_sbc)

  one_replicate <- function(i){
    true_params <- draw_true(priors, i)
    dat_list    <- as.list(simulate_fn(true_params))

    fit <- if (is.null(control)) {
      ulam(model, data = dat_list, chains = chains, cores = cores,
           iter = iter, refresh = refresh)
    } else {
      ulam(model, data = dat_list, chains = chains, cores = cores,
           iter = iter, refresh = refresh, control = control)
    }
    ndiv       <- sum(fit@cstanfit$diagnostic_summary(diagnostics = "divergences")$num_divergent)
    # E-BFMI is a per-chain diagnostic, separate from divergences -- a chain
    # can have zero divergences and still have poor energy exploration
    # (Betancourt's <0.3 rule of thumb), which is what originally motivated
    # this model (MDLS2 got very few divergences but one low-E-BFMI chain).
    n_low_ebfmi <- sum(fit@cstanfit$diagnostic_summary(diagnostics = "ebfmi")$ebfmi < 0.3)

    post <- extract.samples(fit)
    cres <- contrast_fn(post, true_params, means_fn)

    # rank = how many posterior draws fall below the true value; should be
    # uniform across replicates if the model is calibrated.
    extra_ranks <- lapply(extra_estimands, function(f){
      e <- f(post, true_params)
      list(rank = sum(e$post_draws < e$true_value), Ns = length(e$post_draws))
    })

    res <- list(
      rank          = sum(cres$post_contrast < cres$true_contrast),
      Ns            = length(cres$post_contrast),
      true_contrast = cres$true_contrast,
      extra_ranks   = extra_ranks,
      n_divergent   = ndiv,
      n_transitions = (iter %/% 2) * chains,
      n_low_ebfmi   = n_low_ebfmi,
      n_chains      = chains
    )
    cat(i, "done, divergences:", ndiv, ", low-E-BFMI chains:", n_low_ebfmi, "\n")
    res
  }

  # Bounded retry: a rare compile-cache race can fail a replicate's launch.
  safe_replicate <- function(i, max_tries = 3){
    for (attempt in seq_len(max_tries)) {
      res <- tryCatch(one_replicate(i), error = function(e) e)
      if (!inherits(res, "error")) return(res)
    }
    stop("replicate ", i, " failed ", max_tries, "x: ", conditionMessage(res))
  }

  if (n_parallel > 1) {
    message('Compiling once before forking...')
    first <- one_replicate(1)
    rest  <- if (n_sbc > 1) parallel::mclapply(2:n_sbc, safe_replicate, mc.cores = n_parallel) else list()
    c(list(first), rest)
  } else {
    lapply(seq_len(n_sbc), one_replicate)
  }
}

# Summarizes a run_sbc() result: divergence/E-BFMI rates, plus three
# complementary rank-uniformity diagnostics -- one for the main contrast
# and, when run_sbc() was given extra_estimands, one per extra parameter
# (faceted together), since one derived contrast being calibrated doesn't
# mean every parameter is:
#
#   - rank histogram: the direct, intuitive view; catches gross violations
#     at a glance, but its apparent shape depends on bin count/placement.
#   - ECDF difference (Sailynoja, Burkner & Vehtari 2022): checks
#     uniformity everywhere at once via a simultaneous (DKW) 95% band,
#     preferred over a plain ECDF or a single global KS p-value because it
#     stays informative at lower replicate counts and localizes *where*
#     (which rank region) a violation happens.
#   - central-interval coverage: the practically-relevant translation of
#     the same ranks -- "does my reported 50%/89% interval actually cover
#     the truth that often" -- which the SBC package's own docs note
#     should always accompany, not replace, the shape-based checks above,
#     since coverage alone can't say whether a gap is under- or
#     over-confidence in a particular tail.
#
# (Plain ECDF and coverage-difference are skipped: the SBC package's own
# guidance is that ECDF-diff is preferable to plain ECDF "in most cases",
# and coverage's own diagonal reference line already shows the deviation
# direction plain coverage-diff would add.)
save_sbc_report <- function(sbc_out, step, dir = hiermod_out_dir){
  ranks         <- sapply(sbc_out, function(x) x$rank / x$Ns)
  ks_test       <- ks.test(ranks, "punif")
  n_divergent   <- sum(sapply(sbc_out, function(x) x$n_divergent))
  n_transitions <- sum(sapply(sbc_out, function(x) x$n_transitions))
  div_rate      <- 100 * n_divergent / n_transitions
  n_low_ebfmi   <- sum(sapply(sbc_out, function(x) x$n_low_ebfmi))
  n_chains      <- sum(sapply(sbc_out, function(x) x$n_chains))
  ebfmi_rate    <- 100 * n_low_ebfmi / n_chains

  label_text <- paste0(
    "n = ", length(ranks), " replicates \n",
    "KS D = ", round(unname(ks_test$statistic), 3), ",  p = ", round(ks_test$p.value, 3), " \n",
    "divergences: ", format(n_divergent, big.mark = ","), " / ",
    format(n_transitions, big.mark = ","), "  (", round(div_rate, 2), "%) \n",
    "low E-BFMI chains: ", n_low_ebfmi, " / ", n_chains, "  (", round(ebfmi_rate, 2), "%)"
  )

  # Extra per-parameter estimands (see run_sbc()'s extra_estimands) -- empty
  # by default, so every plot below is single-panel for existing callers.
  extra_names <- names(sbc_out[[1]]$extra_ranks)
  extra_ranks <- if (length(extra_names) > 0) setNames(lapply(extra_names, function(nm){
    sapply(sbc_out, function(x) x$extra_ranks[[nm]]$rank / x$extra_ranks[[nm]]$Ns)
  }), extra_names) else NULL
  ranks_list <- c(list(contrast = ranks), extra_ranks)
  ks_tests   <- lapply(ranks_list, function(r) ks.test(r, "punif"))

  # ---- Rank histogram (faceted across all tracked estimands) -----------
  # Bin count chosen as a divisor of n_sim close to 10: the SBC package's
  # own rank-visualization guidance warns a bin count that doesn't divide
  # the number of simulations evenly can visually hide real problems (a
  # bins=17-on-100-simulations example), and more bins trades resolution
  # for less power to catch small-sample violations.
  n_sim      <- length(ranks)
  divisors   <- (1:n_sim)[n_sim %% (1:n_sim) == 0]
  candidates <- divisors[divisors >= 5 & divisors <= 20]
  # Falls back to a plain target of 10 when n_sim has no divisor in a
  # sensible range (e.g. a prime n_sbc) -- better an imperfect bin count
  # than the divisor search degenerating to bins=1.
  bins <- if (length(candidates) > 0) candidates[which.min(abs(candidates - 10))] else 10
  expected_per_bin <- n_sim / bins

  p_hist <- ggplot(
    bind_rows(lapply(names(ranks_list), function(nm) tibble(param = nm, rank = ranks_list[[nm]]))),
    aes(x = rank)) +
    geom_histogram(bins = bins, boundary = 0, fill = "#4E79A7", colour = "white", linewidth = 0.3) +
    geom_hline(yintercept = expected_per_bin, linetype = "dashed", colour = "grey40") +
    facet_wrap(~param) +
    labs(title = paste("SBC rank histogram --", step), subtitle = label_text,
         caption = "Dashed line: count expected under perfect calibration",
         x = "rank", y = "count")
  save_gg("SBC_rank_hist", step, p_hist, width = 8, height = 6, dir = dir)

  # ---- ECDF difference, with a simultaneous (DKW) 95% band --------------
  grid <- seq(0, 1, length.out = 200)
  p_ecdf <- ggplot(
    bind_rows(lapply(names(ranks_list), function(nm){
      u <- ranks_list[[nm]]; eps <- sqrt(log(2 / 0.05) / (2 * length(u)))
      tibble(param = nm, x = grid, diff = sapply(grid, function(x) mean(u <= x) - x), band = eps)
    })),
    aes(x = x)) +
    geom_ribbon(aes(ymin = -band, ymax = band), fill = "grey70", alpha = 0.4) +
    geom_hline(yintercept = 0, colour = "grey40", linetype = "dashed") +
    geom_step(aes(y = diff), colour = "#4E79A7", linewidth = 0.6) +
    facet_wrap(~param) +
    labs(title = "SBC rank ECDF difference (95% simultaneous band)",
         subtitle = "Outside the band anywhere = miscalibrated there, not just on average",
         x = "normalized rank", y = "ECDF(rank) - rank")
  save_gg("SBC_ecdf_diff", step, p_ecdf, width = 8, height = 6, dir = dir)

  # ---- Central-interval coverage -----------------------------------------
  # Same normalized ranks, reframed as "does the width-w central interval
  # actually contain the truth w of the time" -- band reuses the ECDF-diff
  # band width as an approximate, slightly conservative visual guide (not
  # a separately-derived exact coverage-band).
  cov_grid <- seq(0, 1, length.out = 100)
  p_cov <- ggplot(
    bind_rows(lapply(names(ranks_list), function(nm){
      u <- ranks_list[[nm]]; eps <- sqrt(log(2 / 0.05) / (2 * length(u)))
      tibble(param = nm, width = cov_grid,
             coverage = sapply(cov_grid, function(w) mean(abs(u - 0.5) < w / 2)), band = eps)
    })),
    aes(x = width)) +
    geom_ribbon(aes(ymin = pmax(width - band, 0), ymax = pmin(width + band, 1)),
                fill = "grey70", alpha = 0.4) +
    geom_abline(slope = 1, intercept = 0, colour = "grey40", linetype = "dashed") +
    geom_line(aes(y = coverage), colour = "#4E79A7", linewidth = 0.6) +
    facet_wrap(~param) +
    labs(title = "SBC central-interval coverage",
         subtitle = "Above the diagonal = wider/more conservative than nominal; below = narrower/overconfident",
         x = "nominal central interval width", y = "empirical coverage")
  save_gg("SBC_coverage", step, p_cov, width = 8, height = 6, dir = dir)

  invisible(list(ks_test = ks_test, n_divergent = n_divergent,
                 n_transitions = n_transitions, n_low_ebfmi = n_low_ebfmi,
                 n_chains = n_chains, ranks = ranks,
                 extra_ranks = extra_ranks, ks_tests = ks_tests))
}
