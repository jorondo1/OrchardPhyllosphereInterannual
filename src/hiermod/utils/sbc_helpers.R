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
#   n_parallel  : replicates run concurrently via parallel::mclapply()
#                 (safe with cmdstanr's subprocess-based chains)
# Returns a list of per-replicate results, for save_sbc_report().
run_sbc <- function(
    model_fit, simulate_fn, means_fn, contrast_fn = contrast_from_means,
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
    ndiv <- sum(fit@cstanfit$diagnostic_summary(diagnostics = "divergences")$num_divergent)

    post <- extract.samples(fit)
    cres <- contrast_fn(post, true_params, means_fn)

    # rank = how many posterior draws fall below the true value; should be
    # uniform across replicates if the model is calibrated.
    res <- list(
      rank          = sum(cres$post_contrast < cres$true_contrast),
      Ns            = length(cres$post_contrast),
      true_contrast = cres$true_contrast,
      n_divergent   = ndiv,
      n_transitions = (iter %/% 2) * chains
    )
    cat(i, "done, divergences:", ndiv, "\n")
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

# Summarizes a run_sbc() result: KS test for rank uniformity, divergence
# rate (as a % of sampling transitions, not a raw count), and a
# rank-histogram plot.
save_sbc_report <- function(sbc_out, step, dir = hiermod_out_dir){
  ranks         <- sapply(sbc_out, function(x) x$rank / x$Ns)
  ks_test       <- ks.test(ranks, "punif")
  n_divergent   <- sum(sapply(sbc_out, function(x) x$n_divergent))
  n_transitions <- sum(sapply(sbc_out, function(x) x$n_transitions))
  div_rate      <- 100 * n_divergent / n_transitions
  expected_per_bin <- length(ranks) / 10

  label_text <- paste0(
    "n = ", length(ranks), " replicates \n",
    "KS D = ", round(unname(ks_test$statistic), 3), ",  p = ", round(ks_test$p.value, 3), " \n",
    "divergences: ", format(n_divergent, big.mark = ","), " / ",
    format(n_transitions, big.mark = ","), "  (", round(div_rate, 2), "%)"
  )

  p <- ggplot(tibble(rank = ranks), aes(x = rank)) +
    geom_histogram(bins = 10, boundary = 0, fill = "#4E79A7", colour = "white", linewidth = 0.3) +
    geom_hline(yintercept = expected_per_bin, linetype = "dashed", colour = "grey40") +
    annotate("label", x = Inf, y = Inf, hjust = 1.05, vjust = 1.3, label = label_text,
             family = "mono", size = 3, fill = "white", label.size = 0.3, colour = "grey20") +
    labs(title = paste("Simulation-based calibration rank histogram", step),
         caption = "Dashed line: count expected under perfect calibration",
         x = "rank", y = "count")

  save_gg("SBC_rank_hist", step, p, width = 7, height = 5, dir = dir)
  invisible(list(ks_test = ks_test, n_divergent = n_divergent,
                 n_transitions = n_transitions, ranks = ranks))
}
