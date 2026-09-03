#  simulation-based calibration harness. Runs fit/simulate/rank loop.


# One prior-draws replicate -> the i-th "true" parameter set, for whatever
# names extract.prior() returned (an indexed parameter like loga is a
# matrix, row i is its slice; a scalar like sigma_loc is a plain vector,
# element i is its slice) -- same shape logic post_full() already handles,
# so this needs no per-model rewrite: every script's Prior predictive
# check/SBC section can call it as-is on its own `prior <- extract.prior(...)`.
draw_true <- function(priors, i){
  purrr::imap(priors, function(x, name){
    x <- as.matrix(x)
    if (ncol(x) == 1) x[i, 1] else x[i, ]
  })
}

 
# SBC harness
# 99% Claude Code
# Model-specific pieces (simulate_fn/means_fn)  supplied as argument
#   model_fit   : a compiled ulam fit of the model to validate -- supplies
#                 both the model spec (re-fit each replicate) and the prior
#                 draws (extract.prior()), so the two can't drift apart.
#   simulate_fn : function(true_params) -> data.frame (via sim_div_*()), with
#                 columns already matching the model's data arguments.
#   means_fn    : function(post) -> list(mean = cbind(v1, v2), ...) -- same
#                 convention postcounts()/post_full() use for the real fit.
#   contrast_fn : function(post, true_params, means_fn) -> list(post_contrast,
#                 true_contrast) -- defaults to contrast_from_means(), which
#                 assumes `mean` is a 2-column matrix and tests column2-column1
#                 (models 1-3's single Mg contrast). A model whose `mean` has
#                 more columns (model 4's 4-cell MDLS) MUST supply its own --
#                 contrast_from_means() would silently grab the wrong pair of
#                 columns instead of erroring, so this isn't optional there.
#   control     : passed to ulam(); match the real fit's control= or SBC will
#                 see spurious divergences from a weaker sampler setting alone.
#   n_parallel  : replicates run at once via parallel::mclapply() (default 1
#                 = current sequential behaviour, unchanged). Safe with
#                 cmdstanr (each chain is an external subprocess, not
#                 in-process compiled code the way rstan is -- forking is a
#                 known hazard for the latter, not this). Budget
#                 n_parallel * cores against your core count: replicates are
#                 fully independent fits, so this is the lever that actually
#                 uses idle cores, unlike within-chain `threads` (reduce_sum),
#                 which only pays off with large per-replicate N -- worth
#                 more for the one real-data fit than for a small SBC
#                 replicate. cat()'s per-replicate progress line will arrive
#                 out of numeric order once n_parallel > 1 (forked stdout).
#
# Returns a list of per-iteration results (rank, Ns, true_contrast, n_divergent)
# -- pass to summarize_sbc() for a pass/fail read.

run_sbc <- function(model_fit, simulate_fn, means_fn, contrast_fn = contrast_from_means,
                    n_sbc = 8, iter = 1000, chains = 2, cores = chains, refresh = 0,
                    control = NULL, n_parallel = 1){

  model  <- model_fit@formula
  message('Extracting parameters from priors...')
  priors <- extract.prior(model_fit, n = n_sbc)

  one_replicate <- function(i){
    # 1. Simulate a dataset from one set of prior parameters
    true_params <- draw_true(priors, i)
    dat_list    <- as.list(simulate_fn(true_params))

    # 2. Fit the model to the fake data
    fit <- if (is.null(control)) {
      ulam(model, data = dat_list, chains = chains, cores = cores,
           iter = iter, refresh = refresh)
    } else {
      ulam(model, data = dat_list, chains = chains, cores = cores,
           iter = iter, refresh = refresh, control = control)
    }
    # cmdstanr backend (@cstanfit), not rstan (@stanfit)
    ndiv <- sum(fit@cstanfit$diagnostic_summary(diagnostics = "divergences")$num_divergent)

    # 3. Recover the fitted posterior and compare it against the true value
    # that actually generated this replicate's data.
    post <- extract.samples(fit)
    cres <- contrast_fn(post, true_params, means_fn)

    # 4. Rank = how many posterior draws fall below the true value. Across
    # n_sbc replicates this rank should be uniformly distributed if the model
    # is calibrated -- that's what summarize_sbc()'s ks.test checks.
    res <- list(
      rank          = sum(cres$post_contrast < cres$true_contrast),
      Ns            = length(cres$post_contrast),
      true_contrast = cres$true_contrast,
      n_divergent   = ndiv
    )
    cat(i, "done, divergences:", ndiv, "\n")
    res
  }

  # Retries a replicate on failure;  the compile-cache race below is rare
  # once primed (~1 in 5 *batches* still hit it in testing, not 1 in 5
  # replicates) but not eliminated; a bounded retry is cheap insurance
  # against a transient subprocess-launch race silently shrinking n_sbc.
  safe_replicate <- function(i, max_tries = 3){
    for (attempt in seq_len(max_tries)) {
      res <- tryCatch(one_replicate(i), error = function(e) e)
      if (!inherits(res, "error")) return(res)
    }
    stop("replicate ", i, " failed ", max_tries, "x: ", conditionMessage(res))
  }

  if (n_parallel > 1) {
    # ulam() recompiles fresh per replicate (new data -> same Stan code),
    # caching the binary to a tempdir() path all forked children inherit
    # unchanged. 
    message('Compiling once before forking...')
    first <- one_replicate(1)
    rest  <- if (n_sbc > 1) parallel::mclapply(2:n_sbc, safe_replicate, mc.cores = n_parallel) else list()
    c(list(first), rest)
  } else {
    lapply(seq_len(n_sbc), one_replicate)
  }
}

# Quick numeric read on a run_sbc() result: rank uniformity (should be
# roughly uniform -> high p-value) and total divergences (should be ~0). Fast
# smoke test with n_sbc in the 8-15 range, not a substitute for a full
# ~20-100-replicate check on a model you're about to commit to.
summarize_sbc <- function(sbc_out){
  ranks <- sapply(sbc_out, function(x) x$rank / x$Ns)
  list(
    ks_test     = ks.test(ranks, "punif"),
    n_divergent = sum(sapply(sbc_out, function(x) x$n_divergent)),
    ranks       = ranks
  )
}

# Saves a summarize_sbc() result the same way save_report()/save_gg() do for
# a regular fit (see hiermod_utils_savers.R) -- a text report (KS test +
# total divergences) plus the rank histogram as a PDF, so an SBC run leaves
# an archived record instead of only living in console scrollback.
save_sbc_report <- function(sbc_out, step, dir = hiermod_out_dir){
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(dir, paste0("SBC_report_", step, ".txt"))
  con <- file(path, open = "wt")
  sink(con)
  on.exit({ sink(); close(con) })
  print(sbc_out$ks_test)
  cat("\ntotal divergences:", sbc_out$n_divergent, "\n")
  invisible(path)
  message("Saved to ", path)

  save_pdf("SBC_rank_hist", step, function()
    hist(sbc_out$ranks, main = paste("SBC rank histogram --", step), xlab = "rank"),
    dir = dir)
}
