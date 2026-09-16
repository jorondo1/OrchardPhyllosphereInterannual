# Simulation-based calibration via SBC package (https://hyunjimoon.github.io/SBC/),
#
# Typical calibration-script usage (see any Models/*_model.R for the
# model-specific dq_X derived_quantities() this expects, and any numbered
# calibration script for the full pattern in context):
#
#       sbc_gen <- make_sbc_generator(
#         fit_sim, simulate_from_priors_X,
#         keep = c("loga", "sigma"), gen_cols = c("Dv", "Mg"),
#         extra_globals = "sim_div_X")
#
#       sbc_X <- run_sbc_pipeline(
#         sbc_gen$generator, sbc_gen$globals, n_sbc = 100,
#         model = model, model_id = model_id, n_iter = 10000,
#         hiermod_out_dir = hiermod_out_dir, dquants = dq_X,
#         control = list(adapt_delta = 0.99))
#
#       plot_sbc_diagnostics(sbc_X, model_id, n_sbc = 100)
#
#       save_sbc_health_report(
#         model_id, sbc_X, n_sbc = 100, n_iter = 10000,
#         variables = c("loga[1]", "loga[2]", "median_contrast"),
#         hiermod_out_dir = hiermod_out_dir)

# Pulls the i-th prior draw out of extract.prior()'s output as a plain
# named list, for feeding to a model's simulate_from_priors_X() as one
# "true" parameter set (an indexed parameter comes back as a matrix row, a
# scalar as one vector element).
draw_true <- function(priors, i){
  purrr::imap(priors, function(x, name){
    x <- as.matrix(x)
    if (ncol(x) == 1) x[i, 1] else x[i, ]
  })
}

# ---- SBC_backend_ulam: lets compute_SBC() fit rethinking::ulam() models directly ----

SBC_backend_ulam <- function(model, ...){
  structure(list(model = model, args = list(...)), class = "SBC_backend_ulam")
}

# chains is always set to whatever `cores` compute_SBC() hands us for this
# fit, never taken from backend$args

SBC_fit.SBC_backend_ulam <- function(backend, generated, cores){
  do.call(
    rethinking::ulam,
    c(
      list(
        flist = backend$model, 
        data = generated, 
        chains = cores, cores = cores), 
      backend$args)
    )
}

SBC_fit_to_draws_matrix.ulam <- function(fit){
  posterior::as_draws_matrix(fit@cstanfit$draws())
}

SBC_fit_to_diagnostics.ulam <- function(fit, fit_output, fit_messages, fit_warnings){
  # quiet = TRUE: this call itself would otherwise print the same
  # divergence/treedepth text a second time on top of cmdstanr's own
  # automatic post-sample message)
  ds <- fit@cstanfit$diagnostic_summary(diagnostics = c("divergences", "treedepth", "ebfmi"), quiet = TRUE)
  data.frame(n_divergent = sum(ds$num_divergent), n_max_treedepth = sum(ds$num_max_treedepth),
             n_low_ebfmi = sum(ds$ebfmi < 0.3))
}

# ---- Generator factory ------------------------------------------------------

# Builds a generate_datasets()-ready generator closure for one model,
# replacing every script's own copy-pasted generate_one_X(): on each call,
# draws one fresh prior sample from `fit`, simulates a dataset via
# `simulate_fn`, and packages both the way the SBC package expects.
#   fit           : the recovery ulam fit supplying the model formula + priors
#                   (extract.prior(fit, ...) is called fresh inside the
#                   closure on every replicate)
#   simulate_fn   : function(true_params, ...) -> data.frame, e.g.
#                   simulate_from_priors_MDL
#   keep          : names to retain from draw_true()'s output, e.g.
#                   c("loga", "sigma", "sigma_loc")
#   gen_cols      : columns of simulate_fn()'s output to hand to ulam(),
#                   e.g. c("Dv", "Mg", "Lo")
#   extra_globals : character vector of any OTHER top-level object names
#                   simulate_fn's own call chain depends on beyond itself
#                   (e.g. "sim_div_MDL"; a two-level chain like MDLSYv's
#                   needs both "sim_div_MDLSY" and "true_sigma_from_ls")
#   ...           : extra fixed arguments forwarded to simulate_fn on every
#                   call (e.g. shift = 1)
# Returns list(generator = <closure>, globals = <named list of VALUES>) --
# pass $generator to run_sbc_pipeline() as-is and $globals straight through
# as its future.globals. Values (not name-strings) are used deliberately:
# future.apply's future.globals accepts a named list of values directly,
# which avoids relying on future's own static code-analysis to trace names
# into this closure (already confirmed NOT to work automatically -- it
# doesn't see past generate_datasets()'s own internal call layer inside the
# SBC package -- and fresh multisession workers don't have rethinking/
# tidyverse attached either, hence the library() calls inside the closure).
make_sbc_generator <- function(fit, simulate_fn, keep, gen_cols,
                               extra_globals = character(0), ...){
  extra_args <- list(...)
  
  generator <- function(){
    library(rethinking); library(tidyverse) # future::multisession workers start fresh (not auto-attached)
    true_params <- suppressMessages(suppressWarnings(
      draw_true(extract.prior(fit, n = 1, refresh = 0), 1)))[keep]
    dat <- do.call(simulate_fn, c(list(true_params), extra_args))
    list(variables = true_params, generated = as.list(dat[, gen_cols]))
  }
  
  globals <- c(
    list(draw_true = draw_true, fit = fit, simulate_fn = simulate_fn, extra_args = extra_args,
         keep = keep, gen_cols = gen_cols),
    mget(extra_globals, envir = .GlobalEnv))
  
  list(generator = generator, globals = globals)
}

# ---- Pipeline: generate + fit ------------------------------------------------

# Generates n_sbc SBC replicate datasets in parallel 
run_sbc_pipeline <- function(generator, globals, n_sbc, model, model_id, n_iter,
                             hiermod_out_dir, dquants = NULL, refresh = 0, ...){
  future::plan(future::multisession)
  
  datasets <- SBC::generate_datasets(
    SBC::SBC_generator_function(
      generator, future.chunk.size = SBC::default_chunk_size(n_sbc),
      future.globals = globals),
    n_sbc)
  
  backend <- SBC_backend_ulam(model, iter = n_iter, refresh = refresh, ...)
  
  SBC::compute_SBC(
    datasets, backend, dquants = dquants,
    #cache_mode = "results", cache_location = file.path(hiermod_out_dir, paste0("sbc_cache_", model_id)),
    globals = c("SBC_fit.SBC_backend_ulam", "SBC_fit_to_draws_matrix.ulam", "SBC_fit_to_diagnostics.ulam"))
}

# ---- Diagnostic plots --------------------------------------------------------

# Rank-histogram / ECDF-diff / coverage plots for one compute_SBC() result,
# saved via save_gg() under this project's "<model_id>_<n_sbc>sbc_iter"
# naming convention. Returns the three ggplot objects invisibly.
#   - rank histogram: the direct, intuitive view; catches gross violations
#     at a glance, but its apparent shape depends on bin count/placement.
#   - ECDF difference (Sailynoja, Burkner & Vehtari 2022): checks
#     uniformity everywhere at once via a simultaneous (DKW) band, stays
#     informative at lower replicate counts and localizes *where* a
#     violation happens.
#   - central-interval coverage: "does my reported interval actually cover
#     the truth that often" -- should accompany, not replace, the
#     shape-based checks above.
plot_sbc_diagnostics <- function(sbc_result, model_id, n_sbc, width = 9, height = 7){
  sbc_step <- paste0(model_id, "_", n_sbc, "sbc_iter")
  p_rank   <- SBC::plot_rank_hist(sbc_result)
  p_ecdf   <- SBC::plot_ecdf_diff(sbc_result)
  p_cover  <- SBC::plot_coverage(sbc_result)
  save_gg("SBC_rank_hist", sbc_step, p_rank,  width = width, height = height)
  save_gg("SBC_ecdf_diff", sbc_step, p_ecdf,  width = width, height = height)
  save_gg("SBC_coverage",  sbc_step, p_cover, width = width, height = height)
  invisible(list(rank_hist = p_rank, ecdf_diff = p_ecdf, coverage = p_cover))
}

# ---- SBC health report -------------------------------------------------------
# One plain-text file per model -- the PRIMARY artifact to check per run
# (the PDFs above are secondary/visual backup) -- so calibration health
# (sampling diagnostics + rank-fraction calibration for the variables that
# actually matter) can be tracked/diffed/grepped across models and reruns
# without re-opening R or re-reading rank-histogram PDFs. Written after
# compute_SBC(); relies only on fields every SBC_backend_ulam-based result
# has: $stats, $default_diagnostics (built into every SBC_results object
# regardless of backend), and $backend_diagnostics with this backend's own
# n_divergent/n_max_treedepth/n_low_ebfmi columns (SBC_fit_to_diagnostics.ulam()
# above).
#
# variables: character vector of tracked stats (e.g. c("loga[1]", "loga[2]"))
# to check rank-fraction calibration for -- typically the ones under live
# investigation, not necessarily every variable in the model.
save_sbc_health_report <- function(model_id, sbc_result, n_sbc, n_iter, variables,
                                   hiermod_out_dir, alpha = 0.01){
  dd <- sbc_result$default_diagnostics
  bd <- sbc_result$backend_diagnostics
  
  se     <- 1 / sqrt(12 * n_sbc) # SE of a mean rank-fraction under perfect calibration
  z_crit <- qnorm(1 - alpha / 2)
  
  rank_summary <- sbc_result$stats |>
    dplyr::filter(variable %in% variables) |>
    dplyr::group_by(variable) |>
    dplyr::summarise(
      mean_rank_frac   = mean(rank / max_rank),
      median_rank_frac = median(rank / max_rank),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      z    = (mean_rank_frac - 0.5) / se,
      flag = ifelse(abs(z) > z_crit, "MISCALIBRATED", "ok")
    )
  
  lines <- c(
    sprintf("SBC health report -- %s", model_id),
    sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    sprintf("n_sbc = %d, n_iter = %d", n_sbc, n_iter),
    "",
    "-- Backend diagnostics (sampling health) --",
    sprintf("  divergences:     %d total", sum(bd$n_divergent)),
    sprintf("  max treedepth:   %d total (worst single replicate: %d)",
            sum(bd$n_max_treedepth), max(bd$n_max_treedepth)),
    sprintf("  low E-BFMI fits: %d total", sum(bd$n_low_ebfmi)),
    sprintf("  max Rhat:        %.3f (%d/%d fits > 1.01)",
            max(dd$max_rhat, na.rm = TRUE), sum(dd$max_rhat > 1.01, na.rm = TRUE), nrow(dd)),
    sprintf("  min ESS bulk:    %.0f", min(dd$min_ess_bulk, na.rm = TRUE)),
    sprintf("  min ESS tail:    %.0f", min(dd$min_ess_tail, na.rm = TRUE)),
    "",
    sprintf("-- Rank-fraction calibration (expected 0.5, SE = %.4f at n=%d) --", se, n_sbc),
    sprintf("  %-20s %10s %10s %8s  %s", "variable", "mean_frac", "med_frac", "z", "flag")
  )
  
  for (i in seq_len(nrow(rank_summary))) {
    r <- rank_summary[i, ]
    lines <- c(lines, sprintf("  %-20s %10.3f %10.3f %8.2f  %s",
                              r$variable, r$mean_rank_frac, r$median_rank_frac, r$z, r$flag))
  }
  
  overall <- if (any(rank_summary$flag == "MISCALIBRATED")) "MISCALIBRATED" else "OK"
  lines <- c(lines, "", sprintf("Overall: %s", overall))
  
  out_path <- file.path(hiermod_out_dir, paste0("sbc_health_", model_id, "_",n_sbc,"iter.txt"))
  writeLines(lines, out_path)
  message("SBC health report written to ", out_path)
  invisible(rank_summary)
}
