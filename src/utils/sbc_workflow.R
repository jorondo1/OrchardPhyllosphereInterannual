# Simulation-based calibration via the SBC package (https://hyunjimoon.github.io/SBC/)
# - dq_X estimands: defined in each Models/*_model.R
# Typical usage:
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

# Contrast recovery -------------------------

# Posterior vs. true May gap / July gap / seasonal change (median-based)
contrast_recovery <- function(fit, means_fn, may_conv, may_org,
                              july_conv_shift, july_org_shift, shift = 0){
  pf <- post_full(fit, means_fn, shift = shift)
  m  <- pf$median
  
  estimands <- estimand_rows(list(
    "Median May gap (Organic - Conventional)"    = m$median_3 - m$median_1,
    "Median July gap (Organic - Conventional)"   = m$median_4 - m$median_2,
    "Seasonal change in median gap (July - May)" = (m$median_4 - m$median_2) - (m$median_3 - m$median_1)
  ))
  
  may_gap  <- may_org - may_conv
  july_gap <- may_org * exp(july_conv_shift + july_org_shift) - may_conv * exp(july_conv_shift)
  
  true_estimands <- tribble(
    ~statistic,                                    ~value,
    "Median May gap (Organic - Conventional)",     may_gap,
    "Median July gap (Organic - Conventional)",    july_gap,
    "Seasonal change in median gap (July - May)",  july_gap - may_gap,
  )
  
  list(estimands = estimands, true_estimands = true_estimands)
}


# i-th prior draw from extract.prior() as a named list ("true" parameter set)
# - vector parameters -> matrix row; scalars -> single value
draw_true <- function(priors, i){
  purrr::imap(priors, function(x, name){
    x <- as.matrix(x)
    if (ncol(x) == 1) x[i, 1] else x[i, ]
  })
}

# ---- SBC_backend_ulam: lets compute_SBC() fit rethinking::ulam() models ----
# S3 methods below are what the SBC package calls for a custom backend

SBC_backend_ulam <- function(model, ...){
  structure(list(model = model, args = list(...)), class = "SBC_backend_ulam")
}

# chains = cores given by compute_SBC(), never from backend$args

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

# Posterior draws as a draws_matrix
SBC_fit_to_draws_matrix.ulam <- function(fit){
  posterior::as_draws_matrix(fit@cstanfit$draws())
}

# Per-fit sampling diagnostics: divergences, max treedepth, low E-BFMI
SBC_fit_to_diagnostics.ulam <- function(fit, fit_output, fit_messages, fit_warnings){
  # quiet: avoid printing diagnostics twice
  ds <- fit@cstanfit$diagnostic_summary(diagnostics = c("divergences", "treedepth", "ebfmi"), quiet = TRUE)
  data.frame(n_divergent = sum(ds$num_divergent), n_max_treedepth = sum(ds$num_max_treedepth),
             n_low_ebfmi = sum(ds$ebfmi < 0.3))
}

# ---- Generator factory ------------------------------------------------------

# SBC dataset generator for one model: each call draws one prior sample
# from `fit` and simulates a dataset with it
#   fit           : recovery ulam fit (supplies formula + priors)
#   simulate_fn   : function(true_params, ...) -> data.frame
#   keep          : parameters to track, e.g. c("loga", "sigma")
#   gen_cols      : simulated columns passed to ulam(), e.g. c("Dv", "Mg")
#   extra_globals : other object names simulate_fn needs, e.g. "sim_div_X"
#   ...           : fixed args for simulate_fn (e.g. shift = 1)
# Returns list(generator, globals)
# - globals: object VALUES for parallel workers (future can't find them by itself)
make_sbc_generator <- function(fit, simulate_fn, keep, gen_cols,
                               extra_globals = character(0), ...){
  extra_args <- list(...)
  
  generator <- function(){
    library(rethinking); library(tidyverse) # parallel workers start with no packages
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

# Generate n_sbc datasets in parallel, then fit each (compute_SBC)
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

# Save SBC plots for one compute_SBC() result
# - rank histogram: quick view of gross violations
# - ECDF difference (Sailynoja et al. 2022): simultaneous band, shows where it fails
# - coverage: do central intervals cover the truth as often as they should?
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
# Plain-text SBC summary per model (main calibration record; tracked in git)
# - sampling health: divergences, treedepth, E-BFMI, Rhat, ESS
# - rank calibration: mean rank fraction vs 0.5, z-score, flag at `alpha`
# - variables: which tracked stats to check, e.g. c("loga[1]", "loga[2]")
save_sbc_health_report <- function(model_id, sbc_result, n_sbc, n_iter, variables,
                                   hiermod_out_dir, alpha = 0.01){
  dd <- sbc_result$default_diagnostics
  bd <- sbc_result$backend_diagnostics
  
  se     <- 1 / sqrt(12 * n_sbc) # SE of a uniform mean rank fraction
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
