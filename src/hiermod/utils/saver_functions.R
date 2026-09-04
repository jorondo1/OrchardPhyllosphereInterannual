
# Saves a fitted ulam object as .rds, so raw draws survive past the R
# session without a refit (e.g. for a later independent check).
#
# ulam()'s default cmdstanr backend keeps the underlying fit (attr(fit,
# "cstanfit"), an R6 CmdStanMCMC) lazily backed by Stan's own temp-dir CSV
# output -- gone the moment this R session exits. A plain saveRDS() looks
# fine here but silently breaks extract.samples()/precis() in a later
# session ("File does not exist: .../ulam_cmdstanr_*.csv"). Forcing these
# accessor methods first caches everything inside the R6 object itself
# before saving -- exactly what cmdstanr's own (otherwise-inapplicable,
# since it'd only save the bare cstanfit, not the whole ulam wrapper)
# CmdStanFit$save_object() does internally.
save_fit <- function(name, step, fit, dir = hiermod_out_dir){
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  cs <- attr(fit, "cstanfit")
  if (inherits(cs, "R6")) {
    invisible(cs$draws())
    try(invisible(cs$sampler_diagnostics()), silent = TRUE)
    try(invisible(cs$init()), silent = TRUE)
    try(invisible(cs$profiles()), silent = TRUE)
  }
  path <- file.path(dir, paste0(name, "_", step, ".rds"))
  saveRDS(fit, path, compress = 'xz')
  invisible(path)
  message("Saved to ", path)
}

# For ggplot objects. type = "png" for point-heavy plots (e.g.
# bayesplot::mcmc_pairs() on tens of thousands of draws) -- a vector PDF
# stores every point as its own object and can balloon to tens of MB; a
# rasterized PNG is a flat pixel grid regardless of point count, typically
# orders of magnitude smaller for the same plot. dpi only matters for png.
save_gg <- function(name, step, plot = ggplot2::last_plot(), width = 10, height = 8,
                    dir = hiermod_out_dir, type = "pdf", dpi = 150){
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(dir, paste0(name, "_", step, ".", type))
  ggsave(path, plot = plot, width = width, height = height, dpi = dpi)
  invisible(path)
  message("Saved to ", path)
}

# For base-graphics / rethinking plots (dens(), hist(), pairs(), traceplot()):
# pass the plotting code as a zero-arg function, e.g.
#   save_pdf("fit_pairs", "MDL", function() pairs(fit_MDL_sim, pars = c("sigma_loc","b[1]","b[2]")))
save_pdf <- function(name, step, plot_call, width = 10, height = 8,
                     dir = hiermod_out_dir){
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(dir, paste0(name, "_", step, ".pdf"))
  pdf(path, width = width, height = height)
  plot_call()
  dev.off()
  invisible(path)
  message("Saved to ", path)
}

# For plain-text numeric reports: writes to <hiermod_out_dir>/<name>_<step>.txt
# post_counts, model, and recovery are all optional; omit post_counts for
# sim-recovery fits where you just want the precis() table, pass model (the
# alist(), e.g. model_MDLb) to record exactly which model spec produced this
# fit alongside the numbers, e.g.
#   save_report("fit_summary", "MDLb", fit_MDLb, post_counts, model_MDLb)
#   save_report("sim_summary", "MDLb", fit_MDLb_sim, model = model_MDLb)
# post_counts: the LONG statistic/group/value tibble from compute_contrasts()/
# estimand_panels()/estimand_rows() (NOT post_full()'s own raw list output --
# that's a named list of per-parameter tibbles with no statistic/group/value
# shape, report_contrasts_full() can't do anything with it) -- or, for
# models 1-2's older postcounts()/plot_contrast_density() pipeline, a
# data.frame with a `contrast` column.
# recovery: a check_recovery() tibble (or several bind_rows()'d together,
# e.g. fixed effects + sigma-by-cell), for a Parameter recovery section.
# model_name: an optional display name for the header, purely cosmetic --
# deparse(substitute(model)) always printed the literal string "model" (every
# call site's variable is named that), which said nothing about which model
# it actually was. Pass a nickname (e.g. "The Bare Bones") to replace it.
save_report <- function(name, step, fit, post_counts = NULL, model = NULL,
                        recovery = NULL, depth = 2, dir = hiermod_out_dir,
                        model_name = NULL){
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(dir, paste0(name, "_", step, ".txt"))
  con <- file(path, open = "wt")
  sink(con)
  on.exit({ sink(); close(con) })

  if (!is.null(model)) {
    header <- if (!is.null(model_name)) model_name else deparse(substitute(model))
    cat("==== model:", header, "====\n\n")
    print(model)
    cat("\n\n")
  }

  cat("==== precis:", deparse(substitute(fit)), "====\n\n")
  print(round(precis(fit, depth = depth), 3))

  if (!is.null(recovery)) {
    cat("\n\n==== parameter recovery ====\n")
    print(as.data.frame(recovery %>% mutate(across(where(is.numeric), ~round(.x, 3)))))
    cat("covered:", sum(recovery$covered), "/", nrow(recovery), "\n")
  }

  if (!is.null(post_counts)) {
    cat("\n\n==== posterior contrast ====\n")
    if ("statistic" %in% names(post_counts)) {
      # long statistic/group/value tibble (compute_contrasts()/estimand_panels()/
      # estimand_rows()): mean/median/sigma/estimand contrasts, all at once
      report_contrasts_full(post_counts) %>%
        purrr::pwalk(function(statistic, mean, median, PI89_lower, PI89_upper, HPDI_lower, HPDI_upper, ...){
          cat(as.character(statistic), "-- mean:", round(mean,3), " median:", round(median,3),
              " 89% PI: [", round(PI89_lower,3), ",", round(PI89_upper,3), "]",
              " 89% HPDI: [", round(HPDI_lower,3), ",", round(HPDI_upper,3), "]\n")
        })
    } else {
      cat("mean:  ", mean(post_counts$contrast), "\n")
      cat("median:", median(post_counts$contrast), "\n")
      cat("89% PI:\n")
      print(PI(post_counts$contrast))
    }
  }

  invisible(path)
  message("Saved to ", path)

}

# Analyses reports

# post-predictive check stats
# 


