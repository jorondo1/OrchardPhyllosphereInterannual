# saver_functions.R -- consistent, per-fit output helpers (fits, plots,
# text reports) so every model script writes to the same place the same way.

# Saves a ulam fit to .rds. Forces cmdstanr's lazy accessors first, since a
# plain saveRDS() otherwise loses the fit once its temp Stan CSV output is
# gone (breaks extract.samples()/precis() in a later session).
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

# Saves a ggplot to file. Use type = "png" for point-heavy plots (e.g.
# mcmc_pairs() on tens of thousands of draws): a vector PDF of that many
# points balloons to tens of MB, a rasterized PNG doesn't.
save_gg <- function(name, step, plot = ggplot2::last_plot(), width = 10, height = 8,
                    dir = hiermod_out_dir, type = "pdf", dpi = 150){
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(dir, paste0(name, "_", step, ".", type))
  ggsave(path, plot = plot, width = width, height = height, dpi = dpi)
  invisible(path)
  message("Saved to ", path)
}

# Saves base-graphics output (dens(), hist(), pairs(), traceplot()) to a
# PDF; pass the plotting code as a zero-arg function, e.g.
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

# Writes one plain-text report (model spec, precis() table, parameter
# recovery, posterior contrasts) per fit. All three optional pieces are
# assembled in one place so nothing needs re-deriving to read the numbers
# back later. post_counts: the long statistic/group/value tibble from
# compute_contrasts()/estimand_panels()/estimand_rows() (not post_full()'s
# raw list -- wrong shape); recovery: a check_recovery() tibble. Header
# name resolution: model_name= override > attr(model, "name") (set once
# per model in its own model file, e.g.
# attr(model_MD_16S, "name") <- "Samwise the Steadfast" -- an attribute,
# not a list element, since ulam() iterates every element of its `model`
# argument expecting a formula quote() and errors on anything else, same
# reason fits carry cstanfit as attr(fit, "cstanfit") rather than a slot)
# > deparse(substitute(model)) (always just prints "model" -- every call
# site names its variable that -- so only a fallback for a model object
# that hasn't been given a name yet). model_name= stays available for the
# rare script that reports two models at once.
save_report <- function(name, step, fit, post_counts = NULL, model = NULL,
                        recovery = NULL, depth = 2, dir = hiermod_out_dir,
                        model_name = NULL){
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(dir, paste0(name, "_", step, ".txt"))
  con <- file(path, open = "wt")
  sink(con)
  on.exit({ sink(); close(con) })

  if (!is.null(model)) {
    header <- if (!is.null(model_name)) model_name
              else if (!is.null(attr(model, "name"))) attr(model, "name")
              else deparse(substitute(model))
    cat("==== model:", header, "====\n\n")
    print(model)
    cat("\n\n")
  }

  cat("==== precis:", deparse(substitute(fit)), "====\n\n")
  precis_fit <- precis(fit, depth = depth)
  # ess_bulk as a fraction of total post-warmup draws -- can exceed 1 under
  # negative autocorrelation (a real, documented feature of rank-normalized
  # ESS, not a bug: HMC sometimes samples weakly-identified parameters
  # slightly anti-correlated, which is *more* efficient than i.i.d.).
  # Values well below ~0.1-0.2 are the ones worth a second look.
  precis_fit$ess_ratio <- precis_fit$ess_bulk / NROW(extract.samples(fit)[[1]])
  print(round(precis_fit, 3))

  if (!is.null(recovery)) {
    cat("\n\n==== parameter recovery ====\n")
    print(as.data.frame(recovery %>% mutate(across(where(is.numeric), ~round(.x, 3)))))
    cat("covered:", sum(recovery$covered), "/", nrow(recovery), "\n")
  }

  if (!is.null(post_counts)) {
    cat("\n\n==== posterior contrast ====\n")
    if ("statistic" %in% names(post_counts)) {
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

# Comprehensive posterior summary as a styled HTML table (kableExtra), one
# row per statistic/group -- mean, median, 89% PI, 89% HPDI. Companion to
# save_report(): that one is model+precis only (a fit's own diagnostic
# record); this is the "results report" -- every interpretable posterior a
# reader might want a number for, not just the headline contrast. pc_full:
# a compute_contrasts()/estimand_panels()-shaped statistic/group/value
# tibble (bind_rows() together whatever the calling script already built,
# e.g. compute_contrasts(pf, keep=...) + pc_estimands_means +
# pc_estimands_medians -- see any *_16S_analysis.R script for the pattern).
save_posterior_kable <- function(name, step, pc_full, dir = hiermod_out_dir, caption = NULL){
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(dir, paste0(name, "_", step, ".html"))

  report_contrasts_full(pc_full) %>%
    mutate(across(where(is.numeric), ~round(.x, 3))) %>%
    kableExtra::kable("html", caption = caption %||% paste("Posterior summary:", step)) %>%
    kableExtra::kable_styling(full_width = FALSE) %>%
    kableExtra::save_kable(file = path)

  invisible(path)
  message("Saved to ", path)
}
