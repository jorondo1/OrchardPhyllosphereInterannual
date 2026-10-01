# saver_functions.R -- output helpers (fits, plots, reports), same paths/naming everywhere

# Save a ulam fit to .rds
# - loads cmdstanr's lazy fields first; otherwise lost with the temp CSVs
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

# Save a ggplot; type = "png" for point-heavy plots (PDFs get huge)
save_gg <- function(name, step, plot = ggplot2::last_plot(), width = 10, height = 8,
                    dir = hiermod_out_dir, type = "pdf", dpi = 150, prefix = hiermod_marker){
  
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(dir, paste0(prefix, "_", name, "_", step, ".", type))
  ggsave(path, plot = plot, width = width, height = height, dpi = dpi)
  invisible(path)
  message("Saved to ", path)
}

# Save base-graphics output to PDF; plotting code passed as a zero-arg function, e.g.
#   save_pdf("fit_pairs", "MDL", function() pairs(fit_MDL_sim, pars = c("sigma_loc","b[1]","b[2]")))
save_pdf <- function(name, step, plot_call, width = 10, height = 8,
                     dir = hiermod_out_dir, prefix = hiermod_marker){
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(dir, paste0(prefix, "_", name, "_", step, ".pdf"))
  pdf(path, width = width, height = height)
  plot_call()
  dev.off()
  invisible(path)
  message("Saved to ", path)
}

# Plain-text fit report: model spec, precis(), optional recovery + contrasts
# - post_counts: long statistic/group/value tibble (not post_full()'s list)
# - recovery: check_recovery() tibble
# - header name: model_name > attr(model, "name") > variable name
#   (name is an attribute: ulam() errors on non-formula list elements)
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
  # ess_bulk / total draws; > 1 possible (anti-correlated draws); < ~0.1-0.2 worth a look
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

# HTML posterior summary table: median, 89% HPDI, pd per statistic x group
# - "fold" statistics also get an inverted row
# - 3 decimals for variance-partition rows (or all rows with three_digits = TRUE)
save_posterior_kable <- function(
    name, step, pc_full, 
    dir = hiermod_out_dir, prefix = hiermod_marker,
    caption = NULL, all_groups = TRUE, three_digits = FALSE){
  
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(dir, paste0(prefix, name, "_", step, ".html"))
  
  base_tbl <- report_contrasts_full(pc_full, all_groups = all_groups) %>% 
    dplyr::select(statistic, group, median, HPDI_lower, HPDI_upper, pd)
  
  with_inverse <- base_tbl %>%
    bind_rows(
      filter(., str_detect(statistic, "fold")) %>%
        mutate(
          median = 1 / median,
          # 1/x swaps the bounds; pd unchanged
          HPDI_lower_inv = 1 / HPDI_upper,
          HPDI_upper     = 1 / HPDI_lower,
          HPDI_lower     = HPDI_lower_inv,
          statistic = paste(statistic, "(inverse)")
        ) %>%
        dplyr::select(-HPDI_lower_inv)
    )
  
  is_vp <- three_digits | str_detect(with_inverse$statistic, "Variance partition")
  
  with_inverse %>%
    mutate(across(
      where(is.numeric),
      ~ ifelse(
        is_vp, 
        sprintf("%.3f", .x), 
        ifelse( # 2 decimals below 1, else 1
          (abs(.x) < 1 & .x != 0),sprintf("%.2f", .x), sprintf("%.1f", .x)))
    )) %>%
    kableExtra::kable("html", caption = caption %||% paste0(
      "Posterior summary: ", step,
      ". pd = probability of direction (fraction of the posterior on the median's side of 0).")) %>%
    kableExtra::kable_styling(
      full_width = FALSE,
      bootstrap_options = c('striped', 'hover')) %>%
    kableExtra::save_kable(file = path)
  
  invisible(path)
  message("Saved to ", path)
}
