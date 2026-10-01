# predictive_checks.R -- prior and posterior predictive check plots
# - prior_pred data (simulate_from_priors() over draws) is built per model

# ---- Prior predictive check -------------------------------------------------

# Prior predictive "spaghetti": one density per prior draw
# - model (optional): prints the alist() priors on the plot
# - observed (optional): overlays one dataset's density for comparison
prior_predictive_spaghetti <- function(
    prior_pred, model = NULL, value_col = "Dv", draw_col = "draw",
    group_col = NULL, n_sample = 100, upper_q = 0.99,
    observed = NULL, title = "Prior predictive check"){

  fmt <- function(x) format(signif(x, 3), big.mark = ",", scientific = FALSE)

  xlim_upper <- quantile(prior_pred[[value_col]], upper_q, na.rm = TRUE)
  all_draw_ids <- unique(prior_pred[[draw_col]])

  # % of replicates with >= 1 value beyond the cutoff (few extreme draws vs. thin spread)
  prop_extreme <- prior_pred %>%
    group_by(.data[[draw_col]]) %>%
    summarise(extreme = any(.data[[value_col]] > xlim_upper), .groups = "drop") %>%
    pull(extreme) %>% mean() %>% {100 * .} %>% signif(3)

  bulk <- prior_pred[prior_pred[[value_col]] <= xlim_upper, ]
  shared_bw <- bw.nrd0(bulk[[value_col]])

  bulk_draw_ids <- unique(bulk[[draw_col]])
  sample_ids <- sample(bulk_draw_ids, min(n_sample, length(bulk_draw_ids)))
  bulk_sample <- bulk[bulk[[draw_col]] %in% sample_ids, ]

  stat_label <- paste0(
    "mean: ", fmt(mean(prior_pred[[value_col]])), "  \n",
    "median: ", fmt(median(prior_pred[[value_col]])), "  \n",
    "3rd quartile: ", fmt(quantile(prior_pred[[value_col]], 0.75)), "  \n",
    "sd: ", fmt(sd(prior_pred[[value_col]])), "  "
  )

  p <- ggplot(bulk_sample, aes(x = .data[[value_col]], group = .data[[draw_col]])) +
    geom_density(bw = shared_bw, alpha = 0.15, linewidth = 0.2, colour = "steelblue") +
    coord_cartesian(xlim = c(0, xlim_upper)) +
    theme_light() +
    labs(title = title,
         subtitle = paste0(
           length(sample_ids), " of ", length(all_draw_ids), " replicates shown -- ",
           prop_extreme, "% of all replicates have a point beyond the ",
           format(100 * upper_q), "th percentile (", fmt(xlim_upper), "; not shown)")) +
    annotate("text", x = Inf, y = Inf, hjust = 1, vjust = 1.2, size = 3,
             colour = "grey30", label = stat_label, family = "mono")

  if (!is.null(observed)) {
    # same trim and bandwidth as the ensemble
    observed_bulk <- observed[observed <= xlim_upper]
    p <- p + geom_density(
      data = tibble(!!value_col := observed_bulk),
      aes(x = .data[[value_col]]), inherit.aes = FALSE,
      bw = shared_bw, colour = "black", linewidth = 0.9) +
      labs(caption = "Black line: density of the actual simulated dataset (observed=).")
  }

  if (!is.null(model)) {
    model_text <- paste(
      sapply(model, function(x) paste(trimws(deparse(x, width.cutoff = 80)), collapse = "\n")),
      collapse = "\n")
    p <- p + annotate("text", x = Inf, y = 0, hjust = 1, vjust = -0.2, size = 3,
                       family = "mono", colour = "grey20", label = model_text)
  }

  if (!is.null(group_col)) p <- p + facet_wrap(ggplot2::vars(.data[[group_col]]))
  p
}

# ---- Prior predictive tail diagnosis ----------------------------------------

# Which prior parameters drive the prior predictive tail?
# - candidates: named list of prior draws (length n_prior, extract.prior() order)
# - reports each candidate's median in extreme vs. normal draws
# - and its correlation with log(max simulated value): the better driver signal
diagnose_extreme_tail <- function(prior_pred, candidates, value_col = "Dv", upper_q = 0.99){
  n_prior <- length(candidates[[1]])
  stopifnot(all(lengths(candidates) == n_prior))

  max_by_draw <- prior_pred %>%
    group_by(draw) %>%
    summarise(max_val = max(.data[[value_col]]), .groups = "drop") %>%
    mutate(draw = as.integer(draw))

  xlim_upper <- quantile(prior_pred[[value_col]], upper_q, na.rm = TRUE)
  extreme_draws <- max_by_draw$draw[max_by_draw$max_val > xlim_upper]

  cand_tbl <- as_tibble(candidates)
  cand_tbl$draw    <- seq_len(n_prior)
  cand_tbl$extreme <- cand_tbl$draw %in% extreme_draws

  by_group <- cand_tbl %>%
    group_by(extreme) %>%
    summarise(n = n(), across(all_of(names(candidates)), median), .groups = "drop")

  ordered <- cand_tbl[match(max_by_draw$draw, cand_tbl$draw), names(candidates), drop = FALSE]
  correlations <- sort(
    sapply(ordered, function(x) cor(log(max_by_draw$max_val), x)),
    decreasing = TRUE)

  cat(sprintf("Extreme draws (max %s > %.0f%% quantile = %s): %d / %d\n",
              value_col, 100 * upper_q, format(signif(xlim_upper, 4), big.mark = ","),
              length(extreme_draws), n_prior))
  print(by_group)
  cat("\nCorrelation with log(max ", value_col, ") -- biggest driver first:\n", sep = "")
  print(round(correlations, 3))

  invisible(list(xlim_upper = xlim_upper, extreme_draws = extreme_draws,
                 by_group = by_group, correlations = correlations))
}

# ---- Posterior predictive checks (bayesplot) -------------------------------

# sim() draws x obs matrix -> long tibble
postpred_as_long_tibble <- function(post_pred) {
  post_pred %>%
    as_tibble(.name_repair = "minimal") %>%
    set_names(seq_len(ncol(.))) %>%
    mutate(sim = row_number()) %>%
    pivot_longer(-sim, names_to = "obs", values_to = "Dv_sim") %>%
    mutate(obs = as.integer(obs))
}

# Posterior predictive density overlay by group (bayesplot wrapper)
plot_ppc_overlay <- function(fit, dat, group, n = 50, xlim = NULL){
  yrep <- sim(fit, dat, n = n)
  p <- bayesplot::ppc_dens_overlay_grouped(dat$Dv, yrep, group = group)
  if (!is.null(xlim)) p <- p + coord_cartesian(xlim = xlim)
  p
}

# Two-group contrast statistic for ppc_stat(): FUN(group 2) - FUN(group 1)
# - force(): fixes FUN/group at creation (closure)
contrast_stat <- function(FUN, group){
  force(FUN); force(group)
  function(y) FUN(y[group == 2]) - FUN(y[group == 1])
}

# PPC of one contrast_stat()
plot_ppc_contrast_stat <- function(fit, dat, group, FUN, stat_name, group_labels, n = 1000){
  yrep <- sim(fit, dat, n = n)
  bayesplot::ppc_stat(dat$Dv, yrep, stat = contrast_stat(FUN, group)) +
    labs(title = paste0("PPC: ", stat_name, " contrast (", group_labels[2], " - ", group_labels[1], ")"))
}

# PPC of the May gap, July gap and seasonal change in gap (Organic - Conventional)
# - dat: $Dv, $Mg, $Mo (1/2-coded; Mg 2 = Organic)
plot_ppc_season_contrast_stats <- function(fit, dat, n = 1000){
  yrep <- sim(fit, dat, n = n)

  may_gap_stat  <- function(y) median(y[dat$Mg==2 & dat$Mo==1]) - median(y[dat$Mg==1 & dat$Mo==1])
  july_gap_stat <- function(y) median(y[dat$Mg==2 & dat$Mo==2]) - median(y[dat$Mg==1 & dat$Mo==2])
  change_stat   <- function(y) july_gap_stat(y) - may_gap_stat(y)

  p_may    <- bayesplot::ppc_stat(dat$Dv, yrep, stat = may_gap_stat) +
    labs(title = "PPC: May gap (Organic - Conventional)")
  p_july   <- bayesplot::ppc_stat(dat$Dv, yrep, stat = july_gap_stat) +
    labs(title = "PPC: July gap (Organic - Conventional)")
  p_change <- bayesplot::ppc_stat(dat$Dv, yrep, stat = change_stat) +
    labs(title = "PPC: Seasonal change in gap")

  p_may / p_july / p_change + plot_layout(guides = 'collect')
}

# ---- Model comparison --------------------------------------------------

# PSIS comparison of two ulam fits (both fitted with log_lik = TRUE)
# - warns on Pareto k > 0.7 (unreliable; common with few obs per tree)
psis_compare <- function(fit_a, fit_b){
  for (f in list(fit_a, fit_b)) {
    k <- suppressWarnings(PSIS(f, pointwise = TRUE)$k)
    n_bad <- sum(k > 0.7, na.rm = TRUE)
    if (n_bad > 0) {
      cat(sprintf("Pareto k > 0.7 for %d/%d observations -- PSIS comparison may not be trustworthy.\n",
                  n_bad, length(k)))
    }
  }
  rethinking::compare(fit_a, fit_b, func = PSIS)
}
