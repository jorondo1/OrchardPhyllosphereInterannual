# predictive_checks.R -- prior predictive (spaghetti) and posterior
# predictive (density overlay, contrast test-statistic) plotting helpers.
# Building the long-format prior_pred data frame itself (simulate_from_priors()
# looped with map_dfr(..., .id = "draw")) is done per-model.

# ---- Prior predictive check -------------------------------------------------

# Overlays many prior-predictive-simulated densities (one per prior draw),
# to check whether the priors imply plausible data before ever fitting to
# anything real. Not for posterior predictive checks (those use bayesplot,
# below). model (optional): prints the alist() on the plot, to document
# which priors produced it. observed (optional): overlays one real/simulated
# dataset's own density for comparison against the prior ensemble.
prior_predictive_spaghetti <- function(
    prior_pred, model = NULL, value_col = "Dv", draw_col = "draw",
    group_col = NULL, n_sample = 100, upper_q = 0.99,
    observed = NULL, title = "Prior predictive check"){

  fmt <- function(x) format(signif(x, 3), big.mark = ",", scientific = FALSE)

  xlim_upper <- quantile(prior_pred[[value_col]], upper_q, na.rm = TRUE)
  all_draw_ids <- unique(prior_pred[[draw_col]])

  # % of replicates (not raw points) with >=1 value beyond the cutoff --
  # shows whether the tail is a few extreme replicates or spread thin
  # across everyone.
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
    # same xlim_upper trim as the ensemble, so bandwidth/shape stay comparable
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

# Which prior parameter(s) explain the extreme tail of a prior predictive
# check? Every calibration script used to hand-roll this (and one -- ITS's
# own 2.2_MDL_ITS_calibration.R ppc1 stage -- had `xlim_upper_ppc1 <- 0.99`
# instead of `quantile(..., 0.99)`, silently flagging almost every draw as
# "extreme").
#
# `prior_pred` is the long-format tibble built by looping
# simulate_from_priors_*() with map_dfr(..., .id = "draw"); `candidates` is
# a named list of length-n_prior vectors, one per prior parameter suspected
# of driving the tail (e.g. list(sigma_loc = extracted_prior$sigma_loc,
# ...)), each indexed 1:n_prior in the same draw order
# extract.prior()/draw_true() produced them -- not resorted to match
# prior_pred's own (possibly lexically-sorted) draw labels.
#
# Reports two things: (1) each candidate's median among "extreme" vs.
# "normal" draws (the eyeball check every script already did by hand), and
# (2) each candidate's correlation with log(max simulated value) -- added
# because medians can look different between groups even when the
# relationship is weak (MDLS2v: sigma_loc/sigma_tr's medians differed by
# group but correlated only ~0.2; loga -- not in the original group
# breakdown at all -- correlated ~0.7 and was the real driver).
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

# Reshapes sim()'s wide draws-by-observation matrix into a long tibble.
postpred_as_long_tibble <- function(post_pred) {
  post_pred %>%
    as_tibble(.name_repair = "minimal") %>%
    set_names(seq_len(ncol(.))) %>%
    mutate(sim = row_number()) %>%
    pivot_longer(-sim, names_to = "obs", values_to = "Dv_sim") %>%
    mutate(obs = as.integer(obs))
}

# Posterior predictive density overlay, split by group. Thin wrapper
# around bayesplot::ppc_dens_overlay_grouped() since the sim()+labeling
# step repeats identically across models 1-3.
plot_ppc_overlay <- function(fit, dat, group, n = 50, xlim = NULL){
  yrep <- sim(fit, dat, n = n)
  p <- bayesplot::ppc_dens_overlay_grouped(dat$Dv, yrep, group = group)
  if (!is.null(xlim)) p <- p + coord_cartesian(xlim = xlim)
  p
}

# Builds a two-group contrast statistic for ppc_stat() (FUN on group==2
# minus FUN on group==1). Only fits a clean 1/2 split; a model with more
# groups (e.g. Mg x Mo cells) needs its own stat function written directly.
contrast_stat <- function(FUN, group){
  force(FUN); force(group)
  function(y) FUN(y[group == 2]) - FUN(y[group == 1])
}

# sim() + ppc_stat() + a title, for one contrast_stat().
plot_ppc_contrast_stat <- function(fit, dat, group, FUN, stat_name, group_labels, n = 1000){
  yrep <- sim(fit, dat, n = n)
  bayesplot::ppc_stat(dat$Dv, yrep, stat = contrast_stat(FUN, group)) +
    labs(title = paste0("PPC: ", stat_name, " contrast (", group_labels[2], " - ", group_labels[1], ")"))
}

# The May gap / July gap / seasonal-change PPC used from Model 4 onward, as
# three stacked ppc_stat() panels. dat needs $Dv/$Mg/$Mo (1/2-coded);
# Mg==2 is assumed to be the "treatment" side (Organic).
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

  p_may / p_july / p_change
}

# ---- Model comparison --------------------------------------------------

# PSIS model comparison between two ulam fits (both need log_lik = TRUE at
# fit time), with a Pareto-k sanity check first -- k > 0.7 means that
# observation's importance-sampling estimate is unreliable, and if a large
# share of observations are flagged, the comparison table itself (weights,
# dPSIS) isn't trustworthy, not just that one point. Common in hierarchical
# models with few observations per group (e.g. Tree here).
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
