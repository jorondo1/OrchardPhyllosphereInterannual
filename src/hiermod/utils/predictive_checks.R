# predictive_checks.R 

# - prior predictive (spaghetti) and posterior predictive (overlay density, 
# contrast test-statistic) plotting helpers.

# Building the long-format prior_pred data frame itself (draw_true_X()/
# simulate_X(), looped with map_dfr(..., .id = "draw")) done per-model.

# ---- Prior predictive check -------------------------------------------------

# Many prior draws -> many simulated datasets -> one density line per
# replicate, overlaid. 
# 
# NOT for posterior predictive checks: those use bayesplot instead (below).
#
# model (optional): pass the alist() itself to print it on the plot, to
# document which priors produced it. Also prints mean/median/sd computed on ALL
# replicates (n_prior), not just the n_sample plotted.
#
# observed (optional): a numeric vector -- e.g. dat_sim$Dv, the one actual
# simulated dataset a Parameter recovery fit was trained on -- overlaid as
# its own density line, so the prior predictive envelope can be eyeballed
# against the one draw that's actually in play, not just the ensemble.
#
#   prior_pred_MDL <- map_dfr(seq_len(n_prior), function(i) simulate_MDL(draw_true_MDL(i)), .id = "draw")
#   prior_predictive_spaghetti(prior_pred_MDL, model = model_MDL, title = "Prior predictive check -- model_MDL")

prior_predictive_spaghetti <- function(
    prior_pred, model = NULL, value_col = "Dv", draw_col = "draw",
    group_col = NULL, n_sample = 100, upper_q = 0.99,
    observed = NULL, title = "Prior predictive check"){

  fmt <- function(x) format(signif(x, 3), big.mark = ",", scientific = FALSE)

  xlim_upper <- quantile(prior_pred[[value_col]], upper_q, na.rm = TRUE)
  all_draw_ids <- unique(prior_pred[[draw_col]])

  # % of REPLICATES (not raw points) with >=1 value beyond the cutoff --
  # asks how concentrated the tail is (a handful of replicates going way out
  # vs. everyone contributing a little), not the ~constant raw-point fraction
  # (1 - upper_q) that a quantile trivially guarantees.
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
    # same xlim_upper trim as the spaghetti ensemble (bulk) -- keeps the two
    # densities' bandwidth/shape comparable, one outlier shouldn't distort it.
    # inherit.aes = FALSE: the base plot's group = draw_col aes doesn't apply
    # here, this is one single vector, not per-replicate draws.
    observed_bulk <- observed[observed <= xlim_upper]
    p <- p + geom_density(
      data = tibble(!!value_col := observed_bulk),
      aes(x = .data[[value_col]]), inherit.aes = FALSE,
      bw = shared_bw, colour = "black", linewidth = 0.9) +
      labs(caption = "Black line: density of the actual simulated dataset (observed=).")
  }

  if (!is.null(model)) {  #trimws removes large spaces created by breaks
    model_text <- paste(sapply(model, function(x) paste(trimws(deparse(x, width.cutoff = 40)), collapse = " ")), collapse = "\n", " ")
    p <- p + annotate("text", x = Inf, y = 0, hjust = 1, vjust = -0.9, size = 3,
                       family = "mono", colour = "grey20", label = model_text)
  }

  if (!is.null(group_col)) p <- p + facet_wrap(ggplot2::vars(.data[[group_col]]))
  p
}

# ---- Posterior predictive checks (bayesplot) -------------------------------

# convert sim() matrix output to a long tibble
postpred_as_long_tibble <- function(post_pred) {
  post_pred %>%
    as_tibble(.name_repair = "minimal") %>%
    set_names(seq_len(ncol(.))) %>%
    mutate(sim = row_number()) %>%
    pivot_longer(-sim, names_to = "obs", values_to = "Dv_sim") %>%
    mutate(obs = as.integer(obs))
}

# Density overlay split by a grouping factor -- thin wrapper around
# bayesplot::ppc_dens_overlay_grouped() since the sim()+group-label step is
# identical everywhere it's used (models 1-3). group: a label factor/vector
# for dat$Dv's rows, e.g. idx$Mg$to_label(dat$Mg).
plot_ppc_overlay <- function(fit, dat, group, n = 50, xlim = NULL){
  yrep <- sim(fit, dat, n = n)
  p <- bayesplot::ppc_dens_overlay_grouped(dat$Dv, yrep, group = group)
  if (!is.null(xlim)) p <- p + coord_cartesian(xlim = xlim)
  p
}

# Factory for a two-group contrast test statistic, for bayesplot::ppc_stat().
# group must be coded 1/2 (e.g. dat$Mg): returns FUN(group==2) - FUN(group==1).
# Only fits a clean binary split -- once a model's groups aren't a single 1/2
# index (model 4's Mg x Mo cells), write the per-cell/per-contrast stat
# function directly in that script instead.
contrast_stat <- function(FUN, group){
  force(FUN); force(group)
  function(y) FUN(y[group == 2]) - FUN(y[group == 1])
}

# sim() + ppc_stat() + title, for a contrast_stat(). group_labels: the two
# group names in 1/2 order (e.g. idx$Mg$levels), for the title.
plot_ppc_contrast_stat <- function(fit, dat, group, FUN, stat_name, group_labels, n = 1000){
  yrep <- sim(fit, dat, n = n)
  bayesplot::ppc_stat(dat$Dv, yrep, stat = contrast_stat(FUN, group)) +
    labs(title = paste0("PPC: ", stat_name, " contrast (", group_labels[2], " - ", group_labels[1], ")"))
}

# Same idea as plot_ppc_contrast_stat(), but for the Mg x Mo (May/July gap +
# seasonal change) design used from Model 4 onward -- 3 stacked ppc_stat()
# panels via patchwork. dat needs $Dv/$Mg/$Mo (Mg,Mo coded 1/2 as everywhere
# else); Mg==2 assumed to be the "treatment" side (Organic).
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

  p_may / p_july / p_change   # patchwork stack
}
