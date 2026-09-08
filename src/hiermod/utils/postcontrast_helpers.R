# postcontrast_helpers.R -- turn posterior draws (via each model's own
# means_X()) into contrast tables and plots. One consistent shape so every
# model, from a simple 2-group difference to a multi-estimand model, can
# reuse the same report/plot functions.

# ---- Single-contrast helpers (mean only) ---------------------------------

# Prints mean/median/PI/HPDI of a contrast vector. Quick console check.
report_postcount_stats <- function(post_counts) {
  message("Mean: ", round(mean(post_counts$contrast),2))
  message("Median: ", round(median(post_counts$contrast),2))
  pint <- PI(post_counts$contrast)
  message("PI: 5%: ", round(pint[1],2), "; 94%: ", round(pint[2],2))
  hpdi <- HPDI(post_counts$contrast)
  message("HPDI: 5%: ", round(hpdi[1],2), "; 94%: ", round(hpdi[2],2))
}

# Extracts a fit's posterior, runs means_fn() for a 2-group estimand, and
# returns group1/group2/their difference as a data.frame. Models 1-2's
# simpler two-group pipeline.
postcounts <- function(fit, means_fn){
  post <- extract.samples(fit)
  m <- means_fn(post)
  post_counts <- data.frame(mean1 = m$mean[,1], mean2 = m$mean[,2]) %>%
    mutate(contrast = mean2 - mean1)
  report_postcount_stats(post_counts)
  post_counts
}

# SBC counterpart to postcounts(): compares means_fn()'s posterior contrast
# against the true (simulating) contrast instead of just reporting it.
contrast_from_means <- function(post, true_params, means_fn){
  m <- means_fn(post)
  mean_true <- lognormal_mean(true_params$loga, true_params$sigma^2)
  list(post_contrast = m$mean[,2] - m$mean[,1], true_contrast = mean_true[2] - mean_true[1])
}

# Density plot of one or more posterior columns (e.g. mean1/mean2/contrast),
# cropped to a quantile range with a caption noting how much was cropped.
# Keeps extreme tail draws from stretching the plotted axis into uselessness.
plot_contrast_density <- function(post_counts, quant = 1, group_name){

  quant <- if (length(quant) == 1) c(0, quant) else sort(quant)
  lower_q <- quant[1]; upper_q <- quant[2]

  post_upper <- post_counts %>%
    lapply(function(x) quantile(x, probs = upper_q)) %>%
    as_vector() %>% max()

  post_lower <- post_counts %>%
    lapply(function(x) quantile(x, probs = lower_q)) %>%
    as_vector() %>% min()

  prop_dropped <- post_counts %>%
    sapply(., function(x) sum(x > post_upper | x < post_lower) / nrow(.)) %>%
    mean() %>% {100*.} %>% signif(digits = 3)

  p <- post_counts %>%
    pivot_longer(cols = everything(), names_to = group_name) %>%
    ggplot(aes(x = value, group = !!sym(group_name))) +
    geom_density(aes(colour = !!sym(group_name), fill = !!sym(group_name)), alpha = 0.5, linewidth = 0.2) +
    coord_cartesian(xlim = c(floor(post_lower), ceiling(post_upper)))

  if(prop_dropped>0){
    dropped_above <- post_counts %>% lapply(function(x) x[x > post_upper]) %>% unlist(use.names = FALSE)
    dropped_below <- post_counts %>% lapply(function(x) x[x < post_lower]) %>% unlist(use.names = FALSE)

    fmt <- function(x) if (abs(x) >= 1e5) formatC(x, format = "e", digits = 2)
    else format(round(x), big.mark = ",", scientific = FALSE)

    range_bits <- c(
      if (length(dropped_below) > 0) paste0("below: [", fmt(min(dropped_below)), "-", fmt(max(dropped_below)), "]"),
      if (length(dropped_above) > 0) paste0("above: ", fmt(min(dropped_above)), "-", fmt(max(dropped_above)), "]")
    )

    extent_text <- if (lower_q == 0) {
      paste0("exceed the ", format(100*upper_q), "th percentile")
    } else if (upper_q == 1) {
      paste0("fall below the ", format(100*lower_q), "th percentile")
    } else {
      paste0("fall outside the ", format(100*lower_q), "th-", format(100*upper_q), "th percentile range")
    }

    p <- p + labs(caption = paste0(
      "Approximately ", prop_dropped, "% of samples ", extent_text,
      " (", paste(range_bits, collapse = "; "), ") and are not shown on this plot."))
  }
  return(p)
}

# ---- Full posterior, model-agnostic list -----------------------------------

# Every model parameter's raw posterior draws, plus whatever means_fn()
# derives from them, as a named list of small tibbles. The common input
# every contrast/variance-panel function below builds on.
post_full <- function(fit, means_fn = NULL, ...){
  post <- extract.samples(fit)
  if (!is.null(means_fn)) post <- c(post, means_fn(post, ...))

  purrr::imap(post, function(x, name){
    x <- as.matrix(x)
    colnames(x) <- if (ncol(x) == 1) name else paste0(name, "_", seq_len(ncol(x)))
    as_tibble(x)
  })
}

# Reshapes post_full() output into one long statistic/group/value tibble,
# adding a Contrast row for any 2-category parameter (or a single Population
# row for a scalar, or one row per column for >2 categories with no single
# well-defined contrast). Needed so contrast_plot_panels() has one consistent
# shape regardless of how many parameters/categories a model has.
compute_contrasts <- function(pf, keep = NULL, labels = NULL, group_levels = c("1", "2")){
  if (is.null(keep))   keep   <- names(pf)
  if (is.null(labels)) labels <- character(0)

  stat_tibble <- function(name){
    x <- as.matrix(pf[[name]])
    n_cat <- ncol(x)

    if (n_cat == 2){
      bind_rows(
        tibble(statistic = name, group = group_levels[1], value = x[,1]),
        tibble(statistic = name, group = group_levels[2], value = x[,2]),
        tibble(statistic = name, group = "Contrast",       value = x[,2] - x[,1])
      )
    } else if (n_cat == 1){
      tibble(statistic = name, group = "Population", value = x[,1])
    } else {
      map_dfr(seq_len(n_cat), function(i)
        tibble(statistic = name, group = colnames(x)[i], value = x[,i]))
    }
  }

  long <- map_dfr(keep, stat_tibble)

  display <- function(nm) unname(ifelse(nm %in% names(labels), labels[nm], nm))
  level_order <- unique(c(intersect(names(labels), keep), setdiff(keep, names(labels))))

  long %>%
    mutate(statistic = display(statistic)) %>%
    mutate(statistic = factor(statistic, levels = display(level_order)))
}

# ---- Management-Month contrasts (model 4 onward) -----------------

# Wraps a named list of posterior vectors (one per estimand) into
# compute_contrasts()'s statistic/group/value shape, for estimands that are
# already a single computed quantity (e.g. a gap-in-gaps), not a raw
# group1/group2 difference.
estimand_rows <- function(named_values){
  purrr::imap_dfr(named_values, function(v, stat)
    tibble(statistic = factor(stat, levels = names(named_values)),
           group = "Contrast", value = v))
}

# Combines paired group1/group2 estimands (e.g. May/July) with
# estimand_rows()-style Contrast-only ones (e.g. a seasonal change) into one
# ready-to-plot tibble, panel order set automatically. Saves repeating the
# same bind_rows()+factor() boilerplate in every model script.
estimand_panels <- function(pairs, extra = NULL, group_levels = c("1", "2")){
  paired <- purrr::imap_dfr(pairs, function(v, name) bind_rows(
    tibble(statistic = name, group = group_levels[1], value = v[[1]]),
    tibble(statistic = name, group = group_levels[2], value = v[[2]]),
    tibble(statistic = name, group = "Contrast",       value = v[[2]] - v[[1]])
  ))
  full <- if (is.null(extra)) paired else bind_rows(paired, estimand_rows(extra))
  full %<>% mutate(statistic = factor(statistic, levels = c(names(pairs), names(extra))))
  message('Statistic factor levels :')
  message(cat(levels(full$statistic), sep = "\n"))
  return(full)
}

# One summary row (mean/median/89% PI/HPDI) per statistic, from a
# compute_contrasts()-shaped tibble. Feeds both console/report text and
# contrast_plot_panels()'s in-panel labels.
report_contrasts_full <- function(pc_full){
  filtered <- pc_full %>% filter(group %in% c("Contrast", "Population"))
  if (nrow(filtered) == 0) {
    return(tibble(statistic = factor(character(0), levels = levels(pc_full$statistic)),
                   group = character(0), mean = numeric(0), median = numeric(0),
                   PI89_lower = numeric(0), PI89_upper = numeric(0),
                   HPDI_lower = numeric(0), HPDI_upper = numeric(0)))
  }
  filtered %>%
    group_by(statistic, group) %>%
    summarise(
      mean   = mean(value),
      median = median(value),
      PI89_lower = PI(value)[1],
      PI89_upper = PI(value)[2],
      HPDI_lower = HPDI(value)[1],
      HPDI_upper = HPDI(value)[2],
      .groups = "drop"
    )
}

# One row summary (median + 89% PI) for a single vector of draws.
cell_summary <- function(x) tibble(median = median(x), lo89 = PI(x)[1], hi89 = PI(x)[2])

# Prints report_contrasts_full()'s summary rows to the console.
message_contrasts_full <- function(pc_full){
  report_contrasts_full(pc_full) %>%
    purrr::pwalk(function(statistic, mean, median, PI89_lower, PI89_upper, HPDI_lower, HPDI_upper, ...){
      message(statistic, " -- mean: ", round(mean,2),
              "; median: ", round(median,2),
              "; 89% PI: [", round(PI89_lower,2), ", ", round(PI89_upper,2), "]",
              "; 89% HPDI: [", round(HPDI_lower,2), ", ", round(HPDI_upper,2), "]")
    })
}

# One density panel per statistic, both groups plus their Contrast
# overlaid, with a median/PI/HPDI label and an optional true-value
# reference line. The main posterior-vs-truth plot used across every
# model's validation script.
contrast_plot_panels <- function(
    pc_full, quant, group_pal, scales = "free",
    true_vals = NULL, true_vals_label = "True value",
    legend_title = "Posteriors (population mean/median)"){

  if(length(quant)!=2){
      stop("quant is not a two-value numeric vector.")
  }
  if(sum(quant<=1)!=2 | sum(quant>=0)!=2) {
    stop("quant values must be in [0,1]; lower and upper desired quantiles, e.g. c(0.005, 0.995)")
    }

  trimmed <- pc_full %>%
    group_by(statistic) %>%
    filter(value >= quantile(value, quant[1]), value <= quantile(value, quant[2])) %>%
    ungroup()

  prop_dropped <- 100 * (1 - nrow(trimmed) / nrow(pc_full))

  labels <- report_contrasts_full(pc_full) %>%
    mutate(label = paste0(
      group, " median: ", round(median,2), " ",
      "\n89% PI: [", round(PI89_lower,2), ", ", round(PI89_upper,2), "]",
      "\n89% HPDI: [", round(HPDI_lower,2), ", ", round(HPDI_upper,2), "]"))

  refactor_statistic <- function(df) df %>% mutate(statistic = factor(statistic, levels = levels(pc_full$statistic)))

  linetype_vals <- setNames("dashed", true_vals_label)

  trimmed %>%
    ggplot(aes(x = value, fill = group, colour = group)) +
    geom_density(alpha = 0.5, linewidth = 0.2) +
    geom_vline(xintercept = 0, colour = "grey50") +
    { if (!is.null(true_vals)) {
        tv <- refactor_statistic(true_vals)
        if ("group" %in% names(tv)) {
          geom_vline(data = tv, aes(xintercept = value, colour = group, linetype = true_vals_label),
                     linewidth = 0.7, show.legend = c(colour = FALSE, linetype = TRUE))
        } else {
          geom_vline(data = tv, aes(xintercept = value, linetype = true_vals_label),
                     colour = "grey20", linewidth = 0.7, show.legend = c(linetype = TRUE))
        }
      }
      } +
    geom_text(data = labels, aes(label = label), x = Inf, y = Inf,
              hjust = 1.05, vjust = 1.3, size = 2.8, colour = "grey20",
              inherit.aes = FALSE, family = "mono") +
    facet_wrap(~statistic, scales = scales, ncol = 1) +
    scale_fill_manual(values = group_pal) +
    scale_colour_manual(values = group_pal) +
    scale_linetype_manual(name = NULL, values = linetype_vals) +
    theme(legend.position = "bottom") +
    labs(x = "Species diversity", y = NULL,
         fill = legend_title, colour = legend_title,
         caption = paste0(
           "Each panel cropped independently to its own ", 100*quant[1], "th-", 100*quant[2],
           "th percentile range. ~", signif(prop_dropped, 3), "% of draws overall fall outside ",
           "their panel's range and are not shown."))
}

# Same idea as contrast_plot_panels(), but one plot per statistic with its
# own legend, stacked via patchwork -- for panels with many unrelated
# groups (e.g. per-Location, per-Year, per-Cultivar), where a single shared
# legend would be unreadable. palette: one combined named vector covering
# every group across every statistic; each panel looks up only its own subset.
variance_component_panels <- function(pc_full, quant, palette, sd_stats = character(0)){
  if (length(quant) != 2) stop("quant is not a two-value numeric vector.")
  if (sum(quant <= 1) != 2 | sum(quant >= 0) != 2) {
    stop("quant values must be in [0,1]; lower and upper desired quantiles, e.g. c(0.005, 0.995)")
  }

  trimmed <- pc_full %>%
    group_by(statistic) %>%
    filter(value >= quantile(value, quant[1]), value <= quantile(value, quant[2])) %>%
    ungroup()

  prop_dropped <- 100 * (1 - nrow(trimmed) / nrow(pc_full))

  panels <- trimmed %>%
    group_split(statistic) %>%
    purrr::map(function(df){
      stat_name <- as.character(df$statistic[[1]])
      pal_here  <- palette[intersect(names(palette), unique(df$group))]
      p <- df %>%
        ggplot(aes(x = value, fill = group, colour = group)) +
        geom_density(alpha = 0.5, linewidth = 0.2) +
        geom_vline(xintercept = 0, colour = "grey50") +
        scale_fill_manual(values = pal_here) +
        scale_colour_manual(values = pal_here) +
        theme(legend.position = "right") +
        labs(x = NULL, y = NULL, title = stat_name, fill = NULL, colour = NULL)

      # Dashed per-group mean line: a non-centered product (e.g. b[Lo]*
      # sigma_loc) is often skewed, so its peak sits off from its mean.
      # Skipped for SD panels (sigma_loc, etc.), which have no group-specific
      # mean to distinguish from the peak.
      if (!stat_name %in% sd_stats) {
        means_here <- df %>% group_by(group) %>% summarise(m = mean(value), .groups = "drop")
        p <- p + geom_vline(data = means_here, aes(xintercept = m, colour = group),
                             linetype = "dashed", linewidth = 0.5, show.legend = FALSE)
      }
      p
    })

  patchwork::wrap_plots(panels, ncol = 1) +
    patchwork::plot_annotation(caption = paste0(
      "Each panel cropped independently to its own ", 100*quant[1], "th-", 100*quant[2],
      "th percentile range. ~", signif(prop_dropped, 3), "% of draws overall fall outside ",
      "their panel's range and are not shown."))
}

# ---- Parameter-recovery contrast tables (model 4 onward) -------------------

# Posterior + true-value tables for the May gap / July gap / seasonal-change
# estimand set used by every model since Model 4. Bundles what every
# validation script's "Contrast recovery" section was retyping by hand.
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
