# postcontrast_helpers.R

# turn posterior draws (via each model's own means_X() function) into contrast 
# tables and plots. Handles 2-group diffs and multi-estimand models.

# SBC contrast_fn default for the simple 2-group case (Models 2-3): compares
# means_fn()'s posterior contrast against the true (simulating) contrast.
# Models 4+ supply their own contrast_fn instead, since their estimand isn't
# a single group1/group2 difference.
contrast_from_means <- function(post, true_params, means_fn){
  m <- means_fn(post)
  mean_true <- lognormal_mean(true_params$loga, true_params$sigma^2)
  list(post_contrast = m$mean[,2] - m$mean[,1], true_contrast = mean_true[2] - mean_true[1])
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
    # scale_linetype_manual(name = NULL, values = linetype_vals) +
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
  
  labels <- report_contrasts_full(pc_full) %>%
    mutate(label = paste0(
      group, " median: ", round(median,2), " ",
      "\n89% PI: [", round(PI89_lower,2), ", ", round(PI89_upper,2), "]",
      "\n89% HPDI: [", round(HPDI_lower,2), ", ", round(HPDI_upper,2), "]"))
  
  panels <- trimmed %>%
    group_split(statistic) %>%
    purrr::map(function(df){
      stat_name <- as.character(df$statistic[[1]])
      pal_here  <- palette[intersect(names(palette), unique(df$group))]
      # variance_component_panels() builds one standalone ggplot per
      # statistic (no shared facet_wrap the way contrast_plot_panels() has,
      # which is what lets ITS single geom_text(data=labels) auto-route by
      # facet) -- so `labels` (spanning every statistic) has to be filtered
      # down to this panel's own statistic before plotting, or every panel
      # ends up with every other panel's text stacked on top of its own.
      labels_here <- labels %>% filter(statistic == stat_name)

      p <- df %>%
        ggplot(aes(x = value, fill = group, colour = group)) +
        geom_density(alpha = 0.5, linewidth = 0.2) +
        geom_vline(xintercept = 0, colour = "grey50") +
        geom_text(data = labels_here, aes(label = label), x = Inf, y = Inf,
                  hjust = 1.05, vjust = 1.3, size = 2.8, colour = "grey20",
                  inherit.aes = FALSE, family = "mono") +
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
  
  final_plot <- patchwork::wrap_plots(panels, ncol = 1)
  
  if(sum(quant)<1) {
    final_plot <- final_plot +
      patchwork::plot_annotation(caption = paste0(
        "Each panel cropped independently to its own ", 100*quant[1], "th-", 100*quant[2],
        "th percentile range. ~", signif(prop_dropped, 3), "% of draws overall fall outside ",
        "their panel's range and are not shown."))
  }
  return(final_plot)
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

# Spearman correlation + significance for named pairs of posterior draws,
# with human-readable labels -- formalizes the ad hoc cor() collinearity
# checks used since Model 6 into one table, ready for knitr::kable().
# vars: named list of posterior vectors (caller resolves any indexing,
# e.g. yr_2022 = post$yr[,1]). pairs: list of 2-element name vectors into
# vars. labels: optional name -> display-string map.
posterior_cor_table <- function(vars, pairs, labels = character(0)){
  display <- function(nm) unname(ifelse(nm %in% names(labels), labels[nm], nm))
  
  purrr::map_dfr(pairs, function(p){
    ct <- suppressWarnings(cor.test(vars[[p[1]]], vars[[p[2]]], method = "spearman"))
    tibble(
      `Variable 1` = display(p[1]),
      `Variable 2` = display(p[2]),
      rho = unname(ct$estimate),
      p_value = ct$p.value
    )
  }) %>%
    mutate(signif = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE ~ ""
    ))
}
