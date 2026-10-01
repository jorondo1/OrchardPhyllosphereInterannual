# ---- Full posterior, model-agnostic -----------------------------------

# All posterior draws + means_fn() derived quantities, as a named list of tibbles
post_full <- function(fit, means_fn = NULL, ...){
  post <- extract.samples(fit)
  if (!is.null(means_fn)) post <- c(post, means_fn(post, ...))

  purrr::imap(post, function(x, name){
    x <- as.matrix(x)
    colnames(x) <- if (ncol(x) == 1) name else paste0(name, "_", seq_len(ncol(x)))
    as_tibble(x)
  })
}

# post_full() -> long statistic/group/value tibble
# - 2 categories: both + "Contrast" (2nd - 1st)
# - scalar: one "Population" row
# - >2 categories: one row per column, no contrast
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

  # optional display labels; factor order = labelled first, then the rest
  display <- function(nm) unname(ifelse(nm %in% names(labels), labels[nm], nm))
  level_order <- unique(c(intersect(names(labels), keep), setdiff(keep, names(labels))))
  long %>%
    mutate(statistic = display(statistic)) %>%
    mutate(statistic = factor(statistic, levels = display(level_order)))
}

# ---- Management-Month interaction contrasts -----------------

# Named list of posterior vectors -> "Contrast"-only rows
estimand_rows <- function(named_values){
  purrr::imap_dfr(named_values, function(v, stat)
    tibble(statistic = factor(stat, levels = names(named_values)),
           group = "Contrast", value = v))
}

# Paired estimands (e.g. May/July, group 1 vs 2) + their Contrast, plus
# optional Contrast-only extras (e.g. fold differences)
# - pairs: named list of list(group1_draws, group2_draws)
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

# Summary per statistic x group: mean, median, 89% PI, 89% HPDI, pd
# - pd (probability of direction, Makowski et al. 2019): posterior mass on the median's side of 0
# - default: Contrast/Population rows only; all_groups = TRUE for every row
report_contrasts_full <- function(pc_full, labels = TRUE, all_groups = FALSE){
  filtered <- if (all_groups) pc_full else pc_full %>% filter(group %in% c("Contrast", "Population"))

  # Empty input: return an empty table directly
  # - summarise() probes HPDI() on 0 rows, which errors
  # - same columns as below, so aes(label = label) still resolves
  if (nrow(filtered) == 0) {
    empty <- tibble(
      statistic = factor(character(0), levels = levels(pc_full$statistic)),
      group = character(0), mean = numeric(0), median = numeric(0),
      HPDI_lower = numeric(0), HPDI_upper = numeric(0),
      PI89_lower = numeric(0), PI89_upper = numeric(0),
      pd = numeric(0))
    return(if (labels) mutate(empty, label = character(0)) else empty)
  }

  filtered %>%
    group_by(statistic, group) %>%
    summarise(
      mean   = mean(value),
      median = median(value),
      HPDI_lower = HPDI(value)[1],
      HPDI_upper = HPDI(value)[2],
      PI89_lower = PI(value)[1],
      PI89_upper = PI(value)[2],
      pd = max(mean(value > 0), mean(value < 0)),
      .groups = "drop"
    ) %>% {
      if(labels) {
        mutate(., label = paste0(
          group, " median: ", round(median,2), " ",
          "\n89% PI: [", round(PI89_lower,2), ", ", round(PI89_upper,2), "]",
          "\n89% HPDI: [", round(HPDI_lower,2), ", ", round(HPDI_upper,2), "]"))
      } else {.}
    }

}

# Model-building report plot: one density panel per statistic
# - groups + Contrast overlaid, median/PI/HPDI label
# - optional true-value lines (calibration)
# - ratio_stats: fold-change statistics, drawn in a separate bottom panel
contrast_plot_panels <- function(
    pc_full, quant, group_pal, scales = "free",
    true_vals = NULL, true_vals_label = "True value",
    legend_title = "Posteriors (population mean/median)",
    ratio_stats = character(0), ratio_pal = NULL){

  if(length(quant)!=2){
    stop("quant is not a two-value numeric vector.")
  }
  if(sum(quant<=1)!=2 | sum(quant>=0)!=2) {
    stop("quant values must be in [0,1]; lower and upper desired quantiles, e.g. c(0.005, 0.995)")
  }

  # crop each statistic to its own quantile range
  trimmed <- pc_full %>%
    group_by(statistic) %>%
    filter(value >= quantile(value, quant[1]), value <= quantile(value, quant[2])) %>%
    ungroup()

  prop_dropped <- 100 * (1 - nrow(trimmed) / nrow(pc_full))

  labels <- report_contrasts_full(pc_full)

  refactor_statistic <- function(df) df %>% mutate(statistic = factor(statistic, levels = levels(pc_full$statistic)))

  main_levels  <- setdiff(levels(pc_full$statistic), ratio_stats)
  ratio_levels <- intersect(levels(pc_full$statistic), ratio_stats)

  # one legend across main + ratio panels
  main_groups <- trimmed %>% filter(statistic %in% main_levels) %>% pull(group) %>% unique()
  combined_pal <- c(group_pal[intersect(names(group_pal), main_groups)], ratio_pal)

  main_plot <- trimmed %>%
    filter(statistic %in% main_levels) %>%
    ggplot(aes(x = value, fill = group, colour = group)) +
    geom_density(alpha = 0.5, linewidth = 0.2) +
    geom_vline(xintercept = 0, colour = "grey50") +

    # true values: per group if a group column is given, else one line
    { if (!is.null(true_vals)) {
      tv <- refactor_statistic(true_vals)
      if ("group" %in% names(tv)) {
        geom_vline(
          data = tv,
          aes(xintercept = value, colour = group, linetype = true_vals_label),
          linewidth = 0.7, show.legend = c(colour = FALSE, linetype = TRUE))
      } else {
        geom_vline(
          data = tv,
          aes(xintercept = value, linetype = true_vals_label),
          colour = "grey20", linewidth = 0.7, show.legend = c(linetype = TRUE))
      }}} +
    geom_text(data = labels %>% filter(statistic %in% main_levels),
              aes(label = label), x = Inf, y = Inf,
              hjust = 1.05, vjust = 1.3, size = 2.8, colour = "grey20",
              inherit.aes = FALSE, family = "mono") +
    facet_wrap(~statistic, scales = scales, ncol = 1) +
    scale_fill_manual(values = combined_pal, limits = names(combined_pal), guide = "none") +
    scale_colour_manual(values = combined_pal, limits = names(combined_pal), guide = "none") +
    labs(x = "Effective number of ASVs", y = NULL)

  caption <- paste0(
    "Each panel cropped independently to its own ", 100*quant[1], "th-", 100*quant[2],
    "th percentile range. ~", signif(prop_dropped, 3), "% of draws overall fall outside ",
    "their panel's range and are not shown.")

  if (length(ratio_levels) == 0) return(main_plot + labs(caption = caption))

  # ratio panel annotation: one line per statistic, median + 89% HPDI
  ratio_label_text <- report_contrasts_full(pc_full, labels = FALSE) %>%
    filter(statistic %in% ratio_levels) %>%
    arrange(statistic) %>%
    mutate(line = paste0(
      statistic, " posterior median: ", round(median, 2),
      " (89% HPDI [", round(HPDI_lower, 2), ", ", round(HPDI_upper, 2), "])")) %>%
    pull(line) %>%
    paste(collapse = "\n")

  # Dummy rows for the main panel's groups
  # - otherwise their keys in the collected legend are blank (no data here)
  ratio_dummy <- tibble(group = main_groups, x = 1, y = 0)

  ratio_plot <- trimmed %>%
    filter(statistic %in% ratio_levels) %>%
    mutate(group = as.character(statistic)) %>%
    ggplot(aes(x = value, fill = group, colour = group)) +
    geom_density(alpha = 0.5, linewidth = 0.2) +
    geom_area(
      data = ratio_dummy, aes(x = x, y = y, fill = group, colour = group),
      alpha = 0.5, linewidth = 0.2, inherit.aes = FALSE) +
    geom_vline(xintercept = 1, colour = "grey50") +
    annotate("text", x = Inf, y = Inf, label = ratio_label_text,
             hjust = 1.05, vjust = 1.3, size = 2.8, colour = "grey20", family = "mono") +
    scale_fill_manual(values = combined_pal, limits = names(combined_pal)) +
    scale_colour_manual(values = combined_pal, limits = names(combined_pal)) +
    theme(legend.position = "bottom") +
    labs(x = "Fold change", y = NULL, fill = legend_title, colour = legend_title)

  patchwork::wrap_plots(list(main_plot, ratio_plot), ncol = 1,
                        heights = c(length(main_levels), 1)) +
    patchwork::plot_annotation(caption = caption)
}

# Like contrast_plot_panels(), but one plot per statistic, each with its own legend
# - palette: one named vector for all groups; each panel uses its own subset
# - sd_stats: SD panels, no mean line
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

  labels <- report_contrasts_full(pc_full)

  panels <- trimmed %>%
    group_split(statistic) %>%
    purrr::map(function(df){
      stat_name <- as.character(df$statistic[[1]])
      pal_here  <- palette[intersect(names(palette), unique(df$group))]

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

      # dashed per-group mean line (non-centered products are often skewed)
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

# ---- Parameter-recovery contrast tables (Mg x Mo models) -------------------

# Mean- and median-based Mg x Mo estimands from post_full()
# - pf$mean / pf$median: 4 columns in means_*() order (1/3 = May pair, 2/4 = July pair)
# - adds May/July fold differences (col 1 / col 3, col 2 / col 4)
build_pc_estimands <- function(pf, group_levels = c("1", "2")){
  one <- function(x, label){
    pair_names <- paste(c("May", "July"), label)
    estimand_panels(
      pairs = setNames(list(list(x[[1]], x[[3]]), list(x[[2]], x[[4]])), pair_names),
      extra = setNames(list(x[[1]] / x[[3]], x[[2]] / x[[4]]), names(Fold_change_palette)),
      group_levels = group_levels
    )
  }
  list(means = one(pf$mean, "mean"), medians = one(pf$median, "median"))
}

# ---- Variance partition ------------------------------------
# Inputs for both partition functions:
# - terms: named list, each a (draws x n_obs) matrix of that term's contribution
#   to the linear predictor
# - residual_var: residual variance per draw

# "By margin" shortcut: var(full) - var(full minus term), per term
# - no refit: coefficients fixed at their joint posterior values
# - NOT a Type III SS; correlated terms can come out negative
# - superseded by variance_partition_lmg() (default in model files)
variance_partition_panels <- function(terms, residual_var){
  message('Warning: This is a shortcut, not true type III ANOVA SS !!!! ')
  nm   <- names(terms)
  full <- Reduce(`+`, terms)
  explained <- apply(full, 1, var)
  total <- explained + residual_var

  marg_vals <- purrr::imap(terms, function(mat, name) explained - apply(full - mat, 1, var))

  bind_rows(
    purrr::imap_dfr(marg_vals, ~ tibble(group = .y, value = .x / total)),
    tibble(group = "Residual", value = residual_var / total)
  ) %>%
    mutate(
      statistic = "Variance partition (By margin)",
      group = factor(group, levels = c(nm, "Residual"))
    )
}

# Effect-coded Management x Season split (Gelman 2005 "batches")
# - input: combined Mg x Mo contribution (loga[Mg] + gamma*(Mo-1)), draws x n_obs
# - output: Management, Season main effects (level mean - grand mean) + interaction (remainder)
# - orthogonal when cell counts are proportional; independent of reference levels
# - grand mean dropped (constant per draw, no variance)
mgmo_effect_terms <- function(mgmo, Mg, Mo){
  grand <- rowMeans(mgmo)
  level_means <- function(g){
    lv <- sort(unique(g))
    M  <- sapply(lv, function(l) rowMeans(mgmo[, g == l, drop = FALSE]))
    M[, base::match(g, lv), drop = FALSE] - grand
  }
  mg_eff <- level_means(Mg)
  mo_eff <- level_means(Mo)
  list(
    "Management"          = mg_eff,
    "Season"              = mo_eff,
    "Management x Season" = mgmo - grand - mg_eff - mo_eff
  )
}

# LMG/Shapley partition (Lindeman et al. 1980; Groemping 2007), per draw
# - sensitivity::lmg(): R2 gain of each term, averaged over all entry orders
# - full R2 = 1 (linear predictor = sum of terms) -> shares of explained variance,
#   rescaled to shares of total (explained + residual)
# - shares >= 0; 2^K regressions per draw, so run on n_draws evenly spaced draws
variance_partition_lmg <- function(terms, residual_var, n_draws = 1000){
  if (!requireNamespace("sensitivity", quietly = TRUE)) stop("Needs the 'sensitivity' package.")
  nm <- names(terms)
  K  <- length(nm)
  draws <- unique(round(seq(1, length(residual_var), length.out = n_draws)))

  res <- t(vapply(draws, function(s){
    X <- as.data.frame(sapply(terms, function(m) m[s, ]))
    y <- rowSums(X)
    # exact fit by construction -> lm() "perfect fit" warnings, muffled
    lmg_vals <- suppressWarnings(sensitivity::lmg(X, y))$lmg[, 1]

    explained <- var(y)
    c(lmg_vals * explained, residual_var[s]) / (explained + residual_var[s])
  }, numeric(K + 1)))
  colnames(res) <- c(nm, "Residual")

  purrr::imap_dfr(as_tibble(res), ~ tibble(group = .y, value = .x)) %>%
    mutate(
      statistic = "Variance partition (Shapley/LMG)",
      group = factor(group, levels = c(nm, "Residual"))
    )
}
