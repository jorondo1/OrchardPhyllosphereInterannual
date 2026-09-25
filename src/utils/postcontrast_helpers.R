# ---- Full posterior, model-agnostic -----------------------------------

# Every model parameter's raw posterior draws, plus whatever means_fn()
# derives from them, as a named list of small tibbles.
post_full <- function(fit, means_fn = NULL, ...){
  post <- extract.samples(fit)
  if (!is.null(means_fn)) post <- c(post, means_fn(post, ...))
  
  purrr::imap(post, function(x, name){
    x <- as.matrix(x)
    colnames(x) <- if (ncol(x) == 1) name else paste0(name, "_", seq_len(ncol(x)))
    as_tibble(x)
  })
}

# Reshapes post_full() out into  long statistic/group/value tibble,
# add 'contrast' row for any 2-category parameter (or a single Population
# row for a scalar, or one row per column for >2 categories with no single
# well-defined contrast). 
# Used by contrast_plot_panels() 

compute_contrasts <- function(pf, keep = NULL, labels = NULL, group_levels = c("1", "2")){
  if (is.null(keep))   keep   <- names(pf)
  if (is.null(labels)) labels <- character(0)
  
  stat_tibble <- function(name){
    x <- as.matrix(pf[[name]])
    n_cat <- ncol(x)
    
    # Compute contrast when 2 categories
    if (n_cat == 2){ 
      bind_rows(
        tibble(statistic = name, group = group_levels[1], value = x[,1]),
        tibble(statistic = name, group = group_levels[2], value = x[,2]),
        tibble(statistic = name, group = "Contrast",       value = x[,2] - x[,1])
      )
    } else if (n_cat == 1){ # otherwise it's population level
      tibble(statistic = name, group = "Population", value = x[,1])
    } else { # if multiple categories, one per, but no contrast
      map_dfr(seq_len(n_cat), function(i)
        tibble(statistic = name, group = colnames(x)[i], value = x[,i]))
    }
  }
  
  # only keep specified vars
  long <- map_dfr(keep, stat_tibble)
  
  display <- function(nm) unname(ifelse(nm %in% names(labels), labels[nm], nm))
  level_order <- unique(c(intersect(names(labels), keep), setdiff(keep, names(labels))))
  # Add stat 
  long %>%
    mutate(statistic = display(statistic)) %>%
    mutate(statistic = factor(statistic, levels = display(level_order)))
}

# ---- Management-Month interaction contrasts -----------------

# Helper:
# Wraps a named list of posterior vectors (one per estimand) into
# compute_contrasts()'s statistic/group/value shape
estimand_rows <- function(named_values){
  purrr::imap_dfr(named_values, function(v, stat)
    tibble(statistic = factor(stat, levels = names(named_values)),
           group = "Contrast", value = v))
}

# Combines paired group1/group2 estimands (e.g. May/July) Contrast-only ones 
# (e.g. a seasonal change) into tibble

# Pairs is a list of named lists of two long posterior tibble columns 
# e.g. two columns of a tibble from the output of post_full

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

# One summary row (mean/median/89% PI/HPDI/pd) per statistic, from a
# compute_contrasts() output  tibble.

# pd ("probability of direction", Makowski et al. 2019): fraction of
# posterior mass on the same side of 0 as the median 

report_contrasts_full <- function(pc_full, labels = TRUE, all_groups = FALSE){
  filtered <- if (all_groups) pc_full else pc_full %>% filter(group %in% c("Contrast", "Population"))
  
  # Debug with Claude Code::::::::::::::
  # summarise() still probes these aggregation expressions on a degenerate
  # 0/1-row input to infer output types even when `filtered` has 0 matching
  # rows (e.g. every variance_component_panels() call whose groups are
  # custom labels like Year/Cultivar/Location, not literally "Contrast"/
  # "Population") . HPDI() throws ("obj must have nsamp > 1") on that probe
  # input, so this has to short-circuit before summarise() ever runs, not
  # rely on 0-row grouping to no-op safely (confirmed: it doesn't).
  if (nrow(filtered) == 0) {
    return(tibble(
      statistic = factor(character(0), levels = levels(pc_full$statistic)),
      group = character(0), mean = numeric(0), median = numeric(0),
      HPDI_lower = numeric(0), HPDI_upper = numeric(0),
      PI89_lower = numeric(0), PI89_upper = numeric(0),
      pd = numeric(0)))
  }
  # /debug :::::::::::: thank you CLaude
  
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

# Standardised plot for reporting in model building:

# One density panel per statistic, both groups plus their Contrast
# overlaid, with a median/PI/HPDI label and an optional true-value
# reference line (for predictive checks). 

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
  
  trimmed <- pc_full %>%
    group_by(statistic) %>%
    filter(value >= quantile(value, quant[1]), value <= quantile(value, quant[2])) %>%
    ungroup()
  
  prop_dropped <- 100 * (1 - nrow(trimmed) / nrow(pc_full))
  
  # Stats label for reporting;
  labels <- report_contrasts_full(pc_full)
  
  refactor_statistic <- function(df) df %>% mutate(statistic = factor(statistic, levels = levels(pc_full$statistic)))
  
  main_levels  <- setdiff(levels(pc_full$statistic), ratio_stats)
  ratio_levels <- intersect(levels(pc_full$statistic), ratio_stats)
  
  # Unified legend across main_plot + ratio_plot (e.g. Contrast/Conventional/
  # Organic + May/July fold difference, 5 entries)
  main_groups <- trimmed %>% filter(statistic %in% main_levels) %>% pull(group) %>% unique()
  combined_pal <- c(group_pal[intersect(names(group_pal), main_groups)], ratio_pal)
  
  # Main plot
  main_plot <- trimmed %>%
    filter(statistic %in% main_levels) %>%
    ggplot(aes(x = value, fill = group, colour = group)) +
    geom_density(alpha = 0.5, linewidth = 0.2) +
    geom_vline(xintercept = 0, colour = "grey50") +
    
    # Add true valus if provided:
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
    # contrast labels 
    geom_text(data = labels %>% filter(statistic %in% main_levels), 
              aes(label = label), x = Inf, y = Inf,
              hjust = 1.05, vjust = 1.3, size = 2.8, colour = "grey20",
              inherit.aes = FALSE, family = "mono") +
    facet_wrap(~statistic, scales = scales, ncol = 1) +
    scale_fill_manual(values = combined_pal, limits = names(combined_pal), guide = "none") +
    scale_colour_manual(values = combined_pal, limits = names(combined_pal), guide = "none") +
    labs(x = "Effective number of ASVs", y = NULL)
  
  # Custom caption if dropping thin tails with quant parameter
  caption <- paste0(
    "Each panel cropped independently to its own ", 100*quant[1], "th-", 100*quant[2],
    "th percentile range. ~", signif(prop_dropped, 3), "% of draws overall fall outside ",
    "their panel's range and are not shown.")
  
  if (length(ratio_levels) == 0) return(main_plot + labs(caption = caption))
  
  # Ratio/fold-change statistics: Median and 89% HPDI
  
  ratio_label_text <- report_contrasts_full(pc_full, labels = FALSE) %>%
    filter(statistic %in% ratio_levels) %>%
    arrange(statistic) %>%
    # one-liner stats per statistic 
    mutate(line = paste0(
      statistic, " posterior median: ", round(median, 2),
      " (89% HPDI [", round(HPDI_lower, 2), ", ", round(HPDI_upper, 2), "])")) %>%
    pull(line) %>%
    paste(collapse = "\n")
  
  # ratio_plot (last panel in the stack) 
  
  # Claude debug :::::::::::::: 
  # collect legends in a patchwork where each plot shares the levels, but 
  # the input data has 0 lines for some of these levels; will produce a 
  # blank/uncoloured legend, even when using explicit scale_* statements
  # with the same palettes.
  ratio_dummy <- tibble(group = main_groups, x = 1, y = 0)
  
  ratio_plot <- trimmed %>%
    filter(statistic %in% ratio_levels) %>%
    mutate(group = as.character(statistic)) %>%
    ggplot(aes(x = value, fill = group, colour = group)) +
    geom_density(alpha = 0.5, linewidth = 0.2) +
    # this is the fix: 
    # add a geom_area with same aes, but with other plot's data:
    geom_area(
      data = ratio_dummy, aes(x = x, y = y, fill = group, colour = group),
      alpha = 0.5, linewidth = 0.2, inherit.aes = FALSE) +
    #/debug :::::::: thanks claude
    geom_vline(xintercept = 1, colour = "grey50") +
    annotate("text", x = Inf, y = Inf, label = ratio_label_text,
             hjust = 1.05, vjust = 1.3, size = 2.8, colour = "grey20", family = "mono") +
    scale_fill_manual(values = combined_pal, limits = names(combined_pal)) +
    scale_colour_manual(values = combined_pal, limits = names(combined_pal)) +
    theme(legend.position = "bottom") +
    labs(x = "Fold change", y = NULL, fill = legend_title, colour = legend_title)
  
  # robust patchwork
  patchwork::wrap_plots(list(main_plot, ratio_plot), ncol = 1,
                        heights = c(length(main_levels), 1)) +
    patchwork::plot_annotation(caption = caption)
}

# Same idea as contrast_plot_panels(), but one plot per statistic with its
# own legend, stacked via patchwork. palette: one combined named vector covering
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
  
  labels <- report_contrasts_full(pc_full) 
  
  panels <- trimmed %>%
    # one standalone ggplot per statistic:
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
      
      # Dashed per-group mean line: a non-centered product (e.g. b[Lo]*
      # sigma_loc) is often skewed.Skipped for SD panels.
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


# One call replaces the copy-pasted means/medians estimand_panels() pair in
# every *.4_*_analysis.R script. pf: post_full() output (needs $mean/$median,
# each a 4-column tibble in means_*()'s own cbind order: cols 1/3 = the May
# pair, cols 2/4 = the July pair, for whatever two groups means_*() cbinds).
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

# terms: named list, each a (draws x n_obs) matrix of that term's
# per-observation contribution to the linear predictor. 
#
# residual_var: length-draws vector (the noise term)
#
# Reports a marginal ("by margin": var(full) - var(full minus this term),
# order-free) share for every term, plus Residual. NB this is NOT the same
# guarantee as a real ANOVA/PERMANOVA Type III SS! 

# !!!!! THIS IS A SHORTCUT

# CLAUDE explains::::::::::::::
# every term's coefficients are already fixed at their joint
# full-model posterior values, and we just zero out one term's own
# contribution rather than re-estimating anything. So a term whose raw
# (not mean-centered) per-observation values are strongly (anti-)correlated
# with the rest of the linear predictor can come out negative -- this isn't
# a bug, it's a real signature of non-orthogonal/confounded term coding
# under this no-refit shortcut, and it's exactly why Management x Season's
# main effects + interaction get combined back into one term (see
# variance_partition_MDSTYCV()'s own comment) rather than split three ways.
# A real Shapley/LMG (order-averaged) decomposition would be closer to the
# ANOVA convention -- deferred for now.
# /thanks claude :::::::::::::::

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
