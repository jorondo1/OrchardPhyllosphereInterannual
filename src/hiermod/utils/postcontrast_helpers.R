# postcontrast_helpers.R

# posterior contrast reporting/plotting. Each model script defines its own
# means_X(post) -> list(mean = cbind(v1, v2), ...) (see lognormal_mean()'s
# note in hiermod_core.R for why total_var stays inline per model); one
# 2-column matrix per derived statistic, column order matching Mg (or
# whichever 2-level grouping the model uses).
#   - Single-contrast helpers (postcounts()/plot_contrast_density()): models
#     1-2's simulation/SBC pipeline, mean only.
#   - Full-posterior helpers (post_full()/compute_contrasts()): model 3,
#     any number of statistics (mean, median, sigma, plus scalars like
#     sigma_g) in one panel set, built from 2-category means_X() output.
#   - Hand-built multi-estimand helpers (post_full()/estimand_rows()): model
#     4 onward, whenever the estimand is plural and doesn't reduce to a
#     single group1/group2/Contrast shape (e.g. MDLS's May/July/Season-change
#     gaps) -- which cells combine into which estimand is written by hand
#     per script, estimand_rows() just does the generic stacking step.


# ---- Single-contrast helpers (mean only) ---------------------------------

# Everything below is the generic part, identical regardless of how many
# variance components means_X() sums over.

# Reporter function; simple verbose for posterior counts extractor/converter below
report_postcount_stats <- function(post_counts) {
  message("Mean: ", round(mean(post_counts$contrast),2))
  message("Median: ", round(median(post_counts$contrast),2))
  pint <- PI(post_counts$contrast)
  message("PI: 5%: ", round(pint[1],2), "; 94%: ", round(pint[2],2))
  hpdi <- HPDI(post_counts$contrast)
  message("HPDI: 5%: ", round(hpdi[1],2), "; 94%: ", round(hpdi[2],2))
}

# Real-fit backtransform + report: extract.samples(), run the model's own
# means_fn(post) -> list(mean = cbind(v1, v2), ...), package into a
# data.frame with `contrast`.
postcounts <- function(fit, means_fn){
  post <- extract.samples(fit)
  m <- means_fn(post)                                      # -> mean[,1]/mean[,2]
  post_counts <- data.frame(mean1 = m$mean[,1], mean2 = m$mean[,2]) %>%
    mutate(contrast = mean2 - mean1)                        # group2 - group1
  report_postcount_stats(post_counts)
  post_counts
}

# SBC version: same means_fn(post), but compares against true_params instead
# of reporting. true_params$loga/$sigma are always the observation-level
# prior draws -- population-level components (sigma_loc, ...) enter the
# simulator as actual random offsets (see each model's sim_div_*()), not this
# analytic backtransform, so the true mean only needs the single-component
# form of lognormal_mean().

contrast_from_means <- function(post, true_params, means_fn){
  m <- means_fn(post)
  mean_true <- lognormal_mean(true_params$loga, true_params$sigma^2)
  list(post_contrast = m$mean[,2] - m$mean[,1], true_contrast = mean_true[2] - mean_true[1])
}

# Density-overlay for a data.frame of named columns (e.g. mean1/mean2/contrast).
plot_contrast_density <- function(post_counts, quant = 1, group_name){
  
  # quant: either a single upper cutoff (0.99: crops only the right tail) or
  # c(lower, upper) (e.g. c(0.005, 0.995): crops both tails, for quantities
  # like `contrast` that can go negative).
  quant <- if (length(quant) == 1) c(0, quant) else sort(quant)
  lower_q <- quant[1]; upper_q <- quant[2]
  
  post_upper <- post_counts %>%
    lapply(function(x) quantile(x, probs = upper_q)) %>%  # per-column upper quantile
    as_vector() %>% max()                                 # widest column wins

  post_lower <- post_counts %>%
    lapply(function(x) quantile(x, probs = lower_q)) %>%  # per-column lower quantile
    as_vector() %>% min()                                 # widest column wins

  prop_dropped <- post_counts %>%
    sapply(., function(x) sum(x > post_upper | x < post_lower) / nrow(.)) %>%  # % cropped per column
    mean() %>% {100*.} %>% signif(digits = 3)                                  # average across columns, as a %
  
  p <- post_counts %>%
    pivot_longer(cols = everything(), names_to = group_name) %>%
    ggplot(aes(x = value, group = !!sym(group_name))) +
    geom_density(aes(colour = !!sym(group_name), fill = !!sym(group_name)), alpha = 0.5, linewidth = 0.2) +
    coord_cartesian(xlim = c(floor(post_lower), ceiling(post_upper)))
  
  if(prop_dropped>0){
    # range of the dropped (off-plot) tail(s) -- lets a reader gauge how
    # extreme "not shown" actually is, not just how much of it there is.
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

# Every model parameter's raw posterior draws, plus whatever it derives
# from them using a backtransform function if provided. 
# Returned as a named list, one small tibble per parameter, i-e 1 column 
# for a scalar like sigma_g, 2+ for an indexed/derived one

post_full <- function(fit, means_fn = NULL, ...){
  post <- extract.samples(fit)                              # raw draws
  if (!is.null(means_fn)) post <- c(post, means_fn(post, ...))    # + derived quantities

  purrr::imap(post, function(x, name){                       # one tibble per name
    x <- as.matrix(x)
    colnames(x) <- if (ncol(x) == 1) name else paste0(name, "_", seq_len(ncol(x)))
    as_tibble(x)
  })
}

# Long tibble (statistic/group/value), ready for contrast_plot_panels(). For
# each name in `keep`: a 2-column entry (2 categories, e.g. mean's Mg
# columns) becomes group1/group2/Contrast rows; a 1-column/scalar entry
# (e.g. sigma_g) becomes a single Contrast row -- it has no pair to take a
# difference against, so it just rides along as its own posterior. A >2
# column entry (e.g. g[Lo] across 5 Locations) has no single well-defined
# contrast either -- each category rides along under its own column name
# (e.g. "g_1".."g_5"), with no Contrast row, purely for visual comparison.
#
# labels renames the raw names for display (e.g. mean -> "Mean diversity")
# and, via its own order, sets the panel display order; names absent from
# labels keep their raw name and follow, in `keep`'s order. keep = NULL
# (default) keeps every parameter post_full() returned.
compute_contrasts <- function(pf, keep = NULL, labels = NULL, group_levels = c("1", "2")){
  if (is.null(keep))   keep   <- names(pf)          # everything post_full() returned
  if (is.null(labels)) labels <- character(0)       # keeps %in%/indexing below well-defined

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
      tibble(statistic = name, group = "Population", value = x[,1])   # population-level, not a group difference
    } else {
      # no pairwise contrast for >2 categories -- one row block per column,
      # group = that column's own name, no Contrast row.
      map_dfr(seq_len(n_cat), function(i)
        tibble(statistic = name, group = colnames(x)[i], value = x[,i]))
    }
  }

  long <- map_dfr(keep, stat_tibble)                 # one statistic at a time, then stack

  display <- function(nm) unname(ifelse(nm %in% names(labels), labels[nm], nm))
  level_order <- unique(c(intersect(names(labels), keep), setdiff(keep, names(labels))))

  long %>%
    mutate(statistic = display(statistic)) %>%                          # rename where labelled
    mutate(statistic = factor(statistic, levels = display(level_order)))  # + panel order
}


# ---- Management-Month contrasts (model 4 onward) -----------------

# Model 4+'s estimands aren't a single group1/group2/Contrast shape anymore
# (e.g. MDLS's May gap, July gap, Seasonal change).

# The followming handles a generic step: each element of `named_values`
# is one estimand's full posterior, so this just wraps each into a 
# Contrast row and stacks them, the same  statistic/group/value shape 
# report_contrasts_full()/contrast_plot_panels() (below) already expect, unmodified.

estimand_rows <- function(named_values){
  purrr::imap_dfr(named_values, function(v, stat)
    tibble(statistic = factor(stat, levels = names(named_values)),
           group = "Contrast", value = v))
}

# Combines gap_tibble()-style paired panels (raw group1/group2 + their
# Contrast, e.g. May/July) with estimand_rows()-style Contrast-only panels
# (e.g. a second-order "change in gap") into one ready-to-plot tibble, with
# panel order set automatically from names(pairs) then names(extra) --
# no manual mutate(statistic = factor(...)) needed at the call site.
#   pairs: named list of list(group1_values, group2_values), one per panel
#          that should show both raw groups plus their Contrast.
#   extra: named list of already-computed contrast vectors (as passed to
#          estimand_rows()), for panels with no raw groups to show --
#          appended after `pairs`, in the given order.
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

# One row per statistic, summarizing its Contrast draws (or, for a scalar
# parameter with no group to contrast, its Population draws -- same one
# summary row per statistic either way). `group` is kept in the output (not
# just `statistic`) so callers -- contrast_plot_panels()'s label text below --
# can tell which of the two it actually is; every statistic here only ever
# has one or the other, never both, so this doesn't multiply rows.
#
# A pc_full built entirely from >2-column entries (e.g. Location/Year effect
# panels -- every group its own named category, no Contrast/Population row
# anywhere) filters down to 0 rows here; dplyr's summarise() still runs a
# type-inference pass on the empty `value` before checking there are 0
# groups, and HPDI(numeric(0)) errors on that -- hence the early return.
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

# One row summary (median + 89% PI) for a single vector of draws -- used to
# build a per-cell interaction_df (median diversity by Management x Month).
cell_summary <- function(x) tibble(median = median(x), lo89 = PI(x)[1], hi89 = PI(x)[2])

message_contrasts_full <- function(pc_full){
  report_contrasts_full(pc_full) %>%
    purrr::pwalk(function(statistic, mean, median, PI89_lower, PI89_upper, HPDI_lower, HPDI_upper, ...){
      message(statistic, " -- mean: ", round(mean,2),
              "; median: ", round(median,2),
              "; 89% PI: [", round(PI89_lower,2), ", ", round(PI89_upper,2), "]",
              "; 89% HPDI: [", round(HPDI_lower,2), ", ", round(HPDI_upper,2), "]")
    })
}

# One panel per statistic, each showing both group posteriors plus their
# Contrast overlaid on a shared axis, median + 89% PI annotated top-right.
# Contrast uses a neutral colour (group_pal["Contrast"]), not an RGB blend --
# a blend has no established reading as "the delta between these two."
#
# true_vals (optional): one reference vline per statistic/group (a
# Parameter-recovery true value, or -- Model 5 scripts -- a raw observed
# mean/median). Its `statistic` column is re-factored to pc_full's own level
# order first: geom_vline()'s data has no factor levels of its own, so an
# unfactored `statistic` on a second layer makes facet_wrap() silently fall
# back to alphabetical panel order (a general ggplot2 gotcha, not specific
# to true_vals). Optional `group` column colours each line via group_pal; a
# panel missing one group's row just skips that line, doesn't touch the
# shared Posteriors legend. No `group` column -> single neutral colour.
#
# true_vals_label names the always-shown linetype legend entry for this
# layer ("True value" / "Observed") -- its own aesthetic, so it's never
# affected by which groups have a line in a given panel.
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
    group_by(statistic) %>%                            # each panel crops its own tails...
    filter(value >= quantile(value, quant[1]), value <= quantile(value, quant[2])) %>%
    ungroup()

  prop_dropped <- 100 * (1 - nrow(trimmed) / nrow(pc_full))  # overall % of draws cropped out

  labels <- report_contrasts_full(pc_full) %>%          # one summary row per statistic...
    mutate(label = paste0(                              # ...formatted into the panel annotation text
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

# One ggplot PER statistic, each with its own right-side legend, stacked
# vertically via patchwork -- for panels with many unrelated colours per
# statistic (e.g. one per Location, a different set per Year, another for
# Cultivar): contrast_plot_panels()'s single shared bottom legend turns into
# an unreadably long combined list once there are this many categories
# across panels. Each sub-plot already has its own independent scale (no
# facet_wrap involved), so there's no `scales=` argument to pass here --
# every panel is free by construction.
#
# palette: ONE combined named vector covering every group across every
# statistic (e.g. c(idx$Lo$palette_n(), idx$Yr$palette_n(), sigma_loc = ...));
# each panel looks up only the subset of names it actually uses.
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

      # Non-SD panels here are mostly non-centered products (e.g. b[Lo]*
      # sigma_loc, yr[Yr]*sigma_yr): a symmetric z-score times a positive,
      # right-skewed sigma posterior yields a skewed *derived* posterior
      # whose peak (mode) sits off-center from its mean -- a thin dashed
      # mean line per group makes that distinction visible instead of
      # leaving the reader to read central tendency off the peak, which
      # understates whichever direction the skew leans. Not added for SD
      # panels (sigma_loc/sigma_tr/sigma_yr etc.) -- those are direct
      # positive-scale parameters, not group-varying derived effects, so
      # there's no group-specific mean to distinguish from the peak.
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
