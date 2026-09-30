# Alpha diversity figuress

pacman::p_load(tidyverse, patchwork, ggridges)

source('src/hiermod/0_SETUP.R') # also loads postcontrast_helpers.R, saver_functions.R

## Setup -----------------------------------------------

# Variance partition
varpart <- bind_rows(
  readRDS("out/hiermod/16S_7_tree_full_MDSTYCV/varpart_MDSTYCV.rds") %>%
    mutate(Kingdom = 'Bacteria'),
  readRDS("out/hiermod/ITS_7_tree_full_MDSTYCV/varpart_MDSTYCV.rds")%>%
    mutate(Kingdom = 'Fungi')
)
# Posterior parameters
effects <- bind_rows(
  readRDS("out/hiermod/16S_7_tree_full_MDSTYCV/effect_panels_MDSTYCV.rds") %>%
    mutate(Kingdom = 'Bacteria'),
  readRDS("out/hiermod/ITS_7_tree_full_MDSTYCV/effect_panels_MDSTYCV.rds") %>%
    mutate(Kingdom = 'Fungi')
)

## Shared series colours --------------------------------------------------------
# One legend for every panel: the two kingdoms, plus the Bacteria - Fungi
# contrast (posterior-parameter panels only).
series_pal <- c( Contrast = "#B82D2C", Bacteria = "#2C7FB8", Fungi = "#B8B62C")

scales_series <- list(
  scale_fill_manual(values = series_pal, limits = names(series_pal)),
  scale_colour_manual(values = series_pal, limits = names(series_pal))
)

ridge_panel <- function(df, xlab, scale = 0.8, show_legend = NA, labels = waiver()){
  df %>%
    ggplot(aes(x = value, y = group, fill = series, colour = series,
               group = interaction(group, series), height = after_stat(ndensity))) +
    ggridges::geom_density_ridges(
      stat = "density", # per-ridge bandwidth (default stat uses one joint bandwidth)
      alpha = 0.8, linewidth = 0.2,
      scale = scale, rel_min_height = 0.01,
      # TRUE keeps keys for levels absent from this panel's data (e.g. Contrast in A)
      show.legend = show_legend) +
    geom_vline(xintercept = 0, colour = "grey50", linetype = "dashed") +
    # ridges rise `scale` row-units above their baseline and nothing is drawn
    # below it: pad the top by the ridge height, the bottom only slightly
    scale_y_discrete(expand = expansion(add = c(0.1, scale + 0.05)), labels = labels) +
    scales_series +
    labs(x = xlab, y = NULL, fill = NULL, colour = NULL) +
    theme(panel.grid.minor = element_blank())
}

## Posterior parameters: Bacteria, Fungi and Bacteria - Fungi -------------------
# Fits are independent, so draws are paired by index to get the between-kingdom
# contrast per parameter level. Caveat: the Read count contrast compares different
# variables (seq_depth_z is each marker's own z-scored log read count).
# statistic is character in some panels, factor in others -- harmonised first.

params_long <- effects %>%
  mutate(pc_full = purrr::map(pc_full, \(x) mutate(x, statistic = as.character(statistic)))) %>%
  select(Kingdom, pc_full) %>%
  unnest(pc_full) %>%
  mutate(
    statistic = str_remove(statistic, fixed(" (log scale)")),
    statistic = recode(statistic, sigma = "Residual SD", sigma_tr = "Between-tree SD"),
    group     = recode(
      group, 
      "Seq. depth" = "Read count", 
      "Contrast" = "Organic - Conventional",
      "Population" = "Tree"))

params_paired <- params_long %>%
  group_by(Kingdom, statistic, group) %>%
  mutate(draw = row_number()) %>%
  ungroup() %>%
  pivot_wider(names_from = Kingdom, values_from = value) %>%
  mutate(Contrast = Bacteria - Fungi) %>%
  pivot_longer(c(Bacteria, Fungi, Contrast), names_to = "series")

params <- params_paired %>%
  # plot only: each ridge's own 0.5% tails dropped (shapes differ a lot, esp. the SDs)
  group_by(statistic, group, series) %>%
  filter(value <= quantile(value, 0.995), value >= quantile(value, 0.005)) %>%
  ungroup() %>%
  mutate(group  = factor(group, levels = rev(unique(params_long$group))), # first level on top
         series = factor(series, levels = names(series_pal)))

vc_stats <- c("Residual SD", "Between-tree SD")

# B. 
p_effects <- params %>%
  filter(!statistic %in% vc_stats) %>%
  ridge_panel("Log-fold differences / effect of a 1-SD deviation on diversity") +
  theme(legend.position = "none")

# C.
p_vc <- params %>%
  filter(statistic %in% vc_stats & group != 'Organic - Conventional') %>%
  ridge_panel("Standard deviation (log scale)",
              # sigma with the level as subscript (plotmath)
              labels = \(x) as.expression(lapply(x, \(v) bquote(sigma[.(v)])))) +
  theme(legend.position = "none")

## A. Variance partition ------------------------------------------------------------

# One ridge row per variance component, Bacteria and Fungi overlaid on it
p_varpart <- varpart %>%
  mutate(series = factor(Kingdom, levels = names(series_pal)),
         group  = fct_recode(group, "Read count" = "Reads count")) %>%
  ridge_panel("Fraction of total variance (bayesian R-squared)", scale = 0.8,
              show_legend = TRUE) +
  # the figure's only legend: boxed, in panel A's empty bottom-right corner
  theme(legend.position = "inside",
        legend.position.inside = c(0.9, 0.4),
        legend.justification = c(1, 0),
        legend.background = element_rect(colour = "black", linewidth = 0.3, fill = "white"))

## Combined -----------------------------------------------------------------------

p_supp <- p_varpart / p_effects / p_vc +
  plot_layout(heights = c(10, 12, 3)) +
  plot_annotation(tag_levels = "A") &
  theme(axis.title.x = element_text(size = 8)); p_supp

ggsave(plot = p_supp, filename = "out/manuscript/2_alpha_var.pdf", bg = "white",
       width = 2400, height = 3200, units = 'px', dpi = 300)

## Summary tables ----------------------------------------------------------------
# Same format as the *_results_report_MDSTYCV.html reports (save_posterior_kable()),
# from the untrimmed draws: one table for the variance partition, one for the
# posterior parameters.

# Explained share (Bayesian R2) per draw = 1 - that draw's residual share
varpart_tbl <- bind_rows(
  varpart %>%
    transmute(statistic = fct_recode(group, "Read count" = "Reads count"),
              group = Kingdom, value),
  varpart %>%
    filter(group == "Residual") %>%
    transmute(statistic = "Explained (Bayesian R2)", group = Kingdom, value = 1 - value)
) %>%
  mutate(statistic = factor(statistic, levels = unique(statistic)),
         group     = factor(group, levels = names(series_pal)))

save_posterior_kable(
  "2_alpha_supp_report_MDSTYCV", "varpart", varpart_tbl,
  dir = "out/manuscript", prefix = "", three_digits = TRUE,
  caption = paste(
    "Variance partition (Bayesian R2): posterior median and 89% HPDI of each term's",
    "share of total variance, per kingdom. Explained = total share of all model terms",
    "(Bayesian R2 = 1 - Residual). pd = probability of direction."))

params_tbl <- params_paired %>%
  transmute(statistic = group, group =  factor(series, levels = names(series_pal)), value) 
  # transmute(statistic, group,# = paste0(statistic, ": ", group),
  #          # statistic = factor(statistic, levels = unique(statistic)),
  #           group = factor(series, levels = names(series_pal)), value)

save_posterior_kable(
  "2_alpha_supp_report_MDSTYCV", "parameters", params_tbl,
  dir = "out/manuscript", prefix = "",
  caption = paste(
    "Posterior parameters summary statistics. Median and 89% HPDI (log scale). Contrast = Bacteria - Fungi, from draws of the two independent fits paired by index (approximate: shared samples not modelled). Read count is each marker's own z-scored log read count, so its contrast compares different variables. pd = probability of direction"))
