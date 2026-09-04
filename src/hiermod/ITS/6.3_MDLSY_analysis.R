# 6.3_MDLSY_analysis.R -- MODEL 6 (MDLSY): posterior contrast, variance
# components, and the Bayesian R2 / variance-partition (VPC) view, run
# against the fit saved by 6.2_MDLSY_validation.R -- no refit needed.

source('src/hiermod/ITS/0_SETUP.R')
source('src/hiermod/ITS/6.1_MDLSY_model.R') # model, means_MDLSY(), variance_partition_MDLSY()
hiermod_out_dir <- "out/hiermod/ITS_6_lognormal_MDLSY"

fit_MDLSY <- readRDS(file.path(hiermod_out_dir, "fit_MDLSY.rds"))

dat <- list(
  Dv = div$Hill_1,
  Mg = idx$Mg$to_index(div$Management),
  Lo = idx$Lo$to_index(div$Location),
  Mo = idx$Mo$to_index(div$Time),
  Yr = idx$Yr$to_index(div$Year),
  Tr = idx$Tr$to_index(div$Tree_id),
  Cv = idx$Cv$to_index(div$Cultivar),
  deg_h_z = div$deg_h_z,
  precip_72h_z = div$precip_72h_z
)
dat$cell <- (dat$Mg - 1) * 2 + dat$Mo

post <- extract.samples(fit_MDLSY)
pf   <- post_full(fit_MDLSY, means_MDLSY)
m    <- pf$mean
md   <- pf$median

## Management x Season contrasts -----------------------------------------------
# Reported net of weather (deg_h_z = precip_72h_z = 0, this sample's own
# average weather -- means_MDLSY()'s default) and at reference Cultivar
# (cv[Cv] not in the estimand at all, same as Year pre-Model-6 -- see
# 6.1_MDLSY_model.R).

pc_estimands_means <- estimand_panels(
  pairs = list(`May mean` = list(m$mean_1, m$mean_3),
               `July mean` = list(m$mean_2, m$mean_4)),
  extra = list("Change in management contrast, from May to July" =
                 (m$mean_4 - m$mean_2) - (m$mean_3 - m$mean_1)),
  group_levels = c("Conventional", "Organic")
)

# fit_MDLSY@formula (not the freshly-sourced `model`) -- the real fit was
# built from model_vbc (6.2's variance-budget-calibrated priors), which
# this script never sees; the fit's own compiled formula is the only copy
# guaranteed to match what was actually fit.
save_report("fit_summary", "MDLSY", fit_MDLSY, pc_estimands_means, fit_MDLSY@formula, model_name = "The Weatherman")

p_contrast_mean <- contrast_plot_panels(
  pc_estimands_means, quant = c(0.005, 0.995), scales = 'free_y',
  group_pal = Management_palette,
  legend_title = "Posteriors (population means)"); p_contrast_mean

save_gg("fit_contrast_mean", "MDLSY", p_contrast_mean)

pc_estimands_medians <- estimand_panels(
  pairs = list(`May median` = list(md$median_1, md$median_3),
               `July median`= list(md$median_2, md$median_4)),
  extra = list("Change in management contrast, from May to July" =
                 (md$median_4 - md$median_2) - (md$median_3 - md$median_1)),
  group_levels = c("Conventional", "Organic")
)

p_contrast_median <- contrast_plot_panels(
  pc_estimands_medians, quant = c(0.005,.995), scales = 'free_y',
  group_pal = Management_palette,
  legend_title = "Posteriors (population medians)"); p_contrast_median
save_gg("fit_contrast_median", "MDLSY", p_contrast_median)

## Variance components -----------------------------------------------------------
# Same per-panel-legend treatment as Model 5's (5.3_MDLS2_analysis.R),
# extended: Year now goes through sigma_yr (properly pooled as of this
# model, was raw/additive pre-Model-6), plus a Cultivar panel (raw cv[Cv],
# stays plainly additive -- unpooled, no sigma_cv) and a weather-slope panel
# (b_deg/b_precip, each a single scalar posterior).

stat_levels <- c("Location effect (log scale)",
                  "Year effect (log scale)",
                  "Cultivar effect (log scale)",
                  "Population SD: Location / Tree / Year",
                  "Weather slope (per SD)")

pf_re <- list(
  `Location effect (log scale)` = as_tibble(as.matrix(pf$b) * pf$sigma_loc[[1]]) %>%
    setNames(as.character(idx$Lo$to_label_n(seq_along(idx$Lo$levels)))),
  `Year effect (log scale)` = as_tibble(as.matrix(pf$yr) * pf$sigma_yr[[1]]) %>%
    setNames(as.character(idx$Yr$to_label_n(seq_along(idx$Yr$levels)))),
  `Cultivar effect (log scale)` = as_tibble(pf$cv) %>%
    setNames(as.character(idx$Cv$to_label_n(seq_along(idx$Cv$levels))))
)

sigma_re <- bind_rows(
  tibble(statistic = stat_levels[4], group = "sigma_loc", value = pf$sigma_loc$sigma_loc),
  tibble(statistic = stat_levels[4], group = "sigma_tr",  value = pf$sigma_tr$sigma_tr),
  tibble(statistic = stat_levels[4], group = "sigma_yr",  value = pf$sigma_yr$sigma_yr)
)

weather_re <- bind_rows(
  tibble(statistic = stat_levels[5], group = "b_deg",    value = pf$b_deg$b_deg),
  tibble(statistic = stat_levels[5], group = "b_precip", value = pf$b_precip$b_precip)
)

re_pal <- c(
  idx$Lo$palette_n(),
  idx$Yr$palette_n(),
  idx$Cv$palette_n(),
  sigma_loc = Management_palette[["Population"]],
  sigma_tr  = "grey40",
  sigma_yr  = "#8C6D31",
  b_deg     = "#4E79A7",
  b_precip  = "#E15759"
)

pc_random_effects <- bind_rows(
  compute_contrasts(pf_re, keep = names(pf_re)), sigma_re, weather_re) %>%
  mutate(statistic = factor(statistic, levels = stat_levels))

p_random_effects <- variance_component_panels(
  pc_random_effects, quant = c(0.005, 0.995), palette = re_pal); p_random_effects

save_gg("fit_variance_components", "MDLSY", p_random_effects, width = 8, height = 15)

## Variance partition (VPC) / Bayesian R2 -----------------------------------

# variance_partition_MDLSY() (6.1_MDLSY_model.R): 5 components -- Explained
# (fixed effects, i.e. Bayesian R2), Location, Tree, Year, Residual -- each
# a share of total variance, summing to 1 per posterior draw (a genuine
# composition). Stacked bar of posterior medians for the compact headline
# view; ridge densities above, x-axis aligned to the bar, for the full
# posterior shape per component -- R2 and VPC are the same decomposition,
# just two views of it.

vp <- variance_partition_MDLSY(post, dat)

vp_levels <- c("Explained (fixed effects)", "Location", "Tree", "Year", "Residual")
vp <- vp %>% mutate(group = factor(group, levels = vp_levels))

vp_pal <- c(
  `Explained (fixed effects)` = "#499894",
  Location = "#F28E2B",
  Tree     = "#B07AA1",
  Year     = "#8C6D31",
  Residual = "grey50"
)

vp_medians <- vp %>% group_by(group) %>% summarise(median = median(value), .groups = "drop")

p_vp_bar <- vp_medians %>%
  ggplot(aes(x = median, y = "Variance partition", fill = group)) +
  geom_col(position = "stack", width = 0.6) +
  scale_fill_manual(values = vp_pal) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(x = "Share of total variance (posterior medians)", y = NULL, fill = NULL)

p_vp_ridge <- vp %>%
  ggplot(aes(x = value, y = group, fill = group)) +
  ggridges::geom_density_ridges(alpha = 0.7, colour = "white") +
  scale_fill_manual(values = vp_pal) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  guides(fill = "none") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  labs(x = NULL, y = NULL)

p_variance_partition <- p_vp_ridge / p_vp_bar +
  patchwork::plot_layout(heights = c(4, 1)); p_variance_partition

save_gg("fit_variance_partition", "MDLSY", p_variance_partition, width = 8, height = 8)
