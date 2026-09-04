# 5b.3_MDLS2_shifted_analysis.R -- MODEL 5, SHIFTED: posterior contrast and
# variance components, run against the fit saved by
# 5b.2_MDLS2_shifted_validation.R -- no refit needed.

source('src/hiermod/ITS/0_SETUP.R')
source('src/hiermod/ITS/5.1_MDLS2_model.R') # model, means_MDLS2()
hiermod_out_dir <- "out/hiermod/ITS_5_lognormal_MDLS2_shifted"

fit_MDLS2_shifted <- readRDS(file.path(hiermod_out_dir, "fit_MDLS2_shifted.rds"))

pf <- post_full(fit_MDLS2_shifted, means_MDLS2, shift = 1)
m  <- pf$mean
md <- pf$median

## Management x Season contrasts -----------------------------------------------

pc_estimands_means <- estimand_panels(
  pairs = list(`May mean` = list(m$mean_1, m$mean_3),
               `July mean` = list(m$mean_2, m$mean_4)),
  extra = list("Seasonal change in mean diversity difference between management, from May to July" =
                 (m$mean_4 - m$mean_2) - (m$mean_3 - m$mean_1)),
  group_levels = c("Conventional", "Organic")
); pc_estimands_means

# fit_MDLS2_shifted@formula (not the freshly-sourced `model`) -- the real
# fit was built from model_vbc (5b.2's variance-budget-calibrated priors),
# which this script never sees; the fit's own compiled formula is the only
# copy guaranteed to match what was actually fit.
save_report("fit_summary", "MDLS2_shifted", fit_MDLS2_shifted, pc_estimands_means, fit_MDLS2_shifted@formula, model_name = "The Floor Raiser")

p_contrast_mean <- contrast_plot_panels(
  pc_estimands_means, quant = c(0.005, 0.995), scales = 'free_y',
  group_pal = Management_palette,
  legend_title = "Posteriors (population means)"); p_contrast_mean

save_gg("fit_contrast_mean", "MDLS2_shifted", p_contrast_mean, width = 10, height =6)

pc_estimands_medians <- estimand_panels(
  pairs = list(`May median` = list(md$median_1, md$median_3),
               `July median`= list(md$median_2, md$median_4)),
  extra = list("Seasonal change in mean diversity difference between management, from May to July" =
                 (md$median_4 - md$median_2) - (md$median_3 - md$median_1)),
  group_levels = c("Conventional", "Organic")
)

p_contrast_median <- contrast_plot_panels(
  pc_estimands_medians, quant = c(0.005,.995), scales = 'free_y',
  group_pal = Management_palette,
  legend_title = "Posteriors (population medians)"); p_contrast_median
save_gg("fit_contrast_median", "MDLS2_shifted", p_contrast_median, width = 10, height =6)

## Variance components -----------------------------------------------------------

# Same as plain Model 5's (5.3_MDLS2_analysis.R): b[Lo], sigma_loc,
# sigma_tr, yr[Yr] are untouched by shift (shift only enters means_MDLS2()'s
# mean/median backtransform).

stat_levels <- c("Location effect (log scale)",
                  "Location vs Tree effect SD (population)",
                  "Year effect (log scale)")

pf_re <- list(
  `Location effect (log scale)` = as_tibble(as.matrix(pf$b) * pf$sigma_loc[[1]]) %>%
    setNames(as.character(idx$Lo$to_label_n(seq_along(idx$Lo$levels)))),
  `Year effect (log scale)` = as_tibble(pf$yr) %>%
    setNames(as.character(idx$Yr$to_label_n(seq_along(idx$Yr$levels))))
)

sigma_re <- bind_rows(
  tibble(statistic = stat_levels[2], group = "sigma_loc", value = pf$sigma_loc$sigma_loc),
  tibble(statistic = stat_levels[2], group = "sigma_tr",  value = pf$sigma_tr$sigma_tr)
)

re_pal <- c(
  idx$Lo$palette_n(),
  idx$Yr$palette_n(),
  sigma_loc = Management_palette[["Population"]],
  sigma_tr  = "grey40"
)

pc_random_effects <- bind_rows(compute_contrasts(pf_re, keep = names(pf_re)), sigma_re) %>%
  mutate(statistic = factor(statistic, levels = stat_levels))

p_random_effects <- variance_component_panels(
  pc_random_effects, quant = c(0.001, 0.999), palette = re_pal); p_random_effects

save_gg("fit_location_year_effects", "MDLS2_shifted", p_random_effects, width = 8, height = 10)
