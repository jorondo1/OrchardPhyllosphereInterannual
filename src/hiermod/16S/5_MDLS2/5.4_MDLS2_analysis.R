# MODEL 5 (MDLS2), 16S, SHIFTED: posterior contrast and variance components,
# run against the saved fit.

hiermod_marker <- "16S"
source('src/hiermod/0_SETUP.R')
source('src/hiermod/Models/MDLS2_model.R') # means_MDLS2()
hiermod_out_dir <- "out/hiermod/16S_5_lognormal_MDLS2_shifted"

fit_MDLS2 <- readRDS(file.path(hiermod_out_dir, "fit_MDLS2_shifted.rds"))
dat <- readRDS(file.path(hiermod_out_dir, "dat_MDLS2_shifted.rds"))
model_ppc1 <- readRDS(file.path(hiermod_out_dir, "model_ppc1.rds"))

pf <- post_full(fit_MDLS2, means_MDLS2, shift = 1)
m  <- pf$mean
md <- pf$median

## Management x Season contrasts -----------------------------------------------

pc_estimands_means <- estimand_panels(
  pairs = list(`May mean` = list(m$mean_1, m$mean_3),
               `July mean` = list(m$mean_2, m$mean_4)),
  extra = list("Change in management contrast, from May to July" =
                 (m$mean_4 - m$mean_2) - (m$mean_3 - m$mean_1)),
  group_levels = c("Conventional", "Organic")
)

save_report("fit_summary", "MDLS2_shifted", fit_MDLS2, pc_estimands_means, model_ppc1, model_name = "The Splitter")

p_contrast_mean <- contrast_plot_panels(
  pc_estimands_means, quant = c(0.005, 0.995), scales = 'free_y',
  group_pal = Management_palette,
  legend_title = "Posteriors (population means)"); p_contrast_mean

save_gg("fit_contrast_mean", "MDLS2_shifted", p_contrast_mean)

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
save_gg("fit_contrast_median", "MDLS2_shifted", p_contrast_median, width = 10, height = 6)

## Variance components -----------------------------------------------------------
# Same treatment as ITS Model 5's (see ITS's 5.4_MDLS2_analysis.R /
# 5b.4_MDLS2_shifted_analysis.R): b[Lo], sigma_loc, sigma_tr, yr[Yr] are
# untouched by shift (shift only enters means_MDLS2()'s mean/median
# backtransform).

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
  pc_random_effects, quant = c(0.005, 0.995), palette = re_pal,
  sd_stats = stat_levels[2]); p_random_effects

save_gg("fit_location_year_effects", "MDLS2_shifted", p_random_effects, width = 8, height = 10)
