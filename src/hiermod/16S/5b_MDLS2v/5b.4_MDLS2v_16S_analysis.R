# MODEL 5v (MDLS2v), 16S, SHIFTED: posterior contrast and variance
# components, run against the saved fit.

hiermod_marker <- "16S"
source('src/hiermod/0_SETUP.R')
source('src/hiermod/Models/MDLS2v_model.R') # model_MDLS2v_16S, means_MDLS2v()
hiermod_out_dir <- "out/hiermod/16S_5b_lognormal_MDLS2v_shifted"

fit_MDLS2v <- readRDS(file.path(hiermod_out_dir, "fit_MDLS2v_shifted.rds"))
dat <- readRDS(file.path(hiermod_out_dir, "dat_MDLS2v_shifted.rds"))

pf <- post_full(fit_MDLS2v, means_MDLS2v, shift = 1)
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

save_report("fit_summary", model_id, fit_MDLS2v, pc_estimands_means, model_MDLS2v_16S, model_name = "The Structured Splitter")

p_contrast_mean <- contrast_plot_panels(
  pc_estimands_means, quant = c(0.005, 0.995), scales = 'free_y',
  group_pal = Management_palette,
  legend_title = "Posteriors (population means)"); p_contrast_mean

save_gg("fit_contrast_mean", model_id, p_contrast_mean)

pc_estimands_medians <- estimand_panels(
  pairs = list(`May median` = list(md$median_1, md$median_3),
               `July median`= list(md$median_2, md$median_4)),
  extra = list("Change in management contrast, from May to July" =
                 (md$median_4 - md$median_2) - (md$median_3 - md$median_1)),
  group_levels = c("Conventional", "Organic")
)

p_contrast_median <- contrast_plot_panels(
  pc_estimands_medians, quant = c(0.01,.995), scales = 'free_y',
  group_pal = Management_palette,
  legend_title = "Posteriors (population medians)"); p_contrast_median
save_gg("fit_contrast_median", model_id, p_contrast_median, width = 10, height = 6)

## Variance structure -----------------------------------------------------------
# The thing this model reports that Model 5 can't: the 4 cell sigmas
# decomposed into an interpretable baseline + Management effect + Season
# effect + interaction, all on the log scale.

pf_ls <- list(
  ls0     = pf$ls0,
  ls_Mg   = pf$ls_Mg,
  ls_Mo   = pf$ls_Mo,
  ls_MgMo = pf$ls_MgMo
)
pc_ls <- bind_rows(lapply(names(pf_ls), function(nm){
  tibble(statistic = "Variance structure (log scale)", group = nm, value = pf_ls[[nm]][[1]])
}))

p_variance_structure <- pc_ls %>%
  ggplot(aes(x = value, y = group, fill = group)) +
  ggridges::geom_density_ridges(alpha = 0.7, colour = "white") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  guides(fill = "none") +
  labs(x = "Log-sigma effect", y = NULL,
       title = "Model 5v: residual SD structure (baseline + Mg + Mo + Mg:Mo)")
save_gg("fit_variance_structure", model_id, p_variance_structure)

## Random effects -----------------------------------------------------------

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

save_gg("fit_location_year_effects", model_id, p_random_effects, width = 8, height = 10)
