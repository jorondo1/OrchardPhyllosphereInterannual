# MODEL 2 (MDS), ITS, SHIFTED: posterior contrast, run against the saved
# fit. Mirrors 2.4_MDS_16S_analysis.R. No random-effect variance-component
# section (MDS has none yet) -- but sigma[Mg] itself is still a modeled
# variance component worth its own panel.

hiermod_marker <- "ITS"
source('src/hiermod/0_SETUP.R')
source('src/hiermod/Models/MDS_model.R') # model_MDS_ITS, means_MDS()
hiermod_out_dir <- "out/hiermod/ITS_2_interaction_MDS"

fit_MDS <- readRDS(file.path(hiermod_out_dir, "fit_MDS.rds"))
dat_MDS <- readRDS(file.path(hiermod_out_dir, "dat_MDS.rds"))

pf <- post_full(fit_MDS, means_MDS, shift = 1)
m  <- pf$mean
md <- pf$median

## Management x Season contrasts -----------------------------------------------

pc_estimands_means <- estimand_panels(
  pairs = list(`May mean`  = list(m$mean_1, m$mean_3),
               `July mean` = list(m$mean_2, m$mean_4)),
  extra = setNames(
    list(m$mean_1 / m$mean_3, m$mean_2 / m$mean_4),
    names(Fold_change_palette)
  ),
  group_levels = idx$Mg$levels
)

save_report("fit_summary", model_id_MDS, fit_MDS, model = model_MDS_ITS)

p_contrast_mean <- contrast_plot_panels(
  pc_estimands_means, quant = c(0.001, 0.999), scales = 'free_y',
  group_pal = Management_palette,
  legend_title = "Posteriors (population means)",
  ratio_stats = names(Fold_change_palette), ratio_pal = Fold_change_palette); p_contrast_mean

save_gg("fit_contrast_mean", model_id_MDS, p_contrast_mean)

pc_estimands_medians <- estimand_panels(
  pairs = list(`May median`  = list(md$median_1, md$median_3),
               `July median` = list(md$median_2, md$median_4)),
  extra = setNames(
    list(md$median_1 / md$median_3, md$median_2 / md$median_4),
    names(Fold_change_palette)
  ),
  group_levels = idx$Mg$levels
)

p_contrast_median <- contrast_plot_panels(
  pc_estimands_medians, quant = c(0.001, 0.999), scales = 'free_y',
  group_pal = Management_palette,
  legend_title = "Posteriors (population medians)",
  ratio_stats = names(Fold_change_palette), ratio_pal = Fold_change_palette); p_contrast_median

save_gg("fit_contrast_median", model_id_MDS, p_contrast_median)

## Comprehensive posterior summary (results report) -----------------------------

pc_all <- bind_rows(
  compute_contrasts(pf, keep = setdiff(names(pf), c("tr", "mean", "median")), group_levels = idx$Mg$levels),
  pc_estimands_means, pc_estimands_medians)
save_posterior_kable("results_report", model_id_MDS, pc_all)

## Residual variance (sigma[Mg]) -------------------------------------------------

pc_sigma <- compute_contrasts(pf, keep = "sigma", group_levels = idx$Mg$levels)

p_sigma <- contrast_plot_panels(
  pc_sigma, quant = c(0.005, 0.995), group_pal = Management_palette,
  legend_title = "Posteriors (residual SD, log scale)") +
  labs(x = "sigma[Mg]"); p_sigma

save_gg("fit_sigma_posterior", model_id_MDS, p_sigma)
