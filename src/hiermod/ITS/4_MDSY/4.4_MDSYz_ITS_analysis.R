# MODEL 4 (MDSYz), ITS: posterior contrasts from the saved fit (Year panel in 4.3)

hiermod_marker <- "ITS"
source('src/hiermod/0_SETUP.R')
source('src/hiermod/Models/MDSYz_model.R') # model_MDSYz_ITS, means_MDSYz()
hiermod_out_dir <- "out/hiermod/ITS_4_year_MDSY"

fit_MDSYz <- readRDS(file.path(hiermod_out_dir, "fit_MDSYz.rds"))
dat_MDSYz <- readRDS(file.path(hiermod_out_dir, "dat_MDSYz.rds"))

pf <- post_full(fit_MDSYz, means_MDSYz, shift = 1)

## Management x Season contrasts -----------------------------------------------

pc_estimands <- build_pc_estimands(pf, group_levels = idx$Mg$levels)
pc_estimands_means   <- pc_estimands$means
pc_estimands_medians <- pc_estimands$medians

save_report("fit_summary", model_id_MDSYz, fit_MDSYz, model = model_MDSYz_ITS)

p_contrast_mean <- contrast_plot_panels(
  pc_estimands_means, quant = c(0.001, 0.999), scales = 'free_y',
  group_pal = Management_palette,
  legend_title = "Posteriors (population means)",
  ratio_stats = names(Fold_change_palette), ratio_pal = Fold_change_palette); p_contrast_mean

save_gg("fit_contrast_mean", model_id_MDSYz, p_contrast_mean)


p_contrast_median <- contrast_plot_panels(
  pc_estimands_medians, quant =  c(0.001, 0.999), scales = 'free_y',
  group_pal = Management_palette,
  legend_title = "Posteriors (population medians)",
  ratio_stats = names(Fold_change_palette), ratio_pal = Fold_change_palette); p_contrast_median

save_gg("fit_contrast_median", model_id_MDSYz, p_contrast_median)

## Comprehensive posterior summary (results report) -----------------------------
# All interpretable posteriors (excl. raw per-tree draws and generic mean/median columns)

pc_all <- bind_rows(
  compute_contrasts(pf, keep = setdiff(names(pf), c("tr", "mean", "median")), group_levels = idx$Mg$levels),
  pc_estimands_means, pc_estimands_medians)
save_posterior_kable("results_report", model_id_MDSYz, pc_all)

## Residual variance (sigma[Mg]) -------------------------------------------------

pc_sigma <- compute_contrasts(pf, keep = "sigma", group_levels = idx$Mg$levels)

p_sigma <- contrast_plot_panels(
  pc_sigma, quant = c(0.005, 0.995), group_pal = Management_palette,
  legend_title = "Posteriors (residual SD, log scale)") +
  labs(x = "sigma[Mg]"); p_sigma

save_gg("fit_sigma_posterior", model_id_MDSYz, p_sigma)
