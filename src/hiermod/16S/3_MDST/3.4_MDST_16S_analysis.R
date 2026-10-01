# MODEL 3 (MDST), 16S: posterior contrasts from the saved fit
# - same reporting as 2.4 (side-by-side comparison); variance components in 3.3

hiermod_marker <- "16S"
source('src/hiermod/0_SETUP.R')
source('src/hiermod/Models/MDST_model.R') # model_MDST_16S, means_MDST()
hiermod_out_dir <- "out/hiermod/16S_3_tree_MDST"

fit_MDST <- readRDS(file.path(hiermod_out_dir, "fit_MDST.rds"))
dat_MDST <- readRDS(file.path(hiermod_out_dir, "dat_MDST.rds"))

pf <- post_full(fit_MDST, means_MDST, shift = 1)

## Management x Season contrasts -----------------------------------------------

pc_estimands <- build_pc_estimands(pf, group_levels = idx$Mg$levels)
pc_estimands_means   <- pc_estimands$means
pc_estimands_medians <- pc_estimands$medians

save_report("fit_summary", model_id_MDST, fit_MDST, model = model_MDST_16S)

p_contrast_mean <- contrast_plot_panels(
  pc_estimands_means, quant = c(0.005, 0.995), scales = 'free_y',
  group_pal = Management_palette,
  legend_title = "Posteriors (population means)",
  ratio_stats = names(Fold_change_palette), ratio_pal = Fold_change_palette); p_contrast_mean

save_gg("fit_contrast_mean", model_id_MDST, p_contrast_mean)


p_contrast_median <- contrast_plot_panels(
  pc_estimands_medians, quant = c(0.01, 0.995), scales = 'free_y',
  group_pal = Management_palette,
  legend_title = "Posteriors (population medians)",
  ratio_stats = names(Fold_change_palette), ratio_pal = Fold_change_palette); p_contrast_median

save_gg("fit_contrast_median", model_id_MDST, p_contrast_median)

## Comprehensive posterior summary (results report) -----------------------------
# All interpretable posteriors (excl. raw per-tree draws and generic mean/median columns)

pc_all <- bind_rows(
  compute_contrasts(pf, keep = setdiff(names(pf), c("tr", "mean", "median")), group_levels = idx$Mg$levels),
  pc_estimands_means, pc_estimands_medians)
save_posterior_kable("results_report", model_id_MDST, pc_all)

