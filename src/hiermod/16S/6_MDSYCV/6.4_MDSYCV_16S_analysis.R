# MODEL 6 (MDSYCV), 16S: posterior contrasts from the saved fit (effect panels in 6.3)

hiermod_marker <- "16S"
source('src/hiermod/0_SETUP.R')
source('src/hiermod/Models/MDSYCV_model.R') # model_MDSYCV_16S, means_MDSYCV()
hiermod_out_dir <- "out/hiermod/16S_6_cultivar_MDSYCV"

fit_MDSYCV <- readRDS(file.path(hiermod_out_dir, "fit_MDSYCV.rds"))
dat_MDSYCV <- readRDS(file.path(hiermod_out_dir, "dat_MDSYCV.rds"))

pf <- post_full(fit_MDSYCV, means_MDSYCV, shift = 1)

## Management x Season contrasts -----------------------------------------------

pc_estimands <- build_pc_estimands(pf, group_levels = idx$Mg$levels)
pc_estimands_means   <- pc_estimands$means
pc_estimands_medians <- pc_estimands$medians

save_report("fit_summary", model_id_MDSYCV, fit_MDSYCV, model = model_MDSYCV_16S)

p_contrast_mean <- contrast_plot_panels(
  pc_estimands_means, quant = c(0.001, 0.999), scales = 'free_y',
  group_pal = Management_palette,
  legend_title = "Posteriors (population means)",
  ratio_stats = names(Fold_change_palette), ratio_pal = Fold_change_palette); p_contrast_mean

save_gg("fit_contrast_mean", model_id_MDSYCV, p_contrast_mean)


p_contrast_median <- contrast_plot_panels(
  pc_estimands_medians, quant = c(0.001, 0.999), scales = 'free_y',
  group_pal = Management_palette,
  legend_title = "Posteriors (population medians)",
  ratio_stats = names(Fold_change_palette), ratio_pal = Fold_change_palette); p_contrast_median

save_gg("fit_contrast_median", model_id_MDSYCV, p_contrast_median)

## Comprehensive posterior summary (results report) -----------------------------
# All interpretable posteriors (excl. raw per-tree draws and generic mean/median columns)

pc_all <- bind_rows(
  compute_contrasts(pf, keep = setdiff(names(pf), c("tr", "mean", "median")), group_levels = idx$Mg$levels),
  pc_estimands_means, pc_estimands_medians)
save_posterior_kable("results_report", model_id_MDSYCV, pc_all)

## Residual variance (sigma[Mg]) -------------------------------------------------
# Does cultivar shrink the organic vs conventional sigma gap seen in MDSYC?

pc_sigma <- compute_contrasts(pf, keep = "sigma", group_levels = idx$Mg$levels)

p_sigma <- contrast_plot_panels(
  pc_sigma, quant = c(0, 1), group_pal = Management_palette,
  legend_title = "Posteriors (residual SD, log scale)") +
  labs(x = "sigma[Mg]"); p_sigma

save_gg("fit_sigma_posterior", model_id_MDSYCV, p_sigma)

sigma_contrast_MDSYCV <- pc_sigma %>% dplyr::filter(group == "Contrast") %>% dplyr::pull(value)
cat(sprintf(
  "sigma[Mg] contrast (Organic - Conventional), MDSYCV: median %.3f, 89%% PI [%.3f, %.3f]\n",
  median(sigma_contrast_MDSYCV), PI(sigma_contrast_MDSYCV)[1], PI(sigma_contrast_MDSYCV)[2]))
cat("Compare against MDSYC's own [0.16, 0.48] (out/hiermod/16S_5_covariates_MDSYC/fit_summary_MDSYC.txt) --\n")
cat("if this interval shrank meaningfully toward 0, Cultivar is explaining part of that residual asymmetry.\n")

## Variance partition -----------------------------------------------------------
# Needs raw extract.samples() matrices, not post_full() tibbles

post_raw <- extract.samples(fit_MDSYCV)
pc_varpart <- variance_partition_MDSYCV(post_raw, dat_MDSYCV)

p_varpart <- pc_varpart %>%
  ggplot(aes(x = value, fill = group, colour = group)) +
  geom_density(alpha = 0.5, linewidth = 0.2) +
  scale_fill_manual(values = Variance_partition_palette) +
  scale_colour_manual(values = Variance_partition_palette) +
  labs(x = "Fraction of total variance", y = NULL, fill = NULL, colour = NULL,
       title = "Variance partition (Bayesian R2)"); p_varpart

save_gg("fit_variance_partition", model_id_MDSYCV, p_varpart, width = 8, height = 4)
