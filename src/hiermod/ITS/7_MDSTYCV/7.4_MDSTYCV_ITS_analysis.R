# MODEL 7 (MDSTYCV, "Saruman the Fool"), ITS, SHIFTED: posterior contrast,
# run against the saved fit. Mirrors 7.4_MDSTYCV_16S_analysis.R -- Year/
# covariate/cultivar/Tree effect panels are already covered in 7.3's own
# fit script.

hiermod_marker <- "ITS"
source('src/hiermod/0_SETUP.R')
source('src/hiermod/Models/MDSTYCV_model.R') # model_MDSTYCV_ITS, means_MDSTYCV(), variance_partition_MDSTYCV()
hiermod_out_dir <- "out/hiermod/ITS_7_tree_full_MDSTYCV"

fit_MDSTYCV <- readRDS(file.path(hiermod_out_dir, "fit_MDSTYCV.rds"))
dat_MDSTYCV <- readRDS(file.path(hiermod_out_dir, "dat_MDSTYCV.rds"))

set.seed(123); pf <- post_full(fit_MDSTYCV, means_MDSTYCV, shift = 1)
m  <- pf$mean
md <- pf$median

## Posterior estimands

pc_estimands         <- build_pc_estimands(pf, group_levels = idx$Mg$levels)
saveRDS(pc_estimands, file.path(hiermod_out_dir, "pc_MDSTYCV.rds"))

pc_estimands_means   <- pc_estimands$means
pc_estimands_medians <- pc_estimands$medians


## Variance partition -----------------------------------------------------------
# variance_partition_MDSTYCV() needs raw extract.samples(). Marginal
# ("by margin") decomposition only -- see variance_partition_MDSTYCV()'s own
# comment (MDSTYCV_model.R) for why Management x Season stays one combined
# term rather than splitting into main effects + interaction.

post_raw <- extract.samples(fit_MDSTYCV)
pc_varpart <- variance_partition_MDSTYCV(post_raw, dat_MDSTYCV)

p_varpart <- pc_varpart %>%
  ggplot(aes(x = value, y = group, fill = group, height = after_stat(ndensity))) +
  ggridges::geom_density_ridges(stat = "density", alpha = 0.7, colour = "white",
                                 scale = 1.5, rel_min_height = 0.01) +
  geom_vline(xintercept = 0, colour = "grey50", linetype = "dashed") +
  scale_fill_manual(values = Variance_partition_palette) +
  labs(x = "Fraction of total variance", y = NULL, fill = NULL,
       title = "Variance partition (Bayesian R2, by margin)") +
  theme(legend.position = "none"); p_varpart

save_gg("fit_variance_partition", model_id_MDSTYCV, p_varpart, width = 8, height = 4)

## Comprehensive posterior summary (results report) -----------------------------
# Every interpretable posterior's own mean/median/89% PI/HPDI/pd, not just
# the headline contrast -- excludes per-tree tr[Tr] raw draws (too many,
# not individually interpretable) and the generic "mean"/"median" 4-column
# entries (already covered, more legibly, by pc_estimands_means/medians'
# own May/July/Conventional/Organic labels).

pc_all <- bind_rows(
  compute_contrasts(pf, keep = setdiff(names(pf), c("tr", "mean", "median")), group_levels = idx$Mg$levels),
  pc_estimands_means, pc_estimands_medians, pc_varpart)
save_posterior_kable("ITS_results_report", model_id_MDSTYCV, pc_all)

## Management x Season contrasts -----------------------------------------------

save_report("fit_summary", model_id_MDSTYCV, fit_MDSTYCV, model = model_MDSTYCV_ITS)

p_contrast_mean <- contrast_plot_panels(
  pc_estimands_means, quant = c(0.001, 0.999), scales = 'free_y',
  group_pal = Management_palette,
  legend_title = "Posteriors (population means)",
  ratio_stats = names(Fold_change_palette), ratio_pal = Fold_change_palette); p_contrast_mean

save_gg("fit_contrast_mean", model_id_MDSTYCV, p_contrast_mean)

p_contrast_median <- contrast_plot_panels(
  pc_estimands_medians, quant = c(0.001, 0.999), scales = 'free_y',
  group_pal = Management_palette,
  legend_title = "Posteriors (population medians)",
  ratio_stats = names(Fold_change_palette), ratio_pal = Fold_change_palette); p_contrast_median

save_gg("fit_contrast_median", model_id_MDSTYCV, p_contrast_median)

## Residual variance (sigma[Mg]) -------------------------------------------------

pc_sigma <- compute_contrasts(pf, keep = "sigma", group_levels = idx$Mg$levels)

p_sigma <- contrast_plot_panels(
  pc_sigma, quant = c(0, 1), group_pal = Management_palette,
  legend_title = "Posteriors (residual SD, log scale)") +
  labs(x = "sigma[Mg]"); p_sigma

save_gg("fit_sigma_posterior", model_id_MDSTYCV, p_sigma)

sigma_contrast_MDSTYCV <- pc_sigma %>% dplyr::filter(group == "Contrast") %>% dplyr::pull(value)
cat(sprintf(
  "sigma[Mg] contrast (Organic - Conventional), MDSTYCV: median %.3f, 89%% PI [%.3f, %.3f]\n",
  median(sigma_contrast_MDSTYCV), PI(sigma_contrast_MDSTYCV)[1], PI(sigma_contrast_MDSTYCV)[2]))
