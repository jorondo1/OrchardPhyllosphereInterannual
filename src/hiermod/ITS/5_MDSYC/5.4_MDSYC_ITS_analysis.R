# MODEL 5 (MDSYC, "Radagast the Grower"), ITS, SHIFTED: posterior contrast,
# run against the saved fit. Mirrors 5.4_MDSYC_16S_analysis.R -- Year/
# covariate effect panels are already covered in 5.3's own fit script.

hiermod_marker <- "ITS"
source('src/hiermod/0_SETUP.R')
source('src/hiermod/Models/MDSYC_model.R') # model_MDSYC_ITS, means_MDSYC(), variance_partition_MDSYC()
hiermod_out_dir <- "out/hiermod/ITS_5_covariates_MDSYC"

fit_MDSYC <- readRDS(file.path(hiermod_out_dir, "fit_MDSYC.rds"))
dat_MDSYC <- readRDS(file.path(hiermod_out_dir, "dat_MDSYC.rds"))

pf <- post_full(fit_MDSYC, means_MDSYC, shift = 1)
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

save_report("fit_summary", model_id_MDSYC, fit_MDSYC, model = model_MDSYC_ITS)

p_contrast_mean <- contrast_plot_panels(
  pc_estimands_means, quant = c(0.001, 0.999), scales = 'free_y',
  group_pal = Management_palette,
  legend_title = "Posteriors (population means)",
  ratio_stats = names(Fold_change_palette), ratio_pal = Fold_change_palette); p_contrast_mean

save_gg("fit_contrast_mean", model_id_MDSYC, p_contrast_mean)

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

save_gg("fit_contrast_median", model_id_MDSYC, p_contrast_median)

## Comprehensive posterior summary (results report) -----------------------------
# Every interpretable posterior's own mean/median/89% PI/HPDI, not just the
# headline contrast -- excludes per-tree tr[Tr] raw draws (too many, not
# individually interpretable) and the generic "mean"/"median" 4-column
# entries (already covered, more legibly, by pc_estimands_means/medians'
# own May/July/Conventional/Organic labels).

pc_all <- bind_rows(
  compute_contrasts(pf, keep = setdiff(names(pf), c("tr", "mean", "median")), group_levels = idx$Mg$levels),
  pc_estimands_means, pc_estimands_medians)
save_posterior_kable("results_report", model_id_MDSYC, pc_all)

## Residual variance (sigma[Mg]) -------------------------------------------------

pc_sigma <- compute_contrasts(pf, keep = "sigma", group_levels = idx$Mg$levels)

p_sigma <- contrast_plot_panels(
  pc_sigma, quant = c(0, 1), group_pal = Management_palette,
  legend_title = "Posteriors (residual SD, log scale)") +
  labs(x = "sigma[Mg]"); p_sigma

save_gg("fit_sigma_posterior", model_id_MDSYC, p_sigma)

## Variance partition -----------------------------------------------------------
# variance_partition_MDSYC() needs raw extract.samples() (plain matrices,
# for the outer()/%*% arithmetic inside it) -- not `pf`, which post_full()
# wraps into tibbles.

post_raw <- extract.samples(fit_MDSYC)
pc_varpart <- variance_partition_MDSYC(post_raw, dat_MDSYC)

p_varpart <- pc_varpart %>%
  ggplot(aes(x = value, fill = group, colour = group)) +
  geom_density(alpha = 0.5, linewidth = 0.2) +
  scale_fill_manual(values = Variance_partition_palette) +
  scale_colour_manual(values = Variance_partition_palette) +
  labs(x = "Fraction of total variance", y = NULL, fill = NULL, colour = NULL,
       title = "Variance partition (Bayesian R2)"); p_varpart

save_gg("fit_variance_partition", model_id_MDSYC, p_varpart, width = 8, height = 4)
