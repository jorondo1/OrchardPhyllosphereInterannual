# MODEL 6 (MDSYCV, "Bombadil the Eldest"), ITS, SHIFTED: posterior
# contrast, run against the saved fit. Mirrors 6.4_MDSYCV_16S_analysis.R --
# Year/covariate/cultivar effect panels are already covered in 6.3's own
# fit script.

hiermod_marker <- "ITS"
source('src/hiermod/0_SETUP.R')
source('src/hiermod/Models/MDSYCV_model.R') # model_MDSYCV_ITS, means_MDSYCV(), variance_partition_MDSYCV()
hiermod_out_dir <- "out/hiermod/ITS_6_cultivar_MDSYCV"

fit_MDSYCV <- readRDS(file.path(hiermod_out_dir, "fit_MDSYCV.rds"))
dat_MDSYCV <- readRDS(file.path(hiermod_out_dir, "dat_MDSYCV.rds"))

pf <- post_full(fit_MDSYCV, means_MDSYCV, shift = 1)
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

save_report("fit_summary", model_id_MDSYCV, fit_MDSYCV, model = model_MDSYCV_ITS)

p_contrast_mean <- contrast_plot_panels(
  pc_estimands_means, quant = c(0.001, 0.999), scales = 'free_y',
  group_pal = Management_palette,
  legend_title = "Posteriors (population means)",
  ratio_stats = names(Fold_change_palette), ratio_pal = Fold_change_palette); p_contrast_mean

save_gg("fit_contrast_mean", model_id_MDSYCV, p_contrast_mean)

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

save_gg("fit_contrast_median", model_id_MDSYCV, p_contrast_median)

## Comprehensive posterior summary (results report) -----------------------------
# Every interpretable posterior's own mean/median/89% PI/HPDI/pd, not just
# the headline contrast -- excludes per-tree tr[Tr] raw draws (n/a here,
# Tree isn't in this model yet) and the generic "mean"/"median" 4-column
# entries (already covered, more legibly, by pc_estimands_means/medians'
# own May/July/Conventional/Organic labels).

pc_all <- bind_rows(
  compute_contrasts(pf, keep = setdiff(names(pf), c("tr", "mean", "median")), group_levels = idx$Mg$levels),
  pc_estimands_means, pc_estimands_medians)
save_posterior_kable("results_report", model_id_MDSYCV, pc_all)

## Residual variance (sigma[Mg]) -------------------------------------------------
# Does accounting for Cultivar composition shrink the Organic-vs-
# Conventional sigma gap MDSYC found?

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
cat("Compare against MDSYC's own numbers (out/hiermod/ITS_5_covariates_MDSYC/fit_summary_MDSYC.txt) --\n")
cat("if this interval shrank meaningfully toward 0, Cultivar is explaining part of that residual asymmetry.\n")

## Variance partition -----------------------------------------------------------
# variance_partition_MDSYCV() needs raw extract.samples() (plain matrices,
# for the outer()/%*% arithmetic inside it) -- not `pf`, which post_full()
# wraps into tibbles.

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
