# MODEL 5 (MDSYC, "Radagast the Grower"), 16S, SHIFTED: posterior contrast,
# run against the saved fit. Year/covariate effect panels are already
# covered in 5.3's own fit script -- this one mirrors 2.4/3.4/4.4's
# Management x Season contrast reporting so every model in the family
# stays directly comparable. Reported at the default (z=0) reference level
# for all three covariates, i.e. this sample's own average weather/seq
# depth (see means_MDSYC()'s own covariate_offset argument if a different
# reference level is ever wanted).

hiermod_marker <- "16S"
source('src/hiermod/0_SETUP.R')
source('src/hiermod/Models/MDSYC_model.R') # model_MDSYC_16S, means_MDSYC()
hiermod_out_dir <- "out/hiermod/16S_5_covariates_MDSYC"

fit_MDSYC <- readRDS(file.path(hiermod_out_dir, "fit_MDSYC.rds"))
dat_MDSYC <- readRDS(file.path(hiermod_out_dir, "dat_MDSYC.rds"))

pf <- post_full(fit_MDSYC, means_MDSYC, shift = 1)
m  <- pf$mean
md <- pf$median

## Management x Season contrasts -----------------------------------------------

pc_estimands_means <- estimand_panels(
  pairs = list(`May mean`  = list(m$mean_1, m$mean_3),
               `July mean` = list(m$mean_2, m$mean_4)),
  extra = list(
    "May fold change (Conventional / Organic)"  = m$mean_1 / m$mean_3,
    "July fold change (Conventional / Organic)" = m$mean_2 / m$mean_4,
    "Contrast between folds (May / July)"       = (m$mean_1 / m$mean_3) / (m$mean_2 / m$mean_4)
  ),
  group_levels = c("Conventional", "Organic")
)

save_report("fit_summary", model_id_MDSYC, fit_MDSYC, model = model_MDSYC_16S)

p_contrast_mean <- contrast_plot_panels(
  pc_estimands_means, quant = c(0.001, 0.999), scales = 'free_y',
  group_pal = Management_palette,
  legend_title = "Posteriors (population means)",
  ratio_stats = names(Fold_change_palette), ratio_pal = Fold_change_palette); p_contrast_mean

save_gg("fit_contrast_mean", model_id_MDSYC, p_contrast_mean)

pc_estimands_medians <- estimand_panels(
  pairs = list(`May median`  = list(md$median_1, md$median_3),
               `July median` = list(md$median_2, md$median_4)),
  extra = list(
    "May fold change (Conventional / Organic)"  = md$median_1 / md$median_3,
    "July fold change (Conventional / Organic)" = md$median_2 / md$median_4,
    "Contrast between folds (May / July)"       = (md$median_1 / md$median_3) / (md$median_2 / md$median_4)
  ),
  group_levels = c("Conventional", "Organic")
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
