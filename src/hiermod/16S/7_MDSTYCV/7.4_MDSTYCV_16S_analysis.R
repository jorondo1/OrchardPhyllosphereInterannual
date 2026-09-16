# MODEL 7 (MDSTYCV, "Saruman the Fool"), 16S, SHIFTED: posterior contrast,
# run against the saved fit. Year/covariate/cultivar/Tree effect panels are
# already covered in 7.3's own fit script -- this one mirrors
# 2.4/3.4/4.4/5.4/6.4's Management x Season contrast reporting so every
# model in the family stays directly comparable.

hiermod_marker <- "16S"
source('src/hiermod/0_SETUP.R')
source('src/hiermod/Models/MDSTYCV_model.R') # model_MDSTYCV_16S, means_MDSTYCV(), variance_partition_MDSTYCV()
hiermod_out_dir <- "out/hiermod/16S_7_tree_full_MDSTYCV"

fit_MDSTYCV <- readRDS(file.path(hiermod_out_dir, "fit_MDSTYCV.rds"))
dat_MDSTYCV <- readRDS(file.path(hiermod_out_dir, "dat_MDSTYCV.rds"))

pf <- post_full(fit_MDSTYCV, means_MDSTYCV, shift = 1)
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

save_report("fit_summary", model_id_MDSTYCV, fit_MDSTYCV, model = model_MDSTYCV_16S)

p_contrast_mean <- contrast_plot_panels(
  pc_estimands_means, quant = c(0.001, 0.999), scales = 'free_y',
  group_pal = Management_palette,
  legend_title = "Posteriors (population means)",
  ratio_stats = names(Fold_change_palette), ratio_pal = Fold_change_palette); p_contrast_mean

save_gg("fit_contrast_mean", model_id_MDSTYCV, p_contrast_mean)

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

save_gg("fit_contrast_median", model_id_MDSTYCV, p_contrast_median)

## Comprehensive posterior summary (results report) -----------------------------
# Every interpretable posterior's own mean/median/89% PI/HPDI, not just the
# headline contrast -- excludes per-tree tr[Tr] raw draws (too many, not
# individually interpretable) and the generic "mean"/"median" 4-column
# entries (already covered, more legibly, by pc_estimands_means/medians'
# own May/July/Conventional/Organic labels).

pc_all <- bind_rows(
  compute_contrasts(pf, keep = setdiff(names(pf), c("tr", "mean", "median")), group_levels = idx$Mg$levels),
  pc_estimands_means, pc_estimands_medians)
save_posterior_kable("results_report", model_id_MDSTYCV, pc_all)

## Residual variance (sigma[Mg]) -------------------------------------------------
# Does the still-open Organic-vs-Conventional sigma asymmetry (unexplained
# by Cultivar in Model 6: 0.285 vs MDSYC's own 0.31) change once Tree is
# estimated separately? Possible either direction -- Tree could absorb some
# of what looked like residual noise, or leave it untouched.

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
cat("Compare against MDSYCV's own 0.285 [0.127, 0.457] and MDSYC's own 0.31 [0.16, 0.48] --\n")
cat("if this interval shrank further, Tree is explaining part of what looked like residual noise.\n")

## Variance partition -----------------------------------------------------------
# variance_partition_MDSTYCV() needs raw extract.samples() (plain matrices,
# for the outer()/%*% arithmetic inside it) -- not `pf`, which post_full()
# wraps into tibbles. Six groups now: Management x Season/Year/Covariates/
# Cultivar (sequential fixed-effect decomposition) + Tree + Residual.

post_raw <- extract.samples(fit_MDSTYCV)
pc_varpart <- variance_partition_MDSTYCV(post_raw, dat_MDSTYCV)

p_varpart <- pc_varpart %>%
  ggplot(aes(x = value, fill = group, colour = group)) +
  geom_density(alpha = 0.5, linewidth = 0.2) +
  scale_fill_manual(values = Variance_partition_palette) +
  scale_colour_manual(values = Variance_partition_palette) +
  labs(x = "Fraction of total variance", y = NULL, fill = NULL, colour = NULL,
       title = "Variance partition (Bayesian R2)"); p_varpart

save_gg("fit_variance_partition", model_id_MDSTYCV, p_varpart, width = 8, height = 4)
