# MODEL 7 (MDSTYCL, "Saruman the Fool"), 16S, SHIFTED: posterior contrast,
# run against the saved fit. Year/covariate/cultivar/Tree effect panels are
# already covered in 7.3's own fit script -- this one mirrors
# 2.4/3.4/4.4/5.4/6.4's Management x Season contrast reporting so every
# model in the family stays directly comparable.

hiermod_marker <- "16S"
source('src/hiermod/0_SETUP.R')
source('src/hiermod/Models/MDSTYCL_model.R') # model_MDSTYCL_16S, means_MDSTYCL(), variance_partition_MDSTYCL()
hiermod_out_dir <- "out/hiermod/16S_9_location_MDSTYCL"

fit_MDSTYCL <- readRDS(file.path(hiermod_out_dir, "fit_MDSTYCL.rds"))
dat_MDSTYCL <- readRDS(file.path(hiermod_out_dir, "dat_MDSTYCL.rds"))

pf <- post_full(fit_MDSTYCL, means_MDSTYCL, shift = 1)
m  <- pf$mean
md <- pf$median

# Posterior estimands
pc_estimands         <- build_pc_estimands(pf, group_levels = idx$Mg$levels)
pc_estimands_means   <- pc_estimands$means
pc_estimands_medians <- pc_estimands$medians

## Variance partition -----------------------------------------------------------
# Needs raw extract.samples() (plain matrices), not `pf`. Shapley/LMG shares
# with Management / Season / interaction effect-coded -- see the model file's
# variance_partition_*() comment and doc/R2_methods.txt.

post_raw <- extract.samples(fit_MDSTYCL)
pc_varpart <- variance_partition_MDSTYCL(post_raw, dat_MDSTYCL)

p_varpart <- pc_varpart %>%
  ggplot(aes(x = value, y = group, fill = group, height = after_stat(ndensity))) +
  ggridges::geom_density_ridges(stat = "density", alpha = 0.7, colour = "white",
                                 scale = 1.5, rel_min_height = 0.01) +
  geom_vline(xintercept = 0, colour = "grey50", linetype = "dashed") +
  scale_fill_manual(values = Variance_partition_palette) +
  labs(x = "Fraction of total variance", y = NULL, fill = NULL,
       title = "Variance partition (Bayesian R2)") +
  theme(legend.position = "none"); p_varpart

save_gg("fit_variance_partition", model_id_MDSTYCL, p_varpart, width = 8, height = 4)


## Comprehensive posterior summary  -----------------------------
pc_all <- bind_rows(
  compute_contrasts(
    pf, keep = setdiff(names(pf), c("tr", "mean", "median")), 
    group_levels = idx$Mg$levels),
  pc_estimands_means, pc_estimands_medians, pc_varpart)
save_posterior_kable("results_report", model_id_MDSTYCL, pc_all)

## Management x Season contrasts -----------------------------------------------

p_contrast_mean <- contrast_plot_panels(
  pc_estimands_means, quant = c(0.001, 0.999), scales = 'free_y',
  group_pal = Management_palette,
  legend_title = "Posteriors (population means)",
  ratio_stats = names(Fold_change_palette), ratio_pal = Fold_change_palette); p_contrast_mean

save_gg("fit_contrast_mean", model_id_MDSTYCL, p_contrast_mean)

p_contrast_median <- contrast_plot_panels(
  pc_estimands_medians, quant = c(0.001, 0.999), scales = 'free_y',
  group_pal = Management_palette,
  legend_title = "Posteriors (population medians)",
  ratio_stats = names(Fold_change_palette), ratio_pal = Fold_change_palette); p_contrast_median

save_gg("fit_contrast_median", model_id_MDSTYCL, p_contrast_median)

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

save_gg("fit_sigma_posterior", model_id_MDSTYCL, p_sigma)

sigma_contrast_MDSTYCL <- pc_sigma %>% dplyr::filter(group == "Contrast") %>% dplyr::pull(value)
cat(sprintf(
  "sigma[Mg] contrast (Organic - Conventional), MDSTYCL: median %.3f, 89%% PI [%.3f, %.3f]\n",
  median(sigma_contrast_MDSTYCL), PI(sigma_contrast_MDSTYCL)[1], PI(sigma_contrast_MDSTYCL)[2]))
cat("Compare against MDSYCV's own 0.285 [0.127, 0.457] and MDSYC's own 0.31 [0.16, 0.48] --\n")
cat("if this interval shrank further, Location is explaining part of what looked like residual noise.\n")
