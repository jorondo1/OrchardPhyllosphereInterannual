# MODEL 6 (MDSYCV, "Bombadil the Eldest"), 16S, SHIFTED: posterior contrast,
# run against the saved fit. Year/covariate/cultivar effect panels are
# already covered in 6.3's own fit script -- this one mirrors
# 2.4/3.4/4.4/5.4's Management x Season contrast reporting so every model
# in the family stays directly comparable.

hiermod_marker <- "16S"
source('src/hiermod/0_SETUP.R')
source('src/hiermod/Models/MDSYCV_model.R') # model_MDSYCV_16S, means_MDSYCV()
hiermod_out_dir <- "out/hiermod/16S_6_cultivar_MDSYCV"

fit_MDSYCV <- readRDS(file.path(hiermod_out_dir, "fit_MDSYCV.rds"))
dat_MDSYCV <- readRDS(file.path(hiermod_out_dir, "dat_MDSYCV.rds"))

pf <- post_full(fit_MDSYCV, means_MDSYCV, shift = 1)
m  <- pf$mean
md <- pf$median

## Management x Season contrasts -----------------------------------------------

pc_estimands_means <- estimand_panels(
  pairs = list(`May mean`  = list(m$mean_1, m$mean_3),
               `July mean` = list(m$mean_2, m$mean_4)),
  extra = list("Change in management contrast, from May to July" =
                 (m$mean_4 - m$mean_2) - (m$mean_3 - m$mean_1)),
  group_levels = c("Conventional", "Organic")
)

save_report("fit_summary", model_id_MDSYCV, fit_MDSYCV, pc_estimands_means, model_MDSYCV_16S, model_name = "Bombadil the Eldest")

p_contrast_mean <- contrast_plot_panels(
  pc_estimands_means, quant = c(0.001, 0.999), scales = 'free_y',
  group_pal = Management_palette,
  legend_title = "Posteriors (population means)"); p_contrast_mean

save_gg("fit_contrast_mean", model_id_MDSYCV, p_contrast_mean)

pc_estimands_medians <- estimand_panels(
  pairs = list(`May median`  = list(md$median_1, md$median_3),
               `July median` = list(md$median_2, md$median_4)),
  extra = list("Change in management contrast, from May to July" =
                 (md$median_4 - md$median_2) - (md$median_3 - md$median_1)),
  group_levels = c("Conventional", "Organic")
)

p_contrast_median <- contrast_plot_panels(
  pc_estimands_medians, quant = c(0.001, 0.999), scales = 'free_y',
  group_pal = Management_palette,
  legend_title = "Posteriors (population medians)"); p_contrast_median

save_gg("fit_contrast_median", model_id_MDSYCV, p_contrast_median)

## Residual variance (sigma[Mg]) -------------------------------------------------
# The actual question Cultivar was added to help answer: does accounting
# for cultivar composition shrink the Organic-vs-Conventional sigma gap
# MDSYC found (89% contrast [0.16, 0.48], excluding 0)?

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
# variance_partition_MDSYCV() needs raw extract.samples() (plain matrices,
# for the outer()/%*% arithmetic inside it) -- not `pf`, which post_full()
# wraps into tibbles.

post_raw <- extract.samples(fit_MDSYCV)
pc_varpart <- variance_partition_MDSYCV(post_raw, dat_MDSYCV)

p_varpart <- pc_varpart %>%
  ggplot(aes(x = value, fill = group, colour = group)) +
  geom_density(alpha = 0.5, linewidth = 0.2) +
  scale_fill_manual(values = c("Explained (fixed effects)" = "#4C72B0", "Residual" = "grey50")) +
  scale_colour_manual(values = c("Explained (fixed effects)" = "#4C72B0", "Residual" = "grey50")) +
  labs(x = "Fraction of total variance", y = NULL, fill = NULL, colour = NULL,
       title = "Variance partition (Bayesian R2)"); p_varpart

save_gg("fit_variance_partition", model_id_MDSYCV, p_varpart, width = 8, height = 4)
