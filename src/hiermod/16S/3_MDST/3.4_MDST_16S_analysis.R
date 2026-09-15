# MODEL 3 (MDST, "Treebeard the Skeptic"), 16S, SHIFTED: posterior contrast,
# run against the saved fit. Variance components (sigma[Mg], sigma_tr) are
# already covered in 3.3's own fit script -- this one just mirrors 2.4's
# Management x Season contrast reporting so the two models stay directly
# comparable side by side.

hiermod_marker <- "16S"
source('src/hiermod/0_SETUP.R')
source('src/hiermod/Models/MDST_model.R') # model_MDST_16S, means_MDST()
hiermod_out_dir <- "out/hiermod/16S_3_tree_MDST"

fit_MDST <- readRDS(file.path(hiermod_out_dir, "fit_MDST.rds"))
dat_MDST <- readRDS(file.path(hiermod_out_dir, "dat_MDST.rds"))

pf <- post_full(fit_MDST, means_MDST, shift = 1)
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

save_report("fit_summary", model_id_MDST, fit_MDST, pc_estimands_means, model_MDST_16S, model_name = "Treebeard the Skeptic")

p_contrast_mean <- contrast_plot_panels(
  pc_estimands_means, quant = c(0.005, 0.995), scales = 'free_y',
  group_pal = Management_palette,
  legend_title = "Posteriors (population means)"); p_contrast_mean

save_gg("fit_contrast_mean", model_id_MDST, p_contrast_mean)

pc_estimands_medians <- estimand_panels(
  pairs = list(`May median`  = list(md$median_1, md$median_3),
               `July median` = list(md$median_2, md$median_4)),
  extra = list("Change in management contrast, from May to July" =
                 (md$median_4 - md$median_2) - (md$median_3 - md$median_1)),
  group_levels = c("Conventional", "Organic")
)

p_contrast_median <- contrast_plot_panels(
  pc_estimands_medians, quant = c(0.01, 0.995), scales = 'free_y',
  group_pal = Management_palette,
  legend_title = "Posteriors (population medians)"); p_contrast_median

save_gg("fit_contrast_median", model_id_MDST, p_contrast_median)

