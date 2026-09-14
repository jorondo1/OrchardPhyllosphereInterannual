# MODEL 4 (MDLS): posterior contrast, run against the saved fit.

hiermod_marker <- "ITS"
source('src/hiermod/0_SETUP.R')
source('src/hiermod/Models/MDLS_model.R') # model_MDLS_ITS, means_MDLS()
model <- model_MDLS_ITS
hiermod_out_dir <- "out/hiermod/ITS_4_lognormal_MDLS"

fit_MDLS <- readRDS(file.path(hiermod_out_dir, "fit_MDLS.rds"))

## Posterior contrast ---------------------------------------------------------

# because of the means_MDLS output structure:
pf <- post_full(fit_MDLS, means_MDLS)
m  <- pf$mean
md <- pf$median

# estimand_panels() (postcontrast_helpers.R): one call builds the
# May/July raw+Contrast panels plus the Contrast-only "change in gap"
# panel, with statistic factor levels set from pairs/extra order --
# supersedes the old local gap_tibble()+estimand_rows()+bind_rows() trio.

pc_estimands_means <- estimand_panels(
  pairs = list(`May mean` = list(m$mean_1, m$mean_3),
               `July mean` = list(m$mean_2, m$mean_4)),
  extra = list("Change in management contrast, from May to July" =
                 (m$mean_4 - m$mean_2) - (m$mean_3 - m$mean_1)),
  group_levels = c("Conventional", "Organic")
)

save_report("fit_summary", model_id, fit_MDLS, pc_estimands_means, model, model_name = "The Season Ticket")

p_contrast_mean <- contrast_plot_panels(
  pc_estimands_means, quant = c(0.005, 0.995), scales = 'free_y',
  group_pal = Management_palette,
  legend_title = "Posteriors (population means)"); p_contrast_mean

save_gg("fit_contrast_mean", model_id, p_contrast_mean)

pc_estimands_medians <- estimand_panels(
  pairs = list(`May median` = list(md$median_1, md$median_3),
               `July median`= list(md$median_2, md$median_4)),
  extra = list("Change in management contrast, from May to July" =
                 (md$median_4 - md$median_2) - (md$median_3 - md$median_1)),
  group_levels = c("Conventional", "Organic")
)

p_contrast_median <- contrast_plot_panels(
  pc_estimands_medians, quant = c(0.005,.995), scales = 'free_y',
  group_pal = Management_palette,
  legend_title = "Posteriors (population medians)"); p_contrast_median
save_gg("fit_contrast_median", model_id, p_contrast_median)
