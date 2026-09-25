# MODEL 4 (MDSYz, "Elrond the Ageless"), 16S, SHIFTED: posterior contrast,
# run against the saved fit. Year's own effect panel is already covered in
# 4.3's own fit script -- this one mirrors 2.4/3.4's Management x Season
# contrast reporting so every model in the family stays directly comparable.

hiermod_marker <- "16S"
source('src/hiermod/0_SETUP.R')
source('src/hiermod/Models/MDSYz_model.R') # model_MDSYz_16S, means_MDSYz()
hiermod_out_dir <- "out/hiermod/16S_4_year_MDSY"

fit_MDSYz <- readRDS(file.path(hiermod_out_dir, "fit_MDSYz.rds"))
dat_MDSYz <- readRDS(file.path(hiermod_out_dir, "dat_MDSYz.rds"))

pf <- post_full(fit_MDSYz, means_MDSYz, shift = 1)
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

save_report("fit_summary", model_id_MDSYz, fit_MDSYz, model = model_MDSYz_16S)

p_contrast_mean <- contrast_plot_panels(
  pc_estimands_means, quant = c(0.001, 0.999), scales = 'free_y',
  group_pal = Management_palette,
  legend_title = "Posteriors (population means)",
  ratio_stats = names(Fold_change_palette), ratio_pal = Fold_change_palette); p_contrast_mean

save_gg("fit_contrast_mean", model_id_MDSYz, p_contrast_mean)

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
  pc_estimands_medians, quant =  c(0.001, 0.999), scales = 'free_y',
  group_pal = Management_palette,
  legend_title = "Posteriors (population medians)",
  ratio_stats = names(Fold_change_palette), ratio_pal = Fold_change_palette); p_contrast_median

save_gg("fit_contrast_median", model_id_MDSYz, p_contrast_median)

## Comprehensive posterior summary (results report) -----------------------------
# Every interpretable posterior's own mean/median/89% PI/HPDI, not just the
# headline contrast -- excludes per-tree tr[Tr] raw draws (too many, not
# individually interpretable) and the generic "mean"/"median" 4-column
# entries (already covered, more legibly, by pc_estimands_means/medians'
# own May/July/Conventional/Organic labels).

pc_all <- bind_rows(
  compute_contrasts(pf, keep = setdiff(names(pf), c("tr", "mean", "median")), group_levels = idx$Mg$levels),
  pc_estimands_means, pc_estimands_medians)
save_posterior_kable("results_report", model_id_MDSYz, pc_all)

## Residual variance (sigma[Mg]) -------------------------------------------------

pc_sigma <- compute_contrasts(pf, keep = "sigma", group_levels = idx$Mg$levels)

p_sigma <- contrast_plot_panels(
  pc_sigma, quant = c(0.005, 0.995), group_pal = Management_palette,
  legend_title = "Posteriors (residual SD, log scale)") +
  labs(x = "sigma[Mg]"); p_sigma

save_gg("fit_sigma_posterior", model_id_MDSYz, p_sigma)
