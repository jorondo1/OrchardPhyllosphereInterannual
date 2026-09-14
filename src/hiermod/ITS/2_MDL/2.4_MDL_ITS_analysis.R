# MODEL 2B (MDLb): posterior contrast, run against the saved fit.

hiermod_marker <- "ITS"
source('src/hiermod/0_SETUP.R')
source('src/hiermod/Models/MDL_model.R') # means_MDL()
hiermod_out_dir <- "out/hiermod/ITS_2_lognormal_MDL"

fitb <- readRDS(file.path(hiermod_out_dir, "fit_MDLb.rds"))
model_ppc1 <- readRDS(file.path(hiermod_out_dir, "model_ppc1.rds"))

## Posterior contrast ---------------------------------------------------------
# note: mean/PI on a raw-scale lognormal contrast are not robust.
# A handful of large-sigma posterior draws dominate the mean. Prefer
# median() + PI() for reporting.

pf <- post_full(fitb, means_MDL)      # raw draws + means_MDL's mean/median

mdl_labels <- c(
  median = "Median diversity",
  mean   = "Mean diversity",
  sigma  = "Residual SD (log scale)")
# sigma_loc isn't Mg-indexed (no group1/group2 pair) -- see it in
# precis(fitb, depth = 2) / fit_summary_MDLb.txt instead.

pc_full <- compute_contrasts(
  pf, keep = names(mdl_labels), labels = mdl_labels,
  group_levels = idx$Mg$levels)

# Plot posterior predictives and their contrasts
p_contrasts_panel <- contrast_plot_panels(
  quant = c(0.001, 0.999),
  pc_full, group_pal =  Management_palette); p_contrasts_panel

# proportion below zero?
pc_full %>%
  filter(group == "Contrast") %>%
  group_by(statistic) %>%
  summarise(n_neg = sum(value<=0), n = n()) %>%
  mutate(prop_neg = n_neg/n)

save_report("fit_summary", model_id, fitb, pc_full, model_ppc1, model_name = "The Tamed Wildcard")
save_gg("fit_contrasts_panel", model_id, p_contrasts_panel)
