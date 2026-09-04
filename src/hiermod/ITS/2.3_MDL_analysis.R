# 2.3_MDL_analysis.R -- MODEL 2B (MDLb): posterior contrast, run against the
# fit saved by 2.2_MDL_validation.R -- no refit needed.

source('src/hiermod/ITS/0_SETUP.R')
source('src/hiermod/ITS/2.1_MDL_model.R') # model, means_MDL()
hiermod_out_dir <- "out/hiermod/ITS_2_lognormal_MDL"

fitb <- readRDS(file.path(hiermod_out_dir, "fit_MDLb.rds"))

dat <- list(
  Dv = div$Hill_1,
  Mg = idx$Mg$to_index(div$Management),
  Lo = idx$Lo$to_index(div$Location)
)

# model_ppc1 (the tightened-prior variant, see 2.2_MDL_validation.R) is what
# fitb actually is -- rebuild it here purely to label the fit_summary report.
model_ppc1 <- model
model_ppc1$prior_sigma = sigma[Mg] ~ dexp(3)
model_ppc1$prior_sigma_loc = sigma_loc ~ dexp(2)

## Posterior contrast ---------------------------------------------------------
# note: mean/PI on a raw-scale lognormal contrast are not robust.
# A handful of large-sigma posterior draws (expected, see prior predictive
# check in 2.2_MDL_validation.R) dominate the mean. Prefer median() + PI()
# for reporting.

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

save_report("fit_summary", "MDLb", fitb, pc_full, model_ppc1, model_name = "The Tamed Wildcard")
save_gg("fit_contrasts_panel", "MDLb", p_contrasts_panel)
