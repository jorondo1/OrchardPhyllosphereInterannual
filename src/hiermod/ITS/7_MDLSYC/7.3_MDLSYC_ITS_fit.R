# MODEL 7 (MDLSYC): three control covariates added on top of Model 6.
# Real fit, real-data collinearity table, and PPC.

hiermod_marker <- "ITS"
source('src/hiermod/0_SETUP.R')
source('src/hiermod/Models/MDLSYC_model.R') # model_MDLSYC_ITS, means_MDLSYC(), sim_div_MDLSYC()
model <- model_MDLSYC_ITS

hiermod_out_dir <- "out/hiermod/ITS_7_lognormal_MDLSYC"

## Model fit ----------------------------------------------------------------

dat <- list(
  Dv = div$Hill_1 - 1, # subtract the floor because Dv must be (0, Inf)-support to match the likelihood
  Mg = idx$Mg$to_index(div$Management),
  Lo = idx$Lo$to_index(div$Location),
  Mo = idx$Mo$to_index(div$Time),
  Yr = idx$Yr$to_index(div$Year),
  Tr = idx$Tr$to_index(div$Tree_id),
  Cv = idx$Cv$to_index(div$Cultivar),
  deg_h_z = div$deg_h_z,
  precip_72h_z = div$precip_72h_z,
  seq_depth_z = div$seq_depth_z
)
dat$cell <- (dat$Mg - 1) * 2 + dat$Mo

fit_MDLSYC <- ulam(
  model,
  data = dat,
  chains = 6, cores = 6, iter = 20000,
  control = list(adapt_delta = 0.99)
)
save_fit("fit", model_id, fit_MDLSYC)
saveRDS(dat, file.path(hiermod_out_dir, "dat_MDLSYC.rds")) # so 7.4 doesn't rebuild it

precis(fit_MDLSYC, depth = 2)

save_pdf("fit_traceplot", model_id, function() traceplot(fit_MDLSYC, n_cols = 6, max_rows = 10))
save_pdf("fit_trankplot", model_id, function() trankplot(fit_MDLSYC, n_cols = 6, max_rows = 10))

# Real-data collinearity table: correlated draws are the expected symptom
# of the correlation noted above, not a red flag on their own. Spearman
# rho + significance stars are a descriptive summary of posterior
# dependence here, not a classical hypothesis test (draws aren't
# independent samples).
post_MDLSYC <- extract.samples(fit_MDLSYC)

param_labels <- c(
  b_deg     = "Degree-hours slope (weather control)",
  b_precip  = "Precipitation slope (weather control)",
  b_seq     = "Sequencing-depth slope (detection-effort control)",
  s_conv    = "Season shift, Conventional (May -> July)",
  gap_shift = "Season x Management interaction",
  yr_2022   = "Year 2022 effect"
)

cor_vars <- list(
  b_deg     = post_MDLSYC$b_deg,
  b_precip  = post_MDLSYC$b_precip,
  b_seq     = post_MDLSYC$b_seq,
  s_conv    = post_MDLSYC$s_conv,
  gap_shift = post_MDLSYC$gap_shift,
  yr_2022   = post_MDLSYC$yr[,1]
)

cor_pairs <- list(
  c("b_deg", "s_conv"), c("b_deg", "gap_shift"), c("b_deg", "yr_2022"),
  c("b_precip", "s_conv"), c("b_seq", "s_conv"), c("b_seq", "yr_2022"),
  c("b_deg", "b_precip"), c("b_deg", "b_seq"), c("b_precip", "b_seq")
)

collinearity_table <- posterior_cor_table(cor_vars, cor_pairs, param_labels)
knitr::kable(collinearity_table, digits = 3,
             caption = "Model 7 (MDLSYC): posterior collinearity, control covariates vs structural parameters")

## Posterior predictive check --------------------------------------------------

### Overall, by Management x Season cell ----

pp_group <- interaction(idx$Mg$to_label(dat$Mg), idx$Mo$to_label(dat$Mo), sep = " ")
(p_postpred <- plot_ppc_overlay(fit_MDLSYC, dat, pp_group, xlim = c(0,150)))
save_gg("postpred_density", model_id, p_postpred)

### Contrast test statistics ----

(p_ppc <- plot_ppc_season_contrast_stats(fit_MDLSYC, dat))
save_gg("postpred_stat", model_id, p_ppc)
