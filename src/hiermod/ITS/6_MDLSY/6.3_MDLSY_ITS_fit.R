# MODEL 6 (MDLSY): Year pooled, Cultivar added, sigma priors
# variance-budget-calibrated (K=3 -> K=4, see 6.2_MDLSY_calibration.R).
# Real fit and PPC.

hiermod_marker <- "ITS"
source('src/hiermod/0_SETUP.R')
source('src/hiermod/Models/MDLSY_model.R') # means_MDLSY(), sim_div_MDLSY()

hiermod_out_dir <- "out/hiermod/ITS_6_lognormal_MDLSY"
model_vbc <- readRDS(file.path(hiermod_out_dir, "model_vbc.rds"))

## Model fit ----------------------------------------------------------------

dat <- list(
  Dv = div$Hill_1 - 1, # subtract the floor because Dv must be (0, Inf)-support to match the likelihood
  Mg = idx$Mg$to_index(div$Management),
  Lo = idx$Lo$to_index(div$Location),
  Mo = idx$Mo$to_index(div$Time),
  Yr = idx$Yr$to_index(div$Year),
  Tr = idx$Tr$to_index(div$Tree_id),
  Cv = idx$Cv$to_index(div$Cultivar)
)
dat$cell <- (dat$Mg - 1) * 2 + dat$Mo

fit_MDLSY <- ulam(
  model_vbc,
  data = dat,
  chains = 6, cores = 6, iter = 20000,
  control = list(adapt_delta = 0.99)
)
save_fit("fit", model_id, fit_MDLSY)
saveRDS(dat, file.path(hiermod_out_dir, "dat_MDLSY.rds"))

precis(fit_MDLSY, depth = 2)

save_pdf("fit_traceplot", model_id, function() traceplot(fit_MDLSY, n_cols = 6, max_rows = 10))
save_pdf("fit_trankplot", model_id, function() trankplot(fit_MDLSY, n_cols = 6, max_rows = 10))


## Posterior predictive check --------------------------------------------------
post_MDLSY <- extract.samples(fit_MDLSY)

### Overall, by Management x Season cell ----

pp_group <- interaction(idx$Mg$to_label(dat$Mg), idx$Mo$to_label(dat$Mo), sep = " ")
(p_postpred <- plot_ppc_overlay(fit_MDLSY, dat, pp_group, xlim = c(0,150)))
save_gg("postpred_density", model_id, p_postpred)

### Contrast test statistics ----

(p_ppc <- plot_ppc_season_contrast_stats(fit_MDLSY, dat))
save_gg("postpred_stat", model_id, p_ppc)
