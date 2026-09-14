# MODEL 4 (MDLS): Management effect split by Season, with Tree random
# effects and Year as a fixed effect. Real fit and PPC.

hiermod_marker <- "ITS"
source('src/hiermod/0_SETUP.R')
source('src/hiermod/Models/MDLS_model.R') # model_MDLS_ITS, means_MDLS(), sim_div_MDLS()
model <- model_MDLS_ITS

hiermod_out_dir <- "out/hiermod/ITS_4_lognormal_MDLS"

## Model fit ----------------------------------------------------------------

dat <- list(
  Dv = div$Hill_1,
  Mg = idx$Mg$to_index(div$Management),
  Lo = idx$Lo$to_index(div$Location),
  Mo = idx$Mo$to_index(div$Time),
  Yr = idx$Yr$to_index(div$Year),
  Tr = idx$Tr$to_index(div$Tree_id)
)

fit_MDLS <- ulam(
  model,
  data = dat,
  chains = 6, cores = 6, iter = 10000,
  control = list(adapt_delta = 0.99)
)
save_fit("fit", model_id, fit_MDLS)

precis(fit_MDLS, depth = 2)

traceplot(fit_MDLS); trankplot(fit_MDLS)
save_pdf("fit_traceplot", model_id, function() traceplot(fit_MDLS))
save_pdf("fit_trankplot", model_id, function() trankplot(fit_MDLS))

## Posterior predictive check --------------------------------------------------

### Overall, by Management x Season cell ----
# ppc_dens_overlay_grouped() (inside plot_ppc_overlay()) isn't limited to a
# binary group, so the 4-cell interaction can be shown in one call.

pp_group <- interaction(idx$Mg$to_label(dat$Mg), idx$Mo$to_label(dat$Mo), sep = " ")
(p_postpred <- plot_ppc_overlay(fit_MDLS, dat, pp_group, xlim = c(0,150)))
save_gg("postpred_density", model_id, p_postpred)

### Contrast test statistics ----
# contrast_stat()/plot_ppc_contrast_stat() only handle a single binary group
# (dat$Mg); with a 4-cell Mg x Mo design, write the stat functions directly.

(p_ppc <- plot_ppc_season_contrast_stats(fit_MDLS, dat))

save_gg("postpred_stat", model_id, p_ppc)
