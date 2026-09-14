# MODEL 5 (MDLS2): residual SD varies by Management AND Season
# (sigma[Mg] -> sigma[cell]). Real fit and PPC.

hiermod_marker <- "ITS"
source('src/hiermod/0_SETUP.R')
source('src/hiermod/Models/MDLS2_model.R') # model_MDLS2_ITS, means_MDLS2(), sim_div_MDLS2()
model <- model_MDLS2_ITS

hiermod_out_dir <- "out/hiermod/ITS_5_lognormal_MDLS2"

## Model fit ----------------------------------------------------------------

dat <- list(
  Dv = div$Hill_1,
  Mg = idx$Mg$to_index(div$Management),
  Lo = idx$Lo$to_index(div$Location),
  Mo = idx$Mo$to_index(div$Time),
  Yr = idx$Yr$to_index(div$Year),
  Tr = idx$Tr$to_index(div$Tree_id)
)
dat$cell <- (dat$Mg - 1) * 2 + dat$Mo   # 1=Conv-May, 2=Conv-July, 3=Org-May, 4=Org-July

fit_MDLS2 <- ulam(
  model,
  data = dat,
  chains = 6, cores = 6, iter = 20000,
  control = list(adapt_delta = 0.99)
)
save_fit("fit", model_id_MDLS2_ITS, fit_MDLS2)

precis(fit_MDLS2, depth = 2)

save_pdf("fit_traceplot", model_id_MDLS2_ITS, function() traceplot(fit_MDLS2))
save_pdf("fit_trankplot", model_id_MDLS2_ITS, function() trankplot(fit_MDLS2))

## Posterior predictive check --------------------------------------------------

### Overall, by Management x Season cell ----

pp_group <- interaction(idx$Mg$to_label(dat$Mg), idx$Mo$to_label(dat$Mo), sep = " ")
(p_postpred <- plot_ppc_overlay(fit_MDLS2, dat, pp_group, xlim = c(0,150)))
save_gg("postpred_density", model_id_MDLS2_ITS, p_postpred)

# Compare against fit_contrast_density_MDLS.pdf (model 4)
# Conventional-July's yrep hug the observed spike much more closely than model 4's did.

### Contrast test statistics ----

(p_ppc <- plot_ppc_season_contrast_stats(fit_MDLS2, dat))

save_gg("postpred_stat", model_id_MDLS2_ITS, p_ppc)
