# MODEL 6v (MDLSYv), 16S: real fit and PPC.

hiermod_marker <- "16S"
source('src/hiermod/0_SETUP.R')
source('src/hiermod/Models/MDLSYv_model.R') # model_MDLSYv_16S, means_MDLSYv()

hiermod_out_dir <- "out/hiermod/16S_6_lognormal_MDLSYv"
model <- model_MDLSYv_16S

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

fit_MDLSYv <- ulam(
  model,
  data = dat,
  chains = 6, cores = 6, iter = 20000,
  control = list(adapt_delta = 0.99)
)
save_fit("fit", model_id, fit_MDLSYv)
saveRDS(dat, file.path(hiermod_out_dir, "dat_MDLSYv.rds"))

precis(fit_MDLSYv, depth = 2)

save_pdf("fit_trankplot", model_id, function() trankplot(fit_MDLSYv, n_cols = 6, max_rows = 10))

## Posterior predictive check --------------------------------------------------
# dat$Dv is on the same (Hill_1 - 1) scale the model was fit to, and sim()
# generates yrep on that same scale.

pp_group <- interaction(idx$Mg$to_label(dat$Mg), idx$Mo$to_label(dat$Mo), sep = " ")
p_postpred <- plot_ppc_overlay(fit_MDLSYv, dat, pp_group, xlim = c(NA, 2000))
save_gg("postpred_density", model_id, p_postpred)

(p_ppc <- plot_ppc_season_contrast_stats(fit_MDLSYv, dat))
save_gg("postpred_stat", model_id, p_ppc)
