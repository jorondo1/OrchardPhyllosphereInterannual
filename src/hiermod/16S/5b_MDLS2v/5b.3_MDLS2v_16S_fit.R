# MODEL 5v (MDLS2v), 16S, SHIFTED (Hill_1 - 1): real fit and PPC.

hiermod_marker <- "16S"
source('src/hiermod/0_SETUP.R')
source('src/hiermod/Models/MDLS2v_model.R') # model_MDLS2v_16S, means_MDLS2v()

hiermod_out_dir <- "out/hiermod/16S_5b_lognormal_MDLS2v_shifted"
model <- model_MDLS2v_16S

## Model fit ----------------------------------------------------------------

dat <- list(
  Dv = div$Hill_1 - 1, # subtract the floor because Dv must be (0, Inf)-support to match the likelihood
  Mg = idx$Mg$to_index(div$Management),
  Lo = idx$Lo$to_index(div$Location),
  Mo = idx$Mo$to_index(div$Time),
  Yr = idx$Yr$to_index(div$Year),
  Tr = idx$Tr$to_index(div$Tree_id)
)

fit_MDLS2v <- ulam(
  model,
  data = dat,
  chains = 6, cores = 6, iter = 20000,
  control = list(adapt_delta = 0.99)
)

save_fit("fit", model_id, fit_MDLS2v)
saveRDS(dat, file.path(hiermod_out_dir, "dat_MDLS2v.rds"))

p <- precis(fit_MDLS2v, depth = 2); p$ess_ratio <- round(p$ess_bulk / NROW(extract.samples(fit_MDLS2v)[[1]]), 2); p

save_pdf(width = 30, height = 50,
  "fit_trankplot", model_id, 
  function() trankplot(fit_MDLS2v, n_cols = 6, max_rows = 25))

## Posterior predictive check --------------------------------------------------

### Overall, by Management x Season cell ----
# dat$Dv is on the same (Hill_1 - 1) scale the model was fit to, and sim()
# generates yrep on that same scale -- no +1 correction needed anywhere in
# this section, since both sides were shifted down by 1 consistently.

pp_group <- interaction(idx$Mg$to_label(dat$Mg), idx$Mo$to_label(dat$Mo), sep = " ")
p_postpred <- plot_ppc_overlay(fit_MDLS2v, dat, pp_group, xlim = c(NA,2000))
save_gg("postpred_density", model_id, p_postpred)

### Contrast test statistics ----

(p_ppc <- plot_ppc_season_contrast_stats(fit_MDLS2v, dat))
save_gg("postpred_stat", model_id, p_ppc)
