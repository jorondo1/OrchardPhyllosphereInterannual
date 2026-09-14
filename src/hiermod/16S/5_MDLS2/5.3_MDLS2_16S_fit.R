# MODEL 5 (MDLS2), 16S, SHIFTED (Hill_1 - 1): real fit and PPC.

hiermod_marker <- "16S"
source('src/hiermod/0_SETUP.R')
source('src/hiermod/Models/MDLS2_model.R') # means_MDLS2(), sim_div_MDLS2()

hiermod_out_dir <- "out/hiermod/16S_5_lognormal_MDLS2_shifted"
model_ppc1 <- readRDS(file.path(hiermod_out_dir, "model_ppc1.rds"))

## Model fit ----------------------------------------------------------------

dat <- list(
  Dv = div$Hill_1 - 1, # subtract the floor because Dv must be (0, Inf)-support to match the likelihood
  Mg = idx$Mg$to_index(div$Management),
  Lo = idx$Lo$to_index(div$Location),
  Mo = idx$Mo$to_index(div$Time),
  Yr = idx$Yr$to_index(div$Year),
  Tr = idx$Tr$to_index(div$Tree_id)
)
dat$cell <- (dat$Mg - 1) * 2 + dat$Mo   # 1=Conv-May, 2=Conv-July, 3=Org-May, 4=Org-July

fit_MDLS2 <- ulam(
  model_ppc1,
  data = dat,
  chains = 6, cores = 6, iter = 20000,
  control = list(adapt_delta = 0.99)
)
save_fit("fit", model_id_MDLS2_16S, fit_MDLS2)
saveRDS(dat, file.path(hiermod_out_dir, "dat_MDLS2_shifted.rds"))

precis(fit_MDLS2, depth = 2)

#save_pdf("fit_traceplot", model_id_MDLS2_16S, function() traceplot(fit_MDLS2, n_cols = 6, max_rows = 10))
save_pdf("fit_trankplot", model_id_MDLS2_16S, function() trankplot(fit_MDLS2, n_cols = 6, max_rows = 10))

## Posterior predictive check --------------------------------------------------

### Overall, by Management x Season cell ----
# dat$Dv is on the same (Hill_1 - 1) scale the model was fit to, and sim()
# generates yrep on that same scale -- no +1 correction needed anywhere in
# this section, since both sides were shifted down by 1 consistently.

pp_group <- interaction(idx$Mg$to_label(dat$Mg), idx$Mo$to_label(dat$Mo), sep = " ")
p_postpred <- plot_ppc_overlay(fit_MDLS2, dat, pp_group, xlim = c(NA,2000)
                                )
save_gg("postpred_density", model_id_MDLS2_16S, p_postpred)

### Contrast test statistics ----

(p_ppc <- plot_ppc_season_contrast_stats(fit_MDLS2, dat))
save_gg("postpred_stat", model_id_MDLS2_16S, p_ppc)
