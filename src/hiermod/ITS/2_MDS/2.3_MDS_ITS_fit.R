# MODEL 2 (MDS), ITS: real fit (Hill_1 - 1) and PPC

hiermod_marker <- "ITS"
source('src/hiermod/0_SETUP.R')
source('src/hiermod/Models/MDS_model.R') # model_MDS_ITS, means_MDS()

hiermod_out_dir <- "out/hiermod/ITS_2_interaction_MDS"

## Model fit ----------------------------------------------------------------

dat_MDS <- list(
  Dv = div$Hill_1 - 1, # subtract the floor because Dv must be (0, Inf)-support to match the likelihood
  Mg = idx$Mg$to_index(div$Management),
  Mo = idx$Mo$to_index(div$Time)
)

fit_MDS <- ulam(
  model_MDS_ITS,
  data = dat_MDS,
  chains = 6, cores = 6, iter = 10000,
  control = list(adapt_delta = 0.99)
)
save_fit("fit", model_id_MDS, fit_MDS)
saveRDS(dat_MDS, file.path(hiermod_out_dir, "dat_MDS.rds"))

precis(fit_MDS, depth = 2)
save_pdf("fit_trankplot", model_id_MDS,
         function() trankplot(fit_MDS, n_cols = 4, max_rows = 10))

## Posterior predictive check --------------------------------------------------

### Overall, by Management x Season cell ----

pp_group <- ppc_group(dat_MDS)
p_postpred <- plot_ppc_overlay(fit_MDS, dat_MDS, pp_group, xlim = c(NA, 100))
save_gg("postpred_density", model_id_MDS, p_postpred)

### Contrast test statistics ----

(p_ppc <- plot_ppc_season_contrast_stats(fit_MDS, dat_MDS))
save_gg("postpred_stat", model_id_MDS, p_ppc)
