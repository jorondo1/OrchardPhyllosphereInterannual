# MODEL 1b (MDv), ITS: real fit and PPC
# - dexp(1) sigma prior kept: calibrated fine for ITS (no override, unlike 16S)

hiermod_marker <- "ITS"
source('src/hiermod/0_SETUP.R')
source('src/hiermod/Models/MD_model.R') # model_MDv_ITS, means_MDv(), sim_div_M()
model_MDv <- model_MDv_ITS

hiermod_out_dir <- "out/hiermod/ITS_1_lognormal_MD"

dat_MD <- list(
  Dv = div$Hill_1,
  Mg = idx$Mg$to_index(div$Management)
)

## MDv -- Allow Management-specific variance (heteroscedasticity) ============

### Model fit ----------------------------------------------------------------

fit_MDv <- ulam(
  model_MDv,
  data = dat_MD,
  chains = 6, cores = 6, iter = 10000
)
save_fit("fit", model_id_MDv, fit_MDv)

precis(fit_MDv, depth = 2)
save_pdf("fit_traceplot", model_id_MDv, function() traceplot(fit_MDv))

### Posterior predictive check --------------------------------------------------

(p_postpred_MDv <- plot_ppc_overlay(fit_MDv, dat_MD, idx$Mg$to_label(dat_MD$Mg), xlim = c(NA, 100)))
save_gg("postpred_density", model_id_MDv, p_postpred_MDv)
