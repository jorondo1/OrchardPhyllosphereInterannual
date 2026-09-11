# MODEL 1 (MD, constant variance) + MODEL 2 (MDv, Management-specific
# variance): real fit and PPC.

hiermod_marker <- "ITS"
source('src/hiermod/0_SETUP.R')
source('src/hiermod/Models/MD_model.R') # model_MD_ITS/model_MDv_ITS, means_MD/means_MDv, sim_div_M()
model_MD  <- model_MD_ITS
model_MDv <- model_MDv_ITS

hiermod_out_dir <- "out/hiermod/ITS_1_lognormal_MD"

## MD -- Mean difference by Management, constant variance ====================

### Model fit ----------------------------------------------------------------

dat_MD <- list(
  Dv = div$Hill_1,
  Mg = idx$Mg$to_index(div$Management)
)

fit_MD <- ulam(
  model_MD,
  data = dat_MD,
  chains = 6, cores = 6, iter = 10000
)
save_fit("fit", "MD", fit_MD)

precis(fit_MD, depth = 2)
traceplot(fit_MD)
save_pdf("fit_traceplot", "MD", function() traceplot(fit_MD))

### Posterior predictive check --------------------------------------------------

p_postpred_MD <- plot_ppc_overlay(fit_MD, dat_MD, idx$Mg$to_label(dat_MD$Mg))
save_gg("postpred_density", "MD", p_postpred_MD)

## MDv -- Allow Management-specific variance (heteroscedasticity) ============

### Model fit ----------------------------------------------------------------

fit_MDv <- ulam(
  model_MDv,
  data = dat_MD,
  chains = 6, cores = 6, iter = 10000
)
save_fit("fit", "MDv", fit_MDv)

precis(fit_MDv, depth = 2)
traceplot(fit_MDv)
save_pdf("fit_traceplot", "MDv", function() traceplot(fit_MDv))

### Posterior predictive check --------------------------------------------------
# Barely affected relative to MD, but the tail is not as heavy.

p_postpred_MDv <- plot_ppc_overlay(fit_MDv, dat_MD, idx$Mg$to_label(dat_MD$Mg))
save_gg("postpred_density", "MDv", p_postpred_MDv)
