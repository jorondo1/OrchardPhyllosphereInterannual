# MODEL 1 (MD, constant variance) + MODEL 2 (MDv, Management-specific
# variance), 16S: real fit and PPC.

hiermod_marker <- "16S"
source('src/hiermod/0_SETUP.R')
source('src/hiermod/Models/MD_model.R') # model_MD_16S/model_MDv_16S, means_MD/means_MDv, sim_div_M()
model_MD  <- model_MD_16S
model_MDv <- model_MDv_16S
# MDv's file default (dexp(2)) is only the calibration investigation's
# starting point, not its answer -- SBC (1.2_MDv_16S_calibration.R) found
# dhalfnorm(0,1) is what actually calibrates, so the real-data fit applies
# it here as a local override, the same way every calibration script did.
model_MDv$prior_sigma <- quote(sigma[Mg] ~ dhalfnorm(0,1))

hiermod_out_dir <- "out/hiermod/16S_1_lognormal_MD"

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
# traceplot(fit_MDv)
save_pdf("fit_traceplot", model_id_MDv, function() traceplot(fit_MDv))

### Posterior predictive check --------------------------------------------------
# Barely affected relative to MD, but the tail is not as heavy.

(p_postpred_MDv <- plot_ppc_overlay(fit_MDv, dat_MD, idx$Mg$to_label(dat_MD$Mg),  xlim = c(NA, 5000)))
save_gg("postpred_density", model_id_MDv, p_postpred_MDv)
