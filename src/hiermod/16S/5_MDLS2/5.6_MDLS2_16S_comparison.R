# MODEL 5 (MDLS2), 16S: model comparison. Needs real posterior log-lik from
# both models, so this only makes sense once a real fit exists -- not
# something to run during calibration. Fits a reduced alternative (sigma[Mg],
# 2 levels, same as Model 2/4's residual structure) to the same real data and
# compares against Model 5's sigma[cell] (4 levels) via PSIS-LOO: does
# splitting residual SD by Season too actually earn its predictive keep?
#
# ulam() defaults to log_lik = FALSE, so neither the saved fit_MDLS2_shifted
# nor a plain reduced fit carries the log-lik PSIS/compare() need -- both
# models are refit here with log_lik = TRUE instead of reusing saved fits.
# iter dropped from the main fit's 20000 to 8000: PSIS/WAIC need far fewer
# draws for a stable estimate than the reported contrasts do.

hiermod_marker <- "16S"
source('src/hiermod/0_SETUP.R')
source('src/hiermod/Models/MDLS2_model.R')

hiermod_out_dir <- "out/hiermod/16S_5_lognormal_MDLS2_shifted"
dat        <- readRDS(file.path(hiermod_out_dir, "dat_MDLS2_shifted.rds"))
model_ppc1 <- readRDS(file.path(hiermod_out_dir, "model_ppc1.rds"))

model_reduced <- model_ppc1
model_reduced$likelihood <- quote(Dv ~ dlnorm(mu, sigma[Mg]))
model_reduced$pr_sigma   <- quote(sigma[Mg] ~ dexp(2.5))

fit_full <- ulam(
  model_ppc1,
  data = dat,
  chains = 6, cores = 6, iter = 8000,
  control = list(adapt_delta = 0.99),
  log_lik = TRUE
)
save_fit("fit", "MDLS2_shifted_full_loglik", fit_full)

fit_reduced <- ulam(
  model_reduced,
  data = dat,
  chains = 6, cores = 6, iter = 8000,
  control = list(adapt_delta = 0.99),
  log_lik = TRUE
)
save_fit("fit", "MDLS2_shifted_reduced_loglik", fit_reduced)

rethinking::compare(fit_full, fit_reduced, func = PSIS)
