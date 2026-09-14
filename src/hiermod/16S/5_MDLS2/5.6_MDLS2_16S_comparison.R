# MODEL 5 (MDLS2), 16S: model comparison. Needs real posterior log-lik from
# both models, so this only makes sense once a real fit exists -- not
# something to run during calibration. Refits both the real model and a
# reduced alternative (sigma[Mg], 2 levels, same as Model 2/4's residual
# structure) at much lower iter than the real fits (log_lik = TRUE, not
# persisted to disk -- log_lik draws are massive) since PSIS needs far fewer
# draws than the reported contrasts do. Tests whether splitting residual SD
# by Season too (sigma[cell]) actually earns its predictive keep.

hiermod_marker <- "16S"
source('src/hiermod/0_SETUP.R')
source('src/hiermod/Models/MDLS2_model.R')

hiermod_out_dir <- "out/hiermod/16S_5_lognormal_MDLS2_shifted"
model_ppc1 <- readRDS(file.path(hiermod_out_dir, "model_ppc1.rds"))

dat <- list(
  Dv = div$Hill_1 - 1,
  Mg = idx$Mg$to_index(div$Management),
  Lo = idx$Lo$to_index(div$Location),
  Mo = idx$Mo$to_index(div$Time),
  Yr = idx$Yr$to_index(div$Year),
  Tr = idx$Tr$to_index(div$Tree_id)
)
dat$cell <- (dat$Mg - 1) * 2 + dat$Mo

model_reduced <- model_ppc1
model_reduced$likelihood <- quote(Dv ~ dlnorm(mu, sigma[Mg]))
model_reduced$pr_sigma   <- quote(sigma[Mg] ~ dexp(2.5))

fit_full_cmp <- ulam(
  model_ppc1,
  data = dat,
  chains = 6, cores = 6, iter = 2000,
  control = list(adapt_delta = 0.99),
  log_lik = TRUE
)

fit_reduced_cmp <- ulam(
  model_reduced,
  data = dat,
  chains = 6, cores = 6, iter = 2000,
  control = list(adapt_delta = 0.99),
  log_lik = TRUE
)

psis_compare(fit_full_cmp, fit_reduced_cmp)

# Dosn't work. Funnel effect is a problem, but not the only one.
# 

# "The gaussian bowl is the best to skate" - McElreath !
# MCMC is made to sample in bowls, not funnels.
