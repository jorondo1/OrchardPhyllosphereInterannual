# MODEL 5 (MDLS2) vs MODEL 4 (MDLS): PSIS model comparison. Refits both at
# much lower iter than the real fits (log_lik = TRUE, not persisted to disk
# -- log_lik draws are massive) since PSIS needs far fewer draws than the
# reported contrasts do. Both fit unshifted Hill_1, so their log-lik is on
# the same scale. Tests whether sigma[cell] (Model 5) beats sigma[Mg]
# (Model 4) predictively.

hiermod_marker <- "ITS"
source('src/hiermod/0_SETUP.R')
source('src/hiermod/Models/MDLS2_model.R') # sources MDLS_model.R too -> model_MDLS_ITS, model_MDLS2_ITS

dat <- list(
  Dv = div$Hill_1,
  Mg = idx$Mg$to_index(div$Management),
  Lo = idx$Lo$to_index(div$Location),
  Mo = idx$Mo$to_index(div$Time),
  Yr = idx$Yr$to_index(div$Year),
  Tr = idx$Tr$to_index(div$Tree_id)
)
dat$cell <- (dat$Mg - 1) * 2 + dat$Mo

fit_MDLS2_cmp <- ulam(
  model_MDLS2_ITS,
  data = dat,
  chains = 6, cores = 6, iter = 2000,
  control = list(adapt_delta = 0.99),
  log_lik = TRUE
)

fit_MDLS_cmp <- ulam(
  model_MDLS_ITS,
  data = dat,
  chains = 6, cores = 6, iter = 2000,
  control = list(adapt_delta = 0.99),
  log_lik = TRUE
)

psis_compare(fit_MDLS2_cmp, fit_MDLS_cmp)
