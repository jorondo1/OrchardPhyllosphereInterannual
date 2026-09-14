# MODEL 6 (MDLSY) vs MODEL 5, shifted: PSIS model comparison. Refits both at
# much lower iter than the real fits (log_lik = TRUE, not persisted to disk
# -- log_lik draws are massive) since PSIS needs far fewer draws than the
# reported contrasts do. Compared against 5b (shifted), not plain 5 -- Model
# 6 onward fits Hill_1 - 1, so only the shifted variant is on the same
# log-lik scale. Tests whether pooling Year + adding Cultivar beats Model
# 5's fixed/unpooled Year predictively.

hiermod_marker <- "ITS"
source('src/hiermod/0_SETUP.R')
source('src/hiermod/Models/MDLSY_model.R') # sources MDLS2_model.R too -> model_MDLS2_ITS, model_MDLSY_ITS

dat <- list(
  Dv = div$Hill_1 - 1,
  Mg = idx$Mg$to_index(div$Management),
  Lo = idx$Lo$to_index(div$Location),
  Mo = idx$Mo$to_index(div$Time),
  Yr = idx$Yr$to_index(div$Year),
  Tr = idx$Tr$to_index(div$Tree_id),
  Cv = idx$Cv$to_index(div$Cultivar)
)
dat$cell <- (dat$Mg - 1) * 2 + dat$Mo

# Model 6's actual fitted spec (variance-budget-calibrated, see
# 6.2_MDLSY_ITS_calibration.R) -- model_MDLSY_ITS alone is the naive,
# pre-VBC spec.
model_vbc <- model_MDLSY_ITS
model_vbc$pr_sigma     <- quote(sigma[cell] ~ dexp(3.46))
model_vbc$pr_sigma_loc <- quote(sigma_loc   ~ dexp(2.31))
model_vbc$pr_sigma_tr  <- quote(sigma_tr    ~ dexp(2.31))
model_vbc$pr_sigma_yr  <- quote(sigma_yr    ~ dexp(2.31))

fit_MDLSY_cmp <- ulam(
  model_vbc,
  data = dat,
  chains = 6, cores = 6, iter = 2000,
  control = list(adapt_delta = 0.99),
  log_lik = TRUE
)

fit_MDLS2_cmp <- ulam(
  model_MDLS2_ITS,
  data = dat,
  chains = 6, cores = 6, iter = 2000,
  control = list(adapt_delta = 0.99),
  log_lik = TRUE
)

psis_compare(fit_MDLSY_cmp, fit_MDLS2_cmp)
