# MDST vs MDS, 16S: PSIS comparison -- predictive value of the Tree effect
# - refits both with log_lik = TRUE (not saved: log_lik draws are huge)
# - iter 2000 gave Pareto k > 0.7 for all points in both models -> raised iter
# - if still bad: rely on seasonal_change interval widths (3.4) and sigma_tr (3.3) instead

hiermod_marker <- "16S"
source('src/hiermod/0_SETUP.R')
source('src/hiermod/Models/MDS_model.R')   # model_MDS_16S
source('src/hiermod/Models/MDST_model.R')  # model_MDST_16S

hiermod_out_dir <- "out/hiermod/16S_3_tree_MDST"

dat_MDST <- list(
  Dv = div$Hill_1 - 1,
  Mg = idx$Mg$to_index(div$Management),
  Mo = idx$Mo$to_index(div$Time),
  Tr = idx$Tr$to_index(div$Tree_id)
)

fit_MDS_cmp <- ulam(
  model_MDS_16S,
  data = dat_MDST, # extra Tr column is simply unused by MDS's own formula
  chains = 6, cores = 6, iter = 10000,
  control = list(adapt_delta = 0.99),
  log_lik = TRUE
)

fit_MDST_cmp <- ulam(
  model_MDST_16S,
  data = dat_MDST,
  chains = 6, cores = 6, iter = 10000,
  control = list(adapt_delta = 0.99),
  log_lik = TRUE
)

psis_compare(fit_MDS_cmp, fit_MDST_cmp)
