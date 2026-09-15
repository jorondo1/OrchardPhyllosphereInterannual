# MODEL 3 (MDST) vs MODEL 2 (MDS), 16S: model comparison. Quantifies the
# predictive value of accounting for repeated measures (Tree), the thing
# this rebuild has been trying to show/justify rather than assume.
#
# Same pattern as the archived 5.6_MDLS2_16S_comparison.R: needs real
# posterior log-lik from both models, so only makes sense once a real fit
# exists, not during calibration. Refits both models fresh at much lower
# iter than the reported fits (log_lik = TRUE, not persisted to disk --
# log_lik draws are massive) since PSIS needs far fewer draws than the
# reported contrasts do -- deliberately NOT reusing fit_MDS.rds/fit_MDST.rds
# (those were fit without log_lik).
#
# Unlike that archived comparison (which had to hand-patch one model into a
# reduced version), MDS and MDST already exist as two clean, independently
# validated model objects -- the "reduced" model here is simply MDS itself.

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
  chains = 6, cores = 6, iter = 2000,
  control = list(adapt_delta = 0.99),
  log_lik = TRUE
)

fit_MDST_cmp <- ulam(
  model_MDST_16S,
  data = dat_MDST,
  chains = 6, cores = 6, iter = 2000,
  control = list(adapt_delta = 0.99),
  log_lik = TRUE
)

psis_compare(fit_MDS_cmp, fit_MDST_cmp)
