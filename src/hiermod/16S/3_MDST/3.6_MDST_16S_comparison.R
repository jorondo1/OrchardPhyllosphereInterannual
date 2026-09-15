# MODEL 3 (MDST) vs MODEL 2 (MDS), 16S: model comparison. Quantifies the
# predictive value of accounting for repeated measures (Tree), the thing
# this rebuild has been trying to show/justify rather than assume.
#
# Same pattern as the archived 5.6_MDLS2_16S_comparison.R: needs real
# posterior log-lik from both models, so only makes sense once a real fit
# exists, not during calibration. Refits both models fresh (log_lik = TRUE,
# not persisted to disk -- log_lik draws are massive) -- deliberately NOT
# reusing fit_MDS.rds/fit_MDST.rds (those were fit without log_lik).
#
# iter was originally 2000 (PSIS usually needs far fewer draws than the
# reported contrasts do) but that reproduced 5.6_MDLS2's own exact failure:
# Pareto k > 0.7 for 242/242 points on BOTH models, some > 1 -- a complete
# PSIS breakdown, not a few influential points. Both models failing
# identically (not just the one with per-tree latents) argues against
# "weakly-identified tree offsets" as the cause and toward plain
# draw-count fragility instead, so bumped to match the real fits' own iter
# before concluding PSIS just doesn't work here the way it didn't in 5.6.
# If Pareto k is still bad at this iter, stop trusting PSIS for this model
# family/sample size -- fall back to comparing MDS vs MDST's own real-data
# interval widths on `seasonal_change` (3.4's own estimand) plus sigma_tr's
# posterior magnitude (3.3) as the evidence for Tree's value instead.
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
