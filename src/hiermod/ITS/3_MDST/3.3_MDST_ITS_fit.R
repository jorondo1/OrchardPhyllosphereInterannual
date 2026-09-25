# MODEL 3 (MDST, "Treebeard the Skeptic"), ITS, SHIFTED (Hill_1 - 1): real
# fit and PPC. Mirrors 3.3_MDST_16S_fit.R.
#
# Calibration caveat: 2.2/3.2's own n_sbc=500 run flagged loga[2] and
# sigma_tr MISCALIBRATED (z=-2.67/+2.64; 0 divergences, so a rank-fraction
# issue, not a sampling failure) -- consistent with why Tree was dropped
# again after this model in the 16S rebuild too (not merged back in until
# Model 7, once Cultivar/Year/Covariates are also in the model to help
# stabilize sigma_tr). Models 1-6 in this family are scaffolding toward
# Model 7, not independently-interpreted results -- 16S's own MDST showed
# comparable SBC fragility (a small-sigma_tr funnel, sparse 1-2 obs/tree)
# yet fit real data cleanly (0 divergences, Rhat=1.000) anyway, so this
# real fit is still worth having as a documented data point, same as 16S's.

hiermod_marker <- "ITS"
source('src/hiermod/0_SETUP.R')
source('src/hiermod/Models/MDST_model.R') # model_MDST_ITS, means_MDST()

hiermod_out_dir <- "out/hiermod/ITS_3_tree_MDST"

## Model fit ----------------------------------------------------------------

dat_MDST <- list(
  Dv = div$Hill_1 - 1, # subtract the floor because Dv must be (0, Inf)-support to match the likelihood
  Mg = idx$Mg$to_index(div$Management),
  Mo = idx$Mo$to_index(div$Time),
  Tr = idx$Tr$to_index(div$Tree_id)
)

fit_MDST <- ulam(
  model_MDST_ITS,
  data = dat_MDST,
  chains = 6, cores = 6, iter = 10000,
  control = list(adapt_delta = 0.99)
)
save_fit("fit", model_id_MDST, fit_MDST)
saveRDS(dat_MDST, file.path(hiermod_out_dir, "dat_MDST.rds"))

precis(fit_MDST, depth = 2)
save_pdf("fit_trankplot", model_id_MDST,
         function() trankplot(fit_MDST, n_cols = 8, max_rows = 30),
         width = 30, height = 50)

# num_divergent/num_max_treedepth/ebfmi per chain directly:
attr(fit_MDST, "cstanfit")$diagnostic_summary(
  diagnostics = c("divergences", "treedepth", "ebfmi"), quiet = TRUE)
# ebfmi should be comfortably >0.3

## Posterior predictive check --------------------------------------------------

### Overall, by Management x Season cell ----

pp_group <- ppc_group(dat_MDST)
p_postpred <- plot_ppc_overlay(fit_MDST, dat_MDST, pp_group, xlim = c(NA, 100))
save_gg("postpred_density", model_id_MDST, p_postpred)

### Contrast test statistics ----

(p_ppc <- plot_ppc_season_contrast_stats(fit_MDST, dat_MDST))
save_gg("postpred_stat", model_id_MDST, p_ppc)

## Added value of Tree: how much tree-to-tree spread is there? -----------------

pf <- post_full(fit_MDST, means_MDST, shift = 1)

pc_sigma_tr <- bind_rows(
  compute_contrasts(pf, keep = "sigma", group_levels = idx$Mg$levels),
  tibble(statistic = "sigma_tr", group = "Population", value = pf$sigma_tr$sigma_tr)
) %>% mutate(statistic = factor(statistic, levels = c("sigma", "sigma_tr")))

p_sigma_tr <- variance_component_panels(
  pc_sigma_tr, quant = c(0, 1),
  palette = Management_palette, # already carries Population's colour
  sd_stats = c("sigma", "sigma_tr")); p_sigma_tr

save_gg("fit_variance_components", model_id_MDST, p_sigma_tr, width = 12, height = 9)
