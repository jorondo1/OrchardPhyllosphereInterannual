# MODEL 6 (MDSYCV, "Bombadil the Eldest"), ITS, SHIFTED (Hill_1 - 1): real
# fit and PPC. Mirrors 6.3_MDSYCV_16S_fit.R. model_MDSYCV_ITS already
# carries its own validated priors (loga[Mg] ~ dnorm(2,2), everything else
# inherited from model_MDSYCV_16S) -- no local override needed, confirmed
# clean by 6.2_MDSYCV_ITS_calibration.R's own n_sbc=500 run (0 divergences,
# "Overall: OK").

hiermod_marker <- "ITS"
source('src/hiermod/0_SETUP.R')
source('src/hiermod/Models/MDSYCV_model.R') # model_MDSYCV_ITS, means_MDSYCV()

hiermod_out_dir <- "out/hiermod/ITS_6_cultivar_MDSYCV"

## Model fit ----------------------------------------------------------------

dat_MDSYCV <- list(
  Dv = div$Hill_1 - 1, # subtract the floor because Dv must be (0, Inf)-support to match the likelihood
  Mg = idx$Mg$to_index(div$Management),
  Mo = idx$Mo$to_index(div$Time),
  Yr = idx$Yr$to_index(div$Year),
  Cv = idx$Cv$to_index(div$Cultivar),
  deg_h_z = div$deg_h_z,
  precip_72h_z = div$precip_72h_z,
  seq_depth_z = div$seq_depth_z
)

fit_MDSYCV <- ulam(
  model_MDSYCV_ITS,
  data = dat_MDSYCV,
  chains = 6, cores = 6, iter = 10000,
  control = list(adapt_delta = 0.99)
)
save_fit("fit", model_id_MDSYCV, fit_MDSYCV)
saveRDS(dat_MDSYCV, file.path(hiermod_out_dir, "dat_MDSYCV.rds"))

precis(fit_MDSYCV, depth = 2)

save_pdf("fit_trankplot", model_id_MDSYCV,
         function() trankplot(fit_MDSYCV, n_cols = 4, max_rows = 10))

## Posterior predictive check --------------------------------------------------

### Overall, by Management x Season cell ----

pp_group <- ppc_group(dat_MDSYCV)
p_postpred <- plot_ppc_overlay(fit_MDSYCV, dat_MDSYCV, pp_group, xlim = c(NA, 100)); p_postpred
save_gg("postpred_density", model_id_MDSYCV, p_postpred)

### Contrast test statistics ----

(p_ppc <- plot_ppc_season_contrast_stats(fit_MDSYCV, dat_MDSYCV))
save_gg("postpred_stat", model_id_MDSYCV, p_ppc)

## Year effect (fixed, not pooled) ---------------------------------------------
# Same as 4.3/5.3's own panel -- yr1/yr2 free, yr3 = -(yr1+yr2) by construction.

pf <- post_full(fit_MDSYCV)
yr3 <- -(pf$yr1$yr1 + pf$yr2$yr2)

pc_year <- bind_rows(
  tibble(group = idx$Yr$levels[1], value = pf$yr1$yr1),
  tibble(group = idx$Yr$levels[2], value = pf$yr2$yr2),
  tibble(group = idx$Yr$levels[3], value = yr3)
) %>% mutate(statistic = "Year effect (log scale)")

p_year <- variance_component_panels(
  pc_year, quant = c(0, 1), palette = idx$Yr$palette); p_year

## Covariate effects (b_deg, b_precip, b_seq) ----------------------------------

pc_covariates <- bind_rows(
  tibble(group = cov_labels[1], value = pf$b_deg$b_deg),
  tibble(group = cov_labels[2], value = pf$b_precip$b_precip),
  tibble(group = cov_labels[3], value = pf$b_seq$b_seq)
) %>% mutate(statistic = "Covariate effects (log scale)")

p_covariates <- variance_component_panels(
  pc_covariates, quant = c(0, 1), palette = cov_pal); p_covariates

## Cultivar effect (fixed, not pooled) -----------------------------------------
# cv_1/cv_3/cv_4/cv_5 free; cv_2 (Liberty) = -(cv_1+cv_3+cv_4+cv_5) by
# construction (Liberty has the most combined observations -- see
# MDSYCV_model.R's own header for why it's the derived slot).

cv_2 <- -(pf$cv_1$cv_1 + pf$cv_3$cv_3 + pf$cv_4$cv_4 + pf$cv_5$cv_5)

pc_cultivar <- bind_rows(
  tibble(group = idx$Cv$levels[1], value = pf$cv_1$cv_1),
  tibble(group = idx$Cv$levels[2], value = cv_2),
  tibble(group = idx$Cv$levels[3], value = pf$cv_3$cv_3),
  tibble(group = idx$Cv$levels[4], value = pf$cv_4$cv_4),
  tibble(group = idx$Cv$levels[5], value = pf$cv_5$cv_5)
) %>% mutate(statistic = "Cultivar effect (log scale)")

p_cultivar <- variance_component_panels(
  pc_cultivar, quant = c(0, 1), palette = idx$Cv$palette); p_cultivar

## Year/Covariate/Cultivar effects, combined -------------------------------------

p_effects <- p_year / p_covariates / p_cultivar
save_gg("fit_effects", model_id_MDSYCV, p_effects, width = 8, height = 14)

## Full pairwise parameter check ------------------------------------------------
# Same rationale/settings as 6.2's own calibration-stage check -- PNG, not
# PDF, n_keep tuned to the standard 1000-draw thinning target for this
# family.

pairs_vars <- c("loga[1]", "loga[2]", "sigma[1]", "sigma[2]", "yr1", "yr2",
                 "cv_1", "cv_3", "cv_4", "cv_5", "b_deg", "b_precip", "b_seq")
p_pairs <- plot_mcmc_pairs(fit_MDSYCV, variables = pairs_vars, n_keep = 1000)
save_gg("fit_mcmc_pairs", model_id_MDSYCV, p_pairs, width = 15, height = 15, type = "png")
