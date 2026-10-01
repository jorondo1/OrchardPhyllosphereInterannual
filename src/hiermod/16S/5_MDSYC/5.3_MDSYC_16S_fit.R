# MODEL 5 (MDSYC, "Radagast the Grower"), 16S: real fit (Hill_1 - 1) and PPC

hiermod_marker <- "16S"
source('src/hiermod/0_SETUP.R')
source('src/hiermod/Models/MDSYC_model.R') # model_MDSYC_16S, means_MDSYC()

hiermod_out_dir <- "out/hiermod/16S_5_covariates_MDSYC"

## Model fit ----------------------------------------------------------------

dat_MDSYC <- list(
  Dv = div$Hill_1 - 1, # subtract the floor because Dv must be (0, Inf)-support to match the likelihood
  Mg = idx$Mg$to_index(div$Management),
  Mo = idx$Mo$to_index(div$Time),
  Yr = idx$Yr$to_index(div$Year),
  deg_h_z = div$deg_h_z,
  precip_72h_z = div$precip_72h_z,
  seq_depth_z = div$seq_depth_z
)

fit_MDSYC <- ulam(
  model_MDSYC_16S,
  data = dat_MDSYC,
  chains = 6, cores = 6, iter = 10000,
  control = list(adapt_delta = 0.99)
)
save_fit("fit", model_id_MDSYC, fit_MDSYC)
saveRDS(dat_MDSYC, file.path(hiermod_out_dir, "dat_MDSYC.rds"))

precis(fit_MDSYC, depth = 2)
save_pdf("fit_trankplot", model_id_MDSYC, 
         function() trankplot(fit_MDSYC, n_cols = 4, max_rows = 10))

## Posterior predictive check --------------------------------------------------

### Overall, by Management x Season cell ----

pp_group <- ppc_group(dat_MDSYC)
p_postpred <- plot_ppc_overlay(fit_MDSYC, dat_MDSYC, pp_group, xlim = c(NA, 2000)); p_postpred
save_gg("postpred_density", model_id_MDSYC, p_postpred)

### Contrast test statistics ----

(p_ppc <- plot_ppc_season_contrast_stats(fit_MDSYC, dat_MDSYC))
save_gg("postpred_stat", model_id_MDSYC, p_ppc)

## Year effect (fixed, not pooled) ---------------------------------------------
# yr1/yr2 free; yr3 = -(yr1 + yr2)

pf <- post_full(fit_MDSYC)
yr3 <- -(pf$yr1$yr1 + pf$yr2$yr2)

pc_year <- bind_rows(
  tibble(group = idx$Yr$levels[1], value = pf$yr1$yr1),
  tibble(group = idx$Yr$levels[2], value = pf$yr2$yr2),
  tibble(group = idx$Yr$levels[3], value = yr3)
) %>% mutate(statistic = "Year effect (log scale)")

p_year <- variance_component_panels(
  pc_year, quant = c(0, 1), palette = idx$Yr$palette); p_year

save_gg("fit_year_effects", model_id_MDSYC, p_year, width = 8, height = 4)

## Covariate effects (b_deg, b_precip, b_seq) ----------------------------------
# New here: covariate slopes, direction and credibility

pc_covariates <- bind_rows(
  tibble(group = cov_labels[1], value = pf$b_deg$b_deg),
  tibble(group = cov_labels[2], value = pf$b_precip$b_precip),
  tibble(group = cov_labels[3], value = pf$b_seq$b_seq)
) %>% mutate(statistic = "Covariate effects (log scale)")

p_covariates <- variance_component_panels(
  pc_covariates, quant = c(0, 1), palette = cov_pal); p_covariates

save_gg("fit_covariate_effects", model_id_MDSYC, p_covariates, width = 8, height = 4)
