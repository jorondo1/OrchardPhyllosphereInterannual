# MODEL 4 (MDSYz), ITS, SHIFTED (Hill_1 - 1): real fit and PPC. Mirrors
# 4.3_MDSYz_16S_fit.R. model_MDSYz_ITS already carries its own validated
# priors (loga[Mg] ~ dnorm(2,2), everything else inherited from
# model_MDSYz_16S) -- no local override needed, confirmed clean by
# 4.2_MDSYz_ITS_calibration.R's own n_sbc=500 run (0 divergences,
# "Overall: OK").

hiermod_marker <- "ITS"
source('src/hiermod/0_SETUP.R')
source('src/hiermod/Models/MDSYz_model.R') # model_MDSYz_ITS, means_MDSYz()

hiermod_out_dir <- "out/hiermod/ITS_4_year_MDSY"

## Model fit ----------------------------------------------------------------

dat_MDSYz <- list(
  Dv = div$Hill_1 - 1, # subtract the floor because Dv must be (0, Inf)-support to match the likelihood
  Mg = idx$Mg$to_index(div$Management),
  Mo = idx$Mo$to_index(div$Time),
  Yr = idx$Yr$to_index(div$Year)
)

fit_MDSYz <- ulam(
  model_MDSYz_ITS,
  data = dat_MDSYz,
  chains = 6, cores = 6, iter = 10000,
  control = list(adapt_delta = 0.99)
)
save_fit("fit", model_id_MDSYz, fit_MDSYz)
saveRDS(dat_MDSYz, file.path(hiermod_out_dir, "dat_MDSYz.rds"))

precis(fit_MDSYz, depth = 2)
save_pdf("fit_trankplot", model_id_MDSYz, function() trankplot(fit_MDSYz, n_cols = 4, max_rows = 10))

## Posterior predictive check --------------------------------------------------

### Overall, by Management x Season cell ----

pp_group <- ppc_group(dat_MDSYz)
p_postpred <- plot_ppc_overlay(fit_MDSYz, dat_MDSYz, pp_group, xlim = c(NA, 100)); p_postpred
save_gg("postpred_density", model_id_MDSYz, p_postpred)

### Contrast test statistics ----

(p_ppc <- plot_ppc_season_contrast_stats(fit_MDSYz, dat_MDSYz))
save_gg("postpred_stat", model_id_MDSYz, p_ppc)

## Year effect (fixed, not pooled) ---------------------------------------------
# yr1/yr2 are the free parameters; yr3 = -(yr1+yr2) by construction
# (sum-to-zero). Shown together as one panel since all three are on the
# same log-scale footing, no separate hyper-SD to report (this is a fixed
# effect, not a variance component).

pf <- post_full(fit_MDSYz)
yr3 <- -(pf$yr1$yr1 + pf$yr2$yr2)

pc_year <- bind_rows(
  tibble(group = idx$Yr$levels[1], value = pf$yr1$yr1),
  tibble(group = idx$Yr$levels[2], value = pf$yr2$yr2),
  tibble(group = idx$Yr$levels[3], value = yr3)
) %>% mutate(statistic = "Year effect (log scale)")

p_year <- variance_component_panels(
  pc_year, quant = c(0, 1), palette = idx$Yr$palette); p_year

save_gg("fit_year_effects", model_id_MDSYz, p_year, width = 8, height = 4)
