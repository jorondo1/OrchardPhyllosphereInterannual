# MODEL 7 (MDSTYCV, "Saruman the Fool"), 16S, SHIFTED (Hill_1 - 1): real fit
# and PPC. model_MDSTYCV_16S already carries its own validated priors, no
# local override needed here.

hiermod_marker <- "16S"
source('src/hiermod/0_SETUP.R')
source('src/hiermod/Models/MDSTYCV_model.R') # model_MDSTYCV_16S, means_MDSTYCV()

hiermod_out_dir <- "out/hiermod/16S_7_tree_full_MDSTYCV"

## Model fit ----------------------------------------------------------------

dat_MDSTYCV <- list(
  Dv = div$Hill_1 - 1, # subtract the floor because Dv must be (0, Inf)-support to match the likelihood
  Mg = idx$Mg$to_index(div$Management),
  Mo = idx$Mo$to_index(div$Time),
  Yr = idx$Yr$to_index(div$Year),
  Cv = idx$Cv$to_index(div$Cultivar),
  Tr = idx$Tr$to_index(div$Tree_id),
  deg_h_z = div$deg_h_z,
  precip_72h_z = div$precip_72h_z,
  seq_depth_z = div$seq_depth_z
)

fit_MDSTYCV <- ulam(
  model_MDSTYCV_16S,
  data = dat_MDSTYCV,
  chains = 6, cores = 6, iter = 10000,
  control = list(adapt_delta = 0.99)
)
save_fit("fit", model_id_MDSTYCV, fit_MDSTYCV)
saveRDS(dat_MDSTYCV, file.path(hiermod_out_dir, "dat_MDSTYCV.rds"))

precis(fit_MDSTYCV, depth = 2)

save_pdf("fit_trankplot", model_id_MDSTYCV,
         function() trankplot(fit_MDSTYCV, n_cols = 6, max_rows = 20),
         width = 24, height = 36)

## Posterior predictive check --------------------------------------------------

### Overall, by Management x Season cell ----

pp_group <- ppc_group(dat_MDSTYCV)
p_postpred <- plot_ppc_overlay(fit_MDSTYCV, dat_MDSTYCV, pp_group, xlim = c(NA, 2000)); p_postpred
save_gg("postpred_density", model_id_MDSTYCV, p_postpred)

### Contrast test statistics ----

(p_ppc <- plot_ppc_season_contrast_stats(fit_MDSTYCV, dat_MDSTYCV))
save_gg("postpred_stat", model_id_MDSTYCV, p_ppc)

## Year effect (fixed, not pooled) ---------------------------------------------
# Same as 4.3/5.3/6.3's own panel -- yr1/yr2 free, yr3 = -(yr1+yr2) by construction.

pf <- post_full(fit_MDSTYCV)
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
# construction (Liberty has the most combined observations).

cv_2 <- -(pf$cv_1$cv_1 + pf$cv_3$cv_3 + pf$cv_4$cv_4 + pf$cv_5$cv_5)

pc_cultivar <- bind_rows(
  tibble(group = idx$Cv$levels[1], value = pf$cv_1$cv_1),
  tibble(group = idx$Cv$levels[2], value = cv_2),
  tibble(group = idx$Cv$levels[3], value = pf$cv_3$cv_3),
  tibble(group = idx$Cv$levels[4], value = pf$cv_4$cv_4),
  tibble(group = idx$Cv$levels[5], value = pf$cv_5$cv_5)
) %>% mutate(statistic = "Cultivar effect (log scale)")

p_cultivar <- variance_component_panels(
  pc_cultivar, quant = c(0, 1), palette = idx$Cv$palette); ap_cultivar

## Tree effect: how much tree-to-tree spread is there in the real fit? --------
# sigma_tr's own posterior magnitude, next to sigma[Mg] for scale -- same
# diagnostic as 3.3_MDST_16S_fit.R's own "added value of Tree" section.
# tr[Tr] itself (~129 levels) isn't plotted individually -- too many for a
# readable panel; precis()/the trankplot above cover those if ever needed.

pc_sigma_tr <- bind_rows(
  compute_contrasts(pf, keep = "sigma", group_levels = idx$Mg$levels),
  tibble(statistic = "sigma_tr", group = "Population", value = pf$sigma_tr$sigma_tr)
) %>% mutate(statistic = factor(statistic, levels = c("sigma", "sigma_tr")))

p_sigma_tr <- variance_component_panels(
  pc_sigma_tr, quant = c(0, 1),
  palette = Management_palette, # already carries Population's colour
  sd_stats = c("sigma", "sigma_tr")); p_sigma_tr

## Year/Covariate/Cultivar/Tree effects, combined ---------------------------------

p_effects <- p_year / p_covariates / p_cultivar / p_sigma_tr
save_gg("fit_effects", model_id_MDSTYCV, p_effects, width = 8, height = 14)

## Full pairwise parameter check ------------------------------------------------
# Same rationale/settings as 7.2's own calibration-stage check -- the
# combination that mattered most there: sigma_tr vs cv_1..cv_4.

pairs_vars <- c("loga[1]", "loga[2]", "sigma[1]", "sigma[2]", "sigma_tr",
                 "yr1", "yr2", "cv_1", "cv_3", "cv_4", "cv_5",
                 "b_deg", "b_precip", "b_seq")
p_pairs <- plot_mcmc_pairs(fit_MDSTYCV,
                           variables = pairs_vars, n_keep = 1000)
save_gg("fit_mcmc_pairs", model_id_MDSTYCV, p_pairs, width = 15, height = 15, type = "png")


# Summary:

save_report("fit_summary", model_id_MDSTYCV, fit_MDSTYCV, model = model_MDSTYCV_16S)

