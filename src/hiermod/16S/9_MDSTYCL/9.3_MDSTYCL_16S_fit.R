# MODEL 9 (MDSTYCL, "Faramir the Judicious"), 16S: real fit (Hill_1 - 1) on the B/D subset, PPC, effect panels

hiermod_marker <- "16S"
source('src/hiermod/0_SETUP.R')
source('src/hiermod/Models/MDSTYCL_model.R') # model_MDSTYCL_16S, means_MDSTYCL()

hiermod_out_dir <- "out/hiermod/16S_9_location_MDSTYCL"

## Data subset: locations B and D (both managements) ---------------------------
div %<>% filter(Location %in% c('B', 'D'))

# Tree index re-compacted to 1..n within the subset
# - idx$Tr uses the global 129-tree index; ulam sizes tr[] by distinct values (~75)
idx_Tr_local <- make_index(div$Tree_id)

## Model fit ----------------------------------------------------------------

dat_MDSTYCL <- list(
  Dv = div$Hill_1 - 1, # subtract the floor because Dv must be (0, Inf)-support to match the likelihood
  Mg = idx$Mg$to_index(div$Management),
  Mo = idx$Mo$to_index(div$Time),
  Yr = idx$Yr$to_index(div$Year),
  Lo = ifelse(idx$Lo$to_index(div$Location) == 2, 1L, 2L), # B=1, D=2
  Tr = idx_Tr_local$to_index(div$Tree_id),
  deg_h_z = div$deg_h_z,
  precip_72h_z = div$precip_72h_z,
  seq_depth_z = div$seq_depth_z
)

fit_MDSTYCL <- ulam(
  model_MDSTYCL_16S,
  data = dat_MDSTYCL,
  chains = 6, cores = 6, iter = 10000,
  control = list(adapt_delta = 0.99)
)

save_fit("fit", model_id_MDSTYCL, fit_MDSTYCL)
saveRDS(dat_MDSTYCL, file.path(hiermod_out_dir, "dat_MDSTYCL.rds"))

precis(fit_MDSTYCL, depth = 2)

save_pdf("fit_trankplot", model_id_MDSTYCL,
         function() trankplot(fit_MDSTYCL, n_cols = 6, max_rows = 30),
         width = 24, height = 36)

## Posterior predictive check --------------------------------------------------

### Overall, by Management x Season cell ----

pp_group <- ppc_group(dat_MDSTYCL)
p_postpred <- plot_ppc_overlay(fit_MDSTYCL, dat_MDSTYCL, pp_group, xlim = c(NA, 2000))
save_gg("postpred_density", model_id_MDSTYCL, p_postpred)

### Contrast test statistics ----

(p_ppc <- plot_ppc_season_contrast_stats(fit_MDSTYCL, dat_MDSTYCL))
save_gg("postpred_stat", model_id_MDSTYCL, p_ppc)

## Year effect (fixed, not pooled) ---------------------------------------------
# yr1/yr2 free; yr3 = -(yr1 + yr2)

pf <- post_full(fit_MDSTYCL)
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

## Location effect (fixed, sum-to-zero) ------------------------------------------
# lo1 = location B's offset; D = -lo1

pc_location <-  tibble(
  group = "B", value = pf$lo1$lo1) %>% mutate(statistic = "Location effect (log scale)")

p_location <- variance_component_panels(
  pc_location, quant = c(0, 1), palette = fill_loc[c("B", "D")]); p_location

## Tree effect: how much tree-to-tree spread is there in the real fit? --------
# sigma_tr next to sigma[Mg] (as in 3.3)

pc_sigma_tr <- bind_rows(
  compute_contrasts(pf, keep = "sigma", group_levels = idx$Mg$levels),
  tibble(statistic = "sigma_tr", group = "Population", value = pf$sigma_tr$sigma_tr)
) %>% mutate(statistic = factor(statistic, levels = c("sigma", "sigma_tr")))

p_sigma_tr <- variance_component_panels(
  pc_sigma_tr, quant = c(0, 1),
  palette = Management_palette, # already carries Population's colour
  sd_stats = c("sigma", "sigma_tr")); p_sigma_tr

## Year/Covariate/Location/Tree effects, combined ---------------------------------

p_effects <- p_year / p_covariates / p_location / p_sigma_tr
save_gg("fit_effects", model_id_MDSTYCL, p_effects, width = 8, height = 14)

## Full pairwise parameter check ------------------------------------------------
# As in 9.2 (key pair: sigma_tr vs lo1)

pairs_vars <- c("loga[1]", "loga[2]", "sigma[1]", "sigma[2]", "sigma_tr",
                 "yr1", "yr2", "lo1",
                 "b_deg", "b_precip", "b_seq")
p_pairs <- plot_mcmc_pairs(fit_MDSTYCL,
                           variables = pairs_vars, n_keep = 1000)
save_gg("fit_mcmc_pairs", model_id_MDSTYCL, p_pairs, width = 15, height = 15, type = "png")


# Summary:

save_report("fit_summary", model_id_MDSTYCL, fit_MDSTYCL, model = model_MDSTYCL_16S)

