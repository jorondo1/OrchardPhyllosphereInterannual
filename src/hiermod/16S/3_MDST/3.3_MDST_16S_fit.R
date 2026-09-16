# MODEL 3 (MDST, "Treebeard the Skeptic"), 16S, SHIFTED (Hill_1 - 1): real
# fit and PPC. model_MDST_16S already carries its own validated priors (see
# MDST_model.R), no local override needed here.

hiermod_marker <- "16S"
source('src/hiermod/0_SETUP.R')
source('src/hiermod/Models/MDST_model.R') # model_MDST_16S, means_MDST()

hiermod_out_dir <- "out/hiermod/16S_3_tree_MDST"

model_MDST_16S_tight <- model_MDST_16S
#model_MDST_16S_tight$pr_sigma_tr <- quote(sigma_tr ~ dlnorm(log(0.15), 0.5))

## Model fit ----------------------------------------------------------------

dat_MDST <- list(
  Dv = div$Hill_1 - 1, # subtract the floor because Dv must be (0, Inf)-support to match the likelihood
  Mg = idx$Mg$to_index(div$Management),
  Mo = idx$Mo$to_index(div$Time),
  Tr = idx$Tr$to_index(div$Tree_id)
)

fit_MDST <- ulam(
  model_MDST_16S,
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
p_postpred <- plot_ppc_overlay(fit_MDST, dat_MDST, pp_group, xlim = c(NA, 2000))
save_gg("postpred_density", model_id_MDST, p_postpred)

### Contrast test statistics ----

(p_ppc <- plot_ppc_season_contrast_stats(fit_MDST, dat_MDST))
save_gg("postpred_stat", model_id_MDST, p_ppc)

## Added value of Tree: how much tree-to-tree spread is there? -----------------
# The direct parameter-level answer to "is Tree worth having" -- sigma_tr's
# own posterior magnitude, next to sigma[Mg] for scale (a Tree effect that's
# small relative to sigma[Mg] earns its keep mainly via the repeated-measures
# correction to uncertainty, not by explaining much variance on its own; a
# large one means real, previously-unmodeled tree-level structure).
# tr[Tr] has 100+ levels -- too many for a readable ridge/density-by-group
# panel, so only its population SD (sigma_tr) gets one; individual tree
# offsets are better read off precis()/the trankplot above if ever needed.

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
