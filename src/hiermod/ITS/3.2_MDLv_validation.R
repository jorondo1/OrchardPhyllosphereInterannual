# 3.2_MDLv_validation.R -- MODEL 3 (Management effect by Location): parameter
# recovery, prior-predictive check, SBC, the real fit, and PPC.

source('~/Repos/orchardPhyllosphere2/src/hiermod/ITS/0_SETUP.R')
source('~/Repos/orchardPhyllosphere2/src/hiermod/ITS/3.1_MDLv_model.R') # model, means_MDLv(), mdlv_labels, sim_div_MLvary(), simulate_from_priors()
hiermod_out_dir <- "out/hiermod/ITS_3_lognormal_MDLv"

## Model specification ---------------------------------------------------------
# Organic-Conventional gap now varies by Location: g[Lo]*sigma_g, same
# non-centered trick as b[Lo]*sigma_loc, applied via (Mg-1) so it only ever
# shifts the Organic prediction. Why this instead of another intercept
# layer: see MODEL_HISTORY.md. loga[2]-loga[1] stays the population-average
# gap; g[Lo] is a Location's own deviation from it, shrunk toward 0.
#
# Caveat: Saint-Benoit/Windsor are Conventional-only/Organic-only, so their
# g[Lo] is entangled with b[Lo] (genuinely unknown, not "gap is zero") --
# Milton/Compton do essentially all the work of identifying sigma_g.

## Parameter recovery -----------------------------------------------------------
# sigma_g = 0.5 on purpose (not 0) -- the point is to check the model can
# find a location-varying gap when one genuinely exists in the data.

# simulate
dat_sim <- sim_div_MLvary(
  Mg = rbern(250)+1,
  Lo = sample(1:5, 250, replace = TRUE),
  loga = log(c(8, 12)), sigma = cv_to_sigma(c(0.5, 0.8)),
  sigma_loc = 0.5,
  sigma_g = 0.5)

# fit
fit_sim <- ulam(
  model,
  data = as.list(dat_sim),
  chains = 4, cores = 4, iter = 4000,
  control = list(adapt_delta = 0.99)
)

#diagnose
precis(fit_sim, depth = 2)

traceplot(fit_sim); trankplot(fit_sim)
save_pdf("sim_traceplot", "MDLv", function() traceplot(fit_sim))
save_pdf("sim_trankplot", "MDLv", function() trankplot(fit_sim))

# draw posterior samples
pf_sim <- post_full(fit_sim, means_MDLv)        # raw draws + mean/median

pc_full_sim <- compute_contrasts(
  pf_sim, keep = NULL, labels = mdlv_labels,
  group_levels = idx$Mg$levels)                 # -> long, panel-ready

# Plot main outcome contrasts
p_sim_contrast <- pc_full_sim %>%
  filter(!statistic %in% c('loga', 'b', 'g')) %>%
  contrast_plot_panels(
    quant = c(0,0.995), group_pal = Management_palette); p_sim_contrast

save_report("sim_summary", "MDLv", fit_sim, pc_full_sim, model, model_name = "The Copycat")
save_gg("sim_contrast_density", "MDLv", p_sim_contrast)

## Prior predictive check -------------------------------------------------------

n_prior <- 1000
prior <- extract.prior(fit_sim, n = n_prior)

# simulate_from_priors() is sourced from 3.1_MDLv_model.R above.
prior_pred <- map_dfr(seq_len(n_prior), function(i){
  simulate_from_priors(draw_true(prior, i))
}, .id = "draw")

summary(prior_pred$Dv) # similar to earlier looking at 3rd quant

p_prior_pc <- prior_predictive_spaghetti(
  prior_pred, upper_q = 0.99, model = model,
  title = "Prior predictive check"); p_prior_pc

save_gg("sim_prior_PC", "MDLv", p_prior_pc)
# if this still looks implausible for Hill_1, tighten sigma_g further before
# moving on -- same read as MDLb's own prior predictive check

## Simulation-based calibration (SBC) --------------------------------------------
# Same as MDLb's -- now also checks sigma_g recovery.

sbc_MDLv <- run_sbc(
  model_fit   = fit_sim,
  means_fn    = means_MDLv,
  simulate_fn = simulate_from_priors,
  n_sbc = 100, iter = 4000, chains = 4, control = list(adapt_delta = 0.99))

sbc_out_MDLv <- summarize_sbc(sbc_MDLv)
save_sbc_report(sbc_out_MDLv, "MDLv")
hist(sbc_out_MDLv$ranks)

## Model fit ----------------------------------------------------------------

dat <- list(
  Dv = div$Hill_1,
  Mg = idx$Mg$to_index(div$Management),
  Lo = idx$Lo$to_index(div$Location)
)

fitb <- ulam(
  model,
  data = dat,
  chains = 6, cores = 6, iter = 10000,
  control = list(adapt_delta = 0.99)
)
save_fit("fit", "MDLv", fitb)

precis(fitb, depth = 2)

traceplot(fitb); trankplot(fitb)
save_pdf("fit_traceplot", "MDLv", function() traceplot(fitb))
save_pdf("fit_trankplot", "MDLv", function() trankplot(fitb))

## Posterior predictive check --------------------------------------------------

pb_postpred <- plot_ppc_overlay(fitb, dat, idx$Mg$to_label(dat$Mg), xlim = c(0,150)); pb_postpred

save_gg("postpred_density", "MDLv", pb_postpred)

### By Location ---------------------------

postpred_long <- sim(fitb, dat, n = 1000) %>% postpred_as_long_tibble()

Loc_labels_n <- tibble(
  Lo   = idx$Lo$to_label(seq_along(idx$Lo$levels)),
  Lo_n = idx$Lo$to_label_n(seq_along(idx$Lo$levels))
)

pp_joined <- tibble(
  obs = seq_len(length(dat$Mg)),
  Mg  = idx$Mg$to_label(dat$Mg),
  Lo  = idx$Lo$to_label(dat$Lo)) %>%
  left_join(postpred_long, by = "obs") %>%
  left_join(Loc_labels_n, by = "Lo")

obs_df <- tibble(
  Dv = dat$Dv,
  Mg = idx$Mg$to_label(dat$Mg),
  Lo = idx$Lo$to_label(dat$Lo)
) %>% left_join(Loc_labels_n, by = "Lo")

xlim_upper_pp <- quantile(pp_joined$Dv_sim, 0.98)

p_postpred_ridges <- pp_joined %>%
  ggplot(aes(x = Dv_sim, y = Lo_n, height = after_stat(density), fill = Mg)) +
  ggridges::geom_density_ridges(
    stat = "density", alpha = 0.55, colour = "white", scale = 0.8, linewidth = 0.3) +
  geom_point(
    data = obs_df, inherit.aes = FALSE,
    shape = 21, stroke = 0.2, size = 2, alpha = 0.6,
    aes(x = Dv, y = as.numeric(Lo_n) - 0.15, fill = Mg),
    position = position_jitter(height = 0.1, width = 0)) +
  scale_fill_manual(values = Management_palette) +
  scale_colour_manual(values = Management_palette) +
  coord_cartesian(xlim = c(0, xlim_upper_pp)) +
  theme_light() +
  guides(colour = 'none') +
  labs(x = "Diversity", y = "Location", fill = "Management",
       caption = "Ridges = posterior predictive density; points = observed data (jittered)"); p_postpred_ridges

save_gg("postpred_ridges", "MDLv", p_postpred_ridges)

## Contrast statistic ------------------------------
# Same idea as MDLb's: does the model reproduce the specific gap we're
# reporting, not just plausible marginal shapes per group?

(p_ppc_median_contrast <- plot_ppc_contrast_stat(fitb, dat, dat$Mg, median, "median", idx$Mg$levels))
(p_ppc_mad_contrast <- plot_ppc_contrast_stat(fitb, dat, dat$Mg, mad, "dispersion (MAD)", idx$Mg$levels))

save_gg("postpred_stat_median_contrast", "MDLv", p_ppc_median_contrast)
save_gg("postpred_stat_mad_contrast", "MDLv", p_ppc_mad_contrast)

# Next step: Tree ID + Season/Year. Deliberately not bundled in here -- see
# write-up. All Dv rows currently get treated as independent even though many
# trees contribute two rows (May + July) within a year; that's the next
# thing to model properly, on its own, not tacked onto this one.
