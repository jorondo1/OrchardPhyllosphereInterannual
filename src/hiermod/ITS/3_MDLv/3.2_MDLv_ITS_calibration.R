# MODEL 3 (MDLv): Organic-Conventional gap now varies by Location
# (g[Lo]*sigma_g) on top of Model 2's Location-varying intercept.
# Parameter recovery, prior-predictive check, and SBC.

hiermod_marker <- "ITS"
source('src/hiermod/0_SETUP.R')
source('src/hiermod/Models/MDLv_model.R') # model_MDLv_ITS, means_MDLv(), mdlv_labels, sim_div_MDLv(), simulate_from_priors_MDLv(), dq_MDLv
model <- model_MDLv_ITS

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
dat_sim <- sim_div_MDLv(
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
save_pdf("sim_traceplot", model_id, function() traceplot(fit_sim))
save_pdf("sim_trankplot", model_id, function() trankplot(fit_sim))

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

save_report("sim_summary", model_id, fit_sim, pc_full_sim, model, model_name = "The Copycat")
save_gg("sim_contrast_density", model_id, p_sim_contrast)

## Prior predictive check -------------------------------------------------------

n_prior <- 1000
prior <- extract.prior(fit_sim, n = n_prior)

# simulate_from_priors_MDLv() is sourced from hiermod/Models/MDLv_model.R above.
prior_pred <- map_dfr(seq_len(n_prior), function(i){
  simulate_from_priors_MDLv(draw_true(prior, i))
}, .id = "draw")

summary(prior_pred$Dv) # similar to earlier looking at 3rd quant

p_prior_pc <- prior_predictive_spaghetti(
  prior_pred, upper_q = 0.99, model = model,
  title = "Prior predictive check"); p_prior_pc

save_gg("sim_prior_PC", model_id, p_prior_pc)
# if this still looks implausible for Hill_1, tighten sigma_g further before
# moving on -- same read as MDLb's own prior predictive check

## Simulation-based calibration (SBC), via the SBC package --------------------
# dq_MDLv (MDLv_model.R) correctly includes sigma_g on the Organic side only
# (matching means_MDLv()'s own total_var) -- the old run_sbc() call here had
# no contrast_fn either, so it silently fell back to contrast_from_means()'s
# default (drops sigma_loc AND sigma_g from the true side entirely).

sbc_gen_MDLv <- make_sbc_generator(
  fit = fit_sim, simulate_fn = simulate_from_priors_MDLv,
  keep = c("loga", "sigma", "sigma_loc", "sigma_g"), gen_cols = c("Dv", "Mg", "Lo"),
  extra_globals = "sim_div_MDLv")

n_sbc  <- 100
n_iter <- 4000

sbc_MDLv <- run_sbc_pipeline(
  generator = sbc_gen_MDLv$generator, globals = sbc_gen_MDLv$globals,
  n_sbc = n_sbc, model = model, model_id = model_id, n_iter = n_iter,
  hiermod_out_dir = hiermod_out_dir, dquants = dq_MDLv,
  control = list(adapt_delta = 0.99))

plot_sbc_diagnostics(sbc_MDLv, model_id, n_sbc)

sbc_MDLv$stats |>
  dplyr::filter(variable %in% c("loga[1]", "loga[2]")) |>
  dplyr::group_by(variable) |>
  dplyr::summarise(mean_rank_frac = mean(rank / max_rank), median_rank_frac = median(rank / max_rank))

save_sbc_health_report(model_id, sbc_MDLv, n_sbc, n_iter,
                        variables = c("loga[1]", "loga[2]", "sigma[1]", "sigma[2]", "sigma_loc", "sigma_g",
                                      "median_contrast", "mean_contrast"),
                        hiermod_out_dir = hiermod_out_dir)
