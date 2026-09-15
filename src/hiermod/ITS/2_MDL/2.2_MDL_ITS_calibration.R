# MODEL 2 (MDL): adds partial pooling across Location (b[Lo]*sigma_loc,
# non-centered) on top of Model 2/MDv's Management-specific variance.
# Calibration: parameter recovery, prior-predictive check, and SBC.

hiermod_marker <- "ITS"
source('src/hiermod/0_SETUP.R')
source('src/hiermod/Models/MDL_model.R') # model_MDL_ITS, means_MDL(), sim_div_MDL(), simulate_from_priors_MDL(), dq_MDL
model <- model_MDL_ITS

hiermod_out_dir <- "out/hiermod/ITS_2_lognormal_MDL"

# MODEL 2 -- Partial pooling across Location (non-centered) =================

# Keeps Management-specific variance from Model 1's MDv.
# Exploratory: priors here turn out to be too loose (see Prior predictive
# check below) -- MODEL 2B is the version actually fit to real data.

## Parameter recovery -----------------------------------------------------------
# Funnel check via pairs(): if divergences cluster where sigma_loc is small,
# that's the classic non-centered-parameterization funnel. Also serves as the
# compiled `ulam` object extract.prior() needs below, in the Prior predictive
# check -- extract.prior() can't run without a compiled fit of the right
# structure, so this has to come first even though conceptually you'd want to
# sanity-check the priors before trusting any fit's recovery.

dat_sim <- sim_div_MDL(
  Mg = rbern(250)+1,
  # 4 locations for simulation; fewer creates a funnel (see below)
  Lo = sample(1:4, 250, replace = TRUE),
  loga = log(c(8, 12)),           # Difference of 4 in mean, raw scale
  sigma = cv_to_sigma(c(0.5, 0.8)))  # keep different variances

fit_sim <- ulam(
  model,
  data = as.list(dat_sim),
  chains = 6, cores = 6, iter = 10000 ) # divergences!

pf_sim <- post_full(fit_sim, means_MDL)
pc_sim <- compute_contrasts(pf_sim, keep = "mean", group_levels = idx$Mg$levels)

# Mean is often apart from median (iterate to see range)

precis(fit_sim, depth = 2 )

p_sim_contrast <- contrast_plot_panels(pc_sim, quant = c(0, 0.999), group_pal = Management_palette); p_sim_contrast

save_report("sim_summary", "MDL", fit_sim, pc_sim, model, model_name = "The Wildcard")
save_gg("sim_contrast_density", "MDL", p_sim_contrast)

# pairs plot: checkf
pairs_plot <- pairs(fit_sim, pars = c("sigma_loc","b[1]","b[2]"))
save_pdf("fit_pairs", "MDL", function()
  pairs(fit_sim, pars = c("sigma_loc","b[1]","b[2]")))
## Prior predictive check -------------------------------------------------------
# How extreme can Dv get under these priors?

# E[Dv] = exp(mu + sigma^2/2) means any prior tail draw of sigma[Mg]
# or sigma_loc amplifies Dv exponentially. Some large sigma draws are
# expected under dexp(1), P(sigma_loc>3) ~ 5%. If this tail looks
# implausible for Hill_1, tighten sigma[Mg]/sigma_loc, e.g. dexp(2)
# or dexp(3)).

n_prior <- 1000
prior <- extract.prior(fit_sim, n = n_prior)

# draw_true() is shared (sbc_workflow.R) -- generic slice of any extract.prior()
# output, no per-model rewrite needed. simulate_from_priors_MDL() is sourced
# from hiermod/Models/MDL_model.R above.

# One simulated dataset per prior draw
prior_pred <- map_dfr(seq_len(n_prior), function(i){
  simulate_from_priors_MDL(draw_true(prior, i))
}, .id = "draw")

summary(prior_pred$Dv)  # raw per-observation draws
# heavy max regardless of prior
# signal of interest: 3rd quantile

p_sim_spaghetti <- prior_predictive_spaghetti(
  prior_pred, upper_q = 0.99, model = model,
  title = "Prior predictive check -- model (dexp(1))"); p_sim_spaghetti

save_gg("sim_prior_PC", "MDL", p_sim_spaghetti)

# Which prior is driving the extreme tail?
diagnose_extreme_tail(
  prior_pred, value_col = "Dv",
  candidates = list(
    sigma_loc = prior$sigma_loc,
    sigma_max = pmax(prior$sigma[,1], prior$sigma[,2])))

# both sigma separate extreme from non-extreme sharply (~1 -> ~2.5);
# we don't check loga because it's in the dnorm space
# E[Dv] = exp(mu + σ²/2) with unbounded priors will always have some values
# fall in the far range of the tail. We can slightly tighten the sigma priors.

# Try tightening sigma priors
model_ppc1 <- model
model_ppc1$prior_sigma = sigma[Mg] ~ dexp(3)
model_ppc1$prior_sigma_loc = sigma_loc ~ dexp(2)

# Minimal, throwaway fit, not a real parameter-recovery fit
fit_ppc1 <- ulam(
  model_ppc1,
  data = as.list(dat_sim),
  chains = 1, cores = 1, iter = 100, messages = FALSE)

# Extract the priors
prior_ppc1 <- extract.prior(fit_ppc1, n = n_prior)

# simulate over each of the n_prior
prior_pred_ppc1 <- map_dfr(seq_len(n_prior), function(i){
  simulate_from_priors_MDL(draw_true(prior_ppc1, i))
}, .id = "draw")

summary(prior_pred_ppc1$Dv)


# Check distribution of parameters between extreme and non-extreme draws
diagnose_extreme_tail(
  prior_pred_ppc1, value_col = "Dv",
  candidates = list(
    sigma_loc = prior_ppc1$sigma_loc,
    sigma_max = pmax(prior_ppc1$sigma[,1], prior_ppc1$sigma[,2])))
# much more balanced

p_sim_spaghetti_ppc1 <- prior_predictive_spaghetti(
  prior_pred_ppc1, upper_q = 0.99,
  model = model_ppc1,
  title = paste("Prior predictive check: sigma tightened to",
                deparse1(model_ppc1$prior_sigma[[3]]))) +
  labs(x = 'Mean Hill diversity'); p_sim_spaghetti_ppc1

save_gg("sim_prior_PC_ppc1", "MDL", p_sim_spaghetti_ppc1)
# More reasonable!
# sd is 4 orders-of-magnitude lower,
# 3rd quartile 25% lower
# mean is nearly 3 OOM lower

## Parameter recovery: updated priors -----------------------------------------------------------

fitb_sim <- ulam(
  model_ppc1,  # <- <- <- Updated model
  data = as.list(dat_sim),
  chains = 4, cores = 4, iter = 4000,
  control = list(adapt_delta = 0.99) ) # smaller step size

precis(fitb_sim, depth = 2 ) # good r_hats

pf_b_sim <- post_full(fitb_sim, means_MDL)
pc_b_sim <- compute_contrasts(pf_b_sim, keep = "mean", group_levels = idx$Mg$levels)
# Iterate, mean sometimes way over target
mean(pf_b_sim$mean$mean_2 - pf_b_sim$mean$mean_1 > 100) # should be very small (~0)

# if divergences cluster where sigma_loc is small, funnel:
pairs_plot <- pairs(fitb_sim, pars = c("sigma_loc","b[1]","b[2]"))
save_pdf("fit_pairs", model_id, function()
  pairs(fitb_sim, pars = c("sigma_loc","b[1]","b[2]")))
# They do, not sure if ok

# Traces look fine:
traceplot(fitb_sim)
trankplot(fitb_sim)
save_pdf("sim_traceplot", model_id, function() traceplot(fitb_sim))
save_pdf("sim_trankplot", model_id, function() trankplot(fitb_sim))

save_report("sim_summary", model_id, fitb_sim, pc_b_sim, model_ppc1, model_name = "The Tamed Wildcard")

pb_sim_contrast <- contrast_plot_panels(pc_b_sim, quant = c(0, 0.995), group_pal = Management_palette) +
  labs(x = 'Estimates for mean Hill number of order 1 and its contrast'); pb_sim_contrast
# Not incredible, but much better (depends on iteration/unstable; let's SBC!)

save_gg("sim_contrast_density", model_id, pb_sim_contrast)

## Simulation-based calibration (SBC), via the SBC package --------------------
# dq_MDL (MDL_model.R) correctly includes sigma_loc on the true-side
# variance -- fixes a gap in the old run_sbc() call here: it had no
# contrast_fn, so it fell back to contrast_from_means()'s default, which
# silently omits sigma_loc from the true side. SBC never actually tested
# calibration of the full model before this migration, just a simplified
# slice of it.

sbc_gen_MDL <- make_sbc_generator(
  fit = fitb_sim, simulate_fn = simulate_from_priors_MDL,
  keep = c("loga", "sigma", "sigma_loc"), gen_cols = c("Dv", "Mg", "Lo"),
  extra_globals = "sim_div_MDL")

n_sbc  <- 100
n_iter <- 10000

sbc_MDL <- run_sbc_pipeline(
  generator = sbc_gen_MDL$generator, globals = sbc_gen_MDL$globals,
  n_sbc = n_sbc, model = model_ppc1, model_id = model_id, n_iter = n_iter,
  hiermod_out_dir = hiermod_out_dir, dquants = dq_MDL,
  control = list(adapt_delta = 0.99))

plot_sbc_diagnostics(sbc_MDL, model_id, n_sbc)

sbc_MDL$stats |>
  dplyr::filter(variable %in% c("loga[1]", "loga[2]")) |>
  dplyr::group_by(variable) |>
  dplyr::summarise(mean_rank_frac = mean(rank / max_rank), median_rank_frac = median(rank / max_rank))

save_sbc_health_report(model_id, sbc_MDL, n_sbc, n_iter,
                        variables = c("loga[1]", "loga[2]", "sigma[1]", "sigma[2]", "sigma_loc",
                                      "median_contrast", "mean_contrast"),
                        hiermod_out_dir = hiermod_out_dir)

saveRDS(model_ppc1, file.path(hiermod_out_dir, "model_ppc1.rds"))
