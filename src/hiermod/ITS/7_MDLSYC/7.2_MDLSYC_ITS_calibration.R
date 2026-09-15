# MODEL 7 (MDLSYC): three standardized control covariates (deg_h_z,
# precip_72h_z, seq_depth_z) added as additive fixed slopes, on top of
# Model 6. Parameter recovery, collinearity-aware recovery check,
# prior-predictive check, and SBC.

hiermod_marker <- "ITS"
source('src/hiermod/0_SETUP.R')
source('src/hiermod/Models/MDLSYC_model.R') # model_MDLSYC_ITS, means_MDLSYC(), variance_partition_MDLSYC(), sim_div_MDLSYC(), contrast_may_gap_MDLSYC(), simulate_from_priors_MDLSYC(), dq_MDLSYC
model <- model_MDLSYC_ITS

hiermod_out_dir <- "out/hiermod/ITS_7_lognormal_MDLSYC"

## Model specification ---------------------------------------------------------
# deg_h_z/precip_72h_z (degree-hours the day before / precipitation 3 days
# before) and seq_depth_z (log sequencing depth). Their coefficients could be
# entangled with s_conv/gap_shift/yr[Yr]*sigma_yr in the real-data posterior:
# Bad for confidence, but will not bias. SBC below tests recoverability
# against synthetic data where the three are independent.

## Parameter recovery -----------------------------------------------------------

may_conv <- 4
may_org <- 7
july_conv_shift <- 0.35
july_org_shift <- 0.2

true_sigma    <- cv_to_sigma(c(0.6, 0.25, 0.8, 0.8)) # conv_May, conv_July, org_May, org_July
true_sigma_yr <- 0.3
true_b_deg    <- 0.15
true_b_precip <- -0.15
true_b_seq    <- -0.4
true_cv       <- c(0.05, -0.05, 0.1, -0.1, 0)

dat_sim <- sim_div_MDLSYC(
  N_samples = 240,
  n_loc     = 4,
  loga      = log(c(may_conv, may_org)),
  s_conv    = july_conv_shift,
  gap_shift = july_org_shift,
  sigma     = true_sigma,
  sigma_yr  = true_sigma_yr,
  b_deg     = true_b_deg,
  b_precip  = true_b_precip,
  b_seq     = true_b_seq,
  cv        = true_cv,
  p_dropout = 0.1,
  shift     = 1
); head(dat_sim)

hist(dat_sim$Dv_shifted)

fit_sim <- ulam(
  model,
  data = as.list(dat_sim),
  chains = 6, cores = 6, iter = 10000,
  control = list(adapt_delta = 0.99)
)

precis(fit_sim, depth = 2)

### Fixed effects recovery -------------

post_sim <- extract.samples(fit_sim)

fixed_recovery <- check_recovery(
  true = c(loga1 = log(may_conv),
           loga2 = log(may_org),
           s_conv = july_conv_shift,
           gap_shift = july_org_shift,
           b_deg = true_b_deg,
           b_precip = true_b_precip,
           b_seq = true_b_seq),
  post_draws = list(
    loga1 = post_sim$loga[,1],
    loga2 = post_sim$loga[,2],
    s_conv = post_sim$s_conv,
    gap_shift = post_sim$gap_shift,
    b_deg = post_sim$b_deg,
    b_precip = post_sim$b_precip,
    b_seq = post_sim$b_seq)
); fixed_recovery

### Sigma / Year / Cultivar recovery -------------

sigma_recovery <- check_recovery(
  true = list(sigma1 = true_sigma[1], sigma2 = true_sigma[2],
              sigma3 = true_sigma[3], sigma4 = true_sigma[4],
              sigma_yr = true_sigma_yr),
  post_draws = list(sigma1 = post_sim$sigma[,1], sigma2 = post_sim$sigma[,2],
                    sigma3 = post_sim$sigma[,3], sigma4 = post_sim$sigma[,4],
                    sigma_yr = post_sim$sigma_yr)
); sigma_recovery

cv_recovery <- check_recovery(
  true = as.list(setNames(true_cv, paste0("cv", 1:5))),
  post_draws = setNames(lapply(1:5, function(i) post_sim$cv[,i]), paste0("cv", 1:5))
); cv_recovery

# all good

### Contrast recovery ------------------------------

cr <- contrast_recovery(fit_sim, means_MDLSYC, may_conv, may_org,
                         july_conv_shift, july_org_shift, shift = 1)

p_sim_contrast <- contrast_plot_panels(
  cr$estimands, quant = c(0, 0.995), scales = 'free_y',
  group_pal = Management_palette,
  true_vals = cr$true_estimands); p_sim_contrast

save_report("sim_summary", model_id, fit_sim, cr$estimands, model,
            recovery = bind_rows(fixed_recovery, sigma_recovery, cv_recovery),
            model_name = "The Weatherman")
save_gg("sim_contrast_density", model_id, p_sim_contrast)

# Extremely good recovery !

## Collinearity-aware parameter recovery --------------------------------------

# Full Claude suggestion! ---
# The real fit shows ~1% divergences and 1/6 chains with E-BFMI < 0.3. SBC
# above draws deg_h_z/precip_72h_z/seq_depth_z independently, so it only
# tests recoverability, not whether the real design's collinearity
# (deg_h ~ Season r=-0.73; Seq_depth targeted here at -0.5 with the
# structural diversity signal) affects the sampler. Calibration
# itself shouldn't be affected by collinearity in a correctly-specified
# model (right?); what's actually at stake is convergence (divergences/E-BFMI),
# which this checks directly across a handful of replicates, cheaper than
# a full SBC re-run. Escalate to a full SBC only if this shows real
# degradation vs. the plain recovery run above.

n_confound_reps <- 15
confound_chains <- 6
confound_iter   <- 10000

confound_diag <- map_dfr(seq_len(n_confound_reps), function(i){
  dat_confound <- sim_div_MDLSYC(
    N_samples = 240, n_loc = 4,
    loga = log(c(may_conv, may_org)), s_conv = july_conv_shift,
    gap_shift = july_org_shift, sigma = true_sigma, sigma_yr = true_sigma_yr,
    b_deg = true_b_deg, b_precip = true_b_precip, b_seq = true_b_seq,
    cv = true_cv, p_dropout = 0.1, shift = 1,
    # New parameters to force colinearity between variables:
    rho_deg_season = -0.73, rho_seq_mu = -0.5
  )

  fit_confound <- ulam(
    model, data = as.list(dat_confound),
    chains = confound_chains, cores = confound_chains, iter = confound_iter,
    control = list(adapt_delta = 0.99)
  )

  diag <- fit_confound@cstanfit$diagnostic_summary(diagnostics = c("divergences", "treedepth", "ebfmi"))
  tibble(
    rep             = i,
    n_divergent     = sum(diag$num_divergent),
    n_max_treedepth = sum(diag$num_max_treedepth),
    min_ebfmi       = min(diag$ebfmi),
    # attenuated below the -0.5 target by residual sigma[cell] noise on
    # top of the structural signal -- expected, see hiermod/Models/MDLSYC_model.R
    realized_cor_seq_Dv = cor(dat_confound$seq_depth_z, dat_confound$Dv)
  )
})

confound_diag
n_transitions <- n_confound_reps * confound_chains * (confound_iter %/% 2)
cat("Divergence rate:", round(100 * sum(confound_diag$n_divergent) / n_transitions, 3), "%\n")
cat("Replicates with any chain E-BFMI < 0.3:", sum(confound_diag$min_ebfmi < 0.3), "/", n_confound_reps, "\n")

# Compare divergence rate / E-BFMI here against the plain fit_sim run
# above -- if similar, the real fit's geometry issue likely isn't (mainly)
# this collinearity; if clearly worse, it's evidence to escalate to a full
# SBC re-run under this confounded simulator (simulate_from_priors_MDLSYC()
# already accepts rho_deg_season/rho_seq_mu for that).

## Prior predictive check -------------------------------------------------------

n_prior <- 1000
extracted_prior <- extract.prior(fit_sim, n = n_prior)

prior_pred <- map_dfr(seq_len(n_prior), function(i){
  simulate_from_priors_MDLSYC(draw_true(extracted_prior, i), shift = 1)
}, .id = "draw")

summary(prior_pred$Dv_shifted)

p_prior_pc <- prior_predictive_spaghetti(
  prior_pred, value_col = "Dv_shifted", upper_q = 0.99, model = model,
  title = "Prior predictive check", observed = dat_sim$Dv_shifted) ; p_prior_pc

save_gg("sim_prior_PC", model_id, p_prior_pc)

## Simulation-based calibration (SBC), via the SBC package --------------------
# No separate variance-budget-calibration pass here: `model`'s sigma priors
# are already Model 6's calibrated rates (inherited from model_MDLSY_ITS),
# since K stays at 4 (b_deg/b_precip/b_seq are additive mu-level fixed
# effects, not new summed variance terms). dq_MDLSYC (MDLSYC_model.R) now
# checks all three estimands (may_gap/july_gap/seasonal_change). This SBC
# draws the three covariates independently (rho_deg_season = rho_seq_mu = 0
# defaults) -- tests recoverability, not calibration under the real
# design's collinearity (see the Collinearity-aware parameter recovery
# section above, and the methodology checklist).

sbc_gen_MDLSYC <- make_sbc_generator(
  fit = fit_sim, simulate_fn = simulate_from_priors_MDLSYC,
  keep = c("loga", "s_conv", "gap_shift", "sigma", "sigma_loc", "sigma_tr", "sigma_yr", "cv",
           "b_deg", "b_precip", "b_seq"),
  gen_cols = c("Dv", "Mg", "Lo", "Mo", "Yr", "Tr", "Cv", "deg_h_z", "precip_72h_z", "seq_depth_z"),
  extra_globals = "sim_div_MDLSYC", shift = 1)

n_sbc  <- 100
n_iter <- 20000

sbc_MDLSYC <- run_sbc_pipeline(
  generator = sbc_gen_MDLSYC$generator, globals = sbc_gen_MDLSYC$globals,
  n_sbc = n_sbc, model = model, model_id = model_id, n_iter = n_iter,
  hiermod_out_dir = hiermod_out_dir, dquants = dq_MDLSYC,
  control = list(adapt_delta = 0.99))

plot_sbc_diagnostics(sbc_MDLSYC, model_id, n_sbc)

save_sbc_health_report(model_id, sbc_MDLSYC, n_sbc, n_iter,
                        variables = c("loga[1]", "loga[2]", "sigma[1]", "sigma[2]", "sigma[3]", "sigma[4]",
                                      "sigma_loc", "sigma_tr", "sigma_yr", "b_deg", "b_precip", "b_seq",
                                      "may_gap", "july_gap", "seasonal_change"),
                        hiermod_out_dir = hiermod_out_dir)
