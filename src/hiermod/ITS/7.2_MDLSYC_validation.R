# Model validation 

source('src/hiermod/ITS/0_SETUP.R')
source('src/hiermod/ITS/7.1_MDLSYC_model.R') # model, means_MDLSYC(), variance_partition_MDLSYC(), sim_div_MDLSYC(), contrast_may_gap_MDLSYC(), simulate_from_priors()
hiermod_out_dir <- "out/hiermod/ITS_7_lognormal_MDLSYC"

## Model specification ---------------------------------------------------------
# Same structure as MODEL 6, plus three standardized control covariates:
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

save_report("sim_summary", "MDLSYC", fit_sim, cr$estimands, model,
            recovery = bind_rows(fixed_recovery, sigma_recovery, cv_recovery),
            model_name = "The Weatherman")
save_gg("sim_contrast_density", "MDLSYC", p_sim_contrast)

# Extremely good recovery !

## Prior predictive check -------------------------------------------------------

n_prior <- 1000
extracted_prior <- extract.prior(fit_sim, n = n_prior)

prior_pred <- map_dfr(seq_len(n_prior), function(i){
  simulate_from_priors(draw_true(extracted_prior, i), shift = 1)
}, .id = "draw")

summary(prior_pred$Dv_shifted)

p_prior_pc <- prior_predictive_spaghetti(
  prior_pred, value_col = "Dv_shifted", upper_q = 0.99, model = model,
  title = "Prior predictive check", observed = dat_sim$Dv_shifted) ; p_prior_pc

save_gg("sim_prior_PC", "MDLSYC", p_prior_pc)

## Simulation-based calibration (SBC) --------------------------------------------
# No separate variance-budget-calibration pass here: `model`'s sigma priors
# are already Model 6's calibrated rates (hardcoded in 7.1_MDLSYC_model.R),
# since K stays at 4 (b_deg/b_precip/b_seq are additive mu-level fixed
# effects, not new summed variance terms).

ncores <- 18
nchains <- 2
n_sbc = 100
n_iter = 20000

sbc_MDLSYC <- run_sbc(
  model_fit   = fit_sim,
  means_fn    = means_MDLSYC,
  contrast_fn = contrast_may_gap_MDLSYC,
  simulate_fn = function(true_params) simulate_from_priors(true_params, shift = 1),
  n_sbc = n_sbc, iter = n_iter,
  n_parallel = ncores/nchains, chains = nchains, cores = nchains,
  control = list(adapt_delta = 0.99))

sbc_out_MDLSYC <- save_sbc_report(sbc_MDLSYC, paste0("MDLSYC_",n_sbc,"iter"))

## Model fit ----------------------------------------------------------------

dat <- list(
  Dv = div$Hill_1 - 1, # subtract the floor because Dv must be (0, Inf)-support to match the likelihood
  Mg = idx$Mg$to_index(div$Management),
  Lo = idx$Lo$to_index(div$Location),
  Mo = idx$Mo$to_index(div$Time),
  Yr = idx$Yr$to_index(div$Year),
  Tr = idx$Tr$to_index(div$Tree_id),
  Cv = idx$Cv$to_index(div$Cultivar),
  deg_h_z = div$deg_h_z,
  precip_72h_z = div$precip_72h_z,
  seq_depth_z = div$seq_depth_z
)
dat$cell <- (dat$Mg - 1) * 2 + dat$Mo

fit_MDLSYC <- ulam(
  model,
  data = dat,
  chains = 6, cores = 6, iter = 20000,
  control = list(adapt_delta = 0.99)
)
save_fit("fit", "MDLSYC", fit_MDLSYC)

precis(fit_MDLSYC, depth = 2) 

save_pdf("fit_traceplot", "MDLSYC", function() traceplot(fit_MDLSYC, n_cols = 6, max_rows = 10))
save_pdf("fit_trankplot", "MDLSYC", function() trankplot(fit_MDLSYC, n_cols = 6, max_rows = 10))

# Real-data collinearity check: correlated draws are the expected symptom
# of the correlation noted above, not a red flag on their own.
post_MDLSYC <- extract.samples(fit_MDLSYC)
cor(post_MDLSYC$b_deg, post_MDLSYC$s_conv)
cor(post_MDLSYC$b_deg, post_MDLSYC$gap_shift)
cor(post_MDLSYC$b_precip, post_MDLSYC$s_conv)
cor(post_MDLSYC$b_deg, post_MDLSYC$yr[,1])
cor(post_MDLSYC$b_seq, post_MDLSYC$s_conv)
cor(post_MDLSYC$b_seq, post_MDLSYC$yr[,1])

## Posterior predictive check --------------------------------------------------

### Overall, by Management x Season cell ----

pp_group <- interaction(idx$Mg$to_label(dat$Mg), idx$Mo$to_label(dat$Mo), sep = " ")
(p_postpred <- plot_ppc_overlay(fit_MDLSYC, dat, pp_group, xlim = c(0,150)))
save_gg("postpred_density", "MDLSYC", p_postpred)

### Contrast test statistics ----

(p_ppc <- plot_ppc_season_contrast_stats(fit_MDLSYC, dat))
save_gg("postpred_stat", "MDLSYC", p_ppc)
