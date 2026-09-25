# MODEL 5 (MDSYC, "Radagast the Grower"), ITS: MDSYz plus three
# standardized control covariates (deg_h_z, precip_72h_z, seq_depth_z) as
# additive fixed slopes. Mirrors 5.2_MDSYC_16S_calibration.R. All three
# covariates are computed once in 0_SETUP.R, identically regardless of
# Kingdom (deg_h_z/precip_72h_z centered WITHIN Season; seq_depth_z --
# ITS's own Fungi-specific sequencing depth column -- centered globally on
# purpose, see 0_SETUP.R's own header).
#
# model_MDSYC_ITS: loga[Mg] ~ dnorm(2,2) is ITS's own scale; b_deg/b_precip/
# b_seq ~ dnorm(0,1) carry over unchanged from model_MDSYC_16S -- additive
# log-scale slopes on standardized covariates, not tied to Hill_1's own
# baseline, and the old ITS lineage's own MDLSYC_model.R already used the
# identical dnorm(0,1) choice.
#
# True b_deg/b_precip/b_seq reused directly from 16S -- the simulator draws
# these covariates independently of Mg/Mo/Yr, so their recoverability
# doesn't depend on Kingdom either.

hiermod_marker <- "ITS"
source('src/hiermod/0_SETUP.R')
source('src/hiermod/Models/MDSYC_model.R') # model_MDSYC_ITS, means_MDSYC(), sim_div_MDSYC(), simulate_from_priors_MDSYC(), dq_MDSYC

model <- model_MDSYC_ITS
model_id <- model_id_MDSYC

hiermod_out_dir <- "out/hiermod/ITS_5_covariates_MDSYC/Calibration"

## Parameter recovery -----------------------------------------------------------

may_conv <- 20
may_org  <- 16
july_conv_shift <- -1.5
july_org_shift  <- 1.3

true_sigma <- cv_to_sigma(c(0.8, 0.6)) # conv, org
true_yr1   <- 0.1
true_yr2   <- 0.3  # implies yr3 = -0.4

true_b_deg    <- 0.2
true_b_precip <- -0.15
true_b_seq    <- 0.3

set.seed(20260916)

dat_sim <- sim_div_MDSYC(
  N_samples = 240,
  loga = log(c(may_conv, may_org)),
  s_conv = july_conv_shift,
  gap_shift = july_org_shift,
  sigma = true_sigma,
  yr1 = true_yr1,
  yr2 = true_yr2,
  b_deg = true_b_deg,
  b_precip = true_b_precip,
  b_seq = true_b_seq,
  shift = 1
)

fit_sim <- ulam(
  model,
  data = as.list(dat_sim),
  chains = 6, cores = 6, iter = 5000,
  control = list(adapt_delta = 0.99))
precis(fit_sim, depth = 2)

post_sim <- extract.samples(fit_sim)

### Fixed effect + sigma recovery ---------

param_recovery <- check_recovery(
  true = list(
    loga1 = log(may_conv), loga2 = log(may_org),
    s_conv = july_conv_shift, gap_shift = july_org_shift,
    sigma1 = true_sigma[1], sigma2 = true_sigma[2],
    yr1 = true_yr1, yr2 = true_yr2,
    b_deg = true_b_deg, b_precip = true_b_precip, b_seq = true_b_seq),
  post_draws = list(
    loga1 = post_sim$loga[,1], loga2 = post_sim$loga[,2],
    s_conv = post_sim$s_conv, gap_shift = post_sim$gap_shift,
    sigma1 = post_sim$sigma[,1], sigma2 = post_sim$sigma[,2],
    yr1 = post_sim$yr1, yr2 = post_sim$yr2,
    b_deg = post_sim$b_deg, b_precip = post_sim$b_precip, b_seq = post_sim$b_seq))

### Contrast recovery -------------

cr <- contrast_recovery(
  fit_sim, means_MDSYC, may_conv, may_org, july_conv_shift, july_org_shift, shift = 1)

p_sim_contrast <- contrast_plot_panels(
  cr$estimands, quant = c(0.001, 0.999), scales = 'free_y',
  group_pal = Management_palette,
  true_vals = cr$true_estimands)

save_report("sim_summary", model_id, recovery = param_recovery, fit_sim, cr$estimands, model)
save_gg("sim_contrast_density", model_id, p_sim_contrast)

save_pdf("sim_trankplot", model_id,
         function() trankplot(fit_sim, max_rows = 50, n_cols = 8),
         width = 16, height = 20)

## Prior predictive check -------------------------------------------------------

n_prior <- 1000
extracted_prior <- extract.prior(fit_sim, n = n_prior)

prior_pred <- map_dfr(seq_len(n_prior), function(i){
  simulate_from_priors_MDSYC(draw_true(extracted_prior, i), shift = 1)
}, .id = "draw")

summary(prior_pred$Dv_shifted) # judge on median/IQR, not mean/SD

p_prior_pc <- prior_predictive_spaghetti(
  prior_pred, value_col = "Dv_shifted", upper_q = 0.99, model = model,
  title = "Prior predictive check", observed = dat_sim$Dv_shifted)

save_gg("sim_prior_PC", model_id, p_prior_pc)

## Full pairwise parameter check ------------------------------------------------

pairs_vars <- c("loga[1]", "loga[2]", "s_conv", "gap_shift", "sigma[1]", "sigma[2]",
                "yr1", "yr2", "b_deg", "b_precip", "b_seq")
p_pairs <- plot_mcmc_pairs(fit_sim, variables = pairs_vars, n_keep = 1000)
save_gg("sim_mcmc_pairs", model_id, p_pairs, width = 15, height = 15, type = "png")

## Simulation-based calibration (SBC), via the SBC package -----------------------

sbc_gen_MDSYC <- make_sbc_generator(
  fit = fit_sim, simulate_fn = simulate_from_priors_MDSYC,
  keep = c("loga", "s_conv", "gap_shift", "sigma", "yr1", "yr2", "b_deg", "b_precip", "b_seq"),
  gen_cols = c("Dv", "Mg", "Mo", "Yr", "deg_h_z", "precip_72h_z", "seq_depth_z"),
  extra_globals = "sim_div_MDSYC", shift = 1)

n_sbc  <- 500
n_iter <- 10000

sbc_MDSYC <- run_sbc_pipeline(
  generator = sbc_gen_MDSYC$generator, globals = sbc_gen_MDSYC$globals,
  n_sbc = n_sbc, model = model, model_id = model_id, n_iter = n_iter,
  hiermod_out_dir = hiermod_out_dir, dquants = dq_MDSYC,
  control = list(adapt_delta = 0.99))

plot_sbc_diagnostics(sbc_MDSYC, model_id, n_sbc)

save_sbc_health_report(
  model_id, sbc_MDSYC, n_sbc, n_iter,
  variables = c("loga[1]", "loga[2]", "s_conv", "gap_shift", "sigma[1]", "sigma[2]",
                "yr1", "yr2", "b_deg", "b_precip", "b_seq",
                "may_gap", "july_gap", "seasonal_change"),
  hiermod_out_dir = hiermod_out_dir)
