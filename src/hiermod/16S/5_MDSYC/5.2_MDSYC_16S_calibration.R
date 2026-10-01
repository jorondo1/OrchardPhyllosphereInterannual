# MODEL 5 (MDSYC), 16S: calibration
# Aim: weather + read-count covariates (see MDSYC_model.R)
# To validate: b_deg, b_precip, b_seq ~ dnorm(0,1)

hiermod_marker <- "16S"
source('src/hiermod/0_SETUP.R')
source('src/hiermod/Models/MDSYC_model.R') 
source('src/utils/sbc_workflow.R') 

model <- model_MDSYC_16S
model_id <- model_id_MDSYC

hiermod_out_dir <- "out/hiermod/16S_5_covariates_MDSYC/Calibration"

## Parameter recovery -----------------------------------------------------------
# Same baseline/gap/year values as MDSYz, plus modest slopes

may_conv <- 180
may_org  <- 120
july_conv_shift <- 0.35
july_org_shift  <- 0.2

true_sigma <- cv_to_sigma(c(0.5, 0.8)) # conv, org
true_yr1   <- 0.1
true_yr2   <- 0.3  # implies yr3 = -0.4

true_b_deg    <- 0.2
true_b_precip <- -0.15
true_b_seq    <- 0.3  # sequencing depth is the one covariate expected to matter most in practice

set.seed(20260911)

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

(param_recovery <- check_recovery(
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
    b_deg = post_sim$b_deg, b_precip = post_sim$b_precip, b_seq = post_sim$b_seq)))

### Contrast recovery -------------
# Covariates at z = 0 cancel: same estimand as MDSYz

cr <- contrast_recovery(
  fit_sim, means_MDSYC, may_conv, may_org, july_conv_shift, july_org_shift, shift = 1)

p_sim_contrast <- contrast_plot_panels(
  cr$estimands, quant = c(0.01, 0.99), scales = 'free_y',
  group_pal = Management_palette,
  true_vals = cr$true_estimands); p_sim_contrast

save_report("sim_summary", model_id, recovery = param_recovery, fit_sim, cr$estimands, model)
save_gg("sim_contrast_density", model_id, p_sim_contrast)

save_pdf("sim_trankplot", model_id,
         function() trankplot(fit_sim, max_rows = 30, n_cols = 5),
         width = 10, height = 14)

## Prior predictive check -------------------------------------------------------

n_prior <- 1000
extracted_prior <- extract.prior(fit_sim, n = n_prior)

prior_pred <- map_dfr(seq_len(n_prior), function(i){
  simulate_from_priors_MDSYC(draw_true(extracted_prior, i), shift = 1)
}, .id = "draw")

summary(prior_pred$Dv_shifted) # judge on median/IQR, not mean/SD

p_prior_pc <- prior_predictive_spaghetti(
  prior_pred, value_col = "Dv_shifted", upper_q = 0.99, model = model,
  title = "Prior predictive check", observed = dat_sim$Dv_shifted); p_prior_pc

save_gg("sim_prior_PC", model_id, p_prior_pc)

## loga/covariate collinearity check -----------------------------------------
# Covariates independent of Mg/Mo/Yr in this simulator: little correlation expected
# - realistic collinearity tested in 6.5

p_funnel <- function(){
  par(mfrow = c(2,2))
  plot(post_sim$b_deg, post_sim$loga[,1],
       xlab = "b_deg", ylab = "loga[1] (Conventional)", pch = 16, col = scales::alpha("black", 0.15))
  plot(post_sim$b_precip, post_sim$loga[,1],
       xlab = "b_precip", ylab = "loga[1] (Conventional)", pch = 16, col = scales::alpha("black", 0.15))
  plot(post_sim$b_seq, post_sim$loga[,1],
       xlab = "b_seq", ylab = "loga[1] (Conventional)", pch = 16, col = scales::alpha("black", 0.15))
  plot(post_sim$b_seq, post_sim$yr1,
       xlab = "b_seq", ylab = "yr1", pch = 16, col = scales::alpha("black", 0.15))
  par(mfrow = c(1,1))
}
save_pdf("loga_sigma_funnel", model_id, p_funnel)

## Simulation-based calibration (SBC), via the SBC package -----------------------

sbc_gen_MDSYC <- make_sbc_generator(
  fit = fit_sim, simulate_fn = simulate_from_priors_MDSYC,
  keep = c("loga", "s_conv", "gap_shift", "sigma", "yr1", "yr2", "b_deg", "b_precip", "b_seq"),
  gen_cols = c("Dv", "Mg", "Mo", "Yr", "deg_h_z", "precip_72h_z", "seq_depth_z"),
  extra_globals = "sim_div_MDSYC", shift = 1)

n_sbc  <- 100
n_iter <- 10000

sbc_MDSYC <- run_sbc_pipeline(
  generator = sbc_gen_MDSYC$generator, globals = sbc_gen_MDSYC$globals,
  n_sbc = n_sbc, model = model, model_id = model_id, n_iter = n_iter,
  hiermod_out_dir = hiermod_out_dir, dquants = dq_MDSYC,
  control = list(adapt_delta = 0.99))

plot_sbc_diagnostics(sbc_MDSYC, model_id, n_sbc)

sbc_MDSYC$stats |>
  dplyr::filter(variable %in% c("loga[1]", "loga[2]", "b_deg", "b_precip", "b_seq")) |>
  dplyr::group_by(variable) |>
  dplyr::summarise(mean_rank_frac = mean(rank / max_rank), median_rank_frac = median(rank / max_rank))

save_sbc_health_report(
  model_id, sbc_MDSYC, n_sbc, n_iter,
  variables = c("loga[1]", "loga[2]", "s_conv", "gap_shift", "sigma[1]", "sigma[2]",
                "yr1", "yr2", "b_deg", "b_precip", "b_seq",
                "may_gap", "july_gap", "seasonal_change"),
  hiermod_out_dir = hiermod_out_dir)

# n = 400 stress test (n = 100 alone not trusted)

n_sbc <- 400

sbc_MDSYC_2 <- run_sbc_pipeline(
  generator = sbc_gen_MDSYC$generator, globals = sbc_gen_MDSYC$globals,
  n_sbc = n_sbc, model = model, model_id = model_id, n_iter = n_iter,
  hiermod_out_dir = hiermod_out_dir, dquants = dq_MDSYC,
  control = list(adapt_delta = 0.99))

plot_sbc_diagnostics(sbc_MDSYC_2, model_id, n_sbc)

sbc_MDSYC_2$stats |>
  dplyr::filter(variable %in% c("loga[1]", "loga[2]", "b_deg", "b_precip", "b_seq")) |>
  dplyr::group_by(variable) |>
  dplyr::summarise(mean_rank_frac = mean(rank / max_rank), median_rank_frac = median(rank / max_rank))

save_sbc_health_report(model_id, sbc_MDSYC_2, n_sbc, n_iter,
                       variables = c("loga[1]", "loga[2]", "s_conv", "gap_shift", "sigma[1]", "sigma[2]",
                                     "yr1", "yr2", "b_deg", "b_precip", "b_seq",
                                     "may_gap", "july_gap", "seasonal_change"),
                       hiermod_out_dir = hiermod_out_dir)
