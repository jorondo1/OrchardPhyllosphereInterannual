# MODEL 6 (MDSYCV, "Bombadil the Eldest"), ITS: MDSYC plus Cultivar as a
# FIXED, sum-to-zero effect (5 levels: Cortland, Liberty, Paulared,
# Honeycrisp, Spartan). Mirrors 6.2_MDSYCV_16S_calibration.R.
#
# Cv index order/derived slot (idx$Cv$levels, Liberty=2 derived) is the
# SAME physical Cultivar factor as 16S -- same trees, same samples, just a
# different sequencing barcode -- so the "data-richest level in the
# derived slot" choice carries over unchanged (Liberty: 53 ITS samples,
# still the most, matching 16S's own count).
#
# model_MDSYCV_ITS: loga[Mg] ~ dnorm(2,2) is ITS's own scale;
# cv_1/cv_3/cv_4/cv_5 ~ dnorm(0,1) carry over unchanged from
# model_MDSYCV_16S. True cv_* values reused directly from 16S -- additive
# log-scale offsets, not tied to Hill_1's own baseline, and the simulator
# draws Cultivar independently of Mg/Mo either way.

hiermod_marker <- "ITS"
source('src/hiermod/0_SETUP.R')
source('src/hiermod/Models/MDSYCV_model.R') # model_MDSYCV_ITS, means_MDSYCV(), sim_div_MDSYCV(), simulate_from_priors_MDSYCV(), dq_MDSYCV

model <- model_MDSYCV_ITS
model_id <- model_id_MDSYCV

hiermod_out_dir <- "out/hiermod/ITS_6_cultivar_MDSYCV/Calibration"

## Parameter recovery -----------------------------------------------------------

may_conv <- 20
may_org  <- 16
july_conv_shift <- -1.5
july_org_shift  <- 1.3

true_sigma <- cv_to_sigma(c(0.8, 0.6)) # conv, org
true_yr1   <- 0.1
true_yr2   <- 0.3  # implies yr3 = -0.4

true_cv_1 <- 0.2   # Cortland
true_cv_3 <- -0.15 # Paulared
true_cv_4 <- 0.1   # Honeycrisp
true_cv_5 <- -0.1  # Spartan
# implied cv_2 (Liberty) = -0.05

true_b_deg    <- 0.2
true_b_precip <- -0.15
true_b_seq    <- 0.3

set.seed(20260916)

dat_sim <- sim_div_MDSYCV(
  N_samples = 240,
  loga = log(c(may_conv, may_org)),
  s_conv = july_conv_shift,
  gap_shift = july_org_shift,
  sigma = true_sigma,
  yr1 = true_yr1,
  yr2 = true_yr2,
  cv_1 = true_cv_1,
  cv_3 = true_cv_3,
  cv_4 = true_cv_4,
  cv_5 = true_cv_5,
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

### Fixed effect + sigma recovery ---------
post_sim <- extract.samples(fit_sim)

(param_recovery <- check_recovery(
  true = list(
    loga1 = log(may_conv), loga2 = log(may_org),
    s_conv = july_conv_shift, gap_shift = july_org_shift,
    sigma1 = true_sigma[1], sigma2 = true_sigma[2],
    yr1 = true_yr1, yr2 = true_yr2,
    cv_1 = true_cv_1, cv_3 = true_cv_3, cv_4 = true_cv_4, cv_5 = true_cv_5,
    b_deg = true_b_deg, b_precip = true_b_precip, b_seq = true_b_seq),
  post_draws = list(
    loga1 = post_sim$loga[,1], loga2 = post_sim$loga[,2],
    s_conv = post_sim$s_conv, gap_shift = post_sim$gap_shift,
    sigma1 = post_sim$sigma[,1], sigma2 = post_sim$sigma[,2],
    yr1 = post_sim$yr1, yr2 = post_sim$yr2,
    cv_1 = post_sim$cv_1, cv_3 = post_sim$cv_3, cv_4 = post_sim$cv_4, cv_5 = post_sim$cv_5,
    b_deg = post_sim$b_deg, b_precip = post_sim$b_precip, b_seq = post_sim$b_seq)))

### Contrast recovery -------------

cr <- contrast_recovery(
  fit_sim, means_MDSYCV, may_conv, may_org, july_conv_shift, july_org_shift, shift = 1)

p_sim_contrast <- contrast_plot_panels(
  cr$estimands, quant = c(0.001, 0.999), scales = 'free_y',
  group_pal = Management_palette,
  true_vals = cr$true_estimands); p_sim_contrast

save_report("sim_summary", model_id, recovery = param_recovery, fit_sim, cr$estimands, model)
save_gg("sim_contrast_density", model_id, p_sim_contrast)

save_pdf("sim_trankplot", model_id,
         function() trankplot(fit_sim, max_rows = 50, n_cols = 8),
         width = 16, height = 20)

## Prior predictive check -------------------------------------------------------

n_prior <- 1000
extracted_prior <- extract.prior(fit_sim, n = n_prior)

prior_pred <- map_dfr(seq_len(n_prior), function(i){
  simulate_from_priors_MDSYCV(draw_true(extracted_prior, i), shift = 1)
}, .id = "draw")

summary(prior_pred$Dv_shifted) # judge on median/IQR, not mean/SD

p_prior_pc <- prior_predictive_spaghetti(
  prior_pred, value_col = "Dv_shifted", upper_q = 0.99, model = model,
  title = "Prior predictive check", observed = dat_sim$Dv_shifted); p_prior_pc

save_gg("sim_prior_PC", model_id, p_prior_pc)

## Full pairwise parameter check ------------------------------------------------

pairs_vars <- c("loga[1]", "loga[2]", "s_conv", "gap_shift", "sigma[1]", "sigma[2]", "yr1", "yr2",
                 "cv_1", "cv_3", "cv_4", "cv_5", "b_deg", "b_precip", "b_seq")
p_pairs <- plot_mcmc_pairs(fit_sim, variables = pairs_vars, n_keep = 1000)
save_gg("sim_mcmc_pairs", model_id, p_pairs, width = 15, height = 15, type = "png")

## Simulation-based calibration (SBC), via the SBC package -----------------------

sbc_gen_MDSYCV <- make_sbc_generator(
  fit = fit_sim, simulate_fn = simulate_from_priors_MDSYCV,
  keep = c("loga", "s_conv", "gap_shift", "sigma", "yr1", "yr2",
           "cv_1", "cv_3", "cv_4", "cv_5", "b_deg", "b_precip", "b_seq"),
  gen_cols = c("Dv", "Mg", "Mo", "Yr", "Cv", "deg_h_z", "precip_72h_z", "seq_depth_z"),
  extra_globals = "sim_div_MDSYCV", shift = 1)

n_sbc  <- 100
n_iter <- 10000

sbc_MDSYCV <- run_sbc_pipeline(
  generator = sbc_gen_MDSYCV$generator, globals = sbc_gen_MDSYCV$globals,
  n_sbc = n_sbc, model = model, model_id = model_id, n_iter = n_iter,
  hiermod_out_dir = hiermod_out_dir, dquants = dq_MDSYCV,
  control = list(adapt_delta = 0.99))

plot_sbc_diagnostics(sbc_MDSYCV, model_id, n_sbc)

sbc_MDSYCV$stats |>
  dplyr::filter(variable %in% c("loga[1]", "loga[2]", "cv_1", "cv_3", "cv_4", "cv_5")) |>
  dplyr::group_by(variable) |>
  dplyr::summarise(mean_rank_frac = mean(rank / max_rank), median_rank_frac = median(rank / max_rank))

save_sbc_health_report(model_id, sbc_MDSYCV, n_sbc, n_iter,
                        variables = c("loga[1]", "loga[2]", "s_conv", "gap_shift", "sigma[1]", "sigma[2]",
                                      "yr1", "yr2", "cv_1", "cv_3", "cv_4", "cv_5",
                                      "b_deg", "b_precip", "b_seq",
                                      "may_gap", "july_gap", "seasonal_change"),
                        hiermod_out_dir = hiermod_out_dir)

# Same n=400 stress test as every model in this rebuild before trusting an
# n=100 "ok".
n_sbc <- 400
sbc_MDSYCV_2 <- run_sbc_pipeline(
  generator = sbc_gen_MDSYCV$generator, globals = sbc_gen_MDSYCV$globals,
  n_sbc = n_sbc, model = model, model_id = model_id, n_iter = n_iter,
  hiermod_out_dir = hiermod_out_dir, dquants = dq_MDSYCV,
  control = list(adapt_delta = 0.99))

plot_sbc_diagnostics(sbc_MDSYCV_2, model_id, n_sbc)

save_sbc_health_report(model_id, sbc_MDSYCV_2, n_sbc, n_iter,
                        variables = c("loga[1]", "loga[2]", "s_conv", "gap_shift", "sigma[1]", "sigma[2]",
                                      "yr1", "yr2", "cv_1", "cv_3", "cv_4", "cv_5",
                                      "b_deg", "b_precip", "b_seq",
                                      "may_gap", "july_gap", "seasonal_change"),
                        hiermod_out_dir = hiermod_out_dir)
