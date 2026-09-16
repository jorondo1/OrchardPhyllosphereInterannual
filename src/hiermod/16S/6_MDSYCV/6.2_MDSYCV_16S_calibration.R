# MODEL 6 (MDSYCV, "Bombadil the Eldest"), 16S: MDSYC plus Cultivar as a
# FIXED, sum-to-zero effect (5 levels).
#
# loga[Mg]/s_conv/gap_shift/sigma[Mg]/yr1/yr2/b_deg/b_precip/b_seq priors are
# MDSYC's own validated answer, hardcoded here as this model's starting
# point. cv_1/cv_3/cv_4/cv_5 ~ dnorm(0,1) is the one new assumption --
# sum-to-zero applied from the start this time (see MDSYCV_model.R header),
# so this is testing the construction itself, not the leak Year already
# taught us to avoid.

hiermod_marker <- "16S"
source('src/hiermod/0_SETUP.R')
source('src/hiermod/Models/MDSYCV_model.R') # model_MDSYCV_16S, means_MDSYCV(), sim_div_MDSYCV(), simulate_from_priors_MDSYCV(), dq_MDSYCV

model <- model_MDSYCV_16S
model_id <- model_id_MDSYCV

hiermod_out_dir <- "out/hiermod/16S_6_cultivar_MDSYCV"

## Parameter recovery -----------------------------------------------------------
# Same baseline/gap/year/covariate values as MDSYC's own calibration, plus
# modest cultivar offsets, so results stay comparable.

may_conv <- 180
may_org  <- 120
july_conv_shift <- 0.35
july_org_shift  <- 0.2

true_sigma <- cv_to_sigma(c(0.5, 0.8)) # conv, org
true_yr1   <- 0.1
true_yr2   <- 0.3  # implies yr3 = -0.4

true_b_deg    <- 0.2
true_b_precip <- -0.15
true_b_seq    <- 0.3

true_cv_1 <- 0.2   # Cortland
true_cv_3 <- -0.15 # Paulared
true_cv_4 <- 0.1   # Honeycrisp
true_cv_5 <- -0.1  # Spartan
# implied cv_2 (Liberty) = -(0.2 - 0.15 + 0.1 - 0.1) = -0.05

set.seed(20260911)

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
    b_deg = true_b_deg, b_precip = true_b_precip, b_seq = true_b_seq,
    cv_1 = true_cv_1, cv_3 = true_cv_3, cv_4 = true_cv_4, cv_5 = true_cv_5),
  post_draws = list(
    loga1 = post_sim$loga[,1], loga2 = post_sim$loga[,2],
    s_conv = post_sim$s_conv, gap_shift = post_sim$gap_shift,
    sigma1 = post_sim$sigma[,1], sigma2 = post_sim$sigma[,2],
    yr1 = post_sim$yr1, yr2 = post_sim$yr2,
    b_deg = post_sim$b_deg, b_precip = post_sim$b_precip, b_seq = post_sim$b_seq,
    cv_1 = post_sim$cv_1, cv_3 = post_sim$cv_3, cv_4 = post_sim$cv_4, cv_5 = post_sim$cv_5)))

### Contrast recovery -------------
# cv_1/cv_3/cv_4/cv_5 cancel at the default (z=0/average-cultivar)
# reference level, same reasoning as Year/covariates -- unaffected recovery
# estimand vs MDSYC.

cr <- contrast_recovery(
  fit_sim, means_MDSYCV, may_conv, may_org, july_conv_shift, july_org_shift, shift = 1)

p_sim_contrast <- contrast_plot_panels(
  cr$estimands, quant = c(0.01, 0.99), scales = 'free_y',
  group_pal = Management_palette,
  true_vals = cr$true_estimands); p_sim_contrast

save_report("sim_summary", model_id, recovery = param_recovery, fit_sim, cr$estimands, model, model_name = "Bombadil the Eldest")
save_gg("sim_contrast_density", model_id, p_sim_contrast)

save_pdf("sim_trankplot", model_id,
         function() trankplot(fit_sim, max_rows = 30, n_cols = 5),
         width = 10, height = 12)

## Prior predictive check -------------------------------------------------------

n_prior <- 1000
extracted_prior <- extract.prior(fit_sim, n = n_prior)

prior_pred <- map_dfr(seq_len(n_prior), function(i){
  simulate_from_priors_MDSYCV(draw_true(extracted_prior, i), shift = 1)
}, .id = "draw")

summary(prior_pred$Dv_shifted) # judge on median/IQR, not mean/SD

p_prior_pc <- prior_predictive_spaghetti(
  prior_pred, value_col = "Dv_shifted", upper_q = 0.99, model = model,
  title = "Prior predictive check", observed = dat_sim$Dv_shifted)

save_gg("sim_prior_PC", model_id, p_prior_pc)

## Full pairwise parameter check -----------------------------------------------
# Replaces the old hand-picked plot() grid (which spot-checked loga[1] only,
# never loga[2], against a few chosen partners -- a panel-count shortcut,
# not a principled choice, and this project has already seen miscalibration
# land on loga[1] in one model and loga[2] in another). mcmc_pairs()
# (plot_mcmc_pairs()/thin_for_pairs(), hiermod_core.R) gives every parameter
# against every other in one grid, with any divergent transitions
# highlighted directly on it. 13x13 is a lot of panels -- if it's too dense
# to read in practice, worth trimming back to a representative subset, but
# starting from the full grid rather than a hand-picked one.

pairs_vars <- c("loga[1]", "loga[2]", "sigma[1]", "sigma[2]", "yr1", "yr2",
                 "cv_1", "cv_3", "cv_4", "cv_5", "b_deg", "b_precip", "b_seq")

p_pairs <- plot_mcmc_pairs(fit_sim, variables = pairs_vars, n_keep = 1000)
save_gg("sim_mcmc_pairs", model_id, p_pairs, width = 15, height = 15, type = "png")

## Simulation-based calibration (SBC), via the SBC package -----------------------

sbc_gen_MDSYCV <- make_sbc_generator(
  fit = fit_sim, simulate_fn = simulate_from_priors_MDSYCV,
  keep = c("loga", "s_conv", "gap_shift", "sigma", "yr1", "yr2",
           "b_deg", "b_precip", "b_seq", "cv_1", "cv_3", "cv_4", "cv_5"),
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
                                      "yr1", "yr2", "b_deg", "b_precip", "b_seq",
                                      "cv_1", "cv_3", "cv_4", "cv_5",
                                      "may_gap", "july_gap", "seasonal_change"),
                        hiermod_out_dir = hiermod_out_dir)

# Same stress test as every model in this rebuild before trusting an "ok"
# result at n=100.

n_sbc <- 400

sbc_MDSYCV_2 <- run_sbc_pipeline(
  generator = sbc_gen_MDSYCV$generator, globals = sbc_gen_MDSYCV$globals,
  n_sbc = n_sbc, model = model, model_id = model_id, n_iter = n_iter,
  hiermod_out_dir = hiermod_out_dir, dquants = dq_MDSYCV,
  control = list(adapt_delta = 0.99))

plot_sbc_diagnostics(sbc_MDSYCV_2, model_id, n_sbc)

sbc_MDSYCV_2$stats |>
  dplyr::filter(variable %in% c("loga[1]", "loga[2]", "cv_1", "cv_3", "cv_4", "cv_5")) |>
  dplyr::group_by(variable) |>
  dplyr::summarise(mean_rank_frac = mean(rank / max_rank), median_rank_frac = median(rank / max_rank))

save_sbc_health_report(model_id, sbc_MDSYCV_2, n_sbc, n_iter,
                        variables = c("loga[1]", "loga[2]", "s_conv", "gap_shift", "sigma[1]", "sigma[2]",
                                      "yr1", "yr2", "b_deg", "b_precip", "b_seq",
                                      "cv_1", "cv_3", "cv_4", "cv_5",
                                      "may_gap", "july_gap", "seasonal_change"),
                        hiermod_out_dir = hiermod_out_dir)
