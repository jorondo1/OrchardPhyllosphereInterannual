# MODEL 4 (MDSYz), 16S: MDSY + sum-to-zero constraint on Year. MDSY's own
# SBC (4.2) found loga[1]/loga[2] severely miscalibrated (biased high) and
# all three yr[] miscalibrated in the opposite direction (biased low) -- a
# one-directional leak, not noise. See MDSYz_model.R's header for the fix
# (2 free scalars yr1/yr2, third level = -(yr1+yr2) by construction) and
# its direct precedent (Location's own MDLS2v -> MDLS2vz history).

hiermod_marker <- "16S"
source('src/hiermod/0_SETUP.R')
source('src/hiermod/Models/MDSYz_model.R') # model_MDSYz_16S, means_MDSYz(), sim_div_MDSYz(), simulate_from_priors_MDSYz(), dq_MDSYz

model <- model_MDSYz_16S
model_id <- model_id_MDSYz

hiermod_out_dir <- "out/hiermod/16S_4_year_MDSY/Calibration"

## Parameter recovery -----------------------------------------------------------
# Same baseline/gap values as MDSY's own calibration, so results stay
# comparable. true_yr1/true_yr2 chosen so the implied yr3 = -(yr1+yr2)
# matches MDSY's own true_yr[3] = -0.2 as closely as possible while keeping
# all three comparable in magnitude (0, 0.3, -0.2 doesn't sum to zero, so
# it can't be reused verbatim under this constraint).

may_conv <- 180
may_org  <- 120
july_conv_shift <- 0.35
july_org_shift  <- 0.6

true_sigma <- cv_to_sigma(c(0.5, 0.8)) # conv, org
true_yr1   <- 0.1
true_yr2   <- 0.3  # implies yr3 = -(0.1+0.3) = -0.4

set.seed(20260911)

dat_sim <- sim_div_MDSYz(
  N_samples = 240,
  loga = log(c(may_conv, may_org)),
  s_conv = july_conv_shift,
  gap_shift = july_org_shift,
  sigma = true_sigma,
  yr1 = true_yr1,
  yr2 = true_yr2,
  shift = 1
); head(dat_sim)

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
    yr1 = true_yr1, yr2 = true_yr2, yr3 = -(true_yr1 + true_yr2)),
  post_draws = list(
    loga1 = post_sim$loga[,1], loga2 = post_sim$loga[,2],
    s_conv = post_sim$s_conv, gap_shift = post_sim$gap_shift,
    sigma1 = post_sim$sigma[,1], sigma2 = post_sim$sigma[,2],
    yr1 = post_sim$yr1, yr2 = post_sim$yr2, yr3 = -(post_sim$yr1 + post_sim$yr2))))
save_report("sim_summary", model_id, recovery = param_recovery, fit_sim, cr$estimands, model)

### Contrast recovery -------------

cr <- contrast_recovery(
  fit_sim, means_MDSYz, may_conv, may_org, july_conv_shift, july_org_shift, shift = 1)

p_sim_contrast <- contrast_plot_panels(
  cr$estimands, quant = c(0.001, 0.999), scales = 'free_y',
  group_pal = Management_palette,
  true_vals = cr$true_estimands); p_sim_contrast

save_gg("sim_contrast_density", model_id, p_sim_contrast)

save_pdf("sim_trankplot", model_id,
         function() trankplot(fit_sim, max_rows = 30, n_cols = 5),
         width = 10, height = 5)

## Prior predictive check -------------------------------------------------------

n_prior <- 1000
extracted_prior <- extract.prior(fit_sim, n = n_prior)

prior_pred <- map_dfr(seq_len(n_prior), function(i){
  simulate_from_priors_MDSYz(draw_true(extracted_prior, i), shift = 1)
}, .id = "draw")

summary(prior_pred$Dv_shifted) # judge on median/IQR, not mean/SD

p_prior_pc <- prior_predictive_spaghetti(
  prior_pred, value_col = "Dv_shifted", upper_q = 0.99, model = model,
  title = "Prior predictive check", observed = dat_sim$Dv_shifted); p_prior_pc

save_gg("sim_prior_PC", model_id, p_prior_pc)

## loga/gamma x yr1/yr2 correlation check ------------------------------------
# The check that actually matters this time -- did the sum-to-zero
# constraint kill the loga/yr collinearity MDSY showed?

p_funnel <- function(){
  par(mfrow = c(2,2))
  plot(post_sim$yr1, post_sim$loga[,1],
       xlab = "yr1", ylab = "loga[1] (Conventional)", pch = 16, col = scales::alpha("black", 0.15))
  plot(post_sim$yr2, post_sim$loga[,1],
       xlab = "yr2", ylab = "loga[1] (Conventional)", pch = 16, col = scales::alpha("black", 0.15))
  plot(post_sim$sigma[,1], post_sim$loga[,1],
       xlab = "sigma[1] (Conventional)", ylab = "loga[1] (Conventional)", pch = 16, col = scales::alpha("black", 0.15))
  plot(post_sim$yr1, post_sim$yr2,
       xlab = "yr1", ylab = "yr2", pch = 16, col = scales::alpha("black", 0.15))
  par(mfrow = c(1,1))
}
save_pdf("loga_sigma_funnel", model_id, p_funnel)

cat("cor(yr1, loga[1]):", cor(post_sim$yr1, post_sim$loga[,1]), "\n")
cat("cor(yr2, loga[1]):", cor(post_sim$yr2, post_sim$loga[,1]), "\n")

## Simulation-based calibration (SBC), via the SBC package -----------------------

sbc_gen_MDSYz <- make_sbc_generator(
  fit = fit_sim, simulate_fn = simulate_from_priors_MDSYz,
  keep = c("loga", "s_conv", "gap_shift", "sigma", "yr1", "yr2"),
  gen_cols = c("Dv", "Mg", "Mo", "Yr"),
  extra_globals = "sim_div_MDSYz", shift = 1)

n_sbc  <- 100
n_iter <- 10000

sbc_MDSYz <- run_sbc_pipeline(
  generator = sbc_gen_MDSYz$generator, globals = sbc_gen_MDSYz$globals,
  n_sbc = n_sbc, model = model, model_id = model_id, n_iter = n_iter,
  hiermod_out_dir = hiermod_out_dir, dquants = dq_MDSYz,
  control = list(adapt_delta = 0.99))

plot_sbc_diagnostics(sbc_MDSYz, model_id, n_sbc)

sbc_MDSYz$stats |>
  dplyr::filter(variable %in% c("loga[1]", "loga[2]", "yr1", "yr2")) |>
  dplyr::group_by(variable) |>
  dplyr::summarise(mean_rank_frac = mean(rank / max_rank), median_rank_frac = median(rank / max_rank))

save_sbc_health_report(
  model_id, sbc_MDSYz, n_sbc, n_iter,
  variables = c("loga[1]", "loga[2]", "s_conv", "gap_shift", "sigma[1]", "sigma[2]",
                "yr1", "yr2", "may_gap", "july_gap", "seasonal_change"),
  hiermod_out_dir = hiermod_out_dir)


# Same stress test as every model in this rebuild before trusting an "ok"
# result at n=100.

n_sbc <- 400

sbc_MDSYz_2 <- run_sbc_pipeline(
  generator = sbc_gen_MDSYz$generator, globals = sbc_gen_MDSYz$globals,
  n_sbc = n_sbc, model = model, model_id = model_id, n_iter = n_iter,
  hiermod_out_dir = hiermod_out_dir, dquants = dq_MDSYz,
  control = list(adapt_delta = 0.99))

plot_sbc_diagnostics(sbc_MDSYz_2, model_id, n_sbc)

sbc_MDSYz_2$stats |>
  dplyr::filter(variable %in% c("loga[1]", "loga[2]", "yr1", "yr2")) |>
  dplyr::group_by(variable) |>
  dplyr::summarise(mean_rank_frac = mean(rank / max_rank), median_rank_frac = median(rank / max_rank))

save_sbc_health_report(model_id, sbc_MDSYz_2, n_sbc, n_iter,
                       variables = c("loga[1]", "loga[2]", "s_conv", "gap_shift", "sigma[1]", "sigma[2]",
                                     "yr1", "yr2", "may_gap", "july_gap", "seasonal_change"),
                       hiermod_out_dir = hiermod_out_dir)
