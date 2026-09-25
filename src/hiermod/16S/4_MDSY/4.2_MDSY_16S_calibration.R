# MODEL 4 (MDSY), 16S: Management x Season interaction plus Year as a FIXED
# effect (3 levels -- see MDSY_model.R header for why fixed, not pooled).
#
# loga[Mg]/s_conv/gap_shift/sigma[Mg] priors are MDS's own validated answer,
# hardcoded here as this model's starting point. yr[Yr] ~ dnorm(0,1) is the
# one new assumption -- and since Year's cardinality is small and fixed (3,
# not scaling with N_samples like Tree's), it's tracked directly in `keep`
# and the health report, unlike tr[Tr].

hiermod_marker <- "16S"
source('src/hiermod/0_SETUP.R')
source('src/hiermod/Models/MDSY_model.R') # model_MDSY_16S, means_MDSY(), sim_div_MDSY(), simulate_from_priors_MDSY(), dq_MDSY

model <- model_MDSY_16S
model_id <- model_id_MDSY

hiermod_out_dir <- "out/hiermod/16S_4_year_MDSY/Calibration"

## Parameter recovery -----------------------------------------------------------
# Same baseline/gap values as MDS/MDST's own calibration, plus year_offset
# (matching MDS2's own convention), so results stay comparable.

may_conv <- 180
may_org  <- 120
july_conv_shift <- 0.35
july_org_shift  <- 0.6

true_sigma <- cv_to_sigma(c(0.5, 0.8)) # conv, org
true_yr    <- c(0, 0.3, -0.2)          # 3 years

set.seed(20260911)

dat_sim <- sim_div_MDSY(
  N_samples = 240,
  loga = log(c(may_conv, may_org)),
  s_conv = july_conv_shift,
  gap_shift = july_org_shift,
  sigma = true_sigma,
  yr = true_yr,
  shift = 1
)

#Fit simulation
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
    yr1 = true_yr[1], yr2 = true_yr[2], yr3 = true_yr[3]),
  post_draws = list(
    loga1 = post_sim$loga[,1], loga2 = post_sim$loga[,2],
    s_conv = post_sim$s_conv, gap_shift = post_sim$gap_shift,
    sigma1 = post_sim$sigma[,1], sigma2 = post_sim$sigma[,2],
    yr1 = post_sim$yr[,1], yr2 = post_sim$yr[,2], yr3 = post_sim$yr[,3])))

### Contrast recovery -------------

cr <- contrast_recovery(
  fit_sim, means_MDSY, may_conv, may_org, july_conv_shift, july_org_shift, shift = 1)

p_sim_contrast <- contrast_plot_panels(
  cr$estimands, quant = c(0.001, 0.999), scales = 'free_y',
  group_pal = Management_palette,
  true_vals = cr$true_estimands); p_sim_contrast

save_report("sim_summary", model_id, recovery = param_recovery, fit_sim, cr$estimands, model)
save_gg("sim_contrast_density", model_id, p_sim_contrast)

save_pdf("sim_trankplot", model_id,
         function() trankplot(fit_sim, max_rows = 30, n_cols = 5),
         width = 10, height = 5)

## Prior predictive check -------------------------------------------------------

n_prior <- 1000
extracted_prior <- extract.prior(fit_sim, n = n_prior)

prior_pred <- map_dfr(seq_len(n_prior), function(i){
  simulate_from_priors_MDSY(draw_true(extracted_prior, i), shift = 1)
}, .id = "draw")

summary(prior_pred$Dv_shifted) # judge on median/IQR, not mean/SD

p_prior_pc <- prior_predictive_spaghetti(
  prior_pred, value_col = "Dv_shifted", upper_q = 0.99, model = model,
  title = "Prior predictive check", observed = dat_sim$Dv_shifted); p_prior_pc

save_gg("sim_prior_PC", model_id, p_prior_pc)

## loga/gamma x sigma[Mg] funnel check ---------------------------------------
# yr[Yr] is a fixed, unpooled location term (no scale parameter of its own),
# so it doesn't carry the same funnel risk sigma_tr did -- worth checking
# anyway since it shares the same likelihood term as loga/s_conv/gap_shift.

p_funnel <- function(){
  par(mfrow = c(2,2))
  plot(post_sim$sigma[,1], post_sim$loga[,1],
       xlab = "sigma[1] (Conventional)", ylab = "loga[1] (Conventional)", pch = 16, col = scales::alpha("black", 0.15))
  plot(post_sim$sigma[,2], post_sim$loga[,2],
       xlab = "sigma[2] (Organic)", ylab = "loga[2] (Organic)", pch = 16, col = scales::alpha("black", 0.15))
  plot(post_sim$yr[,1], post_sim$loga[,1],
       xlab = "yr[1]", ylab = "loga[1] (Conventional)", pch = 16, col = scales::alpha("black", 0.15))
  plot(post_sim$yr[,1], post_sim$s_conv,
       xlab = "yr[1]", ylab = "s_conv", pch = 16, col = scales::alpha("black", 0.15))
  par(mfrow = c(1,1))
}
save_pdf("loga_sigma_funnel", model_id, p_funnel)

## Simulation-based calibration (SBC), via the SBC package -----------------------
# yr[Yr] IS tracked here (unlike tr[Tr]): only 3 levels, fixed cardinality
# regardless of N_samples, so a prior draw of the whole array is exactly
# what generated each replicate's data -- a genuine, fully testable
# fixed-effect recovery, not a nuisance array to exclude.

sbc_gen_MDSY <- make_sbc_generator(
  fit = fit_sim, simulate_fn = simulate_from_priors_MDSY,
  keep = c("loga", "s_conv", "gap_shift", "sigma", "yr"),
  gen_cols = c("Dv", "Mg", "Mo", "Yr"),
  extra_globals = "sim_div_MDSY", shift = 1)

n_sbc  <- 100
n_iter <- 5000

sbc_MDSY <- run_sbc_pipeline(
  generator = sbc_gen_MDSY$generator, globals = sbc_gen_MDSY$globals,
  n_sbc = n_sbc, model = model, model_id = model_id, n_iter = n_iter,
  hiermod_out_dir = hiermod_out_dir, dquants = dq_MDSY,
  control = list(adapt_delta = 0.99))

plot_sbc_diagnostics(sbc_MDSY, model_id, n_sbc)

sbc_MDSY$stats |>
  dplyr::filter(variable %in% c("loga[1]", "loga[2]", "yr[1]", "yr[2]", "yr[3]")) |>
  dplyr::group_by(variable) |>
  dplyr::summarise(mean_rank_frac = mean(rank / max_rank), median_rank_frac = median(rank / max_rank))

save_sbc_health_report(model_id, sbc_MDSY, n_sbc, n_iter,
                        variables = c("loga[1]", "loga[2]", "s_conv", "gap_shift", "sigma[1]", "sigma[2]",
                                      "yr[1]", "yr[2]", "yr[3]", "may_gap", "july_gap", "seasonal_change"),
                        hiermod_out_dir = hiermod_out_dir)
# loga[1]/loga[2] both badly miscalibrated (z=-9.55/-8.97, posterior systematically too high)
# while all three yr are miscalibrated in the opposite direction (z≈+8.6 to +9, posterior 
# systematically too low). Essentially a one-directional leak between the two, 
# not just noise. Collinearity mechanism?? not a funnel: 0 divergences, 
# but 13/100 fits with elevated Rhat (worst 2.09) 

# Solution: zero 