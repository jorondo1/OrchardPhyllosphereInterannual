# MODEL 5 (MDLS2), 16S, SHIFTED (Hill_1 - 1): same structure as ITS Model 5
# (sigma[cell], gamma <- s_conv + gap_shift*(Mg-1), b[Lo], tr[Tr], yr[Yr]
# fixed/unpooled), starting from 16S Model 2's own validated prior
# (loga ~ dnorm(6,2), sigma family ~ dexp(2)) rather than ITS's. Parameter
# recovery, prior-predictive check, and SBC.

hiermod_marker <- "16S"
source('src/hiermod/0_SETUP.R')
source('src/hiermod/Models/MDLS2_model.R') # model_MDLS2_16S, means_MDLS2(), sim_div_MDLS2(), contrast_may_gap_MDLS2(), simulate_from_priors_MDLS2()
model <- model_MDLS2_16S

hiermod_out_dir <- "out/hiermod/16S_5_lognormal_MDLS2_shifted"

## Parameter recovery -----------------------------------------------------------

# Baseline diversity/gap on 16S's own raw scale (same May Conventional/
# Organic values as 2.2_MDL_calibration.R's conv/org), sigma spread picked
# to differ a lot across cells (same pattern as ITS's own Model 5
# calibration), so this tests whether 4 distinguishable cells come back
# distinguishable, not just plausible.

may_conv <- 180         # hill scale, matches 2.2_MDL_calibration.R's conv
may_org  <- 120         # hill scale, matches 2.2_MDL_calibration.R's org
july_conv_shift <- 0.35 # log scale
july_org_shift  <- 0.2  # log scale

true_sigma <- cv_to_sigma(c(0.5, 0.3, 0.8, 0.6)) # conv_May, conv_July, org_May, org_July

set.seed(20260911)

dat_sim <- sim_div_MDLS2(
  N_samples = 240,
  n_loc = 4,
  loga = log(c(may_conv, may_org)),
  s_conv = july_conv_shift,
  gap_shift = july_org_shift,
  sigma = true_sigma,
  year_offset = c(0, 0.3, -0.2),
  p_dropout = 0.1,
  shift = 1
); head(dat_sim)

dat_sim %>% count(Lo, Yr, Mo, Mg) %>% print(n=100)
hist(dat_sim$Dv_shifted, breaks = 30)

fit_sim <- ulam(
  model,
  data = as.list(dat_sim),
  chains = 6, cores = 6, iter = 10000,
  control = list(adapt_delta = 0.99)
)

precis(fit_sim, depth = 2)

post_sim <- extract.samples(fit_sim)

### Fixed effects recovery -------------

fixed_recovery <- check_recovery(
  true = c(loga1 = log(may_conv),
           loga2 = log(may_org),
           s_conv = july_conv_shift,
           gap_shift = july_org_shift),
  post_draws = list(
    loga1 = post_sim$loga[,1],
    loga2 = post_sim$loga[,2],
    s_conv = post_sim$s_conv,
    gap_shift = post_sim$gap_shift)
)

### Sigma recovery -------------

sigma_recovery <- check_recovery(
  true = list(sigma1 = true_sigma[1], sigma2 = true_sigma[2],
              sigma3 = true_sigma[3], sigma4 = true_sigma[4]),
  post_draws = list(sigma1 = post_sim$sigma[,1], sigma2 = post_sim$sigma[,2],
                    sigma3 = post_sim$sigma[,3], sigma4 = post_sim$sigma[,4])
)

### Contrast recovery ------------------------------
# shift = 1: raw mean/median cell values should reflect the floor
# (contrasts wouldn't need it, it cancels, but contrast_recovery()'s own
# cell-level median_1..median_4 columns do).

cr <- contrast_recovery(fit_sim, means_MDLS2, may_conv, may_org, july_conv_shift, july_org_shift, shift = 1)

p_sim_contrast <- contrast_plot_panels(
  cr$estimands, quant = c(0.01, 0.99), scales = 'free_y',
  group_pal = Management_palette,
  true_vals = cr$true_estimands); p_sim_contrast

save_report("sim_summary", "MDLS2_shifted", fit_sim, cr$estimands, model,
            recovery = bind_rows(fixed_recovery, sigma_recovery), model_name = "The Splitter")
save_gg("sim_contrast_density", "MDLS2_shifted", p_sim_contrast)

## Prior predictive check -------------------------------------------------------

n_prior <- 1000
extracted_prior <- extract.prior(fit_sim, n = n_prior)

prior_pred <- map_dfr(seq_len(n_prior), function(i){
  simulate_from_priors_MDLS2(draw_true(extracted_prior, i), shift = 1)
}, .id = "draw")

summary(prior_pred$Dv_shifted)

p_prior_pc <- prior_predictive_spaghetti(
  prior_pred, value_col = "Dv_shifted", upper_q = 0.99, model = model,
  title = "Prior predictive check", observed = dat_sim$Dv_shifted) ; p_prior_pc

save_gg("sim_prior_PC", "MDLS2_shifted", p_prior_pc)

# Which prior is driving the extreme tail?
(xlim_upper <- quantile(prior_pred$Dv_shifted, 0.99))
extreme_draws <- prior_pred %>%
  group_by(draw) %>%
  summarise(max_Dv = max(Dv_shifted)) %>%
  filter(max_Dv > xlim_upper) %>%
  pull(draw) %>% as.integer()

tibble(
  draw       = seq_len(n_prior),
  extreme    = seq_len(n_prior) %in% extreme_draws,
  sigma_loc  = extracted_prior$sigma_loc,
  sigma_tr = extracted_prior$sigma_tr,
  sigma_max  = pmax(extracted_prior$sigma[,1], extracted_prior$sigma[,2])) %>%
  group_by(extreme) %>%
  summarise(across(c(sigma_loc, sigma_max, sigma_tr), median))

# Biggest gap  is with sigma_tr, then sigma_max

### Update sigma priors --------------------
model_ppc1 <- model
model_ppc1$pr_sigma_tr  <- quote(sigma_tr  ~ dexp(2.5)) #tighter
model_ppc1$pr_sigma     <- quote(sigma[cell] ~ dexp(2.5))
model_ppc1$pr_sigma_loc <- quote(sigma_loc ~ dexp(2.5))

fit_sim_ppc1 <- ulam(
  model_ppc1,
  data = as.list(dat_sim),
  chains = 6, cores = 6, iter = 5000 )
precis(fit_sim_ppc1, depth = 2 )

prior_ppc1 <- extract.prior(fit_sim_ppc1, n = n_prior)

# One simulated dataset per prior draw
prior_ppc1_pred <- map_dfr(seq_len(n_prior), function(i){
  simulate_from_priors_MDLS2(draw_true(prior_ppc1, i))
}, .id = "draw")

summary(prior_ppc1_pred$Dv) # raw per-observation draws
# Much better!
# Mean still very high (trillions), but
# sd more reasonable and 4rd quartile also

p_sim_spaghetti_ppc1 <- prior_predictive_spaghetti(
  prior_ppc1_pred, upper_q = 0.99, model = model_ppc1,
  title = "Prior predictive check: adapted priors)"); p_sim_spaghetti_ppc1

save_gg("sim_prior_PC_ppc1", "MDLS2", p_sim_spaghetti_ppc1)


## Parameter recovery: updated priors -----------------------------------------------------------
post_sim_ppc1 <- extract.samples(fit_sim_ppc1)

### Fixed effects recovery -------------

fixed_recovery <- check_recovery(
  true = c(loga1 = log(may_conv),
           loga2 = log(may_org),
           s_conv = july_conv_shift,
           gap_shift = july_org_shift),
  post_draws = list(
    loga1 = post_sim_ppc1$loga[,1],
    loga2 = post_sim_ppc1$loga[,2],
    s_conv = post_sim_ppc1$s_conv,
    gap_shift = post_sim_ppc1$gap_shift)
)

### Sigma recovery -------------

sigma_recovery <- check_recovery(
  true = list(sigma1 = true_sigma[1], sigma2 = true_sigma[2],
              sigma3 = true_sigma[3], sigma4 = true_sigma[4]),
  post_draws = list(sigma1 = post_sim_ppc1$sigma[,1], sigma2 = post_sim_ppc1$sigma[,2],
                    sigma3 = post_sim_ppc1$sigma[,3], sigma4 = post_sim_ppc1$sigma[,4])
)

### Contrast recovery ------------------------------
# shift = 1: raw mean/median cell values should reflect the floor
# (contrasts wouldn't need it, it cancels, but contrast_recovery()'s own
# cell-level median_1..median_4 columns do).

cr <- contrast_recovery(fit_sim_ppc1, means_MDLS2, may_conv, may_org, july_conv_shift, july_org_shift)

p_sim_contrast_ppc1 <- contrast_plot_panels(
  cr$estimands, quant = c(0.01, 0.99), scales = 'free_y',
  group_pal = Management_palette,
  true_vals = cr$true_estimands); p_sim_contrast_ppc1

save_report("sim_summary", "MDLS2_ppc1", fit_sim_ppc1, cr$estimands, model,
            recovery = bind_rows(fixed_recovery, sigma_recovery), model_name = "The Splitter")
save_gg("sim_contrast_density", "MDLS2_ppc1", p_sim_contrast_ppc1)

# Traces look fine:
save_pdf("sim_traceplot", "MDLb_ppc1", function() traceplot(fit_sim_ppc1))
save_pdf("sim_trankplot", "MDLb-ppc1", function() trankplot(fit_sim_ppc1))

## Simulation-based calibration (SBC) --------------------------------------------

sbc_MDLS2_shifted <- run_sbc(
  model_fit   = fit_sim_ppc1,
  means_fn    = means_MDLS2,
  contrast_fn = contrast_may_gap_MDLS2,
  simulate_fn = function(true_params) simulate_from_priors_MDLS2(true_params, shift = 1),
  n_sbc = 30, iter = 15000, n_parallel = 4, chains = 2, cores = 2,
  control = list(adapt_delta = 0.99))

sbc_out_MDLS2_shifted <- save_sbc_report(sbc_MDLS2_shifted, "MDLS2_shifted_30sbc_iter")
hist(sbc_out_MDLS2_shifted$ranks, breaks = 30)

# Export updated model ! -------------
saveRDS(model_ppc1, file.path(hiermod_out_dir, "model_ppc1.rds"))
