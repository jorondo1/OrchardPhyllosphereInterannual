# MODEL 5, SHIFTED: identical structure to Model 5, but the likelihood is
# fit to (Hill_1 - 1) instead of Hill_1, since Hill numbers can't go below 1.
# Parameter recovery, prior-predictive check, and SBC.

hiermod_marker <- "ITS"
source('src/hiermod/0_SETUP.R')
source('src/hiermod/Models/MDLS2_model.R') # model_MDLS2_ITS, means_MDLS2(), sim_div_MDLS2(), contrast_may_gap_MDLS2(), simulate_from_priors_MDLS2()
model <- model_MDLS2_ITS

hiermod_out_dir <- "out/hiermod/ITS_5_lognormal_MDLS2_shifted"

## Model specification ---------------------------------------------------------

# Same as plain Model 5: likelihood stays `Dv ~ dlnorm(mu, sigma[cell])`:
# Dv is meant to be the *unshifted* (0, Inf)-support lognormal draw, and
# that's true whether or not a floor exists in the real world -- the floor
# only has to be added back (or subtracted out of real data) outside the
# likelihood. Concretely, the only two differences from plain Model 5:
#   - sim_div_MDLS2()/simulate_from_priors_MDLS2() called with shift = 1, so
#     dat_sim$Dv_shifted (= 1 + Dv) is also produced, purely for reporting/
#     comparison -- never fed to the likelihood.
#   - the real-data `dat$Dv` is built as `div$Hill_1 - 1` (subtracting the
#     floor to get back to (0, Inf) support), instead of Model 5's raw
#     `div$Hill_1`.
# means_MDLS2()'s `shift = 1` argument (threaded through post_full()'s ...)
# adds the floor back only when backtransforming to a raw mean/median.
# contrasts don't need it since a shared additive constant cancels
# in any group1-group2 difference.

## Parameter recovery -----------------------------------------------------------

# Same as model 5
may_conv <- 4          # hill scale
may_org <- 7           # hill scale
july_conv_shift <- 0.35 # log scale
july_org_shift <- 0.2  # log scale

true_sigma <- cv_to_sigma(c(0.6, 0.25, 0.8, 0.8))

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

# Dv is what the model actually fits; Dv_shifted (= 1 + Dv) is only for
# sanity-checking that the floor looks right; its NOT passed to ulam().
hist(dat_sim$Dv_shifted, breaks = 30)

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
           gap_shift = july_org_shift),
  post_draws = list(
    loga1 = post_sim$loga[,1],
    loga2 = post_sim$loga[,2],
    s_conv = post_sim$s_conv,
    gap_shift = post_sim$gap_shift)
); fixed_recovery

### Sigma recovery -------------

sigma_recovery <- check_recovery(
  true = list(sigma1 = true_sigma[1], sigma2 = true_sigma[2],
              sigma3 = true_sigma[3], sigma4 = true_sigma[4]),
  post_draws = list(sigma1 = post_sim$sigma[,1], sigma2 = post_sim$sigma[,2],
                    sigma3 = post_sim$sigma[,3], sigma4 = post_sim$sigma[,4])
); sigma_recovery # One sigma is a bit off, nothing to panick about

### Contrast recovery ------------------------------
# shift = 1: raw mean/median cell values should reflect the floor
# (contrasts wouldn't need it, it cancels, but contrast_recovery()'s own
# cell-level median_1..median_4 columns do).

cr <- contrast_recovery(fit_sim, means_MDLS2, may_conv, may_org, july_conv_shift, july_org_shift, shift = 1)

p_sim_contrast <- contrast_plot_panels(
  cr$estimands, quant = c(0, 0.995), scales = 'free_y',
  group_pal = Management_palette,
  true_vals = cr$true_estimands); p_sim_contrast

save_report("sim_summary", model_id_MDLS2_ITS_shifted, fit_sim, cr$estimands, model,
            recovery = bind_rows(fixed_recovery, sigma_recovery), model_name = "The Floor Raiser")
save_gg("sim_contrast_density", model_id_MDLS2_ITS_shifted, p_sim_contrast)

## Prior predictive check -------------------------------------------------------

n_prior <- 1000
extracted_prior <- extract.prior(fit_sim, n = n_prior)

# shift = 1 here too, so prior_pred$Dv_shifted matches what real Dv would look like.
prior_pred <- map_dfr(seq_len(n_prior), function(i){
  simulate_from_priors_MDLS2(draw_true(extracted_prior, i), shift = 1)
}, .id = "draw")

summary(prior_pred$Dv_shifted)

p_prior_pc <- prior_predictive_spaghetti(
  prior_pred, value_col = "Dv_shifted", upper_q = 0.99, model = model,
  title = "Prior predictive check", observed = dat_sim$Dv_shifted) ; p_prior_pc

save_gg("sim_prior_PC", model_id_MDLS2_ITS_shifted, p_prior_pc)

## Simulation-based calibration (SBC) --------------------------------------------

# contrast_may_gap_MDLS2() (hiermod/Models/MDLS2_model.R) needs no
# shift-awareness at all -- it's a contrast (org_May - conv_May), and the
# shift cancels. simulate_fn wraps simulate_from_priors_MDLS2() with shift = 1
# baked in, since run_sbc() only ever calls it as simulate_fn(true_params).


ncores <- 30
nchains <- 2
n_sbc = 100
n_iter = 20000

sbc_MDLS2_shifted <- run_sbc(
  model_fit   = fit_sim,
  means_fn    = means_MDLS2,
  contrast_fn = contrast_may_gap_MDLS2,
  simulate_fn = function(true_params) simulate_from_priors_MDLS2(true_params, shift = 1),
  n_sbc = n_sbc, iter = n_iter,
  n_parallel = ncores/nchains, chains = nchains, cores = nchains,
  control = list(adapt_delta = 0.99))

sbc_out_MDLS2_shifted <- save_sbc_report(sbc_MDLS2_shifted, paste0("MDLS2_shifted_",n_sbc,"iter"))
