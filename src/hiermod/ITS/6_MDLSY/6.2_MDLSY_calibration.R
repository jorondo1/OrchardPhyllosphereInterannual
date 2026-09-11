# MODEL 6 (MDLSY): Year pooled into a non-centered random effect (was
# fixed) and Cultivar added as a fixed effect, on top of Model 5. Parameter
# recovery, prior-predictive check, variance-budget calibration (VBC), SBC.

hiermod_marker <- "ITS"
source('src/hiermod/0_SETUP.R')
source('src/hiermod/Models/MDLSY_model.R') # model_MDLSY_ITS, means_MDLSY(), variance_partition_MDLSY(), sim_div_MDLSY(), contrast_may_gap_MDLSY(), simulate_from_priors_MDLSY()
model <- model_MDLSY_ITS

hiermod_out_dir <- "out/hiermod/ITS_6_lognormal_MDLSY"

## Parameter recovery -----------------------------------------------------------

# Same baseline diversity/gap/sigma[cell] values as MODEL 5, plus a
# near-zero Cultivar vector (matches "we're pretty sure they don't matter"
# -- the point is confirming the model finds that too, not picking a value
# it's guaranteed to recover).

may_conv <- 4
may_org <- 7
july_conv_shift <- 0.35
july_org_shift <- 0.2

true_sigma    <- cv_to_sigma(c(0.6, 0.25, 0.8, 0.8)) # conv_May, conv_July, org_May, org_July
true_sigma_yr <- 0.3
true_cv       <- c(0.05, -0.05, 0.1, -0.1, 0)

dat_sim <- sim_div_MDLSY(
  N_samples = 240,
  n_loc = 4,
  loga = log(c(may_conv, may_org)),
  s_conv = july_conv_shift,
  gap_shift = july_org_shift,
  sigma = true_sigma,
  sigma_yr = true_sigma_yr,
  cv = true_cv,
  p_dropout = 0.1,
  shift = 1
); head(dat_sim)

# Dv is what the model actually fits; Dv_shifted (= 1 + Dv) is only for
# sanity-checking that the floor looks right; it's NOT passed to ulam() --
# same convention as 5b.2_MDLS2_shifted_calibration.R.
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
           gap_shift = july_org_shift),
  post_draws = list(
    loga1 = post_sim$loga[,1],
    loga2 = post_sim$loga[,2],
    s_conv = post_sim$s_conv,
    gap_shift = post_sim$gap_shift)
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
); cv_recovery  #All good

### Contrast recovery ------------------------------

cr <- contrast_recovery(fit_sim, means_MDLSY, may_conv, may_org, july_conv_shift, july_org_shift, shift = 1)

p_sim_contrast <- contrast_plot_panels(
  cr$estimands, quant = c(0, 0.995), scales = 'free_y',
  group_pal = Management_palette,
  true_vals = cr$true_estimands); p_sim_contrast

save_report("sim_summary", "MDLSY", fit_sim, cr$estimands, model,
            recovery = bind_rows(fixed_recovery, sigma_recovery, cv_recovery),
            model_name = "The Varietal")
save_gg("sim_contrast_density", "MDLSY", p_sim_contrast)

## Prior predictive check -------------------------------------------------------

n_prior <- 1000
extracted_prior <- extract.prior(fit_sim, n = n_prior)

prior_pred <- map_dfr(seq_len(n_prior), function(i){
  simulate_from_priors_MDLSY(draw_true(extracted_prior, i), shift = 1)
}, .id = "draw")

summary(prior_pred$Dv_shifted)

p_prior_pc <- prior_predictive_spaghetti(
  prior_pred, value_col = "Dv_shifted", upper_q = 0.99, model = model,
  title = "Prior predictive check", observed = dat_sim$Dv_shifted) ; p_prior_pc

save_gg("sim_prior_PC", "MDLSY", p_prior_pc)

## Variance budget calibration ---------------------------------------------------

# total_var in means_MDLSY() sums sigma[cell]+ sigma_loc+sigma_tr+sigma_yr =
# so K=4, up from Model 5's K=3. Implementing R2D2M2-style variance-decomposition
# priors. sigma_yr is new here, given no unscaled rate of its own to inherit from,
# so it starts from the same base rate as its structural peers (sigma_loc/sigma_tr,
# both non-centered population-level SDs.
#
# NOT targeting dexp(1) -- that was the original, admittedly-too-loose
# default from Model 2's own early exploration (see 2.2_MDL_calibration.R),
# already superseded by dexp(2)/dexp(3) once tightened via a real prior-
# predictive check. This anchors to THAT already-validated K=3 reference,
# not back to the original loose one.

K_ref <- 3
K_new <- 4

rate_sigma     <- scale_dexp_rate(3, K_ref, K_new)  # sigma[cell]
rate_sigma_loc <- scale_dexp_rate(2, K_ref, K_new)  # sigma_loc
rate_sigma_tr  <- scale_dexp_rate(2, K_ref, K_new)  # sigma_tr
rate_sigma_yr  <- scale_dexp_rate(2, K_ref, K_new)  # sigma_yr -- new, inherits sigma_loc/sigma_tr's base rate
c(sigma = rate_sigma, sigma_loc = rate_sigma_loc, sigma_tr = rate_sigma_tr, sigma_yr = rate_sigma_yr)
# (~3.46, ~2.31, ~2.31, ~2.31) -- every rate goes up by sqrt(4/3) =~ 1.15,
# holding the expected total variance close to what Model 5 already
# validated instead of letting the 4th summed term inflate it further.

# model_vbc: keeps the original and the calibrated version both inspectable,
# and makes explicit which one the real fit below actually uses (model_vbc,
# not model).
model_vbc <- model
model_vbc$pr_sigma     <- bquote(sigma[cell] ~ dexp(.(rate_sigma)))
model_vbc$pr_sigma_loc <- bquote(sigma_loc   ~ dexp(.(rate_sigma_loc)))
model_vbc$pr_sigma_tr  <- bquote(sigma_tr    ~ dexp(.(rate_sigma_tr)))
model_vbc$pr_sigma_yr  <- bquote(sigma_yr    ~ dexp(.(rate_sigma_yr)))

# Overfitting guard: scale_dexp_rate() never looks at any simulated or real
# data -- it's a closed-form recalculation from K and the already
# prior-predictive-validated K=3 rates, so no data-dependent tuning risk in
# that step itself. What actually guards against over-tightening is
# downstream: SBC's rank-uniformity test below (now pointed at fit_cal) and
# the recovery re-check right after fit_cal (a real red flag would be
# covered=FALSE showing up here that wasn't there for the Parameter
# recovery section above).
#
# Refit with the calibrated priors -- run_sbc() pulls its model spec from
# model_fit@formula directly, not the live `model_vbc` variable, so SBC/the
# real fit below only see this patch if they're pointed at a fit made
# *after* it.
fit_cal <- ulam(
  model_vbc,
  data = as.list(dat_sim),
  chains = 6, cores = 6, iter = 10000,
  control = list(adapt_delta = 0.99)
)
precis(fit_cal, depth = 2)
# Re-run the prior predictive check against fit_cal here if the tail still
# looks implausible -- same pattern as 2.2_MDL_calibration.R's model_ppc1
# iteration, not repeated automatically since "how much is enough" is a
# judgment call, not something to auto-loop.

## 2nd Prior predictive check -------------------------------------------------------
extracted_prior_cal <- extract.prior(fit_cal, n = n_prior)

prior_pred_cal <- map_dfr(seq_len(n_prior), function(i){
  simulate_from_priors_MDLSY(draw_true(extracted_prior_cal, i), shift = 1)
}, .id = "draw")

summary(prior_pred_cal$Dv_shifted)

p_prior_pc_cal <- prior_predictive_spaghetti(
  prior_pred_cal, value_col = "Dv_shifted", upper_q = 0.99, model = model_vbc,
  title = "Prior predictive check", observed = dat_sim$Dv_shifted)

save_gg("sim_prior_PC", "MDLSY_VBCal", p_prior_pc_cal)

# If "extremeness" was spread evenly across replicates (~250 obs each),
# expect roughly 1-(1-0.01)^250= 92% of replicates to touch it just from a
# trivial 1% per-point tail. At 36%, we nowhere near: tail still coming from
# a minority of "unlucky" prior draws, not a universal property of the prior.
# It's worse than Model 5's 19.7%, but still plenty ok. Both mean and sd
# decreased a lot, that's what vbc does.

# Recovery re-check :
post_cal <- extract.samples(fit_cal)
sigma_recovery_cal <- check_recovery(
  true = list(sigma1 = true_sigma[1], sigma2 = true_sigma[2],
              sigma3 = true_sigma[3], sigma4 = true_sigma[4],
              sigma_yr = true_sigma_yr),
  post_draws = list(sigma1 = post_cal$sigma[,1], sigma2 = post_cal$sigma[,2],
                    sigma3 = post_cal$sigma[,3], sigma4 = post_cal$sigma[,4],
                    sigma_yr = post_cal$sigma_yr)
); sigma_recovery_cal

## Simulation-based calibration (SBC) --------------------------------------------

# contrast_may_gap_MDLSY() (hiermod/Models/MDLSY_model.R) -- same reasoning
# as models 4/5's: run_sbc()'s default contrast_fn assumes a 2-column `mean`.
# contrast_may_gap_MDLSY() needs no shift-awareness -- it's a contrast
# (org_May - conv_May), and the shift cancels, same as
# 5b.2_MDLS2_shifted_calibration.R.

ncores <- 18
nchains <- 2
n_sbc = 100
n_iter = 10000

sbc_MDLSY <- run_sbc(
  model_fit   = fit_cal,
  means_fn    = means_MDLSY,
  contrast_fn = contrast_may_gap_MDLSY,
  simulate_fn = function(true_params) simulate_from_priors_MDLSY(true_params, shift = 1),
  n_sbc = n_sbc, iter = n_iter,
  n_parallel = ncores/nchains, chains = nchains, cores = nchains,
  control = list(adapt_delta = 0.99))

sbc_out_MDLSY <- save_sbc_report(sbc_MDLSY, paste0("MDLSY_",n_sbc,"iter"))
hist(sbc_out_MDLSY$ranks, breaks = 30)

saveRDS(model_vbc, file.path(hiermod_out_dir, "model_vbc.rds"))
