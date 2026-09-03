# 6.2_MDLSY_validation.R -- MODEL 6 (MDLSY): parameter recovery,
# prior-predictive check, SBC, the real fit, and PPC.

source('~/Repos/orchardPhyllosphere2/src/hiermod/ITS/0_SETUP.R')
source('~/Repos/orchardPhyllosphere2/src/hiermod/ITS/6.1_MDLSY_model.R') # model, means_MDLSY(), variance_partition_MDLSY(), sim_div_MDLSY(), contrast_may_gap_MDLSY(), simulate_from_priors()
hiermod_out_dir <- "out/hiermod/ITS_6_lognormal_MDLSY"

## Model specification ---------------------------------------------------------
# Same structure as MODEL 5, plus: yr[Yr] pooled (yr[Yr]*sigma_yr, was a
# fixed dnorm(0,1) effect); Cultivar (cv[Cv]) as a simple fixed/unpooled
# effect, expected near-zero but shown explicitly rather than omitted;
# deg_h_z/precip_72h_z (degree-hours the day before / precipitation 3 days
# before, both standardized in 0_SETUP.R) as fixed-effect slopes -- these
# are meant to CONTROL for weather, not describe its effect for its own
# sake, so the headline estimand (May/July gap, seasonal change) is reported
# net of weather (deg_h_z=precip_72h_z=0, i.e. this sample's own average
# weather) -- see means_MDLSY() in 6.1_MDLSY_model.R.
#
# Real, confirmed collinearity, not fixed by this model: deg_h correlates
# -0.73 with Season and varies sharply by Year (940->1060->1559 across
# 2022-2024); precip_72h correlates -0.55 with Season, also varies by Year.
# b_deg/b_precip will be entangled with s_conv/gap_shift/yr[Yr]*sigma_yr in
# the real-data posterior -- expected, not a bug (see the collinearity
# assessment in TODO.md's Discussion notes for the full reasoning: this
# widens uncertainty and correlates parameter estimates, it does not bias
# them, and the SBC below deliberately does NOT reproduce this correlation
# in synthetic data -- it validates that the model CAN recover these
# parameters in principle, not how identifiable they are under the real
# design's collinearity specifically).

## Parameter recovery -----------------------------------------------------------

# Same baseline diversity/gap/sigma[cell] values as MODEL 5, plus modest,
# testable weather slopes and a near-zero Cultivar vector (matches "we're
# pretty sure they don't matter" -- the point is confirming the model finds
# that too, not picking a value it's guaranteed to recover).

may_conv <- 4
may_org <- 7
july_conv_shift <- 0.35
july_org_shift <- 0.2

true_sigma    <- cv_to_sigma(c(0.6, 0.25, 0.8, 0.8)) # conv_May, conv_July, org_May, org_July
true_sigma_yr <- 0.3
true_b_deg    <- 0.15
true_b_precip <- -0.15
true_cv       <- c(0.05, -0.05, 0.1, -0.1, 0)

dat_sim <- sim_div_MDLSY(
  N_samples = 240,
  n_loc = 4,
  loga = log(c(may_conv, may_org)),
  s_conv = july_conv_shift,
  gap_shift = july_org_shift,
  sigma = true_sigma,
  sigma_yr = true_sigma_yr,
  b_deg = true_b_deg,
  b_precip = true_b_precip,
  cv = true_cv,
  p_dropout = 0.1
); head(dat_sim)

hist(dat_sim$Dv)

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
           b_precip = true_b_precip),
  post_draws = list(
    loga1 = post_sim$loga[,1],
    loga2 = post_sim$loga[,2],
    s_conv = post_sim$s_conv,
    gap_shift = post_sim$gap_shift,
    b_deg = post_sim$b_deg,
    b_precip = post_sim$b_precip)
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

### Contrast recovery ------------------------------

pf_sim <- post_full(fit_sim, means_MDLSY)
m_sim  <- pf_sim$median

pc_estimands_sim <- estimand_rows(list(
  "Median May gap (Organic - Conventional)"  = m_sim$median_3 - m_sim$median_1,
  "Median July gap (Organic - Conventional)" = m_sim$median_4 - m_sim$median_2,
  "Seasonal change in median gap (July - May)"= (m_sim$median_4 - m_sim$median_2) - (m_sim$median_3 - m_sim$median_1)
))

may_gap <- may_org-may_conv
july_gap <- may_org*exp(july_conv_shift+july_org_shift)-may_conv*exp(july_conv_shift)

true_estimands <- tribble(
  ~statistic,                                   ~value,
  "Median May gap (Organic - Conventional)",    may_gap,
  "Median July gap (Organic - Conventional)",   july_gap,
  "Seasonal change in median gap (July - May)", july_gap-may_gap,
)

p_sim_contrast <- contrast_plot_panels(
  pc_estimands_sim, quant = c(0, 0.995), scales = 'free_y',
  group_pal = Management_palette,
  true_vals = true_estimands); p_sim_contrast

save_report("sim_summary", "MDLSY", fit_sim, pc_estimands_sim, model,
            recovery = bind_rows(fixed_recovery, sigma_recovery, cv_recovery),
            model_name = "The Weatherman")
save_gg("sim_contrast_density", "MDLSY", p_sim_contrast)

## Prior predictive check -------------------------------------------------------

n_prior <- 1000
extracted_prior <- extract.prior(fit_sim, n = n_prior)

prior_pred <- map_dfr(seq_len(n_prior), function(i){
  simulate_from_priors(draw_true(extracted_prior, i))
}, .id = "draw")

summary(prior_pred$Dv)

p_prior_pc <- prior_predictive_spaghetti(
  prior_pred, upper_q = 0.99, model = model,
  title = "Prior predictive check", observed = dat_sim$Dv) ; p_prior_pc

save_gg("sim_prior_PC", "MDLSY", p_prior_pc)

## Variance budget calibration ---------------------------------------------------

# K genuinely grows here: total_var in means_MDLSY() sums sigma[cell]+
# sigma_loc+sigma_tr+sigma_yr = K=4, up from Model 5's K=3 -- this is where
# scale_dexp_rate() (hiermod_core.R; R2D2M2-style variance-decomposition
# priors, full reasoning + references there) actually changes something,
# not just confirms no-op like 5.2/5b.2 did. K_ref=3 rates are those two
# scripts' own already-prior-predictive-checked ones (sigma[cell]~dexp(3),
# sigma_loc/sigma_tr~dexp(2)); sigma_yr is new here, given no unscaled
# rate of its own to inherit from, so it starts from the same base rate as
# its structural peers (sigma_loc/sigma_tr, both non-centered population-
# level SDs: 2) before the same scaling is applied.
#
# NOT targeting dexp(1) -- that was the original, admittedly-too-loose
# default from Model 2's own early exploration (see 2.2_MDL_validation.R),
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

# Own object, not a mutation of `model` -- same convention as MDLb's
# model_ppc1 (2.2_MDL_validation.R): keeps the original and the calibrated
# version both inspectable, and makes explicit which one the real fit below
# actually uses (model_vbc, not model).
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
# looks implausible -- same pattern as 2.2_MDL_validation.R's model_ppc1
# iteration, not repeated automatically since "how much is enough" is a
# judgment call, not something to auto-loop.

# Recovery re-check -- the overfitting guard itself.
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

# contrast_may_gap_MDLSY() (6.1_MDLSY_model.R) -- same reasoning as models
# 4/5's: run_sbc()'s default contrast_fn assumes a 2-column `mean`.
# simulate_from_priors() draws deg_h_z/precip_72h_z as independent
# rnorm(0,1) -- this validates recoverability in principle, not calibration
# under the real deg_h/precip_72h ~ Season/Year correlation (see the
# collinearity note above and TODO.md).

sbc_MDLSY <- run_sbc(
  model_fit   = fit_cal,
  means_fn    = means_MDLSY,
  contrast_fn = contrast_may_gap_MDLSY,
  simulate_fn = simulate_from_priors,
  n_sbc = 30, iter = 15000, n_parallel = 4, chains = 2, cores = 2,
  control = list(adapt_delta = 0.99))

(sbc_out_MDLSY <- summarize_sbc(sbc_MDLSY))
save_sbc_report(sbc_out_MDLSY, "MDLSY_30sbc_iter")
hist(sbc_out_MDLSY$ranks, breaks = 30)

## Model fit ----------------------------------------------------------------

dat <- list(
  Dv = div$Hill_1,
  Mg = idx$Mg$to_index(div$Management),
  Lo = idx$Lo$to_index(div$Location),
  Mo = idx$Mo$to_index(div$Time),
  Yr = idx$Yr$to_index(div$Year),
  Tr = idx$Tr$to_index(div$Tree_id),
  Cv = idx$Cv$to_index(div$Cultivar),
  deg_h_z = div$deg_h_z,
  precip_72h_z = div$precip_72h_z
)
dat$cell <- (dat$Mg - 1) * 2 + dat$Mo

fit_MDLSY <- ulam(
  model_vbc,
  data = dat,
  chains = 6, cores = 6, iter = 20000,
  control = list(adapt_delta = 0.99)
)
save_fit("fit", "MDLSY", fit_MDLSY)

precis(fit_MDLSY, depth = 2)

save_pdf("fit_traceplot", "MDLSY", function() traceplot(fit_MDLSY, n_cols = 6, max_rows = 10))
save_pdf("fit_trankplot", "MDLSY", function() trankplot(fit_MDLSY, n_cols = 6, max_rows = 10))

# Check the confirmed real-data collinearity's effect on the posterior
# directly: correlated draws are the expected symptom (wider, correlated
# estimates), not a red flag on their own -- see TODO.md if these turn out
# to be extreme.
post_MDLSY <- extract.samples(fit_MDLSY)
cor(post_MDLSY$b_deg, post_MDLSY$s_conv)
cor(post_MDLSY$b_deg, post_MDLSY$gap_shift)
cor(post_MDLSY$b_precip, post_MDLSY$s_conv)
cor(post_MDLSY$b_deg, post_MDLSY$yr[,1])   # deg_h varied sharply by Year (see above) -- 2022's yr[1] is the one to watch

## Posterior predictive check --------------------------------------------------

### Overall, by Management x Season cell ----

pp_group <- interaction(idx$Mg$to_label(dat$Mg), idx$Mo$to_label(dat$Mo), sep = " ")
(p_postpred <- plot_ppc_overlay(fit_MDLSY, dat, pp_group, xlim = c(0,150)))
save_gg("postpred_density", "MDLSY", p_postpred)

### Contrast test statistics ----

(p_ppc <- plot_ppc_season_contrast_stats(fit_MDLSY, dat))
save_gg("postpred_stat", "MDLSY", p_ppc)
