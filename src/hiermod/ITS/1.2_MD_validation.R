# 1.2_MD_validation.R -- MODEL 1 (MD) + MODEL 2 (MDv): parameter recovery,
# the real fit, and PPC. No SBC/prior-PC here -- both models are simple
# enough (no pooling, no hierarchical structure) that a full calibration
# check wasn't judged necessary (see MODEL_HISTORY.md).

source('src/hiermod/ITS/0_SETUP.R')
source('src/hiermod/ITS/1.1_MD_model.R') # model_MD/model_MDv, means_MD/means_MDv, sim_div_M(), postcounts_Model1/2()

hiermod_out_dir <- "out/hiermod/ITS_1_lognormal_MD"

# Explore raw outcome distribution
hist(div$Hill_1, breaks = 30)
save_pdf("hill1_hist", "raw", function() hist(div$Hill_1, breaks = 30))

## MD -- Mean difference by Management, constant variance ====================

### Effect-size sanity check ----------------------------------------------------
# Eyeball whether the assumed group means look like plausible Hill_1 values.
# Informal (hand-picked mean_/cv_, not drawn from priors) -- see MDv's
# Prior predictive check for the real, prior-driven version. Mean diff of 4:
dat_sim_con <- sim_div_M(rep(1,100), mean_ = 8, cv_ = 0.5)
dat_sim_org <- sim_div_M(rep(1,100), mean_ = 12, cv_ = 0.5)

div_range <- c(dat_sim_con$Dv, dat_sim_org$Dv)
dens(dat_sim_con$Dv, lwd =3, xlim = c(floor(min(div_range)),2+ceiling(max(div_range))))
dens(dat_sim_org$Dv, lwd = 3, col =2, add = TRUE)
save_pdf("prior_pred_dens", "MD", function(){
  dens(dat_sim_con$Dv, lwd = 3, xlim = c(floor(min(div_range)), 2+ceiling(max(div_range))))
  dens(dat_sim_org$Dv, lwd = 3, col = 2, add = TRUE)
})

### Parameter recovery -----------------------------------------------------------

dat_sim <- sim_div_M(
  rbern(200)+1,
  mean_ = c(8,12),  # Difference of 4 in mean
  cv_ = c(0.5,0.5))

dat_MD_sim <- list(
  Dv = dat_sim$Dv,
  Mg = dat_sim$Mg
)

fit_MD_sim <- ulam(
  model_MD,
  data = dat_MD_sim,
  chains = 6, cores = 6, iter = 10000
)
precis(fit_MD_sim, depth = 2)

post_counts_MD_sim <- postcounts_Model1(fit_MD_sim)
save_report("sim_summary", "MD", fit_MD_sim, post_counts_MD_sim, model_MD, model_name = "The Bare Bones")

# it's in the vicinity
p_MD_sim_contrast <- plot_contrast_density(post_counts_MD_sim, group_name = 'Posterior mean'); p_MD_sim_contrast
save_gg("sim_contrast_density", "MD", p_MD_sim_contrast)

### Model fit ----------------------------------------------------------------

dat_MD <- list(
  Dv = div$Hill_1,
  Mg = idx$Mg$to_index(div$Management)
)

fit_MD <- ulam(
  model_MD,
  data = dat_MD,
  chains = 6, cores = 6, iter = 10000
)
save_fit("fit", "MD", fit_MD)

precis(fit_MD, depth = 2)
traceplot(fit_MD)
save_pdf("fit_traceplot", "MD", function() traceplot(fit_MD))

### Posterior predictive check --------------------------------------------------

p_postpred_MD <- plot_ppc_overlay(fit_MD, dat_MD, idx$Mg$to_label(dat_MD$Mg))
save_gg("postpred_density", "MD", p_postpred_MD)

## MDv -- Allow Management-specific variance (heteroscedasticity) ============

# Demonstrate the problem: fitting a shared-sigma model when the groups
# actually have unequal variance.

# iterate to look at mean contrast
output_iter <- file.path(hiermod_out_dir, "unequal_variance_iters.txt"); for (i in 1:10) {

  if (i == 1) {
    cat("Iteration\tMean\tMedian\tPI\n", file = output_iter, append = FALSE)
  }
  dat_sim <- sim_div_M(
    rbern(250)+1,
    mean_ = c(8,12),  # same
    cv_ = c(0.5,0.8)) # <- here

  fit_MD_sim <- ulam(
    model_MD,
    data = list(
      Dv = dat_sim$Dv,
      Mg = dat_sim$Mg
    ), chains = 6, cores = 6, iter = 10000)

  post_counts <- postcounts_Model1(fit_MD_sim)

  # Write data row with tabs
  cat(i, "\t",
      mean(post_counts$contrast), "\t",
      median(post_counts$contrast), "\t",
      PI(post_counts$contrast), "\n",
      file = output_iter, append = TRUE)

} # Rarely reaches 4, range 2-4

### Parameter recovery -----------------------------------------------------------

output_iter <- file.path(hiermod_out_dir, "unequal_variance_iters_fixed.txt"); for (i in 1:10) {

  if (i == 1) {
    cat("Iteration\tMean\tMedian\tPI\n", file = output_iter, append = FALSE)
  }

  dat_sim <- sim_div_M(
    rbern(250)+1,
    mean_ = c(8,12),
    cv_ = c(0.5,0.8))

  fit_MD_sim <- ulam(
    model_MDv,
    data = list(
      Dv = dat_sim$Dv,
      Mg = dat_sim$Mg
    ),chains = 2, cores = 2, iter = 1000)

  post_counts <- postcounts_Model2(fit_MD_sim)

  cat(i, "\t",
      mean(post_counts$contrast), "\t",
      median(post_counts$contrast), "\t",
      PI(post_counts$contrast), "\n",
      file = output_iter, append = TRUE)

}  # contrasts hover around 4, range 3-5

fit_MDv_sim <- ulam(
  model_MDv,
  data = list(
    Dv = dat_sim$Dv,
    Mg = dat_sim$Mg
  ), chains = 6, cores = 6, iter = 10000)

precis(fit_MDv_sim, depth = 2 )

post_counts_MDv_sim <- postcounts_Model2(fit_MDv_sim)

# Plot contrast
p_MDv_sim_contrast <- plot_contrast_density(post_counts_MDv_sim, group_name = 'Posterior mean') +
  labs(
    subtitle = "Here, allowing group-specific variances allows the recovery of the true contrast.",
    x = 'Mean Hill number of order 1') +
  theme(legend.position = c(0.8, 0.8)); p_MDv_sim_contrast

save_report("sim_summary", "MDv", fit_MDv_sim, post_counts_MDv_sim, model_MDv, model_name = "The Loose Cannon")
save_gg("sim_contrast_density", "MDv", p_MDv_sim_contrast, width = 8, height = 4)

### Prior predictive check -------------------------------------------------------
# TODO: not yet done for this model. Same pattern as MDLb (prior_predictive_spaghetti(), predictive_checks.R).

### Simulation-based calibration (SBC) --------------------------------------------
# TODO: not yet done for this model. See MDLb (2.2_MDL_validation.R) for the run_sbc()/save_sbc_report() pattern.

### Model fit ----------------------------------------------------------------

fit_MDv <- ulam(
  model_MDv,
  data = dat_MD,
  chains = 6, cores = 6, iter = 10000
)
save_fit("fit", "MDv", fit_MDv)

precis(fit_MDv, depth = 2)
traceplot(fit_MDv)
save_pdf("fit_traceplot", "MDv", function() traceplot(fit_MDv))

### Posterior predictive check --------------------------------------------------
# Barely affected relative to MD, but the tail is not as heavy.

p_postpred_MDv <- plot_ppc_overlay(fit_MDv, dat_MD, idx$Mg$to_label(dat_MD$Mg))
save_gg("postpred_density", "MDv", p_postpred_MDv)
