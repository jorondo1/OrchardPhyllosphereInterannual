# MODEL 1 (MD, constant variance) + MODEL 2 (MDv, Management-specific
# variance), 16S: same models as ITS, on 16S's own data. No formal SBC/
# prior-PC here -- both models are simple enough (no pooling, no
# hierarchical structure) that it wasn't judged necessary.

hiermod_marker <- "16S"
source('src/hiermod/0_SETUP.R')
source('src/hiermod/Models/MD_model.R') # model_MD_16S, model_MDv_16S, sim_div_M(), means_MD(), means_MDv()
model_MD  <- model_MD_16S
model_MDv <- model_MDv_16S

hiermod_out_dir <- "out/hiermod/16S_1_lognormal_MD"

# Explore raw outcome distribution
hist(div$Hill_1, breaks = 30)
mean(div$Hill_1)
save_pdf("hill1_hist", "raw", function() hist(div$Hill_1, breaks = 30))

## MD -- Mean difference by Management, constant variance ====================

### Effect-size sanity check ----------------------------------------------------
# Eyeball whether the assumed group means look like plausible Hill_1 values.
# Informal (hand-picked mean_/cv_, not drawn from priors) -- see MDv's
# Prior predictive check for the real, prior-driven version. Mean diff of 4:
dat_sim_con <- sim_div_M(rep(1,100), mean_ = 180, cv_ = 1)
dat_sim_org <- sim_div_M(rep(1,100), mean_ = 120, cv_ = 1)

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
  mean_ = c(180,120),  # Difference of 60 in mean
  cv_ = c(0.8,0.8))

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

pf_MD_sim <- post_full(fit_MD_sim, means_MD)
pc_MD_sim <- compute_contrasts(pf_MD_sim, keep = "mean", group_levels = idx$Mg$levels)
save_report("sim_summary", "MD", fit_MD_sim, pc_MD_sim, model_MD, model_name = "The Bare Bones")

# it's in the vicinity
p_MD_sim_contrast <- contrast_plot_panels(pc_MD_sim, quant = c(0, 1), group_pal = Management_palette) +
  labs(x = 'Mean Hill number of order 1'); p_MD_sim_contrast
save_gg("sim_contrast_density", "MD", p_MD_sim_contrast)

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
    mean_ = c(180,120),  # same
    cv_ = c(0.5,0.8)) # <- here

  fit_MD_sim <- ulam(
    model_MD,
    data = list(
      Dv = dat_sim$Dv,
      Mg = dat_sim$Mg
    ), chains = 6, cores = 6, iter = 5000)

  pf_iter <- post_full(fit_MD_sim, means_MD)
  contrast_iter <- pf_iter$mean$mean_2 - pf_iter$mean$mean_1

  # Write data row with tabs
  cat(i, "\t",
      mean(contrast_iter), "\t",
      median(contrast_iter), "\t",
      PI(contrast_iter), "\n",
      file = output_iter, append = TRUE)

} # Rarely reaches -60; usually in the -70 to -90 range

### Parameter recovery -----------------------------------------------------------

output_iter <- file.path(hiermod_out_dir, "unequal_variance_iters_fixed.txt"); for (i in 1:10) {

  if (i == 1) {
    cat("Iteration\tMean\tMedian\tPI\n", file = output_iter, append = FALSE)
  }

  dat_sim <- sim_div_M(
    rbern(250)+1,
    mean_ = c(180,120),
    cv_ = c(0.5,0.8))

  fit_MD_sim <- ulam(
    model_MDv,
    data = list(
      Dv = dat_sim$Dv,
      Mg = dat_sim$Mg
    ),chains = 2, cores = 2, iter = 1000)

  pf_iter <- post_full(fit_MD_sim, means_MDv)
  contrast_iter <- pf_iter$mean$mean_2 - pf_iter$mean$mean_1

  cat(i, "\t",
      mean(contrast_iter), "\t",
      median(contrast_iter), "\t",
      PI(contrast_iter), "\n",
      file = output_iter, append = TRUE)

}  # contrasts hover around 4, range 3-5

fit_MDv_sim <- ulam(
  model_MDv,
  data = list(
    Dv = dat_sim$Dv,
    Mg = dat_sim$Mg
  ), chains = 6, cores = 6, iter = 10000)

precis(fit_MDv_sim, depth = 2 )

pf_MDv_sim <- post_full(fit_MDv_sim, means_MDv)
pc_MDv_sim <- compute_contrasts(pf_MDv_sim, keep = "mean", group_levels = idx$Mg$levels)

# Plot contrast
p_MDv_sim_contrast <- contrast_plot_panels(pc_MDv_sim, quant = c(0, 1), group_pal = Management_palette) +
  labs(
    subtitle = "Here, allowing group-specific variances allows the recovery of the true contrast.",
    x = 'Mean Hill number of order 1'); p_MDv_sim_contrast

save_report("sim_summary", "MDv", fit_MDv_sim, pc_MDv_sim, model_MDv, model_name = "The Loose Cannon")
save_gg("sim_contrast_density", "MDv", p_MDv_sim_contrast, width = 8, height = 4)

### Prior predictive check -------------------------------------------------------
# TODO: not yet done for this model. Same pattern as MDLb (prior_predictive_spaghetti(), predictive_checks.R).

### Simulation-based calibration (SBC) --------------------------------------------
# TODO: not yet done for this model. See MDLb (2.3_MDL_fit.R) for the run_sbc()/save_sbc_report() pattern.
