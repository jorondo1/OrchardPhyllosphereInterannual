# MODEL 1 (MD, constant variance) + MODEL 2 (MDv, Management-specific
# variance), 16S: same models as ITS, on 16S's own data.
#
# Formal prior-PC/SBC added for MD specifically: SBC on every later model in
# this family (MDLS2v, MDLS2vz, MDL, MDS2) found loga[]'s posterior
# systematically shrunk toward its prior mean (rank-fraction well below
# 0.5), and that bias persisted identically after removing Location
# entirely (MDS2) -- ruling out any random-effect cardinality issue as the
# cause. MD is the floor test: just loga[Mg] and a single SHARED sigma (no
# per-group variance at all, so there's no sigma[Mg] left to entangle with
# loga group-wise). If the bias shows up even here, it's a base-level
# loga-prior/lognormal-likelihood shrinkage effect, present from the very
# first model onward. loga/sigma priors patched to this investigation's
# established 16S values (dnorm(5,2)/dexp(2)) rather than the original,
# never-validated ITS import (dnorm(2,2)/dexp(1)) -- scoped to this script's
# own `model_MD` only, not 1.3's real-data fit.

hiermod_marker <- "16S"
source('src/hiermod/0_SETUP.R')
source('src/hiermod/Models/MD_model.R') # model_MD_16S, model_MDv_16S, sim_div_M(), sim_div_MD(), means_MD(), means_MDv(), simulate_from_priors_MD()
model_MD  <- model_MD_16S
model_MD$prior_loga  <- quote(loga[Mg] ~ dnorm(5,2))
model_MD$prior_sigma <- quote(sigma ~ dexp(2))
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
save_pdf("prior_pred_dens", model_id_MD, function(){
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
save_report("sim_summary", model_id_MD, fit_MD_sim, pc_MD_sim, model_MD, model_name = "The Bare Bones")

# it's in the vicinity
p_MD_sim_contrast <- contrast_plot_panels(pc_MD_sim, quant = c(0, 1), group_pal = Management_palette) +
  labs(x = 'Mean Hill number of order 1'); p_MD_sim_contrast
save_gg("sim_contrast_density", model_id_MD, p_MD_sim_contrast)

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

save_report("sim_summary", model_id_MDv, fit_MDv_sim, pc_MDv_sim, model_MDv, model_name = "The Loose Cannon")
save_gg("sim_contrast_density", model_id_MDv, p_MDv_sim_contrast, width = 8, height = 4)

## MD -- formal calibration (loga-shrinkage floor test) =======================
# Fresh recovery fit under the patched (5,2)/(2) priors -- fit_MD_sim above
# used model_MD before the patch, and was fit via sim_div_M()'s mean_/cv_
# convention rather than a direct loga/sigma true value, so it isn't reused
# here.

true_conv <- 180
true_org  <- 120
true_sigma_MD <- cv_to_sigma(0.65) # single shared CV, in between MDL's own 0.5/0.8

set.seed(20260911)

dat_sim_MD <- sim_div_MD(Mg = rbern(250) + 1, loga = log(c(true_conv, true_org)), sigma = true_sigma_MD)

fit_MD_cal <- ulam(
  model_MD,
  data = as.list(dat_sim_MD),
  chains = 6, cores = 6, iter = 5000)
precis(fit_MD_cal, depth = 2)

post_MD_cal <- extract.samples(fit_MD_cal)

(param_recovery_MD <- check_recovery(
  true = list(loga1 = log(true_conv), loga2 = log(true_org), sigma = true_sigma_MD),
  post_draws = list(loga1 = post_MD_cal$loga[,1], loga2 = post_MD_cal$loga[,2], sigma = post_MD_cal$sigma)
))

### loga x sigma funnel check ----------------------------------------------
# sigma is a single SHARED scalar here (no per-Mg split at all), the
# cleanest possible check of whether loga trades off against the residual
# scale even with no group-specific variance to entangle with.

p_funnel_MD <- function(){
  par(mfrow = c(1,2))
  plot(post_MD_cal$sigma, post_MD_cal$loga[,1],
       xlab = "sigma (shared)", ylab = "loga[1] (Conventional)", pch = 16, col = scales::alpha("black", 0.15))
  plot(post_MD_cal$sigma, post_MD_cal$loga[,2],
       xlab = "sigma (shared)", ylab = "loga[2] (Organic)", pch = 16, col = scales::alpha("black", 0.15))
  par(mfrow = c(1,1))
}
save_pdf("loga_sigma_funnel", model_id_MD, p_funnel_MD)

### Prior predictive check -------------------------------------------------------

n_prior <- 1000
extracted_prior_MD <- extract.prior(fit_MD_cal, n = n_prior)

prior_pred_MD <- map_dfr(seq_len(n_prior), function(i){
  simulate_from_priors_MD(draw_true(extracted_prior_MD, i))
}, .id = "draw")

summary(prior_pred_MD$Dv) # judge on median/IQR, not mean/SD

p_prior_pc_MD <- prior_predictive_spaghetti(
  prior_pred_MD, value_col = "Dv", upper_q = 0.99, model = model_MD,
  title = "Prior predictive check: loga[Mg] ~ dnorm(5,2)", observed = dat_sim_MD$Dv); p_prior_pc_MD

save_gg("sim_prior_PC", model_id_MD, p_prior_pc_MD)

### Simulation-based calibration (SBC), via the SBC package --------------------
# The floor test: does loga[] still show the ~0.15-0.35 rank-fraction skew
# found in every other model in this family, even with no random effects,
# no Season/Tree/Year, and a single shared (not per-group) sigma?
source('src/utils/sbc_backend_ulam.R')
library(SBC)
future::plan(future::multisession)

generate_one_MD <- function(){
  true_params <- suppressMessages(suppressWarnings(draw_true(extract.prior(fit_MD_cal, n = 1, refresh = 0), 1)))[
    c("loga", "sigma")]
  dat <- simulate_from_priors_MD(true_params)
  list(variables = true_params, generated = as.list(dat[, c("Dv", "Mg")]))
}

dq_MD <- derived_quantities(
  median_contrast = exp(loga[2]) - exp(loga[1]),
  mean_contrast   = exp(loga[2] + sigma^2/2) - exp(loga[1] + sigma^2/2)
)

n_sbc  <- 100
n_iter <- 10000

datasets_path_MD <- file.path(hiermod_out_dir, "sbc_datasets_MD.rds")
if (file.exists(datasets_path_MD)) {
  datasets_MD <- readRDS(datasets_path_MD)
} else {
  datasets_MD <- generate_datasets(SBC_generator_function(generate_one_MD), n_sbc)
  saveRDS(datasets_MD, datasets_path_MD, compress = "xz")
}

backend_MD <- SBC_backend_ulam(model_MD, iter = n_iter, refresh = 0,
                                control = list(adapt_delta = 0.99))

sbc_MD <- compute_SBC(
  datasets_MD, backend_MD, dquants = dq_MD,
  cache_mode = "results", cache_location = file.path(hiermod_out_dir, "sbc_cache_MD"),
  globals = c("SBC_fit.SBC_backend_ulam", "SBC_fit_to_draws_matrix.ulam",
              "SBC_fit_to_diagnostics.ulam"))

p_sbc_rank  <- plot_rank_hist(sbc_MD)
p_sbc_ecdf  <- plot_ecdf_diff(sbc_MD)
p_sbc_cover <- plot_coverage(sbc_MD)

sbc_step <- paste0(model_id_MD, "_", n_sbc, "sbc_iter")
save_gg("SBC_rank_hist", sbc_step, p_sbc_rank, width = 9, height = 7)
save_gg("SBC_ecdf_diff", sbc_step, p_sbc_ecdf, width = 9, height = 7)
save_gg("SBC_coverage",  sbc_step, p_sbc_cover, width = 9, height = 7)

(sbc_diag_summary_MD <- sbc_MD$backend_diagnostics %>%
   dplyr::summarise(total_divergent = sum(n_divergent), total_max_treedepth = sum(n_max_treedepth),
                     total_low_ebfmi = sum(n_low_ebfmi), n_replicates = dplyr::n()))

sbc_MD$stats |>
  dplyr::filter(variable %in% c("loga[1]", "loga[2]")) |>
  dplyr::group_by(variable) |>
  dplyr::summarise(mean_rank_frac = mean(rank / max_rank), median_rank_frac = median(rank / max_rank))
