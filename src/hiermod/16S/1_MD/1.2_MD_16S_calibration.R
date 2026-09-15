# MODEL 1 (MD, constant variance) + MODEL 2 (MDv, Management-specific
# variance), 16S: same models as ITS, on 16S's own data.
#
# Formal prior-PC/SBC added for MD and MDv specifically: SBC on every later
# model in this family (MDLS2v, MDLS2vz, MDL, MDS2) found loga[]'s
# posterior systematically shrunk toward its prior mean (rank-fraction well
# below 0.5), and that bias persisted identically after removing Location
# entirely (MDS2) -- ruling out any random-effect cardinality issue as the
# cause. MD (no per-group sigma, no random effects at all) calibrated
# cleanly (rank-fraction ~0.5, well within noise at n=100) -- ruling out a
# base-level loga-prior/lognormal-likelihood shrinkage effect too. MDv
# isolates the remaining candidate: per-group sigma[Mg], still with no
# random effects. If MDv calibrates cleanly like MD, the cause is
# specifically the random effects (Location/Tree/Year); if it's already
# biased here, it's the loga[Mg]-sigma[Mg] per-group interaction itself.
# loga/sigma priors patched to this investigation's established 16S values
# (dnorm(5,2)/dexp(2)) rather than the original, never-validated ITS import
# (dnorm(2,2)/dexp(1)) -- scoped to this script's own model_MD/model_MDv
# only, not 1.3's real-data fits.

hiermod_marker <- "16S"
source('src/hiermod/0_SETUP.R')
source('src/hiermod/Models/MD_model.R') # model_MD_16S, model_MDv_16S, sim_div_M(), sim_div_MD(), sim_div_MDv(), means_MD(), means_MDv(), simulate_from_priors_MD(), simulate_from_priors_MDv()
model_MD  <- model_MD_16S
model_MD$prior_loga  <- quote(loga[Mg] ~ dnorm(5,2))
model_MD$prior_sigma <- quote(sigma ~ dexp(2))
model_MDv <- model_MDv_16S
model_MDv$prior_loga  <- quote(loga[Mg] ~ dnorm(5,2))
model_MDv$prior_sigma <- quote(sigma[Mg] ~ dexp(2))

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

### Parameter + contrast recovery -------------------------------------------
# One simulated dataset with loga/sigma given directly (sim_div_MD(), not
# sim_div_M()'s mean_/cv_ round-trip), matching every later model's own
# convention -- lets check_recovery() test the raw parameters directly
# against their true values, not just the derived contrast.

true_conv <- 180
true_org  <- 120
true_sigma_MD <- cv_to_sigma(0.65) # single shared CV, in between MDL's own 0.5/0.8

set.seed(20260911)

dat_sim_MD <- sim_div_MD(Mg = rbern(250) + 1, loga = log(c(true_conv, true_org)), sigma = true_sigma_MD)

fit_MD_cal <- ulam(
  model_MD,
  data = as.list(dat_sim_MD),
  chains = 6, cores = 6, iter = 10000)
precis(fit_MD_cal, depth = 2)

post_MD_cal <- extract.samples(fit_MD_cal)

(param_recovery_MD <- check_recovery(
  true = list(loga1 = log(true_conv), loga2 = log(true_org), sigma = true_sigma_MD),
  post_draws = list(loga1 = post_MD_cal$loga[,1], loga2 = post_MD_cal$loga[,2], sigma = post_MD_cal$sigma)
))

### Contrast recovery -------------

true_vals_MD <- tribble(
  ~statistic, ~group,     ~value,
  "median",   "Contrast", true_org - true_conv,
  "mean",     "Contrast", lognormal_mean(log(true_org), true_sigma_MD^2) -
                           lognormal_mean(log(true_conv), true_sigma_MD^2)
)

pf_MD_cal <- post_full(fit_MD_cal, means_MD)
pc_MD_cal <- compute_contrasts(pf_MD_cal, keep = c("mean", "median"), group_levels = idx$Mg$levels)

p_MD_cal_contrast <- contrast_plot_panels(
  pc_MD_cal, quant = c(0.005, 0.99), group_pal = Management_palette,
  true_vals = true_vals_MD); p_MD_cal_contrast

save_report("sim_summary", model_id_MD, recovery = param_recovery_MD, fit_MD_cal, pc_MD_cal, model_MD, model_name = "The Bare Bones")
save_gg("sim_contrast_density", model_id_MD, p_MD_cal_contrast)

## MDv -- Allow Management-specific variance (heteroscedasticity) ============
# Same heteroscedastic-truth demonstration as before (MD's shared sigma
# biases the contrast; MDv's per-group sigma fixes it), modernized to the
# check_recovery()/compute_contrasts() pattern used everywhere else in this
# family -- replaces the old 10x for-loop/text-file approach (repeatability
# is what the formal SBC above now actually tests properly, rather than an
# ad hoc rerun-and-eyeball).

true_sigma_hetero <- cv_to_sigma(c(0.5, 0.8))

dat_sim_hetero <- sim_div_MDv(Mg = rbern(250) + 1, loga = log(c(true_conv, true_org)), sigma = true_sigma_hetero)

true_vals_hetero <- tribble(
  ~statistic, ~group,     ~value,
  "median",   "Contrast", true_org - true_conv,
  "mean",     "Contrast", lognormal_mean(log(true_org), true_sigma_hetero[2]^2) -
                           lognormal_mean(log(true_conv), true_sigma_hetero[1]^2)
)

### MD on heteroscedastic truth -- demonstrates the problem ----

fit_MD_hetero <- ulam(model_MD, data = as.list(dat_sim_hetero), chains = 6, cores = 6, iter = 5000)

pf_MD_hetero <- post_full(fit_MD_hetero, means_MD)
pc_MD_hetero <- compute_contrasts(pf_MD_hetero, keep = c("mean", "median"), group_levels = idx$Mg$levels)

p_MD_hetero_contrast <- contrast_plot_panels(
  pc_MD_hetero, quant = c(0, 1), group_pal = Management_palette, true_vals = true_vals_hetero) +
  labs(subtitle = "Shared sigma under heteroscedastic truth: contrast comes out biased"); p_MD_hetero_contrast

save_gg("sim_contrast_density", "MD_hetero", p_MD_hetero_contrast)

### MDv on the same data -- recovers correctly ----

fit_MDv_sim <- ulam(model_MDv, data = as.list(dat_sim_hetero), chains = 6, cores = 6, iter = 5000)
precis(fit_MDv_sim, depth = 2)

post_MDv_sim <- extract.samples(fit_MDv_sim)

(param_recovery_MDv <- check_recovery(
  true = list(loga1 = log(true_conv), loga2 = log(true_org),
              sigma1 = true_sigma_hetero[1], sigma2 = true_sigma_hetero[2]),
  post_draws = list(loga1 = post_MDv_sim$loga[,1], loga2 = post_MDv_sim$loga[,2],
                     sigma1 = post_MDv_sim$sigma[,1], sigma2 = post_MDv_sim$sigma[,2])
))

pf_MDv_sim <- post_full(fit_MDv_sim, means_MDv)
pc_MDv_sim <- compute_contrasts(pf_MDv_sim, keep = c("mean", "median"), group_levels = idx$Mg$levels)

p_MDv_sim_contrast <- contrast_plot_panels(
  pc_MDv_sim, quant = c(0, 1), group_pal = Management_palette, true_vals = true_vals_hetero) +
  labs(subtitle = "Group-specific variance recovers the true contrast"); p_MDv_sim_contrast

save_report("sim_summary", model_id_MDv, recovery = param_recovery_MDv, fit_MDv_sim, pc_MDv_sim, model_MDv, model_name = "The Loose Cannon")
save_gg("sim_contrast_density", model_id_MDv, p_MDv_sim_contrast, width = 8, height = 4)

## MD -- formal calibration (loga-shrinkage floor test) =======================

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
  library(rethinking); library(tidyverse) # future::multisession workers start fresh -- not auto-attached
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
  # future.chunk.size activates generate_datasets()'s built-in
  # future.apply::future_replicate() path (default is Inf = sequential).
  # future.globals must be named explicitly -- auto-detection (TRUE) can't
  # trace into this closure from inside the SBC package's own internals.
  datasets_MD <- generate_datasets(
    SBC_generator_function(
      generate_one_MD, future.chunk.size = default_chunk_size(n_sbc),
      future.globals = c("draw_true", "fit_MD_cal", "simulate_from_priors_MD", "sim_div_MD")),
    n_sbc)
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

# MD result: 0 divergences, 0 low-EBFMI, only 1 fit barely over the Rhat
# threshold at 1.011, treedepth totals negligible compared to MDS2's
# 145,523. MD is properly calibrated (loga[1]/loga[2] rank-fractions
# 0.508/0.456, both well within n=100 noise).

## MDv -- formal calibration (per-group sigma isolation test) ================
# MD (shared sigma, no random effects) calibrated cleanly above. MDL and
# MDS2 (both per-group sigma[Mg] PLUS random effects) both showed severe
# loga miscalibration. MDv isolates which of those two additions is
# responsible: per-group sigma[Mg], with STILL no random effects at all.
# Reuses fit_MDv_sim from the heteroscedasticity demo above (already fit
# under the patched (5,2)/(2) priors) rather than a fresh calibration fit.

### loga x sigma[Mg] funnel check ----------------------------------------
# Two per-group sigmas to check now, unlike MD's single shared one.

p_funnel_MDv <- function(){
  par(mfrow = c(1,2))
  plot(post_MDv_sim$sigma[,1], post_MDv_sim$loga[,1],
       xlab = "sigma[1] (Conventional)", ylab = "loga[1] (Conventional)", pch = 16, col = scales::alpha("black", 0.15))
  plot(post_MDv_sim$sigma[,2], post_MDv_sim$loga[,2],
       xlab = "sigma[2] (Organic)", ylab = "loga[2] (Organic)", pch = 16, col = scales::alpha("black", 0.15))
  par(mfrow = c(1,1))
}
save_pdf("loga_sigma_funnel", model_id_MDv, p_funnel_MDv)

### Simulation-based calibration (SBC), via the SBC package ----------------
# n_sbc/n_iter reused from MD's own SBC section above.

generate_one_MDv <- function(){
  library(rethinking); library(tidyverse) # future::multisession workers start fresh -- not auto-attached
  true_params <- suppressMessages(suppressWarnings(draw_true(extract.prior(fit_MDv_sim, n = 1, refresh = 0), 1)))[
    c("loga", "sigma")]
  dat <- simulate_from_priors_MDv(true_params)
  list(variables = true_params, generated = as.list(dat[, c("Dv", "Mg")]))
}

dq_MDv <- derived_quantities(
  median_contrast = exp(loga[2]) - exp(loga[1]),
  mean_contrast   = exp(loga[2] + sigma[2]^2/2) - exp(loga[1] + sigma[1]^2/2)
)

datasets_path_MDv <- file.path(hiermod_out_dir, "sbc_datasets_MDv.rds")
if (file.exists(datasets_path_MDv)) {
  datasets_MDv <- readRDS(datasets_path_MDv)
} else {
  datasets_MDv <- generate_datasets(
    SBC_generator_function(
      generate_one_MDv, future.chunk.size = default_chunk_size(n_sbc),
      future.globals = c("draw_true", "fit_MDv_sim", "simulate_from_priors_MDv", "sim_div_MDv")),
    n_sbc)
  saveRDS(datasets_MDv, datasets_path_MDv, compress = "xz")
}

backend_MDv <- SBC_backend_ulam(model_MDv, iter = n_iter, refresh = 0,
                                 control = list(adapt_delta = 0.99))

sbc_MDv <- compute_SBC(
  datasets_MDv, backend_MDv, dquants = dq_MDv,
  cache_mode = "results", cache_location = file.path(hiermod_out_dir, "sbc_cache_MDv"),
  globals = c("SBC_fit.SBC_backend_ulam", "SBC_fit_to_draws_matrix.ulam",
              "SBC_fit_to_diagnostics.ulam"))

p_sbc_rank_MDv  <- plot_rank_hist(sbc_MDv)
p_sbc_ecdf_MDv  <- plot_ecdf_diff(sbc_MDv)
p_sbc_cover_MDv <- plot_coverage(sbc_MDv)

sbc_step_MDv <- paste0(model_id_MDv, "_", n_sbc, "sbc_iter")
save_gg("SBC_rank_hist", sbc_step_MDv, p_sbc_rank_MDv, width = 9, height = 7)
save_gg("SBC_ecdf_diff", sbc_step_MDv, p_sbc_ecdf_MDv, width = 9, height = 7)
save_gg("SBC_coverage",  sbc_step_MDv, p_sbc_cover_MDv, width = 9, height = 7)

(sbc_diag_summary_MDv <- sbc_MDv$backend_diagnostics %>%
   dplyr::summarise(total_divergent = sum(n_divergent), total_max_treedepth = sum(n_max_treedepth),
                     total_low_ebfmi = sum(n_low_ebfmi), n_replicates = dplyr::n()))

sbc_MDv$stats |>
  dplyr::filter(variable %in% c("loga[1]", "loga[2]")) |>
  dplyr::group_by(variable) |>
  dplyr::summarise(mean_rank_frac = mean(rank / max_rank), median_rank_frac = median(rank / max_rank))