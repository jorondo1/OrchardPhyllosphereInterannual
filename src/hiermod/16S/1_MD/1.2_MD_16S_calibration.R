# MODEL 1 (MD, constant variance), 16S: same model as ITS, on 16S's own
# data. MDv (Management-specific variance) is a separate sibling script
# (1.2_MDv_16S_calibration.R) -- split out since each now has its own full
# recovery/prior-PC/SBC investigation rather than sharing one file.
#
# Formal prior-PC/SBC added here specifically as the floor test in a wider
# investigation: SBC on every later model in this family (MDLS2v, MDLS2vz,
# MDL, MDS2) found loga[]'s posterior systematically shrunk toward its
# prior mean (rank-fraction well below 0.5), and that bias persisted
# identically after removing Location entirely (MDS2) -- ruling out any
# random-effect cardinality issue as the cause. MD (no per-group sigma, no
# random effects at all) is the simplest possible version of the shared
# loga[Mg]/lognormal-likelihood structure; if it calibrates cleanly, a
# base-level prior/likelihood shrinkage effect is ruled out too, pointing
# the remaining investigation at per-group sigma[Mg] and/or the random
# effects specifically (see 1.2_MDv_16S_calibration.R for the next isolating
# step). loga/sigma priors patched to this investigation's established 16S
# values (dnorm(5,2)/dexp(2)) rather than the original, never-validated ITS
# import (dnorm(2,2)/dexp(1)) -- scoped to this script's own model_MD only,
# not 1.3's real-data fit.

hiermod_marker <- "16S"
source('src/hiermod/0_SETUP.R')
source('src/hiermod/Models/MD_model.R') 
source('src/utils/sbc_workflow.R') 
model_MD <- model_MD_16S
model_MD$prior_loga  <- quote(loga[Mg] ~ dnorm(5,2))
model_MD$prior_sigma <- quote(sigma ~ dexp(2))

hiermod_out_dir <- "out/hiermod/16S_1_lognormal_MD/Calibration"

save_pdf("hill1_hist", "raw", function() hist(div$Hill_1, breaks = 30))

## MD -- Mean difference by Management, constant variance ====================

### Effect-size sanity check ----------------------------------------------------
# Eyeball whether the assumed group means look like plausible Hill_1 values.
# Informal (hand-picked mean_/cv_, not drawn from priors) -- see the Prior
# predictive check below for the real, prior-driven version. Mean diff of 4:
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

save_report("sim_summary", model_id_MD, recovery = param_recovery_MD, fit_MD_cal, pc_MD_cal, model_MD)
save_gg("sim_contrast_density", model_id_MD, p_MD_cal_contrast)

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
sbc_gen_MD <- make_sbc_generator(
  fit = fit_MD_cal, simulate_fn = simulate_from_priors_MD,
  keep = c("loga", "sigma"), gen_cols = c("Dv", "Mg"),
  extra_globals = "sim_div_MD")

n_sbc  <- 100
n_iter <- 10000

sbc_MD <- run_sbc_pipeline(
  generator = sbc_gen_MD$generator, globals = sbc_gen_MD$globals,
  n_sbc = n_sbc, model = model_MD, model_id = model_id_MD, n_iter = n_iter,
  hiermod_out_dir = hiermod_out_dir, dquants = dq_MD,
  control = list(adapt_delta = 0.99))

plot_sbc_diagnostics(sbc_MD, model_id_MD, n_sbc)

sbc_MD$stats |>
  dplyr::filter(variable %in% c("loga[1]", "loga[2]")) |>
  dplyr::group_by(variable) |>
  dplyr::summarise(mean_rank_frac = mean(rank / max_rank), median_rank_frac = median(rank / max_rank))

save_sbc_health_report(model_id_MD, sbc_MD, n_sbc, n_iter,
                        variables = c("loga[1]", "loga[2]", "sigma", "median_contrast", "mean_contrast"),
                        hiermod_out_dir = hiermod_out_dir)

# MD result: 0 divergences, 0 low-EBFMI, only 1 fit barely over the Rhat
# threshold at 1.011, treedepth totals negligible compared to MDS2's
# 145,523. MD is properly calibrated (loga[1]/loga[2] rank-fractions
# 0.508/0.456, both well within n=100 noise).
