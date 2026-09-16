# MODEL 1 (MDv, Management-specific variance), ITS: mirrors
# 1.2_MDv_16S_calibration.R's own role, on ITS (Fungi) data. Split from MD
# into its own sibling script, same convention as 16S -- each gets its own
# full recovery/prior-PC/SBC investigation.
#
# model_MDv_ITS's file-level prior (sigma[Mg] ~ dexp(1)) is tested here
# AS-IS, not pre-patched to dhalfnorm(0,1) -- the 16S rebuild found MDv's
# own dexp(2) analogue miscalibrated (mode-at-zero shrinkage) and only
# fixed it as dhalfnorm(0,1) after SBC actually showed the problem. Worth
# checking directly for ITS rather than assuming the same fix transfers
# unmodified; if SBC below shows the same pathology, switch pr_sigma to
# dhalfnorm(0,1) here (matching 16S's own resolution) before moving on to
# Model 2.
#
# True values: ITS's own real, season-pooled Hill_1 scale by Management
# (Conventional mean ~11.8, CV ~0.99; Organic mean ~14.9, CV ~0.60) --
# picked to reflect the real asymmetry (Conventional's much higher CV,
# driven by its own huge May->July swing), not a copy of 16S's own numbers.

hiermod_marker <- "ITS"
source('src/hiermod/0_SETUP.R')
source('src/hiermod/Models/MD_model.R') # model_MDv_ITS, sim_div_MDv(), means_MDv(), simulate_from_priors_MDv(), dq_MDv
model_MDv <- model_MDv_ITS

hiermod_out_dir <- "out/hiermod/ITS_1_lognormal_MD/Calibration_MDv"

### Parameter + contrast recovery -------------------------------------------

true_conv     <- 12
true_org      <- 15
true_sigma_MDv <- cv_to_sigma(c(0.9, 0.6)) # Conv, Org -- matches the real CV asymmetry

set.seed(20260916)

dat_sim_MDv <- sim_div_MDv(Mg = rbern(250) + 1, loga = log(c(true_conv, true_org)), sigma = true_sigma_MDv)

fit_MDv_sim <- ulam(
  model_MDv,
  data = as.list(dat_sim_MDv),
  chains = 6, cores = 6, iter = 10000)
precis(fit_MDv_sim, depth = 2)

post_MDv_sim <- extract.samples(fit_MDv_sim)

(param_recovery_MDv <- check_recovery(
  true = list(loga1 = log(true_conv), loga2 = log(true_org),
              sigma1 = true_sigma_MDv[1], sigma2 = true_sigma_MDv[2]),
  post_draws = list(loga1 = post_MDv_sim$loga[,1], loga2 = post_MDv_sim$loga[,2],
                    sigma1 = post_MDv_sim$sigma[,1], sigma2 = post_MDv_sim$sigma[,2])
))

### Contrast recovery -------------

pf_MDv_sim <- post_full(fit_MDv_sim, means_MDv)
pc_MDv_sim <- compute_contrasts(pf_MDv_sim, keep = c("mean", "median"), group_levels = idx$Mg$levels)

p_MDv_sim_contrast <- contrast_plot_panels(
  pc_MDv_sim, quant = c(0, 0.999), group_pal = Management_palette); p_MDv_sim_contrast

save_report("sim_summary", model_id_MDv, recovery = param_recovery_MDv, fit_MDv_sim, pc_MDv_sim, model_MDv)
save_gg("sim_contrast_density", model_id_MDv, p_MDv_sim_contrast, width = 8, height = 4)

### loga x sigma funnel check -------------------------------------------------

p_funnel_MDv <- function(){
  par(mfrow = c(1,2))
  plot(post_MDv_sim$sigma[,1], post_MDv_sim$loga[,1],
       xlab = "sigma[1] (Conventional)", ylab = "loga[1] (Conventional)", pch = 16, col = scales::alpha("black", 0.15))
  plot(post_MDv_sim$sigma[,2], post_MDv_sim$loga[,2],
       xlab = "sigma[2] (Organic)", ylab = "loga[2] (Organic)", pch = 16, col = scales::alpha("black", 0.15))
  par(mfrow = c(1,1))
}
save_pdf("loga_sigma_funnel", model_id_MDv, p_funnel_MDv)

### Prior predictive check -------------------------------------------------------

n_prior <- 1000
extracted_prior_MDv <- extract.prior(fit_MDv_sim, n = n_prior)

prior_pred_MDv <- map_dfr(seq_len(n_prior), function(i){
  simulate_from_priors_MDv(draw_true(extracted_prior_MDv, i))
}, .id = "draw")

summary(prior_pred_MDv$Dv) # judge on median/IQR, not mean/SD

p_prior_pc_MDv <- prior_predictive_spaghetti(
  prior_pred_MDv, value_col = "Dv", upper_q = 0.99, model = model_MDv,
  title = "Prior predictive check: sigma[Mg] ~ dexp(1)", observed = dat_sim_MDv$Dv); p_prior_pc_MDv

save_gg("sim_prior_PC", model_id_MDv, p_prior_pc_MDv)

### Simulation-based calibration (SBC), via the SBC package --------------------

sbc_gen_MDv <- make_sbc_generator(
  fit = fit_MDv_sim, simulate_fn = simulate_from_priors_MDv,
  keep = c("loga", "sigma"), gen_cols = c("Dv", "Mg"),
  extra_globals = "sim_div_MDv")

n_sbc  <- 100
n_iter <- 10000

sbc_MDv <- run_sbc_pipeline(
  generator = sbc_gen_MDv$generator, globals = sbc_gen_MDv$globals,
  n_sbc = n_sbc, model = model_MDv, model_id = model_id_MDv, n_iter = n_iter,
  hiermod_out_dir = hiermod_out_dir, dquants = dq_MDv,
  control = list(adapt_delta = 0.99))

plot_sbc_diagnostics(sbc_MDv, model_id_MDv, n_sbc)

sbc_MDv$stats |>
  dplyr::filter(variable %in% c("loga[1]", "loga[2]", "sigma[1]", "sigma[2]")) |>
  dplyr::group_by(variable) |>
  dplyr::summarise(mean_rank_frac = mean(rank / max_rank), median_rank_frac = median(rank / max_rank))

save_sbc_health_report(
  model_id_MDv, sbc_MDv, n_sbc, n_iter,
  variables = c("loga[1]", "loga[2]", "sigma[1]", "sigma[2]",
                "median_contrast", "mean_contrast"),
  hiermod_out_dir = hiermod_out_dir)

# n_sbc=400 stress test, same discipline as every model in this rebuild
# before trusting an n=100 "ok". If sigma[Mg] shows the same mode-at-zero
# shrinkage 16S's own MDv did, switch model_MDv_ITS$prior_sigma to
# dhalfnorm(0,1) (in MD_model.R) and rerun before moving to Model 2.
n_sbc <- 400
sbc_MDv_2 <- run_sbc_pipeline(
  generator = sbc_gen_MDv$generator, globals = sbc_gen_MDv$globals,
  n_sbc = n_sbc, model = model_MDv, model_id = model_id_MDv, n_iter = n_iter,
  hiermod_out_dir = hiermod_out_dir, dquants = dq_MDv,
  control = list(adapt_delta = 0.99))

plot_sbc_diagnostics(sbc_MDv_2, model_id_MDv, n_sbc)

save_sbc_health_report(model_id_MDv, sbc_MDv_2, n_sbc, n_iter,
                       variables = c("loga[1]", "loga[2]", "sigma[1]", "sigma[2]",
                                     "median_contrast", "mean_contrast"),
                       hiermod_out_dir = hiermod_out_dir)
