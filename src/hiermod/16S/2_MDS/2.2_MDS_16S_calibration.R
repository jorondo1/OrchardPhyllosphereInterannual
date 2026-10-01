# MODEL 2 (MDS), 16S: calibration
# Aim: Management x Season interaction, before any random effect

hiermod_marker <- "16S"
source('src/hiermod/0_SETUP.R')
source('src/hiermod/Models/MDS_model.R')
source('src/utils/sbc_workflow.R') 

model <- model_MDS_16S
model_id <- model_id_MDS

hiermod_out_dir <- "out/hiermod/16S_2_interaction_MDS/Calibration"

## Parameter recovery -----------------------------------------------------------
# Same baseline/gap values throughout the model family (comparability)

may_conv <- 180
may_org  <- 120
july_conv_shift <- 0.35
july_org_shift  <- 0.2

true_sigma <- cv_to_sigma(c(0.5, 0.8)) # conv, org

set.seed(20260911)

dat_sim <- sim_div_MDS(
  N_samples = 240,
  loga = log(c(may_conv, may_org)),
  s_conv = july_conv_shift,
  gap_shift = july_org_shift,
  sigma = true_sigma,
  shift = 1
); head(dat_sim)

fit_sim <- ulam(
  model,
  data = as.list(dat_sim),
  chains = 6, cores = 6, iter = 5000,
  control = list(adapt_delta = 0.99))
precis(fit_sim, depth = 2)

post_sim <- extract.samples(fit_sim)

### Fixed effect + sigma recovery ---------

(param_recovery <- check_recovery(
  true = list(
    loga1 = log(may_conv), loga2 = log(may_org),
    s_conv = july_conv_shift, gap_shift = july_org_shift,
    sigma1 = true_sigma[1], sigma2 = true_sigma[2]),
  post_draws = list(
    loga1 = post_sim$loga[,1], loga2 = post_sim$loga[,2],
    s_conv = post_sim$s_conv, gap_shift = post_sim$gap_shift,
    sigma1 = post_sim$sigma[,1], sigma2 = post_sim$sigma[,2])))

### Contrast recovery -------------

cr <- contrast_recovery(
  fit_sim, means_MDS, may_conv, may_org, july_conv_shift, july_org_shift, shift = 1)

p_sim_contrast <- contrast_plot_panels(
  cr$estimands, quant = c(0.01, 0.99), scales = 'free_y',
  group_pal = Management_palette,
  true_vals = cr$true_estimands); p_sim_contrast

save_report("sim_summary", model_id, recovery = param_recovery, fit_sim, cr$estimands, model)
save_gg("sim_contrast_density", model_id, p_sim_contrast)

save_pdf("sim_trankplot", model_id,
         function() trankplot(fit_sim, max_rows = 30, n_cols = 3),
         width = 10, height = 15)

## Prior predictive check -------------------------------------------------------

n_prior <- 1000
extracted_prior <- extract.prior(fit_sim, n = n_prior)

prior_pred <- map_dfr(seq_len(n_prior), function(i){
  simulate_from_priors_MDS(draw_true(extracted_prior, i), shift = 1)
}, .id = "draw")

summary(prior_pred$Dv_shifted) # judge on median/IQR, not mean/SD

p_prior_pc <- prior_predictive_spaghetti(
  prior_pred, value_col = "Dv_shifted", upper_q = 0.99, model = model,
  title = "Prior predictive check", observed = dat_sim$Dv_shifted); p_prior_pc

save_gg("sim_prior_PC", model_id, p_prior_pc)

## loga/gamma x sigma[Mg] funnel check ---------------------------------------
# New s_conv/gap_shift share the likelihood term: check vs sigma[Mg] too

p_funnel <- function(){
  par(mfrow = c(2,2))
  plot(post_sim$sigma[,1], post_sim$loga[,1],
       xlab = "sigma[1] (Conventional)", ylab = "loga[1] (Conventional)", pch = 16, col = scales::alpha("black", 0.15))
  plot(post_sim$sigma[,2], post_sim$loga[,2],
       xlab = "sigma[2] (Organic)", ylab = "loga[2] (Organic)", pch = 16, col = scales::alpha("black", 0.15))
  plot(post_sim$sigma[,1], post_sim$s_conv,
       xlab = "sigma[1] (Conventional)", ylab = "s_conv", pch = 16, col = scales::alpha("black", 0.15))
  plot(post_sim$sigma[,2], post_sim$gap_shift,
       xlab = "sigma[2] (Organic)", ylab = "gap_shift", pch = 16, col = scales::alpha("black", 0.15))
  par(mfrow = c(1,1))
}
save_pdf("loga_sigma_funnel", model_id, p_funnel)

## Simulation-based calibration (SBC), via the SBC package -----------------------

sbc_gen_MDS <- make_sbc_generator(
  fit = fit_sim, simulate_fn = simulate_from_priors_MDS,
  keep = c("loga", "s_conv", "gap_shift", "sigma"),
  gen_cols = c("Dv", "Mg", "Mo"),
  extra_globals = "sim_div_MDS", shift = 1)

n_sbc  <- 400 # straight up might as well 
n_iter <- 5000

sbc_MDS <- run_sbc_pipeline(
  generator = sbc_gen_MDS$generator, globals = sbc_gen_MDS$globals,
  n_sbc = n_sbc, model = model, model_id = model_id, n_iter = n_iter,
  hiermod_out_dir = hiermod_out_dir, dquants = dq_MDS,
  control = list(adapt_delta = 0.99))

plot_sbc_diagnostics(sbc_MDS, model_id, n_sbc)

sbc_MDS$stats |>
  dplyr::filter(variable %in% c("loga[1]", "loga[2]", "s_conv", "gap_shift")) |>
  dplyr::group_by(variable) |>
  dplyr::summarise(mean_rank_frac = mean(rank / max_rank), median_rank_frac = median(rank / max_rank))

save_sbc_health_report(model_id, sbc_MDS, n_sbc, n_iter,
                        variables = c("loga[1]", "loga[2]", "s_conv", "gap_shift", "sigma[1]", "sigma[2]",
                                      "may_gap", "july_gap", "seasonal_change"),
                        hiermod_out_dir = hiermod_out_dir)
