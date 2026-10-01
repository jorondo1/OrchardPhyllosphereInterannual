# MODEL 3 (MDST, "Treebeard the Skeptic"), ITS: calibration
# - watch for the same sigma_tr fragility as 16S (design property, not bacterial)
# - true values: s_conv/gap_shift as ITS model 2; sigma_tr = 0.3 as 16S

hiermod_marker <- "ITS"
source('src/hiermod/0_SETUP.R')
source('src/hiermod/Models/MDST_model.R') # model_MDST_ITS, means_MDST(), sim_div_MDST(), simulate_from_priors_MDST(), dq_MDST

model <- model_MDST_ITS
model_id <- model_id_MDST

hiermod_out_dir <- "out/hiermod/ITS_3_tree_MDST/Calibration"

## Parameter recovery -----------------------------------------------------------

may_conv <- 20
may_org  <- 16
july_conv_shift <- -1.5
july_org_shift  <- 1.3

true_sigma    <- cv_to_sigma(c(0.8, 0.6)) # conv, org
true_sigma_tr <- 0.3

set.seed(20260916)

dat_sim <- sim_div_MDST(
  N_samples = 240,
  loga = log(c(may_conv, may_org)),
  s_conv = july_conv_shift,
  gap_shift = july_org_shift,
  sigma = true_sigma,
  sigma_tr = true_sigma_tr,
  shift = 1
)

fit_sim <- ulam(
  model,
  data = as.list(dat_sim),
  chains = 6, cores = 6, iter = 5000,
  control = list(adapt_delta = 0.99))
precis(fit_sim, depth = 2)

post_sim <- extract.samples(fit_sim)

### Fixed effect + sigma recovery ---------

param_recovery <- check_recovery(
  true = list(
    loga1 = log(may_conv), loga2 = log(may_org),
    s_conv = july_conv_shift, gap_shift = july_org_shift,
    sigma1 = true_sigma[1], sigma2 = true_sigma[2], sigma_tr = true_sigma_tr),
  post_draws = list(
    loga1 = post_sim$loga[,1], loga2 = post_sim$loga[,2],
    s_conv = post_sim$s_conv, gap_shift = post_sim$gap_shift,
    sigma1 = post_sim$sigma[,1], sigma2 = post_sim$sigma[,2], sigma_tr = post_sim$sigma_tr))

### Contrast recovery -------------

cr <- contrast_recovery(
  fit_sim, means_MDST, may_conv, may_org, july_conv_shift, july_org_shift, shift = 1)

p_sim_contrast <- contrast_plot_panels(
  cr$estimands, quant = c(0.01, 0.99), scales = 'free_y',
  group_pal = Management_palette,
  true_vals = cr$true_estimands)

save_report("sim_summary", model_id, recovery = param_recovery, fit_sim, cr$estimands, model)
save_gg("sim_contrast_density", model_id, p_sim_contrast)

save_pdf("sim_trankplot", model_id,
         function() trankplot(fit_sim, max_rows = 30, n_cols = 5),
         width = 10, height = 15)

## Prior predictive check -------------------------------------------------------

n_prior <- 1000
extracted_prior <- extract.prior(fit_sim, n = n_prior)

prior_pred <- map_dfr(seq_len(n_prior), function(i){
  simulate_from_priors_MDST(draw_true(extracted_prior, i), shift = 1)
}, .id = "draw")

summary(prior_pred$Dv_shifted) # judge on median/IQR, not mean/SD

p_prior_pc <- prior_predictive_spaghetti(
  prior_pred, value_col = "Dv_shifted", upper_q = 0.99, model = model,
  title = "Prior predictive check", observed = dat_sim$Dv_shifted)

save_gg("sim_prior_PC", model_id, p_prior_pc)

## Simulation-based calibration (SBC), via the SBC package -----------------------
# tr[Tr] not tracked; everything else incl. sigma_tr tracked

sbc_gen_MDST <- make_sbc_generator(
  fit = fit_sim, simulate_fn = simulate_from_priors_MDST,
  keep = c("loga", "s_conv", "gap_shift", "sigma", "sigma_tr"),
  gen_cols = c("Dv", "Mg", "Mo", "Tr"),
  extra_globals = "sim_div_MDST", shift = 1)

n_sbc  <- 500
n_iter <- 10000

sbc_MDST <- run_sbc_pipeline(
  generator = sbc_gen_MDST$generator, globals = sbc_gen_MDST$globals,
  n_sbc = n_sbc, model = model, model_id = model_id, n_iter = n_iter,
  hiermod_out_dir = hiermod_out_dir, dquants = dq_MDST,
  control = list(adapt_delta = 0.99))

plot_sbc_diagnostics(sbc_MDST, model_id, n_sbc)

save_sbc_health_report(
  model_id, sbc_MDST, n_sbc, n_iter,
  variables = c("loga[1]", "loga[2]", "s_conv", "gap_shift", "sigma[1]", "sigma[2]",
                "sigma_tr", "may_gap", "july_gap", "seasonal_change"),
  hiermod_out_dir = hiermod_out_dir)
