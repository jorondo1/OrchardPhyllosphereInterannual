# MODEL 5v (MDLS2v), 16S, SHIFTED (Hill_1 - 1): alternative to Model 5 --
# sigma[cell] reparametrized as a log-linear function of Management/Season
# main effects + interaction (ls0/ls_Mg/ls_Mo/ls_MgMo) instead of 4
# independent dexp-distributed values. Parameter recovery, prior-predictive
# check, and SBC.

hiermod_marker <- "16S"
source('src/hiermod/0_SETUP.R')
source('src/hiermod/Models/MDLS2v_model.R') # sim_div_MDLS2(), means_MDLS2v(), simulate_from_priors_MDLS2v(), contrast_may_gap_MDLS2v()

model <- model_MDLS2v_16S

hiermod_out_dir <- "out/hiermod/16S_5b_lognormal_MDLS2v_shifted"

## Parameter recovery -----------------------------------------------------------

# Same baseline diversity/gap as Model 5's own calibration. True ls* values
# chosen to reproduce roughly the same cell-sigma pattern Model 5's real fit
# found (conv_May/org_May much more dispersed than conv_July/org_July),
# so this tests recovery under the same kind of heterogeneity that caused
# Model 5's divergences, not an easy/artificially tame case.
may_conv <- 180
may_org  <- 120
july_conv_shift <- 0.35
july_org_shift  <- 0.2

# choose Sigma parameter values:
true_ls0     <- 0.34   # conv_May sigma = exp(ls0) ~= 1.4
true_ls_Mg   <- 0.17   # org_May sigma = exp(ls0+ls_Mg) ~= 1.66
true_ls_Mo   <- -1.16  # conv_July sigma = exp(ls0+ls_Mo) ~= 0.44
true_ls_MgMo <- 0.14   # org_July sigma = exp(ls0+ls_Mg+ls_Mo+ls_MgMo) ~= 0.6

# Convert Mg-Mo combination into its sigma (each coef is a log addition)
true_sigma <- true_sigma_from_ls(list(ls0 = true_ls0, ls_Mg = true_ls_Mg,
                                      ls_Mo = true_ls_Mo, ls_MgMo = true_ls_MgMo))
true_sigma # conv_May, conv_July, org_May, org_July

set.seed(20260911)

# Simulator:
dat_sim <- sim_div_MDLS2(
  N_samples = 240,
  n_loc = 4,
  loga = log(c(may_conv, may_org)),
  s_conv = july_conv_shift,
  gap_shift = july_org_shift,
  sigma = true_sigma, # <<<- here
  year_offset = c(0, 0.3, -0.2),
  p_dropout = 0.1,
  shift = 1
); head(dat_sim)

dat_sim %>% count(Lo, Yr, Mo, Mg) %>% print(n=100)
hist(dat_sim$Dv_shifted, breaks = 100)

fit_sim <- ulam(
  model,
  data = as.list(dat_sim),
  chains = 6, cores = 6, iter = 10000,
  control = list(adapt_delta = 0.99)
)

precis(fit_sim, depth = 2)

post_sim <- extract.samples(fit_sim)

### Fixed effects recovery -------------

coef_recovery <- check_recovery(
  true = c(
    # True fixed effects
    loga1 = log(may_conv),
    loga2 = log(may_org),
    s_conv = july_conv_shift,
    gap_shift = july_org_shift,
    
    # True Variance-structure 
    ls0 = true_ls0, ls_Mg = true_ls_Mg, 
    ls_Mo = true_ls_Mo, ls_MgMo = true_ls_MgMo),
  
  # Parameter recovery:
  post_draws = list(
    loga1 = post_sim$loga[,1],
    loga2 = post_sim$loga[,2],
    s_conv = post_sim$s_conv,
    gap_shift = post_sim$gap_shift,
    ls0 = post_sim$ls0, ls_Mg = post_sim$ls_Mg,
    ls_Mo = post_sim$ls_Mo, ls_MgMo = post_sim$ls_MgMo)
); coef_recovery

### Contrast recovery ------------------------------

cr <- contrast_recovery(
  fit_sim, means_MDLS2v, may_conv, may_org, july_conv_shift, july_org_shift, shift = 1)

p_sim_contrast <- contrast_plot_panels(
  cr$estimands, quant = c(0.01, 0.99), scales = 'free_y',
  group_pal = Management_palette,
  true_vals = cr$true_estimands); p_sim_contrast

save_report("sim_summary", model_id, fit_sim, cr$estimands, model,
            recovery = coef_recovery, 
            model_name = "The Structured Splitter")
save_gg("sim_contrast_density", model_id, p_sim_contrast)

## Prior predictive check -------------------------------------------------------

n_prior <- 1000
extracted_prior <- extract.prior(fit_sim, n = n_prior)

prior_pred <- map_dfr(seq_len(n_prior), function(i){
  simulate_from_priors_MDLS2v(draw_true(extracted_prior, i), shift = 1)
}, .id = "draw")

summary(prior_pred$Dv_shifted)

p_prior_pc <- prior_predictive_spaghetti(
  prior_pred, value_col = "Dv_shifted", upper_q = 0.99, model = model,
  title = "Prior predictive check", observed = dat_sim$Dv_shifted) ; p_prior_pc

save_gg("sim_prior_PC", model_id, p_prior_pc)

# Hill a very implausible, mean in the 1E+60+++ range is insane

# Which prior is driving the extreme tail? No stored sigma[cell] parameter
# here (unlike MDLS2) -- derive the 4 cell sigmas from ls0/ls_Mg/ls_Mo/ls_MgMo
# via the shared sigma_cell() helper (MDLS2v_model.R) instead.
sigma_max_prior <- pmax(
  sigma_cell(1,1, extracted_prior$ls0, extracted_prior$ls_Mg, extracted_prior$ls_Mo, extracted_prior$ls_MgMo),
  sigma_cell(1,2, extracted_prior$ls0, extracted_prior$ls_Mg, extracted_prior$ls_Mo, extracted_prior$ls_MgMo),
  sigma_cell(2,1, extracted_prior$ls0, extracted_prior$ls_Mg, extracted_prior$ls_Mo, extracted_prior$ls_MgMo),
  sigma_cell(2,2, extracted_prior$ls0, extracted_prior$ls_Mg, extracted_prior$ls_Mo, extracted_prior$ls_MgMo))

diagnose_extreme_tail(
  prior_pred, value_col = "Dv_shifted",
  candidates = list(
    loga_max  = pmax(extracted_prior$loga[,1], extracted_prior$loga[,2]),
    sigma_loc = extracted_prior$sigma_loc,
    sigma_tr  = extracted_prior$sigma_tr,
    sigma_max = sigma_max_prior))

# Loga is a big big driver, highly correlated 

## Simulation-based calibration (SBC) --------------------------------------------

sbc_MDLS2v <- run_sbc(
  model_fit   = fit_sim,
  means_fn    = means_MDLS2v,
  contrast_fn = contrast_may_gap_MDLS2v,
  simulate_fn = function(true_params) simulate_from_priors_MDLS2v(true_params, shift = 1),
  n_sbc = 30, iter = 15000, n_parallel = 4, chains = 2, cores = 2,
  control = list(adapt_delta = 0.99))

sbc_out_MDLS2v <- save_sbc_report(sbc_MDLS2v, paste0(model_id, "_30sbc_iter"))
hist(sbc_out_MDLS2v$ranks, breaks = 30)
