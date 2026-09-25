# MODEL 7 (MDSTYCV, "Saruman the Fool"), ITS: merges Tree (MDST, Model 3)
# into the Year+Covariates+Cultivar branch (MDSYCV, Model 6) -- the last
# model in this family, mirroring 7.2_MDSTYCV_16S_calibration.R. Tree is
# the same deterministic-nesting-in-Cultivar/Location design as 16S (same
# physical trees), so the same "does Cultivar's fixed effect and Tree's
# random effect coexist cleanly" question applies here, including
# sigma_tr's own known fragility (see MDST_model.R's header).
#
# model_MDSTYCV_ITS: loga[Mg] ~ dnorm(2,2) is ITS's own scale; everything
# else carries over unchanged from model_MDSTYCV_16S. True values match
# every earlier ITS calibration script in this family (see
# 2.2_MDS_ITS_calibration.R's header for the real-data May/July rationale).

hiermod_marker <- "ITS"
source('src/hiermod/0_SETUP.R')
source('src/hiermod/Models/MDSTYCV_model.R') # model_MDSTYCV_ITS, means_MDSTYCV(), sim_div_MDSTYCV(), simulate_from_priors_MDSTYCV(), dq_MDSTYCV

model <- model_MDSTYCV_ITS
model_id <- model_id_MDSTYCV

hiermod_out_dir <- "out/hiermod/ITS_7_tree_full_MDSTYCV/Calibration"

## Parameter recovery -----------------------------------------------------------

may_conv <- 20
may_org  <- 16
july_conv_shift <- -1.5
july_org_shift  <- 1.3

true_sigma    <- cv_to_sigma(c(0.8, 0.6)) # conv, org
true_sigma_tr <- 0.3
true_yr1   <- 0.1
true_yr2   <- 0.3  # implies yr3 = -0.4

true_cv_1 <- 0.2   # Cortland
true_cv_3 <- -0.15 # Paulared
true_cv_4 <- 0.1   # Honeycrisp
true_cv_5 <- -0.1  # Spartan
# implied cv_2 (Liberty) = -0.05

true_b_deg    <- 0.2
true_b_precip <- -0.15
true_b_seq    <- 0.3

set.seed(20260916)

dat_sim <- sim_div_MDSTYCV(
  N_samples = 240,
  loga = log(c(may_conv, may_org)),
  s_conv = july_conv_shift,
  gap_shift = july_org_shift,
  sigma = true_sigma,
  yr1 = true_yr1,
  yr2 = true_yr2,
  cv_1 = true_cv_1,
  cv_3 = true_cv_3,
  cv_4 = true_cv_4,
  cv_5 = true_cv_5,
  sigma_tr = true_sigma_tr,
  b_deg = true_b_deg,
  b_precip = true_b_precip,
  b_seq = true_b_seq,
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

(param_recovery <- check_recovery(
  true = list(
    loga1 = log(may_conv), loga2 = log(may_org),
    s_conv = july_conv_shift, gap_shift = july_org_shift,
    sigma1 = true_sigma[1], sigma2 = true_sigma[2], sigma_tr = true_sigma_tr,
    yr1 = true_yr1, yr2 = true_yr2,
    cv_1 = true_cv_1, cv_3 = true_cv_3, cv_4 = true_cv_4, cv_5 = true_cv_5,
    b_deg = true_b_deg, b_precip = true_b_precip, b_seq = true_b_seq),
  post_draws = list(
    loga1 = post_sim$loga[,1], loga2 = post_sim$loga[,2],
    s_conv = post_sim$s_conv, gap_shift = post_sim$gap_shift,
    sigma1 = post_sim$sigma[,1], sigma2 = post_sim$sigma[,2], sigma_tr = post_sim$sigma_tr,
    yr1 = post_sim$yr1, yr2 = post_sim$yr2,
    cv_1 = post_sim$cv_1, cv_3 = post_sim$cv_3, cv_4 = post_sim$cv_4, cv_5 = post_sim$cv_5,
    b_deg = post_sim$b_deg, b_precip = post_sim$b_precip, b_seq = post_sim$b_seq)))

### Contrast recovery -------------

cr <- contrast_recovery(
  fit_sim, means_MDSTYCV, may_conv, may_org, july_conv_shift, july_org_shift, shift = 1)

p_sim_contrast <- contrast_plot_panels(
  cr$estimands, quant = c(0.01, 0.99), scales = 'free_y',
  group_pal = Management_palette,
  true_vals = cr$true_estimands)

save_report("sim_summary", model_id, recovery = param_recovery,
            fit_sim, cr$estimands, model)
save_gg("sim_contrast_density", model_id, p_sim_contrast)

save_pdf("sim_trankplot", model_id,
         function() trankplot(fit_sim, max_rows = 50, n_cols = 10),
         width = 20, height = 30)

## Prior predictive check -------------------------------------------------------

n_prior <- 1000
extracted_prior <- extract.prior(fit_sim, n = n_prior)

prior_pred <- map_dfr(seq_len(n_prior), function(i){
  simulate_from_priors_MDSTYCV(draw_true(extracted_prior, i), shift = 1)
}, .id = "draw")

summary(prior_pred$Dv_shifted) # judge on median/IQR, not mean/SD

p_prior_pc <- prior_predictive_spaghetti(
  prior_pred, value_col = "Dv_shifted", upper_q = 0.99, model = model,
  title = "Prior predictive check", observed = dat_sim$Dv_shifted)

save_gg("sim_prior_PC", model_id, p_prior_pc)

## Full pairwise parameter check ------------------------------------------------
# sigma_tr against cv_1..cv_4 (Tree and Cultivar partition the same trees)
# is the specific new combination to watch, same as 16S's own check.

pairs_vars <- c("loga[1]", "loga[2]", "sigma[1]", "sigma[2]", "sigma_tr",
                "yr1", "yr2", "cv_1", "cv_3", "cv_4", "cv_5",
                "b_deg", "b_precip", "b_seq")
p_pairs <- plot_mcmc_pairs(fit_sim, variables = pairs_vars, n_keep = 1000)
save_gg("mcmc_pairs", model_id, p_pairs, width = 15, height = 15, type = "png")

## Simulation-based calibration (SBC), via the SBC package -----------------------
# tr[Tr] stays out of `keep` (cardinality scales with N_samples, simulator
# draws fresh per-tree offsets each replicate) -- every other parameter,
# including sigma_tr, is tracked directly.

sbc_gen_MDSTYCV <- make_sbc_generator(
  fit = fit_sim, simulate_fn = simulate_from_priors_MDSTYCV,
  keep = c("loga", "s_conv", "gap_shift", "sigma", "sigma_tr", "yr1", "yr2",
           "cv_1", "cv_3", "cv_4", "cv_5", "b_deg", "b_precip", "b_seq"),
  gen_cols = c("Dv", "Mg", "Mo", "Yr", "Cv", "Tr", "deg_h_z", "precip_72h_z", "seq_depth_z"),
  extra_globals = "sim_div_MDSTYCV", shift = 1)

n_sbc  <- 500
n_iter <- 10000

sbc_MDSTYCV <- run_sbc_pipeline(
  generator = sbc_gen_MDSTYCV$generator, globals = sbc_gen_MDSTYCV$globals,
  n_sbc = n_sbc, model = model, model_id = model_id, n_iter = n_iter,
  hiermod_out_dir = hiermod_out_dir, dquants = dq_MDSTYCV,
  control = list(adapt_delta = 0.99))

plot_sbc_diagnostics(sbc_MDSTYCV, model_id, n_sbc)

save_sbc_health_report(
  model_id, sbc_MDSTYCV, n_sbc, n_iter,
  variables = c("loga[1]", "loga[2]", "s_conv", "gap_shift", "sigma[1]", "sigma[2]",
                "sigma_tr", "yr1", "yr2", "cv_1", "cv_3", "cv_4", "cv_5",
                "b_deg", "b_precip", "b_seq",
                "may_gap", "july_gap", "seasonal_change"),
  hiermod_out_dir = hiermod_out_dir)
