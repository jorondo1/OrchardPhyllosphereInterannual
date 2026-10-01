# MODEL 7 (MDSTYCV, "Saruman the Fool"), 16S: calibration
# Aim: add Tree (model 3) to the Year/Covariate/Cultivar model (model 6)
# Question: Tree nested in Cultivar -- does sigma_tr's fragility (MDST) improve, worsen or stay?
# - same sigma_tr prior as MDST, so health is directly comparable

hiermod_marker <- "16S"
source('src/hiermod/0_SETUP.R')
source('src/hiermod/Models/MDSTYCV_model.R')
source('src/utils/sbc_workflow.R') 

model <- model_MDSTYCV_16S
model_id <- model_id_MDSTYCV

hiermod_out_dir <- "out/hiermod/16S_7_tree_full_MDSTYCV/Calibration"

## Parameter recovery -----------------------------------------------------------
# Same values as MDSYCV, plus sigma_tr = 0.3

may_conv <- 180
may_org  <- 120
july_conv_shift <- 0.35
july_org_shift  <- 0.2

true_sigma    <- cv_to_sigma(c(0.5, 0.8)) # conv, org
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

set.seed(20260911)

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
  title = "Prior predictive check", observed = dat_sim$Dv_shifted); p_prior_pc

save_gg("sim_prior_PC", model_id, p_prior_pc)

## Full pairwise parameter check ------------------------------------------------
# New: sigma_tr vs cv_* (Tree and Cultivar over the same trees); loga vs everything

pairs_vars <- c("loga[1]", "loga[2]", "sigma[1]", "sigma[2]", "sigma_tr",
                "yr1", "yr2", "cv_1", "cv_3", "cv_4", "cv_5",
                "b_deg", "b_precip", "b_seq")
p_pairs <- plot_mcmc_pairs(fit_sim, variables = pairs_vars, n_keep = 1000)
save_gg("mcmc_pairs", model_id, p_pairs, width = 15, height = 15, type = "png")

## Simulation-based calibration (SBC), via the SBC package -----------------------
# tr[Tr] not tracked (fresh per-tree offsets); everything else incl. sigma_tr tracked

sbc_gen_MDSTYCV <- make_sbc_generator(
  fit = fit_sim, simulate_fn = simulate_from_priors_MDSTYCV,
  keep = c("loga", "s_conv", "gap_shift", "sigma", "sigma_tr", "yr1", "yr2",
           "cv_1", "cv_3", "cv_4", "cv_5", "b_deg", "b_precip", "b_seq"),
  gen_cols = c("Dv", "Mg", "Mo", "Yr", "Cv", "Tr", "deg_h_z", "precip_72h_z", "seq_depth_z"),
  extra_globals = "sim_div_MDSTYCV", shift = 1)

n_sbc  <- 50
n_iter <- 5000

sbc_MDSTYCV <- run_sbc_pipeline(
  generator = sbc_gen_MDSTYCV$generator, globals = sbc_gen_MDSTYCV$globals,
  n_sbc = n_sbc, model = model, model_id = model_id, n_iter = n_iter,
  hiermod_out_dir = hiermod_out_dir, dquants = dq_MDSTYCV,
  control = list(adapt_delta = 0.99))

plot_sbc_diagnostics(sbc_MDSTYCV, model_id, n_sbc)

sbc_MDSTYCV$stats %>% 
  dplyr::filter(variable %in% c("loga[1]", "loga[2]", "sigma_tr", "cv_1", "cv_3", "cv_4", "cv_5")) %>% 
  dplyr::group_by(variable) %>% 
  dplyr::summarise(mean_rank_frac = mean(rank / max_rank), 
                   median_rank_frac = median(rank / max_rank))

save_sbc_health_report(
  model_id, sbc_MDSTYCV, n_sbc, n_iter,
  variables = c("loga[1]", "loga[2]", "s_conv", "gap_shift", "sigma[1]", "sigma[2]",
                "sigma_tr", "yr1", "yr2", "cv_1", "cv_3", "cv_4", "cv_5",
                "b_deg", "b_precip", "b_seq",
                "may_gap", "july_gap", "seasonal_change"),
  hiermod_out_dir = hiermod_out_dir)

# n = 400 stress test (n = 100 alone not trusted)

n_sbc <- 400
n_iter <- 10000
sbc_MDSTYCV_2 <- run_sbc_pipeline(
  generator = sbc_gen_MDSTYCV$generator, globals = sbc_gen_MDSTYCV$globals,
  n_sbc = n_sbc, model = model, model_id = model_id, n_iter = n_iter,
  hiermod_out_dir = hiermod_out_dir, dquants = dq_MDSTYCV,
  control = list(adapt_delta = 0.99))

plot_sbc_diagnostics(sbc_MDSTYCV_2, model_id, n_sbc)

sbc_MDSTYCV_2$stats %>% 
  dplyr::filter(variable %in% c("loga[1]", "loga[2]", "sigma_tr", "cv_1", "cv_3", "cv_4", "cv_5")) %>% 
  dplyr::group_by(variable) %>% 
  dplyr::summarise(mean_rank_frac = mean(rank / max_rank), 
                   median_rank_frac = median(rank / max_rank))

save_sbc_health_report(
  model_id, sbc_MDSTYCV_2, n_sbc, n_iter,
  variables = c("loga[1]", "loga[2]", "s_conv", "gap_shift", "sigma[1]", "sigma[2]",
                "sigma_tr", "yr1", "yr2", "cv_1", "cv_3", "cv_4", "cv_5",
                "b_deg", "b_precip", "b_seq",
                "may_gap", "july_gap", "seasonal_change"),
  hiermod_out_dir = hiermod_out_dir)

## Compare with MDST : does Cultivar help or hurt sigma_tr? -------
# MDST at n = 400, same sigma_tr prior (sbc_health_MDST_400iter.txt)

mdst_divergences   <- 3218
mdst_pct_bad_rhat  <- 74/400
mdst_max_rhat      <- 2.119

mdstycv_dd <- sbc_MDSTYCV_2$default_diagnostics
mdstycv_bd <- sbc_MDSTYCV_2$backend_diagnostics

cat(sprintf(
  "sigma_tr backend health, MDST (Tree alone) vs MDSTYCV (Tree + Year/Cultivar/Covariates):\n"))
cat(sprintf("  divergences:        %d  ->  %d\n", mdst_divergences, sum(mdstycv_bd$n_divergent)))
cat(sprintf("  %% fits Rhat>1.01:   %.1f%%  ->  %.1f%%\n",
            100*mdst_pct_bad_rhat, 100*mean(mdstycv_dd$max_rhat > 1.01, na.rm = TRUE)))
cat(sprintf("  max Rhat:           %.3f  ->  %.3f\n", mdst_max_rhat, max(mdstycv_dd$max_rhat, na.rm = TRUE)))
# Result vs MDST:
# divergences:        3218  ->  0
# % fits Rhat>1.01:   18.5%  ->  17.8%
# max Rhat:           2.119  ->  2.123
# - cultivar doesn't fix everything, but the divergence drop is a good sign
# - ECDF plots not perfect but overall good; estimands clean

