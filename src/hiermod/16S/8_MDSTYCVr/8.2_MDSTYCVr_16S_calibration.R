# MODEL 8 (MDSTYCVr, "Gimli the Greedy"), 16S: MDSTYCV with Cultivar
# switched from a FIXED sum-to-zero effect to a partially-pooled RANDOM
# effect (cv[Cv]*sigma_cv, non-centered) -- same recipe as Tree's own
# tr[Tr]*sigma_tr.
#
# The actual question this script exists to answer (from discussion, not
# a settled decision): the 5 cultivars here are a specific, deliberate
# choice -- we picked these varieties, we could have picked others -- so
# treating them as an exchangeable sample from a broader cultivar
# population is at least defensible, unlike Year (see
# MDSTYCV_posterior_guide.html section 5). But is 5 levels even enough to
# identify sigma_cv? And since Tree is DETERMINISTICALLY NESTED in
# Cultivar in the real data (129/129 trees map to exactly one cultivar --
# MDSTYCV_model.R's own header), sigma_cv and sigma_tr are now two
# hyperparameters competing for variance at adjacent levels of the SAME
# nesting -- MDST's own sigma_tr already has documented fragility on its
# own (small-sigma_tr funnel, sparse 1-2 obs/tree). This script's pairs-
# check is built specifically around cor(sigma_tr, sigma_cv): if that's
# large, the two are trading off against each other rather than being
# separately identified, which would be the concrete argument against
# this model regardless of what the "5 levels is marginal" heuristic says
# in the abstract.
#
# loga[Mg]/s_conv/gap_shift/sigma[Mg]/yr1/yr2/sigma_tr/b_deg/b_precip/
# b_seq priors are MDSTYCV's own validated answer, hardcoded here as this
# model's starting point. cv[Cv] ~ dnorm(0,1) / sigma_cv ~ dhalfnorm(0,1)
# are the one new assumption to validate -- same prior family as sigma_tr,
# so any difference in behaviour is about the 5-level cardinality/nesting,
# not a different prior choice.

hiermod_marker <- "16S"
source('src/hiermod/0_SETUP.R')
source('src/hiermod/Models/MDSTYCVr_model.R') # model_MDSTYCVr_16S, means_MDSTYCVr(), sim_div_MDSTYCVr(), simulate_from_priors_MDSTYCVr(), dq_MDSTYCVr

model <- model_MDSTYCVr_16S
model_id <- model_id_MDSTYCVr

hiermod_out_dir <- "out/hiermod/16S_8_cultivar_random_MDSTYCVr/Calibration"

## Parameter recovery -----------------------------------------------------------
# Same baseline/gap/year/covariate values as MDSTYCV's own calibration, so
# results stay comparable across the whole family. true_sigma_cv=0.15
# matches the SD implied by MDSTYCV's own fixed cv_1..cv_5 calibration
# values (sd(c(0.2,-0.05,-0.15,0.1,-0.1)) ~ 0.146) -- same effective
# cultivar-to-cultivar spread, just generated as a random draw instead of
# 4 fixed coefficients, so this is an apples-to-apples comparison of the
# two parameterizations, not a different scenario.

may_conv <- 180
may_org  <- 120
july_conv_shift <- 0.35
july_org_shift  <- 0.2

true_sigma    <- cv_to_sigma(c(0.5, 0.8)) # conv, org
true_sigma_tr <- 0.3
true_sigma_cv <- 0.15
true_yr1   <- 0.1
true_yr2   <- 0.3  # implies yr3 = -0.4

true_b_deg    <- 0.2
true_b_precip <- -0.15
true_b_seq    <- 0.3

set.seed(20260917)

dat_sim <- sim_div_MDSTYCVr(
  N_samples = 240,
  loga = log(c(may_conv, may_org)),
  s_conv = july_conv_shift,
  gap_shift = july_org_shift,
  sigma = true_sigma,
  yr1 = true_yr1,
  yr2 = true_yr2,
  sigma_cv = true_sigma_cv,
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
    sigma1 = true_sigma[1], sigma2 = true_sigma[2],
    sigma_tr = true_sigma_tr, sigma_cv = true_sigma_cv,
    yr1 = true_yr1, yr2 = true_yr2,
    b_deg = true_b_deg, b_precip = true_b_precip, b_seq = true_b_seq),
  post_draws = list(
    loga1 = post_sim$loga[,1], loga2 = post_sim$loga[,2],
    s_conv = post_sim$s_conv, gap_shift = post_sim$gap_shift,
    sigma1 = post_sim$sigma[,1], sigma2 = post_sim$sigma[,2],
    sigma_tr = post_sim$sigma_tr, sigma_cv = post_sim$sigma_cv,
    yr1 = post_sim$yr1, yr2 = post_sim$yr2,
    b_deg = post_sim$b_deg, b_precip = post_sim$b_precip, b_seq = post_sim$b_seq)))

### Contrast recovery -------------

cr <- contrast_recovery(
  fit_sim, means_MDSTYCVr, may_conv, may_org, july_conv_shift, july_org_shift, shift = 1)

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
  simulate_from_priors_MDSTYCVr(draw_true(extracted_prior, i), shift = 1)
}, .id = "draw")

summary(prior_pred$Dv_shifted) # judge on median/IQR, not mean/SD

p_prior_pc <- prior_predictive_spaghetti(
  prior_pred, value_col = "Dv_shifted", upper_q = 0.99, model = model,
  title = "Prior predictive check", observed = dat_sim$Dv_shifted); p_prior_pc

save_gg("sim_prior_PC", model_id, p_prior_pc)

## Full pairwise parameter check ------------------------------------------------
# The specific new combination this whole model exists to test: sigma_tr
# against sigma_cv (Tree nested in Cultivar, both now random, for the
# first time). Also loga against everything, given this family's track
# record of small loga leaks showing up in unexpected places.

pairs_vars <- c("loga[1]", "loga[2]", "sigma[1]", "sigma[2]", "sigma_tr", "sigma_cv",
                "yr1", "yr2", "b_deg", "b_precip", "b_seq")
p_pairs <- plot_mcmc_pairs(fit_sim, variables = pairs_vars, n_keep = 1000)
save_gg("mcmc_pairs", model_id, p_pairs, width = 15, height = 15, type = "png")

cat(sprintf("cor(sigma_tr, sigma_cv) in the posterior: %.3f\n",
            cor(post_sim$sigma_tr, post_sim$sigma_cv)))
# A large magnitude here (say, |r| > 0.5) would mean the two variance
# components are trading off against each other rather than being
# separately identified -- the concrete version of the "does nesting
# Cultivar-as-random with Tree actually work" question, not just the
# "5 levels is marginal" heuristic in the abstract.

## Simulation-based calibration (SBC), via the SBC package -----------------------
# tr[Tr]/cv[Cv] stay out of `keep` (cardinality scales with N_samples,
# simulator draws fresh per-tree/per-cultivar offsets each replicate) --
# every hyperparameter, including sigma_tr/sigma_cv, is tracked directly.

sbc_gen_MDSTYCVr <- make_sbc_generator(
  fit = fit_sim, simulate_fn = simulate_from_priors_MDSTYCVr,
  keep = c("loga", "s_conv", "gap_shift", "sigma", "sigma_tr", "sigma_cv",
           "yr1", "yr2", "b_deg", "b_precip", "b_seq"),
  gen_cols = c("Dv", "Mg", "Mo", "Yr", "Cv", "Tr", "deg_h_z", "precip_72h_z", "seq_depth_z"),
  extra_globals = "sim_div_MDSTYCVr", shift = 1)

n_sbc  <- 50
n_iter <- 5000

sbc_MDSTYCVr <- run_sbc_pipeline(
  generator = sbc_gen_MDSTYCVr$generator, globals = sbc_gen_MDSTYCVr$globals,
  n_sbc = n_sbc, model = model, model_id = model_id, n_iter = n_iter,
  hiermod_out_dir = hiermod_out_dir, dquants = dq_MDSTYCVr,
  control = list(adapt_delta = 0.99))

plot_sbc_diagnostics(sbc_MDSTYCVr, model_id, n_sbc)

sbc_MDSTYCVr$stats %>%
  dplyr::filter(variable %in% c("loga[1]", "loga[2]", "sigma_tr", "sigma_cv")) %>%
  dplyr::group_by(variable) %>%
  dplyr::summarise(mean_rank_frac = mean(rank / max_rank),
                   median_rank_frac = median(rank / max_rank))

save_sbc_health_report(
  model_id, sbc_MDSTYCVr, n_sbc, n_iter,
  variables = c("loga[1]", "loga[2]", "s_conv", "gap_shift", "sigma[1]", "sigma[2]",
                "sigma_tr", "sigma_cv", "yr1", "yr2",
                "b_deg", "b_precip", "b_seq",
                "may_gap", "july_gap", "seasonal_change"),
  hiermod_out_dir = hiermod_out_dir)

# Same stress test as every model in this rebuild before trusting an "ok"
# result at n=50. EXPENSIVE (comparable models take ~30min at n_sbc=400 on
# a 6-core laptop) -- run on a cluster/more cores, not as part of a first
# pass.

n_sbc <- 500
n_iter <- 10000
sbc_MDSTYCVr_2 <- run_sbc_pipeline(
  generator = sbc_gen_MDSTYCVr$generator, globals = sbc_gen_MDSTYCVr$globals,
  n_sbc = n_sbc, model = model, model_id = model_id, n_iter = n_iter,
  hiermod_out_dir = hiermod_out_dir, dquants = dq_MDSTYCVr,
  control = list(adapt_delta = 0.99))

plot_sbc_diagnostics(sbc_MDSTYCVr_2, model_id, n_sbc)

sbc_MDSTYCVr_2$stats %>%
  dplyr::filter(variable %in% c("loga[1]", "loga[2]", "sigma_tr", "sigma_cv")) %>%
  dplyr::group_by(variable) %>%
  dplyr::summarise(mean_rank_frac = mean(rank / max_rank),
                   median_rank_frac = median(rank / max_rank))

save_sbc_health_report(
  model_id, sbc_MDSTYCVr_2, n_sbc, n_iter,
  variables = c("loga[1]", "loga[2]", "s_conv", "gap_shift", "sigma[1]", "sigma[2]",
                "sigma_tr", "sigma_cv", "yr1", "yr2",
                "b_deg", "b_precip", "b_seq",
                "may_gap", "july_gap", "seasonal_change"),
  hiermod_out_dir = hiermod_out_dir)

## Compare with MDSTYCV: does a random Cultivar help or hurt sigma_tr? -----------
# MDSTYCV's own standalone numbers (7.2_MDSTYCV_16S_calibration.R,
# sbc_health_MDSTYCV_400iter.txt), same sigma_tr ~ dhalfnorm(0,1) prior --
# n_sbc differs (400 vs this script's 500), close enough for a qualitative
# before/after read, not a formal test. Fill in MDSTYCV's own numbers
# before comparing.

mdstycv_divergences  <- NA # from sbc_health_MDSTYCV_400iter.txt
mdstycv_pct_bad_rhat <- NA
mdstycv_max_rhat     <- NA

mdstycvr_dd <- sbc_MDSTYCVr_2$default_diagnostics
mdstycvr_bd <- sbc_MDSTYCVr_2$backend_diagnostics

cat(sprintf(
  "sigma_tr backend health, MDSTYCV (Cultivar fixed) vs MDSTYCVr (Cultivar random):\n"))
cat(sprintf("  divergences:        %s  ->  %d\n", mdstycv_divergences, sum(mdstycvr_bd$n_divergent)))
cat(sprintf("  %% fits Rhat>1.01:   %s  ->  %.1f%%\n",
            mdstycv_pct_bad_rhat, 100*mean(mdstycvr_dd$max_rhat > 1.01, na.rm = TRUE)))
cat(sprintf("  max Rhat:           %s  ->  %.3f\n", mdstycv_max_rhat, max(mdstycvr_dd$max_rhat, na.rm = TRUE)))
