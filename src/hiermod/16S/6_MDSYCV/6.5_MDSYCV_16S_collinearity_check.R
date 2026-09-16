# MODEL 6 (MDSYCV): collinearity-aware SBC stress test. 

#  real weather/sequencing-depth collinearity: cor(deg_h_z, Mo)=0.73, 
# cor(precip_72h_z, Mo)=0.55, cor(seq_depth_z, Mg)=0.31 (src/check_seq_depth_confound.R,
# src/hiermod/16S/TODO.md). Every calibration script up to now (including
# MDSYCV's own 6.2) drew covariates independently of Mg/Mo -- this checks
# whether b_deg/b_precip/b_seq/s_conv/gap_shift stay identifiable once the
# simulator matches the real data's own collinearity structure, using
# sim_div_MDSYCV()'s new rho_deg_season/rho_precip_season/rho_seq_mg
# arguments (MDSYCV_model.R, template borrowed from the ITS lineage's own
# MDLSYC_model.R, adapted to the specific correlations actually confirmed
# in this data rather than the generic "correlate with mu" version there).

hiermod_marker <- "16S"
source('src/hiermod/0_SETUP.R')
source('src/hiermod/Models/MDSYCV_model.R') # model_MDSYCV_16S, means_MDSYCV(), sim_div_MDSYCV(), simulate_from_priors_MDSYCV(), dq_MDSYCV

model <- model_MDSYCV_16S
model_id <- "MDSYCV_collin"

hiermod_out_dir <- "out/hiermod/16S_6_cultivar_MDSYCV"

## Parameter recovery, under realistic collinearity ----------------------------
# Same baseline/gap/year/cultivar/covariate values as 6.2's own calibration
# (comparable), but deg_h_z/precip_72h_z correlated with Season and
# seq_depth_z correlated with Management at the real data's own measured
# strength -- not independent draws.

may_conv <- 180
may_org  <- 120
july_conv_shift <- 0.35
july_org_shift  <- 0.2

true_sigma <- cv_to_sigma(c(0.5, 0.8)) # conv, org
true_yr1   <- 0.1
true_yr2   <- 0.3  # implies yr3 = -0.4

true_cv_1 <- 0.2   # Cortland
true_cv_3 <- -0.15 # Paulared
true_cv_4 <- 0.1   # Honeycrisp
true_cv_5 <- -0.1  # Spartan

true_b_deg    <- 0.2
true_b_precip <- -0.15
true_b_seq    <- 0.3

# Real, measured correlations (see header) -- not guesses.
rho_deg_season    <- 0.735
rho_precip_season <- 0.547
rho_seq_mg        <- 0.310

set.seed(20260911)

dat_sim <- sim_div_MDSYCV(
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
  b_deg = true_b_deg,
  b_precip = true_b_precip,
  b_seq = true_b_seq,
  shift = 1,
  rho_deg_season = rho_deg_season,
  rho_precip_season = rho_precip_season,
  rho_seq_mg = rho_seq_mg
)

# Sanity check: did the simulator actually achieve the intended collinearity?
cat("Achieved cor(deg_h_z, Mo):", cor(dat_sim$deg_h_z, dat_sim$Mo), "(target", rho_deg_season, ")\n")
cat("Achieved cor(precip_72h_z, Mo):", cor(dat_sim$precip_72h_z, dat_sim$Mo), "(target", rho_precip_season, ")\n")
cat("Achieved cor(seq_depth_z, Mg):", cor(dat_sim$seq_depth_z, dat_sim$Mg), "(target", rho_seq_mg, ")\n")

fit_sim <- ulam(
  model,
  data = as.list(dat_sim),
  chains = 6, cores = 6, iter = 5000,
  control = list(adapt_delta = 0.99))
precis(fit_sim, depth = 2)


### Fixed effect + sigma recovery ---------
post_sim <- extract.samples(fit_sim)

# The parameters most at risk under this specific collinearity: s_conv/
# gap_shift (share Season with deg_h_z/precip_72h_z) and loga/sigma[Mg]
# (share Management with seq_depth_z).

(param_recovery <- check_recovery(
  true = list(
    loga1 = log(may_conv), loga2 = log(may_org),
    s_conv = july_conv_shift, gap_shift = july_org_shift,
    sigma1 = true_sigma[1], sigma2 = true_sigma[2],
    yr1 = true_yr1, yr2 = true_yr2,
    cv_1 = true_cv_1, cv_3 = true_cv_3, cv_4 = true_cv_4, cv_5 = true_cv_5,
    b_deg = true_b_deg, b_precip = true_b_precip, b_seq = true_b_seq),
  post_draws = list(
    loga1 = post_sim$loga[,1], loga2 = post_sim$loga[,2],
    s_conv = post_sim$s_conv, gap_shift = post_sim$gap_shift,
    sigma1 = post_sim$sigma[,1], sigma2 = post_sim$sigma[,2],
    yr1 = post_sim$yr1, yr2 = post_sim$yr2,
    cv_1 = post_sim$cv_1, cv_3 = post_sim$cv_3, cv_4 = post_sim$cv_4, cv_5 = post_sim$cv_5,
    b_deg = post_sim$b_deg, b_precip = post_sim$b_precip, b_seq = post_sim$b_seq)))

## Collinearity check: does the correlation show up in the posterior too? ------
# The direct question -- b_deg vs s_conv/gap_shift, b_seq vs loga/sigma[Mg].

pairs_vars <- c("loga[1]", "loga[2]", "s_conv", "gap_shift", "sigma[1]", "sigma[2]",
                 "b_deg", "b_precip", "b_seq")
p_pairs <- plot_mcmc_pairs(fit_sim, variables = pairs_vars, n_keep = 1000)
save_gg("mcmc_pairs", model_id, p_pairs, width = 15, height = 15, type = "png")

cat("cor(b_deg, s_conv):", cor(post_sim$b_deg, post_sim$s_conv), "\n")
cat("cor(b_deg, gap_shift):", cor(post_sim$b_deg, post_sim$gap_shift), "\n")
cat("cor(b_seq, loga[1]):", cor(post_sim$b_seq, post_sim$loga[,1]), "\n")
cat("cor(b_seq, sigma[1]):", cor(post_sim$b_seq, post_sim$sigma[,1]), "\n")

## Simulation-based calibration (SBC), under realistic collinearity ------------
# extra_globals doesn't need anything new -- rho_* are passed as fixed
# extra args via make_sbc_generator()'s own `...` (identical to how `shift`
# is already forwarded), applied identically on every replicate.

sbc_gen_collin <- make_sbc_generator(
  fit = fit_sim, simulate_fn = simulate_from_priors_MDSYCV,
  keep = c("loga", "s_conv", "gap_shift", "sigma", "yr1", "yr2",
           "cv_1", "cv_3", "cv_4", "cv_5", "b_deg", "b_precip", "b_seq"),
  gen_cols = c("Dv", "Mg", "Mo", "Yr", "Cv", "deg_h_z", "precip_72h_z", "seq_depth_z"),
  extra_globals = "sim_div_MDSYCV", shift = 1,
  rho_deg_season = rho_deg_season, rho_precip_season = rho_precip_season, rho_seq_mg = rho_seq_mg)

n_sbc  <- 100
n_iter <- 10000

sbc_collin <- run_sbc_pipeline(
  generator = sbc_gen_collin$generator, globals = sbc_gen_collin$globals,
  n_sbc = n_sbc, model = model, model_id = model_id, n_iter = n_iter,
  hiermod_out_dir = hiermod_out_dir, dquants = dq_MDSYCV,
  control = list(adapt_delta = 0.99))

plot_sbc_diagnostics(sbc_collin, model_id, n_sbc)

sbc_collin$stats |>
  dplyr::filter(variable %in% c("loga[1]", "loga[2]", "s_conv", "gap_shift", "b_deg", "b_precip", "b_seq")) |>
  dplyr::group_by(variable) |>
  dplyr::summarise(mean_rank_frac = mean(rank / max_rank), median_rank_frac = median(rank / max_rank))

save_sbc_health_report(model_id, sbc_collin, n_sbc, n_iter,
                        variables = c("loga[1]", "loga[2]", "s_conv", "gap_shift", "sigma[1]", "sigma[2]",
                                      "yr1", "yr2", "cv_1", "cv_3", "cv_4", "cv_5",
                                      "b_deg", "b_precip", "b_seq",
                                      "may_gap", "july_gap", "seasonal_change"),
                        hiermod_out_dir = hiermod_out_dir)

# Compare this report against 6.2's own baseline (sbc_health_MDSYCV_400iter.txt,
# independent-covariate draws) -- if b_deg/b_precip/s_conv/gap_shift/b_seq
# newly show up flagged here but not there, that's the realistic-collinearity
# cost this rebuild hasn't had to pay yet. If the health report looks the
# same as the baseline, the real data's own collinearity isn't a practical
# identifiability problem for this model, just a reporting-caveat one
# (already documented in TODO.md/the posterior guides).
