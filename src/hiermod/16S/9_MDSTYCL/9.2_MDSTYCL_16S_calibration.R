# MODEL 9 (MDSTYCL, "Faramir the Judicious"), 16S: MDSTYCV with Cultivar's
# fixed effect replaced by Location's fixed effect (lo1, 2 levels).
#
# Why Location instead of Cultivar-as-random (Model 8, MDSTYCVr): that
# attempt failed badly on real SBC -- 297 divergences, 154/500 (30.8%)
# fits Rhat>1.01, loga[1]/loga[2] severely MISCALIBRATED (z=-11.96/-14.22),
# not just sigma_cv (z=-5.76). Nesting sigma_cv with sigma_tr over the
# same 129 trees leaked into the headline estimand, not just the new
# parameter. A fixed effect doesn't have that funnel/variance-competition
# failure mode, so this returns to the same recipe Cultivar (Model 6) and
# Year already used, just for a different grouping factor.
#
# Why Location is worth testing at all: never in this from-scratch
# family's lineage before, and Tree is deterministically nested in
# Location the same way it's nested in Cultivar (129/129) -- structurally
# the same "does Tree's own random effect coexist with a coarser fixed
# grouping over the same trees" question MDSTYCV already answered cleanly
# for Cultivar, so there's real reason to expect this one to behave.
#
# Why Cultivar is dropped, not kept alongside Location: confirmed by
# direct query of the real data (not assumed) -- Location x Management is
# A=Conventional-only (50 rows), C=Organic-only (51), B and D have both.
# Model 9's real fit will restrict to the B/D subset (141/242 rows) so
# Location isn't confounded with Management the way it is for A/C. But
# within that B/D subset, Location x Cultivar is a HARD alias for 3 of 5
# cultivars: Location D has ONLY Honeycrisp/Spartan (0 Cortland/Liberty/
# Paulared) -- Cortland/Liberty/Paulared are perfectly collinear with
# "Location B" there, not just correlated. Cultivar has to go.
#
# This calibration script itself simulates Lo independently of Mg/Cv --
# it tests whether the MODEL STRUCTURE (2-level Location fixed effect +
# Tree nested in it, alongside Year/Covariates) is recoverable in
# principle, same scope as every calibration script in this family. The
# real B/D subset's own confound is a real-fit-stage (9.3, not yet built)
# concern, not a calibration one.
#
# loga[Mg]/s_conv/gap_shift/sigma[Mg]/yr1/yr2/sigma_tr/b_deg/b_precip/
# b_seq priors are MDSTYCV's own validated answer, hardcoded here as this
# model's starting point. lo1 ~ dnorm(0,1) is the one new assumption to
# validate.

hiermod_marker <- "16S"
source('src/hiermod/0_SETUP.R')
source('src/hiermod/Models/MDSTYCL_model.R') # model_MDSTYCL_16S, means_MDSTYCL(), sim_div_MDSTYCL(), simulate_from_priors_MDSTYCL(), dq_MDSTYCL

model <- model_MDSTYCL_16S
model_id <- model_id_MDSTYCL

hiermod_out_dir <- "out/hiermod/16S_9_location_MDSTYCL/Calibration"

## Parameter recovery -----------------------------------------------------------
# Same baseline/gap/year/covariate/sigma_tr values as MDSTYCV's own
# calibration, so results stay comparable across the whole family.
# true_lo1=0.15 is an arbitrary but modest offset, comparable in magnitude
# to a single cv_* level from Model 8's own calibration.

may_conv <- 180
may_org  <- 120
july_conv_shift <- 0.35
july_org_shift  <- 0.2

true_sigma    <- cv_to_sigma(c(0.5, 0.8)) # conv, org
true_sigma_tr <- 0.3
true_lo1      <- 0.15 # implies lo2 = -0.15
true_yr1   <- 0.1
true_yr2   <- 0.3  # implies yr3 = -0.4

true_b_deg    <- 0.2
true_b_precip <- -0.15
true_b_seq    <- 0.3

set.seed(20260917)

dat_sim <- sim_div_MDSTYCL(
  N_samples = 240,
  loga = log(c(may_conv, may_org)),
  s_conv = july_conv_shift,
  gap_shift = july_org_shift,
  sigma = true_sigma,
  yr1 = true_yr1,
  yr2 = true_yr2,
  lo1 = true_lo1,
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
    yr1 = true_yr1, yr2 = true_yr2, lo1 = true_lo1,
    b_deg = true_b_deg, b_precip = true_b_precip, b_seq = true_b_seq),
  post_draws = list(
    loga1 = post_sim$loga[,1], loga2 = post_sim$loga[,2],
    s_conv = post_sim$s_conv, gap_shift = post_sim$gap_shift,
    sigma1 = post_sim$sigma[,1], sigma2 = post_sim$sigma[,2], sigma_tr = post_sim$sigma_tr,
    yr1 = post_sim$yr1, yr2 = post_sim$yr2, lo1 = post_sim$lo1,
    b_deg = post_sim$b_deg, b_precip = post_sim$b_precip, b_seq = post_sim$b_seq)))

### Contrast recovery -------------

cr <- contrast_recovery(
  fit_sim, means_MDSTYCL, may_conv, may_org, july_conv_shift, july_org_shift, shift = 1)

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
  simulate_from_priors_MDSTYCL(draw_true(extracted_prior, i), shift = 1)
}, .id = "draw")

summary(prior_pred$Dv_shifted) # judge on median/IQR, not mean/SD

p_prior_pc <- prior_predictive_spaghetti(
  prior_pred, value_col = "Dv_shifted", upper_q = 0.99, model = model,
  title = "Prior predictive check", observed = dat_sim$Dv_shifted); p_prior_pc

save_gg("sim_prior_PC", model_id, p_prior_pc)

## Full pairwise parameter check ------------------------------------------------
# sigma_tr against lo1 (Tree nested in Location, the new combination this
# model exists to test) -- Model 8's own failure mode was a large
# sigma_tr/sigma_cv trade-off, so check the analogous pair here too, even
# though lo1 is a mean-level fixed effect, not a variance component (a
# smaller risk in principle, worth confirming rather than assuming).

pairs_vars <- c("loga[1]", "loga[2]", "sigma[1]", "sigma[2]", "sigma_tr", "lo1",
                "yr1", "yr2", "b_deg", "b_precip", "b_seq")
p_pairs <- plot_mcmc_pairs(fit_sim, variables = pairs_vars, n_keep = 1000)
save_gg("mcmc_pairs", model_id, p_pairs, width = 15, height = 15, type = "png")

cat(sprintf("cor(sigma_tr, lo1) in the posterior: %.3f\n",
            cor(post_sim$sigma_tr, post_sim$lo1)))

## Simulation-based calibration (SBC), via the SBC package -----------------------
# tr[Tr] stays out of `keep` (cardinality scales with N_samples, simulator
# draws fresh per-tree offsets each replicate) -- every other parameter,
# including sigma_tr/lo1, is tracked directly.

sbc_gen_MDSTYCL <- make_sbc_generator(
  fit = fit_sim, simulate_fn = simulate_from_priors_MDSTYCL,
  keep = c("loga", "s_conv", "gap_shift", "sigma", "sigma_tr", "yr1", "yr2", "lo1",
           "b_deg", "b_precip", "b_seq"),
  gen_cols = c("Dv", "Mg", "Mo", "Yr", "Lo", "Tr", "deg_h_z", "precip_72h_z", "seq_depth_z"),
  extra_globals = "sim_div_MDSTYCL", shift = 1)

n_sbc  <- 100
n_iter <- 5000

sbc_MDSTYCL <- run_sbc_pipeline(
  generator = sbc_gen_MDSTYCL$generator, globals = sbc_gen_MDSTYCL$globals,
  n_sbc = n_sbc, model = model, model_id = model_id, n_iter = n_iter,
  hiermod_out_dir = hiermod_out_dir, dquants = dq_MDSTYCL,
  control = list(adapt_delta = 0.99))

plot_sbc_diagnostics(sbc_MDSTYCL, model_id, n_sbc)

save_sbc_health_report(
  model_id, sbc_MDSTYCL, n_sbc, n_iter,
  variables = c("loga[1]", "loga[2]", "s_conv", "gap_shift", "sigma[1]", "sigma[2]",
                "sigma_tr", "yr1", "yr2", "lo1",
                "b_deg", "b_precip", "b_seq",
                "may_gap", "july_gap", "seasonal_change"),
  hiermod_out_dir = hiermod_out_dir)

# Same stress test as every model in this rebuild before trusting an "ok"
# result at n=100. EXPENSIVE (comparable models take ~30min at n_sbc=400-
# 500 on a 6-core laptop) -- run on a cluster/more cores, not as part of a
# first pass.

n_sbc <- 500
n_iter <- 10000
sbc_MDSTYCL_2 <- run_sbc_pipeline(
  generator = sbc_gen_MDSTYCL$generator, globals = sbc_gen_MDSTYCL$globals,
  n_sbc = n_sbc, model = model, model_id = model_id, n_iter = n_iter,
  hiermod_out_dir = hiermod_out_dir, dquants = dq_MDSTYCL,
  control = list(adapt_delta = 0.99))

plot_sbc_diagnostics(sbc_MDSTYCL_2, model_id, n_sbc)

save_sbc_health_report(
  model_id, sbc_MDSTYCL_2, n_sbc, n_iter,
  variables = c("loga[1]", "loga[2]", "s_conv", "gap_shift", "sigma[1]", "sigma[2]",
                "sigma_tr", "yr1", "yr2", "lo1",
                "b_deg", "b_precip", "b_seq",
                "may_gap", "july_gap", "seasonal_change"),
  hiermod_out_dir = hiermod_out_dir)
