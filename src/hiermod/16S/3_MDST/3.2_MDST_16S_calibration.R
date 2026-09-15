# MODEL 3 (MDST, "Treebeard the Skeptic"), 16S: Management x Season
# interaction plus a Tree random effect (repeated measures: each tree
# contributes one May row and one July row).
#
# MDS treated a tree's two rows as independent draws, overstating effective N 
# and risks sucking real tree-to-tree variation into sigma[Mg]. Here we quantify
# how much of the variation is actually tree-level heterogeneity before trusting
# a (potentially overconfident) season/management contrast.
#
# loga[Mg]/s_conv/gap_shift/sigma[Mg] priors are MDS's own validated answer,
# hardcoded here as this model's starting point (see MDST_model.R header).

# TO VALIDATE:
# sigma_tr ~ dhalfnorm(0,1) is the one new assmuption SBC will check

hiermod_marker <- "16S"
source('src/hiermod/0_SETUP.R')
source('src/hiermod/Models/MDST_model.R') # model_MDST_16S, means_MDST(), sim_div_MDST(), simulate_from_priors_MDST(), dq_MDST

model <- model_MDST_16S
model_id <- model_id_MDST

hiermod_out_dir <- "out/hiermod/16S_3_tree_MDST"

## Parameter recovery -----------------------------------------------------------
# Same baseline/gap values as MDS/MDS2's own calibration, plus sigma_tr=0.3
# (matching MDS2's own convention), so results stay comparable.

may_conv <- 180
may_org  <- 120
july_conv_shift <- 0.35
july_org_shift  <- 0.2

true_sigma    <- cv_to_sigma(c(0.5, 0.8)) # conv, org
true_sigma_tr <- 0.3

set.seed(20260911)

# Simulator
dat_sim <- sim_div_MDST(
  N_samples = 240,
  loga = log(c(may_conv, may_org)),
  s_conv = july_conv_shift,
  gap_shift = july_org_shift,
  sigma = true_sigma,
  sigma_tr = true_sigma_tr,
  shift = 1
)

#Fit simulation
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
    sigma_tr = true_sigma_tr),
  post_draws = list(
    loga1 = post_sim$loga[,1], loga2 = post_sim$loga[,2],
    s_conv = post_sim$s_conv, gap_shift = post_sim$gap_shift,
    sigma1 = post_sim$sigma[,1], sigma2 = post_sim$sigma[,2],
    sigma_tr = post_sim$sigma_tr)))

### Contrast recovery -------------

cr <- contrast_recovery(
  fit_sim, means_MDST, may_conv, may_org, july_conv_shift, july_org_shift, shift = 1)

p_sim_contrast <- contrast_plot_panels(
  cr$estimands, quant = c(0.01, 0.99), scales = 'free_y',
  group_pal = Management_palette,
  true_vals = cr$true_estimands); p_sim_contrast

save_report("sim_summary", model_id, recovery = param_recovery, fit_sim, cr$estimands, model, model_name = "Treebeard the Skeptic")
save_gg("sim_contrast_density", model_id, p_sim_contrast)

save_pdf("sim_trankplot", model_id,
         function() trankplot(fit_sim, max_rows = 30, n_cols = 5),
         width = 30, height = 50)

## Prior predictive check -------------------------------------------------------

n_prior <- 1000
extracted_prior <- extract.prior(fit_sim, n = n_prior)

prior_pred <- map_dfr(seq_len(n_prior), function(i){
  simulate_from_priors_MDST(draw_true(extracted_prior, i), shift = 1)
}, .id = "draw")

summary(prior_pred$Dv_shifted) # judge on median/IQR, not mean/SD

p_prior_pc <- prior_predictive_spaghetti(
  prior_pred, value_col = "Dv_shifted", upper_q = 0.99, model = model,
  title = "Prior predictive check", observed = dat_sim$Dv_shifted); p_prior_pc

save_gg("sim_prior_PC", model_id, p_prior_pc)

## loga/gamma x sigma[Mg]/sigma_tr funnel check ------------------------------
# sigma_tr is the new scale parameter sharing the same likelihood term as
# loga/s_conv/gap_shift: the exact kind of entanglement risk sigma[Mg]
# already turned out to have in Model 1.

p_funnel <- function(){
  par(mfrow = c(2,3))
  plot(post_sim$sigma[,1], post_sim$loga[,1],
       xlab = "sigma[1] (Conventional)", ylab = "loga[1] (Conventional)", pch = 16, col = scales::alpha("black", 0.15))
  plot(post_sim$sigma[,2], post_sim$loga[,2],
       xlab = "sigma[2] (Organic)", ylab = "loga[2] (Organic)", pch = 16, col = scales::alpha("black", 0.15))
  plot(post_sim$sigma_tr, post_sim$loga[,1],
       xlab = "sigma_tr", ylab = "loga[1] (Conventional)", pch = 16, col = scales::alpha("black", 0.15))
  plot(post_sim$sigma_tr, post_sim$loga[,2],
       xlab = "sigma_tr", ylab = "loga[2] (Organic)", pch = 16, col = scales::alpha("black", 0.15))
  plot(post_sim$sigma_tr, post_sim$s_conv,
       xlab = "sigma_tr", ylab = "s_conv", pch = 16, col = scales::alpha("black", 0.15))
  plot(post_sim$sigma_tr, post_sim$gap_shift,
       xlab = "sigma_tr", ylab = "gap_shift", pch = 16, col = scales::alpha("black", 0.15))
  par(mfrow = c(1,1))
}
save_pdf("loga_sigma_funnel", model_id, p_funnel)

## Simulation-based calibration  -----------------------
# tr[Tr] stays out of `variables`/`keep` because the simulator draws a fresh
# per-tree offset internally at each replicate, it's not a predetermined quantity.
# sigma_tr is tracked.

# sigma_tr doesn't sit outside the likelihood like an independent nuisance 
# parameter;it contributes to mu, scaling each tree's own z-score before the sum
# is passed to the likelihood.

sbc_gen_MDST <- make_sbc_generator(
  fit = fit_sim, simulate_fn = simulate_from_priors_MDST,
  keep = c("loga", "s_conv", "gap_shift", "sigma", "sigma_tr"),
  gen_cols = c("Dv", "Mg", "Mo", "Tr"),
  extra_globals = "sim_div_MDST", shift = 1)

n_sbc  <- 400
n_iter <- 10000

sbc_MDST <- run_sbc_pipeline(
  generator = sbc_gen_MDST$generator,
  globals = sbc_gen_MDST$globals,
  n_sbc = n_sbc,  n_iter = n_iter,
  model = model, model_id = model_id,
  hiermod_out_dir = hiermod_out_dir, dquants = dq_MDST,
  control = list(adapt_delta = 0.99))

plot_sbc_diagnostics(sbc_MDST, model_id, n_sbc)

sbc_MDST$stats |>
  dplyr::filter(variable %in% c("loga[1]", "loga[2]", "s_conv", "gap_shift", "sigma_tr")) |>
  dplyr::group_by(variable) |>
  dplyr::summarise(mean_rank_frac = mean(rank / max_rank), median_rank_frac = median(rank / max_rank))

save_sbc_health_report(
  model_id, sbc_MDST, n_sbc, n_iter,
  variables = c("loga[1]", "loga[2]", "s_conv", "gap_shift", "sigma[1]", "sigma[2]",
                "sigma_tr", "may_gap", "july_gap", "seasonal_change"),
  hiermod_out_dir = hiermod_out_dir)

# file.remove('out/hiermod/16S_3_tree_MDST/sbc_cache_MDST.rds')

# Miscalibrations on gap_shift and loga[2]!
# tighten adapt_delta in case this is a funnel problem:
model_id <- "MDST_999"

sbc_MDST_999 <- run_sbc_pipeline(
  generator = sbc_gen_MDST$generator,
  globals = sbc_gen_MDST$globals,
  n_sbc = n_sbc,  n_iter = n_iter,
  model = model, model_id = model_id,
  hiermod_out_dir = hiermod_out_dir, dquants = dq_MDST,
  control = list(adapt_delta = 0.999))

plot_sbc_diagnostics(sbc_MDST_999, model_id, n_sbc)

sbc_MDST_999$stats |>
  dplyr::filter(variable %in% c("loga[1]", "loga[2]", "s_conv", "gap_shift", "sigma_tr")) |>
  dplyr::group_by(variable) |>
  dplyr::summarise(mean_rank_frac = mean(rank / max_rank), median_rank_frac = median(rank / max_rank))

save_sbc_health_report(
  model_id, sbc_MDST_999, n_sbc, n_iter,
  variables = c("loga[1]", "loga[2]", "s_conv", "gap_shift", "sigma[1]", "sigma[2]",
                "sigma_tr", "may_gap", "july_gap", "seasonal_change"),
  hiermod_out_dir = hiermod_out_dir)

# doesn'T help, now both loga are miscalibrated / overconfident.
# Divergences vanished, but treedepth usage nearly doubled and more fits now have bad Rhat
# Not much we can do about this except drop the random effect
# Let's tighten the sigma_tr prior (we can revert back to 0.99)

## Calibration: tighter, more realistic sigma_tr prior -----------------------
# half-normal's density is highest AT zero, so narrowing its scale only packs
# MORE mass into the exact near-zero region that triggers the funnel --
# log-normal has zero density at zero and can center on the real fit's own
# ~0.15 magnitude instead.

model_MDST_tight <- model
model_MDST_tight$pr_sigma_tr <- quote(sigma_tr ~ dlnorm(log(0.15), 0.5))
model_id <- "MDST_tight"

n_sbc  <- 400

fit_sim_tight <- ulam(
  model_MDST_tight, data = as.list(dat_sim),
  chains = 6, cores = 6, iter = 5000,
  control = list(adapt_delta = 0.99))
precis(fit_sim_tight, depth = 2)

sbc_gen_MDST_tight <- make_sbc_generator(
  fit = fit_sim_tight, simulate_fn = simulate_from_priors_MDST,
  keep = c("loga", "s_conv", "gap_shift", "sigma", "sigma_tr"),
  gen_cols = c("Dv", "Mg", "Mo", "Tr"),
  extra_globals = "sim_div_MDST", shift = 1)

sbc_MDST_tight <- run_sbc_pipeline(
  generator = sbc_gen_MDST_tight$generator,
  globals = sbc_gen_MDST_tight$globals,
  n_sbc = n_sbc, n_iter = n_iter,
  model = model_MDST_tight, model_id = model_id,
  hiermod_out_dir = hiermod_out_dir, dquants = dq_MDST,
  control = list(adapt_delta = 0.99))

plot_sbc_diagnostics(sbc_MDST_tight, model_id, n_sbc)

sbc_MDST_tight$stats |>
  dplyr::filter(variable %in% c("loga[1]", "loga[2]", "s_conv", "gap_shift", "sigma_tr")) |>
  dplyr::group_by(variable) |>
  dplyr::summarise(mean_rank_frac = mean(rank / max_rank), median_rank_frac = median(rank / max_rank))

save_sbc_health_report(
  model_id, sbc_MDST_tight, n_sbc, n_iter,
  variables = c("loga[1]", "loga[2]", "s_conv", "gap_shift", "sigma[1]", "sigma[2]",
                "sigma_tr", "may_gap", "july_gap", "seasonal_change"),
  hiermod_out_dir = hiermod_out_dir)

