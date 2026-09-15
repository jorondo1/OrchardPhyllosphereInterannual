# MODEL 2 (MDv, Management-specific variance), 16S: allows heteroscedasticity
# across Mg groups (sigma[Mg] instead of MD's shared sigma). Sibling script
# to 1.2_MD_16S_calibration.R (split out of what used to be one combined
# file) -- both share the "16S_1_lognormal_MD" output directory and the
# 1.3 real-data fit script (MDv gets its own real-data fit there; see that
# script's own header for the current state of that split).
#
# Two purposes here: (1) demonstrate why MDv exists at all -- MD's shared
# sigma gives a biased contrast under heteroscedastic truth, MDv's
# per-group sigma fixes it; (2) the next isolating step in the wider
# loga-miscalibration investigation (see 1.2_MD_16S_calibration.R's header):
# MD (no per-group sigma, no random effects) calibrated cleanly under SBC.
# MDL and MDS2 (both per-group sigma[Mg] PLUS random effects) both showed
# severe loga miscalibration. MDv isolates which addition is responsible:
# per-group sigma[Mg], with STILL no random effects at all. If MDv
# calibrates cleanly like MD, the cause is specifically the random effects
# (Location/Tree/Year); if it's already biased here, it's the
# loga[Mg]-sigma[Mg] per-group interaction itself. loga/sigma priors
# patched to this investigation's established 16S values (dnorm(5,2)/
# dexp(2)) rather than the original, never-validated ITS import
# (dnorm(2,2)/dexp(1)) -- scoped to this script's own model_MD/model_MDv
# only, not 1.3's real-data fits.

hiermod_marker <- "16S"
source('src/hiermod/0_SETUP.R')
source('src/hiermod/Models/MD_model.R') # model_MD_16S, model_MDv_16S, sim_div_MDv(), means_MD(), means_MDv(), simulate_from_priors_MDv(), dq_MDv
model_MD <- model_MD_16S
model_MD$prior_loga  <- quote(loga[Mg] ~ dnorm(5,2))
model_MD$prior_sigma <- quote(sigma ~ dexp(2))
model_MDv <- model_MDv_16S
model_MDv$prior_loga  <- quote(loga[Mg] ~ dnorm(5,2))
model_MDv$prior_sigma <- quote(sigma[Mg] ~ dexp(2))

hiermod_out_dir <- "out/hiermod/16S_1_lognormal_MD"

## MDv -- Allow Management-specific variance (heteroscedasticity) ============
# One simulated dataset under heteroscedastic truth (cv_ 0.5/0.8), fit
# under BOTH models: MD's shared sigma should give a biased contrast, MDv's
# per-group sigma should recover it -- the actual justification for MDv
# existing as its own model.

true_conv <- 180
true_org  <- 120
true_sigma_hetero <- cv_to_sigma(c(0.5, 0.8))

set.seed(20260911)

dat_sim_hetero <- sim_div_MDv(Mg = rbern(250) + 1, loga = log(c(true_conv, true_org)), sigma = true_sigma_hetero)

true_vals_hetero <- tribble(
  ~statistic, ~group,     ~value,
  "median",   "Contrast", true_org - true_conv,
  "mean",     "Contrast", lognormal_mean(log(true_org), true_sigma_hetero[2]^2) -
                           lognormal_mean(log(true_conv), true_sigma_hetero[1]^2)
)

### MD on heteroscedastic truth -- demonstrates the problem ----

fit_MD_hetero <- ulam(model_MD, data = as.list(dat_sim_hetero), chains = 6, cores = 6, iter = 5000)

pf_MD_hetero <- post_full(fit_MD_hetero, means_MD)
pc_MD_hetero <- compute_contrasts(pf_MD_hetero, keep = c("mean", "median"), group_levels = idx$Mg$levels)

p_MD_hetero_contrast <- contrast_plot_panels(
  pc_MD_hetero, quant = c(0, 1), group_pal = Management_palette, true_vals = true_vals_hetero) +
  labs(subtitle = "Shared sigma under heteroscedastic truth: contrast comes out biased"); p_MD_hetero_contrast

save_gg("sim_contrast_density", "MD_hetero", p_MD_hetero_contrast)

### MDv on the same data -- recovers correctly ----

fit_MDv_sim <- ulam(model_MDv, data = as.list(dat_sim_hetero), chains = 6, cores = 6, iter = 5000)
precis(fit_MDv_sim, depth = 2)

post_MDv_sim <- extract.samples(fit_MDv_sim)

(param_recovery_MDv <- check_recovery(
  true = list(loga1 = log(true_conv), loga2 = log(true_org),
              sigma1 = true_sigma_hetero[1], sigma2 = true_sigma_hetero[2]),
  post_draws = list(loga1 = post_MDv_sim$loga[,1], loga2 = post_MDv_sim$loga[,2],
                     sigma1 = post_MDv_sim$sigma[,1], sigma2 = post_MDv_sim$sigma[,2])
))

pf_MDv_sim <- post_full(fit_MDv_sim, means_MDv)
pc_MDv_sim <- compute_contrasts(pf_MDv_sim, keep = c("mean", "median"), group_levels = idx$Mg$levels)

p_MDv_sim_contrast <- contrast_plot_panels(
  pc_MDv_sim, quant = c(0, 1), group_pal = Management_palette, true_vals = true_vals_hetero) +
  labs(subtitle = "Group-specific variance recovers the true contrast"); p_MDv_sim_contrast

save_report("sim_summary", model_id_MDv, recovery = param_recovery_MDv, fit_MDv_sim, pc_MDv_sim, model_MDv, model_name = "The Loose Cannon")
save_gg("sim_contrast_density", model_id_MDv, p_MDv_sim_contrast, width = 8, height = 4)

## MDv -- formal calibration (per-group sigma isolation test) ================
# Reuses fit_MDv_sim above (already fit under the patched (5,2)/(2) priors)
# rather than a fresh calibration fit.

### loga x sigma[Mg] funnel check ----------------------------------------
# Two per-group sigmas to check, unlike MD's single shared one.

p_funnel_MDv <- function(){
  par(mfrow = c(1,2))
  plot(post_MDv_sim$sigma[,1], post_MDv_sim$loga[,1],
       xlab = "sigma[1] (Conventional)", ylab = "loga[1] (Conventional)", pch = 16, col = scales::alpha("black", 0.15))
  plot(post_MDv_sim$sigma[,2], post_MDv_sim$loga[,2],
       xlab = "sigma[2] (Organic)", ylab = "loga[2] (Organic)", pch = 16, col = scales::alpha("black", 0.15))
  par(mfrow = c(1,1))
}
save_pdf("loga_sigma_funnel", model_id_MDv, p_funnel_MDv)

### Simulation-based calibration (SBC), via the SBC package ----------------

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
  dplyr::filter(variable %in% c("loga[1]", "loga[2]")) |>
  dplyr::group_by(variable) |>
  dplyr::summarise(mean_rank_frac = mean(rank / max_rank), median_rank_frac = median(rank / max_rank))

save_sbc_health_report(model_id_MDv, sbc_MDv, n_sbc, n_iter,
                        variables = c("loga[1]", "loga[2]", "sigma[1]", "sigma[2]", "median_contrast", "mean_contrast"),
                        hiermod_out_dir = hiermod_out_dir)
