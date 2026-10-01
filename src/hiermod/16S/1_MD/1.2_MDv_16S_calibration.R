# MODEL 1b (MDv, Management-specific sigma), 16S: calibration
# Aims:
# - show why MDv: MD's shared sigma biases the contrast under heteroscedastic truth
# - isolate the loga miscalibration: per-group sigma alone, still no random effects
# Priors: 16S starting values (dnorm(5,2), dexp(2)); sigma prior revised below

hiermod_marker <- "16S"
source('src/hiermod/0_SETUP.R')
source('src/hiermod/Models/MD_model.R') 
source('src/utils/sbc_workflow.R') 
model_MDv <- model_MDv_16S
model_MDv$prior_loga  <- quote(loga[Mg] ~ dnorm(5,2))
model_MDv$prior_sigma <- quote(sigma[Mg] ~ dexp(2))

hiermod_out_dir <- "out/hiermod/16S_1_lognormal_MD/Calibration"

## MDv -- Allow Management-specific variance (heteroscedasticity) ============
# One heteroscedastic dataset (CV 0.5 / 0.8), fit with MD and MDv

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

fit_MDv_sim <- ulam(
  model_MDv, data = as.list(dat_sim_hetero), 
  chains = 6, cores = 6, iter = 5000)
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

save_report("sim_summary", model_id_MDv, recovery = param_recovery_MDv, fit_MDv_sim, pc_MDv_sim, model_MDv)
save_gg("sim_contrast_density", model_id_MDv, p_MDv_sim_contrast, width = 8, height = 4)

## MDv -- formal calibration (per-group sigma isolation test) ================
# Reuses fit_MDv_sim (no new calibration fit)

### loga x sigma[Mg] funnel check ----------------------------------------
# Two per-group sigmas

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

save_sbc_health_report(
  model_id_MDv, sbc_MDv, n_sbc, n_iter,
  variables = c("loga[1]", "loga[2]", "sigma[1]", "sigma[2]", "median_contrast", "mean_contrast"),
  hiermod_out_dir = hiermod_out_dir)

# One loga slightly off; maybe chance (see report)

n_sbc  <- 400

sbc_MDv_2 <- run_sbc_pipeline(
  generator = sbc_gen_MDv$generator, globals = sbc_gen_MDv$globals,
  n_sbc = n_sbc, model = model_MDv, model_id = model_id_MDv, n_iter = n_iter,
  hiermod_out_dir = hiermod_out_dir, dquants = dq_MDv,
  control = list(adapt_delta = 0.99))
plot_sbc_diagnostics(sbc_MDv_2, model_id_MDv, n_sbc)

sbc_MDv_2$stats |>
  dplyr::filter(variable %in% c("loga[1]", "loga[2]")) |>
  dplyr::group_by(variable) |>
  dplyr::summarise(mean_rank_frac = mean(rank / max_rank), median_rank_frac = median(rank / max_rank))

save_sbc_health_report(
  model_id_MDv, sbc_MDv_2, n_sbc, n_iter,
  variables = c("loga[1]", "loga[2]", "sigma[1]", "sigma[2]", "median_contrast", "mean_contrast"),
  hiermod_out_dir = hiermod_out_dir)
# n = 400: worse -- mean_contrast and loga[1] miscalibrated
# - sigma prior shape distorts the mean-based estimand, leaks into loga[1]

### Calibratation: dexp(1) --------
model_MDv_dexp1 <- model_MDv
model_MDv_dexp1$prior_sigma <- quote(sigma[Mg] ~ dexp(1))
model_id_MDv <- 'MDv_dexp1'
fit_MDv_sim_dexp1 <- ulam(
  model_MDv_dexp1, data = as.list(dat_sim_hetero), 
  chains = 6, cores = 6, iter = 5000)
precis(fit_MDv_sim, depth = 2)

sbc_gen_MDv_dexp1 <- make_sbc_generator(
  fit = fit_MDv_sim_dexp1, simulate_fn = simulate_from_priors_MDv,
  keep = c("loga", "sigma"), gen_cols = c("Dv", "Mg"),
  extra_globals = "sim_div_MDv")

n_sbc  <- 100
n_iter <- 10000

sbc_MDv_dexp1 <- run_sbc_pipeline(
  generator = sbc_gen_MDv_dexp1$generator, globals = sbc_gen_MDv_dexp1$globals,
  n_sbc = n_sbc, model = model_MDv_dexp1, model_id = model_id_MDv, n_iter = n_iter,
  hiermod_out_dir = hiermod_out_dir, dquants = dq_MDv,
  control = list(adapt_delta = 0.99))

plot_sbc_diagnostics(sbc_MDv_dexp1, model_id_MDv, n_sbc)

sbc_MDv_dexp1$stats |>
  dplyr::filter(variable %in% c("loga[1]", "loga[2]")) |>
  dplyr::group_by(variable) |>
  dplyr::summarise(mean_rank_frac = mean(rank / max_rank), median_rank_frac = median(rank / max_rank))

save_sbc_health_report(
  model_id_MDv, sbc_MDv_dexp1, n_sbc, n_iter,
  variables = c("loga[1]", "loga[2]", "sigma[1]", "sigma[2]", "median_contrast", "mean_contrast"),
  hiermod_out_dir = hiermod_out_dir)
# Much better

### Calibration: half-normal(0,1) --------
# half-normal(0,1) vs dexp(1): similar scale (mean 0.8 vs 1)
# - flatter near 0 and lighter tail -> less pull toward tiny sigmas when
#   several variance terms are added later
model_MDv_halfnorm <- model_MDv
model_MDv_halfnorm$prior_sigma <- quote(sigma[Mg] ~ dhalfnorm(0,1))
model_id_MDv <- 'MDv_halfnorm'
fit_MDv_sim_halfnorm <- ulam(
  model_MDv_halfnorm,
  data = as.list(dat_sim_hetero),
  chains = 6, cores = 6, iter = 5000)
precis(fit_MDv_sim_halfnorm, depth = 2)

sbc_gen_MDv_halfnorm <- make_sbc_generator(
  fit = fit_MDv_sim_halfnorm, simulate_fn = simulate_from_priors_MDv,
  keep = c("loga", "sigma"), gen_cols = c("Dv", "Mg"),
  extra_globals = "sim_div_MDv")

n_sbc  <- 100
n_iter <- 10000

sbc_MDv_halfnorm <- run_sbc_pipeline(
  generator = sbc_gen_MDv_halfnorm$generator, 
  globals = sbc_gen_MDv_halfnorm$globals,
  n_sbc = n_sbc, model = model_MDv_halfnorm, 
  model_id = model_id_MDv, 
  n_iter = n_iter,
  hiermod_out_dir = hiermod_out_dir, 
  dquants = dq_MDv,
  control = list(adapt_delta = 0.99))

plot_sbc_diagnostics(sbc_MDv_halfnorm, model_id_MDv, n_sbc)

sbc_MDv_halfnorm$stats |>
  dplyr::filter(variable %in% c("loga[1]", "loga[2]")) |>
  dplyr::group_by(variable) |>
  dplyr::summarise(mean_rank_frac = mean(rank / max_rank), median_rank_frac = median(rank / max_rank))

save_sbc_health_report(
  model_id_MDv, sbc_MDv_halfnorm, n_sbc, n_iter,
  variables = c("loga[1]", "loga[2]", "sigma[1]", "sigma[2]", "median_contrast", "mean_contrast"),
  hiermod_out_dir = hiermod_out_dir)

### Posterior predictive check -- 3rd quartile ---------------------------
# Q3 depends on sigma (unlike the median): does the upper tail drift?
p_ppc_q3_contrast <- plot_ppc_contrast_stat(
  fit_MDv_sim_halfnorm, dat_sim_hetero, dat_sim_hetero$Mg,
  function(y) quantile(y, 0.75), "3rd quartile", idx$Mg$levels)
p_ppc_q3_contrast

save_gg("ppc_q3_contrast", model_id_MDv, p_ppc_q3_contrast)

# Clean at n = 100, but so was dexp(2); one bad fit (Rhat 1.72) -> n = 400 check

n_sbc <- 400

sbc_MDv_halfnorm_2 <- run_sbc_pipeline(
  generator = sbc_gen_MDv_halfnorm$generator, 
  globals = sbc_gen_MDv_halfnorm$globals,
  n_sbc = n_sbc, n_iter = n_iter,
  model = model_MDv_halfnorm, model_id = model_id_MDv, 
  hiermod_out_dir = hiermod_out_dir, 
  dquants = dq_MDv,
  control = list(adapt_delta = 0.99))

plot_sbc_diagnostics(sbc_MDv_halfnorm_2, model_id_MDv, n_sbc)

sbc_MDv_halfnorm_2$stats |>
  dplyr::filter(variable %in% c("loga[1]", "loga[2]")) |>
  dplyr::group_by(variable) |>
  dplyr::summarise(mean_rank_frac = mean(rank / max_rank), median_rank_frac = median(rank / max_rank))

save_sbc_health_report(
  model_id_MDv, sbc_MDv_halfnorm_2, n_sbc, n_iter,
  variables = c("loga[1]", "loga[2]", "sigma[1]", "sigma[2]", "median_contrast", "mean_contrast"),
  hiermod_out_dir = hiermod_out_dir)
