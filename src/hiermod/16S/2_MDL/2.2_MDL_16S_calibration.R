# MODEL 2 (MDLb), 16S: partial pooling across Location (non-centered
# b[Lo]*sigma_loc). ITS's default loga ~ dnorm(2,2) doesn't make sense here --
# 16S's raw diversity scale is much higher -- so this starts straight from
# loga ~ dnorm(5,2) (same scale already used for MDLS2v/MDLS2vz, so any SBC
# comparison against those is apples-to-apples) rather than re-demonstrating
# the mismatch.

hiermod_marker <- "16S"
source('src/hiermod/0_SETUP.R')
source('src/hiermod/Models/MDL_model.R') # model_MDL_16S, means_MDL(), sim_div_MDL(), simulate_from_priors_MDL()

model <- model_MDL_16S
model$prior_loga     <- quote(loga[Mg] ~ dnorm(5,2))
model$prior_sigma     <- quote(sigma[Mg] ~ dexp(2))
model$prior_sigma_loc <- quote(sigma_loc ~ dexp(2))

hiermod_out_dir <- "out/hiermod/16S_2_lognormal_MDL"

## Parameter recovery -----------------------------------------------------------

conv <- 180
org  <- 120
true_sigma     <- cv_to_sigma(c(0.5, 0.8)) # keep different variances
true_sigma_loc <- 0.5

dat_sim <- sim_div_MDL(
  Mg = rbern(250)+1,
  # 4 locations for simulation; fewer creates a funnel (see below)
  Lo = sample(1:4, 250, replace = TRUE),
  loga = log(c(conv, org)),           # Difference of 60 in median, raw scale
  sigma = true_sigma,
  sigma_loc = true_sigma_loc)

fit_sim <- ulam(
  model,
  data = as.list(dat_sim),
  chains = 6, cores = 6, iter = 5000,
  control = list(adapt_delta = 0.99))
precis(fit_sim, depth = 2)

post_sim <- extract.samples(fit_sim)

### Fixed effect + sigma recovery ---------

check_recovery(
  true = c(loga1 = log(conv), loga2 = log(org)),
  post_draws = list(loga1 = post_sim$loga[,1], loga2 = post_sim$loga[,2])
)

check_recovery(
  true = list(sigma1 = true_sigma[1], sigma2 = true_sigma[2]),
  post_draws = list(sigma1 = post_sim$sigma[,1], sigma2 = post_sim$sigma[,2])
)

### Contrast recovery -------------
# means_MDL() already returns both mean and median. Median needs no
# total_var (median = exp(loga), the raw conv/org values loga was built
# from); mean does (sigma[Mg]^2 + sigma_loc^2, matching means_MDL()'s own
# total_var).

true_median_contrast <- org - conv
true_mean_contrast <- lognormal_mean(log(org), true_sigma[2]^2 + true_sigma_loc^2) -
  lognormal_mean(log(conv), true_sigma[1]^2 + true_sigma_loc^2)

true_vals <- tribble(
  ~statistic, ~group,     ~value,
  "median",   "Contrast", true_median_contrast,
  "mean",     "Contrast", true_mean_contrast
)

pf_sim <- post_full(fit_sim, means_MDL)
pc_sim <- compute_contrasts(pf_sim, keep = c("mean", "median"), group_levels = idx$Mg$levels)

p_sim_contrast <- contrast_plot_panels(
  pc_sim, quant = c(0.005, 0.99), group_pal = Management_palette,
  true_vals = true_vals); p_sim_contrast

save_report("sim_summary", model_id, fit_sim, pc_sim, model, model_name = "The Wildcard")
save_gg("sim_contrast_density", model_id, p_sim_contrast)

save_pdf("sim_traceplot", model_id, function() traceplot(fit_sim))
save_pdf("sim_trankplot", model_id, function() trankplot(fit_sim))

## Prior predictive check ---------------------------------------------------

n_prior <- 10000
prior <- extract.prior(fit_sim, n = n_prior)

# draw_true() is shared (sbc_helpers.R) -- generic slice of any extract.prior()
# output, no per-model rewrite needed.
prior_pred <- map_dfr(seq_len(n_prior), function(i){
  simulate_from_priors_MDL(draw_true(prior, i))
}, .id = "draw")

summary(prior_pred$Dv)  # judge on median/IQR, not mean/SD -- multiplicative
                         # blowups make the raw mean/sd meaningless here

p_sim_spaghetti <- prior_predictive_spaghetti(
  prior_pred, upper_q = 0.99, model = model,
  title = "Prior predictive check: loga[Mg] ~ dnorm(5,2)"); p_sim_spaghetti

save_gg("sim_prior_PC", model_id, p_sim_spaghetti)

## loga x sigma_loc funnel check ---------------------------------------------
# With only n_loc=4, loga and sigma_loc can trade off against each other
# (the classic NCP funnel): motivated by MDLS2v/MDLS2vz's SBC finding that
# loga's bias correlates with sigma_loc. Eyeball whether that geometry is
# already visible here, in the simplest model that includes Location.

p_funnel <- function(){
  par(mfrow = c(1,2))
  plot(post_sim$sigma_loc, post_sim$loga[,1],
       xlab = "sigma_loc", ylab = "loga[1] (Conventional)", pch = 16, col = scales::alpha("black", 0.15))
  plot(post_sim$sigma_loc, post_sim$loga[,2],
       xlab = "sigma_loc", ylab = "loga[2] (Organic)", pch = 16, col = scales::alpha("black", 0.15))
  par(mfrow = c(1,1))
}
save_pdf("loga_sigmaloc_funnel", model_id, p_funnel)

## Simulation-based calibration (SBC), via the SBC package -----------------------
# Checks whether loga[]/sigma_loc's SBC miscalibration found in MDLS2v (and
# still present after the sum-to-zero fix in MDLS2vz) already shows up in
# the simplest model that includes Location at all -- MDL has no
# covariates, no Tree, no Year, just Mg + Lo. If it shows up here too, the
# pathology is intrinsic to a 4-level Location random effect, not
# something later structure (covariates splitting loga, Tree/Year layers)
# amplifies.
# b[Lo] itself stays out of `variables` for the same reason as MDLS2v: the
# simulator draws a fresh per-location offset internally, so a prior draw
# of the array isn't what generated the data.
source('src/utils/sbc_backend_ulam.R')
library(SBC)
future::plan(future::multisession)

generate_one_MDL <- function(){
  true_params <- suppressMessages(suppressWarnings(draw_true(extract.prior(fit_sim, n = 1, refresh = 0), 1)))[
    c("loga", "sigma", "sigma_loc")]
  dat <- simulate_from_priors_MDL(true_params)
  list(variables = true_params, generated = as.list(dat[, c("Dv", "Mg", "Lo")]))
}

# median/mean contrast, matching means_MDL()'s own total_var construction.
dq_MDL <- derived_quantities(
  median_contrast = exp(loga[2]) - exp(loga[1]),
  mean_contrast =
    exp(loga[2] + (sigma[2]^2 + sigma_loc^2) / 2) -
    exp(loga[1] + (sigma[1]^2 + sigma_loc^2) / 2)
)

n_sbc  <- 100
n_iter <- 10000

datasets_path_MDL <- file.path(hiermod_out_dir, "sbc_datasets_MDL.rds")
if (file.exists(datasets_path_MDL)) {
  datasets_MDL <- readRDS(datasets_path_MDL)
} else {
  datasets_MDL <- generate_datasets(SBC_generator_function(generate_one_MDL), n_sbc)
  saveRDS(datasets_MDL, datasets_path_MDL, compress = "xz")
}

backend_MDL <- SBC_backend_ulam(model, iter = n_iter, refresh = 0,
                                 control = list(adapt_delta = 0.99))

sbc_MDL <- compute_SBC(
  datasets_MDL, backend_MDL, dquants = dq_MDL,
  cache_mode = "results", cache_location = file.path(hiermod_out_dir, "sbc_cache_MDL"),
  globals = c("SBC_fit.SBC_backend_ulam", "SBC_fit_to_draws_matrix.ulam",
              "SBC_fit_to_diagnostics.ulam"))

p_sbc_rank  <- plot_rank_hist(sbc_MDL)
p_sbc_ecdf  <- plot_ecdf_diff(sbc_MDL)
p_sbc_cover <- plot_coverage(sbc_MDL)

sbc_step <- paste0(model_id, "_", n_sbc, "sbc_iter")
save_gg("SBC_rank_hist", sbc_step, p_sbc_rank, width = 9, height = 7)
save_gg("SBC_ecdf_diff", sbc_step, p_sbc_ecdf, width = 9, height = 7)
save_gg("SBC_coverage",  sbc_step, p_sbc_cover, width = 9, height = 7)

(sbc_diag_summary <- sbc_MDL$backend_diagnostics %>%
   dplyr::summarise(total_divergent = sum(n_divergent), total_max_treedepth = sum(n_max_treedepth),
                     total_low_ebfmi = sum(n_low_ebfmi), n_replicates = dplyr::n()))

# The actual test: does loga show the same rank-fraction skew here as in
# MDLS2v (~0.16, expected ~0.5)?
sbc_MDL$stats |>
  dplyr::filter(variable %in% c("loga[1]", "loga[2]")) |>
  dplyr::group_by(variable) |>
  dplyr::summarise(mean_rank_frac = mean(rank / max_rank), median_rank_frac = median(rank / max_rank))

# Saved as model_ppc1.rds for compatibility with 2.3_MDL_16S_fit.R's existing
# readRDS() call -- rerun 2.3/2.4 against this if the fit itself needs updating.
saveRDS(model, file.path(hiermod_out_dir, "model_ppc1.rds"))
