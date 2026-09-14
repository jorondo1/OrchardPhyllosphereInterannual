# MODEL 6v (MDLSYv), 16S: pools yr[Yr] into a non-centered random effect
# and adds Cultivar (cv[Cv]), on top of MDLS2v's structured variance.
# Parameter recovery, prior-predictive check, and SBC. No variance-budget
# calibration here -- priors get tuned by hand against the plots below if
# needed, not via automatic dexp-rate rescaling.

hiermod_marker <- "16S"
source('src/hiermod/0_SETUP.R')
source('src/hiermod/Models/MDLSYv_model.R') # model_MDLSYv_16S, means_MDLSYv(), simulate_from_priors_MDLSYv(), contrast_may_gap_MDLSYv()

model <- model_MDLSYv_16S

hiermod_out_dir <- "out/hiermod/16S_6_lognormal_MDLSYv"

## Parameter recovery -----------------------------------------------------------

# Same May/July/gap baseline and true ls* values MDLS2v's own calibration
# validated (reproduces Model 5's real heterogeneity pattern) -- this tests
# recovery of the NEW terms (sigma_yr, cv) under that same, already-known-
# hard heterogeneity, not an easier synthetic case.
may_conv <- 180
may_org  <- 120
july_conv_shift <- 0.35
july_org_shift  <- 0.2

true_ls0     <- 0.34
true_ls_Mg   <- 0.17
true_ls_Mo   <- -1.16
true_ls_MgMo <- 0.14
true_sigma_yr <- 0.3
true_cv <- c(0.05, -0.05, 0.1, -0.1, 0)  # 5 Cultivar levels, same magnitudes ITS's own MDLSY calibration used

true_sigma <- true_sigma_from_ls(list(ls0 = true_ls0, ls_Mg = true_ls_Mg,
                                       ls_Mo = true_ls_Mo, ls_MgMo = true_ls_MgMo))
true_sigma # conv_May, conv_July, org_May, org_July

true_coefs <- c(
  loga1 = log(may_conv), loga2 = log(may_org),
  s_conv = july_conv_shift, gap_shift = july_org_shift,
  ls0 = true_ls0, ls_Mg = true_ls_Mg, ls_Mo = true_ls_Mo, ls_MgMo = true_ls_MgMo,
  sigma_yr = true_sigma_yr)

set.seed(20260911)

dat_sim <- simulate_from_priors_MDLSYv(
  list(loga = log(c(may_conv, may_org)), s_conv = july_conv_shift, gap_shift = july_org_shift,
       ls0 = true_ls0, ls_Mg = true_ls_Mg, ls_Mo = true_ls_Mo, ls_MgMo = true_ls_MgMo,
       sigma_loc = 0.5, sigma_tr = 0.3, sigma_yr = true_sigma_yr, cv = true_cv),
  N_samples = 240, n_loc = 4, p_dropout = 0.1, shift = 1
)
head(dat_sim)
dat_sim %>% count(Lo, Yr, Mo, Mg, Cv) %>% print(n = 100)
hist(dat_sim$Dv_shifted, breaks = 100)

fit_sim <- ulam(
  model,
  data = as.list(dat_sim),
  chains = 6, cores = 6, iter = 5000,
  control = list(adapt_delta = 0.99)
)

precis(fit_sim, depth = 2)

post_sim <- extract.samples(fit_sim)

### Fixed effects + variance-structure recovery -------------

(coef_recovery <- check_recovery(
  true = true_coefs,
  post_draws = list(
    loga1 = post_sim$loga[,1], loga2 = post_sim$loga[,2],
    s_conv = post_sim$s_conv, gap_shift = post_sim$gap_shift,
    ls0 = post_sim$ls0, ls_Mg = post_sim$ls_Mg,
    ls_Mo = post_sim$ls_Mo, ls_MgMo = post_sim$ls_MgMo,
    sigma_yr = post_sim$sigma_yr)
))

### Contrast recovery ------------------------------

cr <- contrast_recovery(
  fit_sim, means_MDLSYv, may_conv, may_org, july_conv_shift, july_org_shift, shift = 1)

p_sim_contrast <- contrast_plot_panels(
  cr$estimands, quant = c(0.01, 0.99), scales = 'free_y',
  group_pal = Management_palette,
  true_vals = cr$true_estimands); p_sim_contrast

save_report("sim_summary", model_id, fit_sim, cr$estimands, model,
            recovery = coef_recovery, model_name = "The Varietal Splitter")
save_gg("sim_contrast_density", model_id, p_sim_contrast)

## Prior predictive check -------------------------------------------------------

n_prior <- 1000
extracted_prior <- extract.prior(fit_sim, n = n_prior)

prior_pred <- map_dfr(seq_len(n_prior), function(i){
  simulate_from_priors_MDLSYv(draw_true(extracted_prior, i), shift = 1)
}, .id = "draw")

summary(prior_pred$Dv_shifted)  # judge on median/IQR, not mean/SD -- see MDLS2v's own calibration notes

p_prior_pc <- prior_predictive_spaghetti(
  prior_pred, value_col = "Dv_shifted", upper_q = 0.99, model = model,
  title = "Prior predictive check", observed = dat_sim$Dv_shifted); p_prior_pc

save_gg("sim_prior_PC", model_id, p_prior_pc)

diagnose_extreme_tail(
  prior_pred, value_col = "Dv_shifted",
  candidates = list(
    loga_max  = pmax(extracted_prior$loga[,1], extracted_prior$loga[,2]),
    sigma_loc = extracted_prior$sigma_loc,
    sigma_tr  = extracted_prior$sigma_tr,
    sigma_yr  = extracted_prior$sigma_yr,
    sigma_max = pmax(
      sigma_cell(1,1, extracted_prior$ls0, extracted_prior$ls_Mg, extracted_prior$ls_Mo, extracted_prior$ls_MgMo),
      sigma_cell(1,2, extracted_prior$ls0, extracted_prior$ls_Mg, extracted_prior$ls_Mo, extracted_prior$ls_MgMo),
      sigma_cell(2,1, extracted_prior$ls0, extracted_prior$ls_Mg, extracted_prior$ls_Mo, extracted_prior$ls_MgMo),
      sigma_cell(2,2, extracted_prior$ls0, extracted_prior$ls_Mg, extracted_prior$ls_Mo, extracted_prior$ls_MgMo))
  ))

## Simulation-based calibration -----------------------
# Same template as 5b.2_MDLS2v_16S_calibration.R

source('src/utils/sbc_backend_ulam.R')
library(SBC)
future::plan(future::multisession)3

generate_one_MDLSYv <- function(){3333
  true_params <- suppressMessages(suppressWarnings(
    draw_true(extract.prior(fit_sim, n = 1, refresh = 0), 1)))[
    c("loga", "s_conv", "gap_shift", "ls0", "ls_Mg", "ls_Mo", "ls_MgMo",
      "sigma_loc", "sigma_tr", "sigma_yr", "cv")]
  dat <- simulate_from_priors_MDLSYv(true_params, shift = 1)
  list(variables = true_params,
       generated = as.list(dat[, c("Dv", "Mg", "Lo", "Mo", "Yr", "Tr", "Cv")]))
}

# Same three contrasts as MDLS2v's own SBC (may_gap/july_gap/seasonal_change),
# with sigma_yr^2 added into every total_var term.
dq_MDLSYv <- derived_quantities(
  may_gap =
    exp(loga[2] + (exp(ls0 + ls_Mg)^2 + sigma_loc^2 + sigma_tr^2 + sigma_yr^2) / 2) -
    exp(loga[1] + (exp(ls0)^2         + sigma_loc^2 + sigma_tr^2 + sigma_yr^2) / 2),
  july_gap =
    exp(loga[2] + s_conv + gap_shift + (exp(ls0 + ls_Mg + ls_Mo + ls_MgMo)^2 + sigma_loc^2 + sigma_tr^2 + sigma_yr^2) / 2) -
    exp(loga[1] + s_conv +             (exp(ls0 + ls_Mo)^2                 + sigma_loc^2 + sigma_tr^2 + sigma_yr^2) / 2),
  seasonal_change =
    (exp(loga[2] + s_conv + gap_shift + (exp(ls0 + ls_Mg + ls_Mo + ls_MgMo)^2 + sigma_loc^2 + sigma_tr^2 + sigma_yr^2) / 2) -
     exp(loga[1] + s_conv +             (exp(ls0 + ls_Mo)^2                 + sigma_loc^2 + sigma_tr^2 + sigma_yr^2) / 2)) -
    (exp(loga[2] + (exp(ls0 + ls_Mg)^2 + sigma_loc^2 + sigma_tr^2 + sigma_yr^2) / 2) -
     exp(loga[1] + (exp(ls0)^2         + sigma_loc^2 + sigma_tr^2 + sigma_yr^2) / 2))
)

n_sbc  <- 100
n_iter <- 10000

# Dataset generation left serial -- see 5b.2's own comment for why
# parallelizing it isn't worth chasing (future workers lack rethinking
# and every custom function the generator's call chain touches).
datasets_MDLSYv <- generate_datasets(SBC_generator_function(generate_one_MDLSYv), n_sbc)

backend_MDLSYv <- SBC_backend_ulam(model, chains = 2, iter = n_iter,
                                    refresh = 0, control = list(adapt_delta = 0.99))

sbc_MDLSYv <- compute_SBC(
  datasets_MDLSYv, backend_MDLSYv, dquants = dq_MDLSYv,
  cache_mode = "results", cache_location = file.path(hiermod_out_dir, "sbc_cache_MDLSYv"),
  globals = c("SBC_fit.SBC_backend_ulam", "SBC_fit_to_draws_matrix.ulam",
              "SBC_fit_to_diagnostics.ulam"))

p_sbc_rank  <- plot_rank_hist(sbc_MDLSYv)
p_sbc_ecdf  <- plot_ecdf_diff(sbc_MDLSYv)
p_sbc_cover <- plot_coverage(sbc_MDLSYv)

sbc_step <- paste0(model_id, "_", n_sbc, "sbc_iter")
save_gg("SBC_rank_hist", sbc_step, p_sbc_rank, width = 9, height = 7)
save_gg("SBC_ecdf_diff", sbc_step, p_sbc_ecdf, width = 9, height = 7)
save_gg("SBC_coverage",  sbc_step, p_sbc_cover, width = 9, height = 7)

(sbc_diag_summary <- sbc_MDLSYv$backend_diagnostics %>%
   dplyr::summarise(total_divergent = sum(n_divergent), total_max_treedepth = sum(n_max_treedepth),
                     total_low_ebfmi = sum(n_low_ebfmi), n_replicates = dplyr::n()))
