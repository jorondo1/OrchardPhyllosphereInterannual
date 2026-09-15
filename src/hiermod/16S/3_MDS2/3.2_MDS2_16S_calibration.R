# MODEL 3 (MDS2), 16S, SHIFTED (Hill_1 - 1): Management x Season, no
# Location (see MDS2_model.R's header for why -- Location A/C being
# single-management makes it unfixably confounded with Management, so 16S
# drops it rather than continuing to patch around low cardinality).
# sigma stays Management-only (no cell-level split yet). Parameter
# recovery, prior-predictive check, and SBC via the SBC package -- same
# template as MDL's own rewritten 2.2 script.

hiermod_marker <- "16S"
source('src/hiermod/0_SETUP.R')
source('src/hiermod/Models/MDS2_model.R') # model_MDS2_16S, means_MDS2(), sim_div_MDS2(), simulate_from_priors_MDS2()

model <- model_MDS2_16S

hiermod_out_dir <- "out/hiermod/16S_3_lognormal_MDS2_shifted"

## Parameter recovery -----------------------------------------------------------
# Same baseline/gap values as MDL's and the old MDLS2 family's own
# calibration, so results stay comparable.

may_conv <- 180
may_org  <- 120
july_conv_shift <- 0.35
july_org_shift  <- 0.2

true_sigma <- cv_to_sigma(c(0.5, 0.8)) # conv, org

set.seed(20260911)

dat_sim <- sim_div_MDS2(
  N_samples = 240,
  loga = log(c(may_conv, may_org)),
  s_conv = july_conv_shift,
  gap_shift = july_org_shift,
  sigma = true_sigma,
  year_offset = c(0, 0.3, -0.2),
  sigma_tr = 0.3,
  shift = 1
); head(dat_sim)

hist(dat_sim$Dv_shifted, breaks = 100)

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
    sigma1 = true_sigma[1], sigma2 = true_sigma[2]),
  post_draws = list(
    loga1 = post_sim$loga[,1], loga2 = post_sim$loga[,2],
    sigma1 = post_sim$sigma[,1], sigma2 = post_sim$sigma[,2])
)
)
### Contrast recovery -------------

cr <- contrast_recovery(
  fit_sim, means_MDS2, may_conv, may_org, july_conv_shift, july_org_shift, shift = 1)

p_sim_contrast <- contrast_plot_panels(
  cr$estimands, quant = c(0.01, 0.99), scales = 'free_y',
  group_pal = Management_palette,
  true_vals = cr$true_estimands); p_sim_contrast

save_report("sim_summary", model_id, recovery = param_recovery, fit_sim, cr$estimands, model, model_name = "The Locationless")
save_gg("sim_contrast_density", model_id, p_sim_contrast)

save_pdf("sim_trankplot", model_id, 
         function() trankplot(fit_sim, max_rows = 30, n_cols = 5),
         width = 30, height = 50)

## Prior predictive check -------------------------------------------------------

n_prior <- 1000
extracted_prior <- extract.prior(fit_sim, n = n_prior)

prior_pred <- map_dfr(seq_len(n_prior), function(i){
  simulate_from_priors_MDS2(draw_true(extracted_prior, i), shift = 1)
}, .id = "draw")

summary(prior_pred$Dv_shifted) # judge on median/IQR, not mean/SD

p_prior_pc <- prior_predictive_spaghetti(
  prior_pred, value_col = "Dv_shifted", upper_q = 0.99, model = model,
  title = "Prior predictive check", observed = dat_sim$Dv_shifted); p_prior_pc

save_gg("sim_prior_PC", model_id, p_prior_pc)

## Simulation-based calibration (SBC), via the SBC package -----------------------
# Same template as MDL's/MDLS2v's own calibration. yr[Yr]/tr[Tr] stay out
# of `variables` for the same reason as those models: the simulator draws
# fresh per-level offsets internally, so a prior draw of the array isn't
# what generated the data.
source('src/utils/sbc_backend_ulam.R')
library(SBC)
future::plan(future::multisession)

generate_one_MDS2 <- function(){
  true_params <- suppressMessages(suppressWarnings(draw_true(extract.prior(fit_sim, n = 1, refresh = 0), 1)))[
    c("loga", "s_conv", "gap_shift", "sigma", "sigma_tr")]
  dat <- simulate_from_priors_MDS2(true_params, shift = 1)
  list(variables = true_params,
       generated = as.list(dat[, c("Dv", "Mg", "Mo", "Yr", "Tr")]))
}

dq_MDS2 <- derived_quantities(
  may_gap =
    exp(loga[2] + (sigma[2]^2 + sigma_tr^2) / 2) -
    exp(loga[1] + (sigma[1]^2 + sigma_tr^2) / 2),
  july_gap =
    exp(loga[2] + s_conv + gap_shift + (sigma[2]^2 + sigma_tr^2) / 2) -
    exp(loga[1] + s_conv +             (sigma[1]^2 + sigma_tr^2) / 2),
  seasonal_change =
    (exp(loga[2] + s_conv + gap_shift + (sigma[2]^2 + sigma_tr^2) / 2) -
       exp(loga[1] + s_conv +             (sigma[1]^2 + sigma_tr^2) / 2)) -
    (exp(loga[2] + (sigma[2]^2 + sigma_tr^2) / 2) -
       exp(loga[1] + (sigma[1]^2 + sigma_tr^2) / 2))
)

n_sbc  <- 100
n_iter <- 10000

datasets_path_MDS2 <- file.path(hiermod_out_dir, "sbc_datasets_MDS2.rds")
if (file.exists(datasets_path_MDS2)) {
  datasets_MDS2 <- readRDS(datasets_path_MDS2)
} else {
  datasets_MDS2 <- generate_datasets(SBC_generator_function(generate_one_MDS2), n_sbc)
  saveRDS(datasets_MDS2, datasets_path_MDS2, compress = "xz")
}

backend_MDS2 <- SBC_backend_ulam(model, iter = n_iter, refresh = 0,
                                 control = list(adapt_delta = 0.99))

sbc_MDS2 <- compute_SBC(
  datasets_MDS2, backend_MDS2, dquants = dq_MDS2,
  cache_mode = "results", cache_location = file.path(hiermod_out_dir, "sbc_cache_MDS2"),
  globals = c("SBC_fit.SBC_backend_ulam", "SBC_fit_to_draws_matrix.ulam",
              "SBC_fit_to_diagnostics.ulam"))

p_sbc_rank  <- plot_rank_hist(sbc_MDS2)
p_sbc_ecdf  <- plot_ecdf_diff(sbc_MDS2)
p_sbc_cover <- plot_coverage(sbc_MDS2)

sbc_step <- paste0(model_id, "_", n_sbc, "sbc_iter")
save_gg("SBC_rank_hist", sbc_step, p_sbc_rank, width = 9, height = 7)
save_gg("SBC_ecdf_diff", sbc_step, p_sbc_ecdf, width = 9, height = 7)
save_gg("SBC_coverage",  sbc_step, p_sbc_cover, width = 9, height = 7)

(sbc_diag_summary <- sbc_MDS2$backend_diagnostics %>%
    dplyr::summarise(total_divergent = sum(n_divergent), total_max_treedepth = sum(n_max_treedepth),
                     total_low_ebfmi = sum(n_low_ebfmi), n_replicates = dplyr::n()))

sbc_MDS2$stats |>
  dplyr::filter(variable %in% c("loga[1]", "loga[2]")) |>
  dplyr::group_by(variable) |>
  dplyr::summarise(mean_rank_frac = mean(rank / max_rank), median_rank_frac = median(rank / max_rank))
