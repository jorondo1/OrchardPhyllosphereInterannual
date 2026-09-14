# MODEL 5c (MDLS2vz), 16S, SHIFTED (Hill_1 - 1): MDLS2v + a sum-to-zero
# constraint on the Location random effect b[Lo] (see MDLS2vz_model.R's
# header for the SBC finding that motivated this). Parameter recovery,
# prior-predictive check, and SBC -- same template as
# 5b.2_MDLS2v_16S_calibration.R.

hiermod_marker <- "16S"
source('src/hiermod/0_SETUP.R')
source('src/hiermod/Models/MDLS2vz_model.R') # model_MDLS2vz_16S, sim_div_MDLS2vz(), simulate_from_priors_MDLS2vz()

model <- model_MDLS2vz_16S

hiermod_out_dir <- "out/hiermod/16S_5c_lognormal_MDLS2vz_shifted"

## Parameter recovery -----------------------------------------------------------

# Same baseline/true values as MDLS2v's own calibration -- this tests
# recovery under the exact same heterogeneity, isolating the effect of the
# sum-to-zero change rather than introducing a new scenario.
may_conv <- 180
may_org  <- 120
july_conv_shift <- 0.35
july_org_shift  <- 0.2

true_ls0     <- 0.34
true_ls_Mg   <- 0.17
true_ls_Mo   <- -1.16
true_ls_MgMo <- 0.14

true_sigma <- true_sigma_from_ls(list(ls0 = true_ls0, ls_Mg = true_ls_Mg,
                                       ls_Mo = true_ls_Mo, ls_MgMo = true_ls_MgMo))
true_sigma # conv_May, conv_July, org_May, org_July

true_coefs <- c(
  loga1 = log(may_conv), loga2 = log(may_org),
  s_conv = july_conv_shift, gap_shift = july_org_shift,
  ls0 = true_ls0, ls_Mg = true_ls_Mg, ls_Mo = true_ls_Mo, ls_MgMo = true_ls_MgMo)

set.seed(20260911)

dat_sim <- sim_div_MDLS2vz(
  N_samples = 240,
  n_loc = 4,
  loga = log(c(may_conv, may_org)),
  s_conv = july_conv_shift,
  gap_shift = july_org_shift,
  sigma = true_sigma,
  year_offset = c(0, 0.3, -0.2),
  p_dropout = 0.1,
  shift = 1
); head(dat_sim)

dat_sim %>% count(Lo, Yr, Mo, Mg) %>% print(n = 100)
hist(dat_sim$Dv_shifted, breaks = 100)

fit_sim <- ulam(
  model,
  data = as.list(dat_sim),
  chains = 6, cores = 6, iter = 5000,
  control = list(adapt_delta = 0.99)
)

precis(fit_sim, depth = 2)

post_sim <- extract.samples(fit_sim)

### Sum-to-zero sanity check -------------
# b1/b2/b3 are the raw scalar parameters (i.i.d. dnorm(0,1) priors);
# b_loc4 = -(b1+b2+b3) is reconstructed here in R post-hoc (ulam()'s own
# gq> mechanism fails even for a trivial scalar-from-scalar derivation;
# see MDLS2vz_model.R's own comment). Sum should be ~0 for every draw;
# confirms the constraint is actually enforced, not just declared. b_loc4
# has 3x the others' prior variance by construction; worth eyeballing
# here, not assuming. ---> is it a problem?

b_loc4 <- -(post_sim$b1 + post_sim$b2 + post_sim$b3)
b_mat <- cbind(post_sim$b1, post_sim$b2, post_sim$b3, b_loc4)
b_sums <- rowSums(b_mat)
cat("sum(b1,b2,b3,b_loc4) per draw : mean:", mean(b_sums), " max abs:", max(abs(b_sums)), "\n")
cat("marginal SD per location [1..4] (expect the 4th ~sqrt(3)x the others):\n")
print(apply(b_mat, 2, sd))

### Fixed effects + variance-structure recovery -------------

(coef_recovery <- check_recovery(
  true = true_coefs,
  post_draws = list(
    loga1 = post_sim$loga[,1], loga2 = post_sim$loga[,2],
    s_conv = post_sim$s_conv, gap_shift = post_sim$gap_shift,
    ls0 = post_sim$ls0, ls_Mg = post_sim$ls_Mg,
    ls_Mo = post_sim$ls_Mo, ls_MgMo = post_sim$ls_MgMo)
))

### Contrast recovery ------------------------------

cr <- contrast_recovery(
  fit_sim, means_MDLS2v, may_conv, may_org, july_conv_shift, july_org_shift, shift = 1)

p_sim_contrast <- contrast_plot_panels(
  cr$estimands, quant = c(0.01, 0.99), scales = 'free_y',
  group_pal = Management_palette,
  true_vals = cr$true_estimands); p_sim_contrast

save_report("sim_summary", model_id, fit_sim, cr$estimands, model,
            recovery = coef_recovery, model_name = "The Zero-Summed Splitter")
save_gg("sim_contrast_density", model_id, p_sim_contrast)

## Prior predictive check -------------------------------------------------------

n_prior <- 1000
extracted_prior <- extract.prior(fit_sim, n = n_prior)

prior_pred <- map_dfr(seq_len(n_prior), function(i){
  simulate_from_priors_MDLS2vz(draw_true(extracted_prior, i), shift = 1)
}, .id = "draw")

summary(prior_pred$Dv_shifted)  # judge on median/IQR, not mean/SD

p_prior_pc <- prior_predictive_spaghetti(
  prior_pred, value_col = "Dv_shifted", upper_q = 0.99, model = model,
  title = "Prior predictive check", observed = dat_sim$Dv_shifted); p_prior_pc

save_gg("sim_prior_PC", model_id, p_prior_pc)

## Simulation-based calibration (SBC), via the SBC package -----------------------
# Same template as 5b.2_MDLS2v_16S_calibration.R. b[Lo]/tr[Tr]/yr[Yr] stay
# out of `variables` for the same reason as MDLS2v: the simulator draws
# fresh per-level offsets internally (now sum-to-zero for Location), so a
# prior draw of the array isn't what generated the data.
source('src/utils/sbc_backend_ulam.R')
library(SBC)
future::plan(future::multisession)

generate_one_MDLS2vz <- function(){
  true_params <- suppressMessages(suppressWarnings(draw_true(extract.prior(fit_sim, n = 1, refresh = 0), 1)))[
    c("loga", "s_conv", "gap_shift", "ls0", "ls_Mg", "ls_Mo", "ls_MgMo", "sigma_loc", "sigma_tr")]
  dat <- simulate_from_priors_MDLS2vz(true_params, shift = 1)
  list(variables = true_params,
       generated = as.list(dat[, c("Dv", "Mg", "Lo", "Mo", "Yr", "Tr")]))
}

# Same three contrasts as MDLS2v's own SBC.
dq_MDLS2vz <- derived_quantities(
  may_gap =
    exp(loga[2] + (exp(ls0 + ls_Mg)^2 + sigma_loc^2 + sigma_tr^2) / 2) -
    exp(loga[1] + (exp(ls0)^2         + sigma_loc^2 + sigma_tr^2) / 2),
  july_gap =
    exp(loga[2] + s_conv + gap_shift + (exp(ls0 + ls_Mg + ls_Mo + ls_MgMo)^2 + sigma_loc^2 + sigma_tr^2) / 2) -
    exp(loga[1] + s_conv +             (exp(ls0 + ls_Mo)^2                 + sigma_loc^2 + sigma_tr^2) / 2),
  seasonal_change =
    (exp(loga[2] + s_conv + gap_shift + (exp(ls0 + ls_Mg + ls_Mo + ls_MgMo)^2 + sigma_loc^2 + sigma_tr^2) / 2) -
     exp(loga[1] + s_conv +             (exp(ls0 + ls_Mo)^2                 + sigma_loc^2 + sigma_tr^2) / 2)) -
    (exp(loga[2] + (exp(ls0 + ls_Mg)^2 + sigma_loc^2 + sigma_tr^2) / 2) -
     exp(loga[1] + (exp(ls0)^2         + sigma_loc^2 + sigma_tr^2) / 2))
)

n_sbc  <- 100
n_iter <- 10000

datasets_path_MDLS2vz <- file.path(hiermod_out_dir, "sbc_datasets_MDLS2vz.rds")
if (file.exists(datasets_path_MDLS2vz)) {
  datasets_MDLS2vz <- readRDS(datasets_path_MDLS2vz)
} else {
  datasets_MDLS2vz <- generate_datasets(SBC_generator_function(generate_one_MDLS2vz), n_sbc)
  saveRDS(datasets_MDLS2vz, datasets_path_MDLS2vz, compress = "xz")
}

backend_MDLS2vz <- SBC_backend_ulam(model, iter = n_iter,
                                     refresh = 0, control = list(adapt_delta = 0.99))

sbc_MDLS2vz <- compute_SBC(
  datasets_MDLS2vz, backend_MDLS2vz, dquants = dq_MDLS2vz,
  cache_mode = "results", cache_location = file.path(hiermod_out_dir, "sbc_cache_MDLS2vz"),
  globals = c("SBC_fit.SBC_backend_ulam", "SBC_fit_to_draws_matrix.ulam",
              "SBC_fit_to_diagnostics.ulam"))

p_sbc_rank  <- plot_rank_hist(sbc_MDLS2vz)
p_sbc_ecdf  <- plot_ecdf_diff(sbc_MDLS2vz)
p_sbc_cover <- plot_coverage(sbc_MDLS2vz)

sbc_step <- paste0(model_id, "_", n_sbc, "sbc_iter")
save_gg("SBC_rank_hist", sbc_step, p_sbc_rank, width = 9, height = 7)
save_gg("SBC_ecdf_diff", sbc_step, p_sbc_ecdf, width = 9, height = 7)
save_gg("SBC_coverage",  sbc_step, p_sbc_cover, width = 9, height = 7)

(sbc_diag_summary <- sbc_MDLS2vz$backend_diagnostics %>%
   dplyr::summarise(total_divergent = sum(n_divergent), total_max_treedepth = sum(n_max_treedepth),
                     total_low_ebfmi = sum(n_low_ebfmi), n_replicates = dplyr::n()))

# Compare loga's calibration specifically against MDLS2v's own result --
# this is the actual test of whether the fix worked.
sbc_MDLS2vz$stats |>
  dplyr::filter(variable %in% c("loga[1]", "loga[2]")) |>
  dplyr::group_by(variable) |>
  dplyr::summarise(mean_rank_frac = mean(rank / max_rank), median_rank_frac = median(rank / max_rank))
