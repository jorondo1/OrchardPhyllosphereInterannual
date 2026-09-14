# MODEL 5v (MDLS2v), 16S, SHIFTED (Hill_1 - 1): alternative to Model 5 --
# sigma[cell] reparametrized as a log-linear function of Management/Season
# main effects + interaction (ls0/ls_Mg/ls_Mo/ls_MgMo) instead of 4
# independent dexp-distributed values. Parameter recovery, prior-predictive
# check, and SBC.

hiermod_marker <- "16S"
source('src/hiermod/0_SETUP.R')
source('src/hiermod/Models/MDLS2v_model.R') # sim_div_MDLS2(), means_MDLS2v(), simulate_from_priors_MDLS2v(), contrast_may_gap_MDLS2v()

model <- model_MDLS2v_16S

hiermod_out_dir <- "out/hiermod/16S_5b_lognormal_MDLS2v_shifted"

## Parameter recovery -----------------------------------------------------------

# Same baseline diversity/gap as Model 5's own calibration. True ls* values
# chosen to reproduce roughly the same cell-sigma pattern Model 5's real fit
# found (conv_May/org_May much more dispersed than conv_July/org_July),
# so this tests recovery under the same kind of heterogeneity that caused
# Model 5's divergences, not an easy/artificially tame case.
may_conv <- 180
may_org  <- 120
july_conv_shift <- 0.35
july_org_shift  <- 0.2

# choose Sigma parameter values:
true_ls0     <- 0.34   # conv_May sigma = exp(ls0) ~= 1.4
true_ls_Mg   <- 0.17   # org_May sigma = exp(ls0+ls_Mg) ~= 1.66
true_ls_Mo   <- -1.16  # conv_July sigma = exp(ls0+ls_Mo) ~= 0.44
true_ls_MgMo <- 0.14   # org_July sigma = exp(ls0+ls_Mg+ls_Mo+ls_MgMo) ~= 0.6

# Convert Mg-Mo combination into its sigma (each coef is a log addition)
true_sigma <- true_sigma_from_ls(list(ls0 = true_ls0, ls_Mg = true_ls_Mg,
                                      ls_Mo = true_ls_Mo, ls_MgMo = true_ls_MgMo))
true_sigma # conv_May, conv_July, org_May, org_July

# TRUE EFFECT LIST
true_coefs <- c(
  # True fixed effects
  loga1 = log(may_conv),
  loga2 = log(may_org),
  s_conv = july_conv_shift,
  gap_shift = july_org_shift,
  
  # True Variance-structure 
  ls0 = true_ls0, ls_Mg = true_ls_Mg, 
  ls_Mo = true_ls_Mo, ls_MgMo = true_ls_MgMo)

set.seed(20260911)

# Simulator:
dat_sim <- sim_div_MDLS2(
  N_samples = 240,
  n_loc = 4,
  loga = log(c(may_conv, may_org)),
  s_conv = july_conv_shift,
  gap_shift = july_org_shift,
  sigma = true_sigma, # <<<- here
  year_offset = c(0, 0.3, -0.2),
  p_dropout = 0.1,
  shift = 1
); head(dat_sim)

dat_sim %>% count(Lo, Yr, Mo, Mg) %>% print(n=100)
hist(dat_sim$Dv_shifted, breaks = 100)

fit_sim <- ulam(
  model,
  data = as.list(dat_sim),
  chains = 6, cores = 6, iter = 5000,
  control = list(adapt_delta = 0.99)
)

precis(fit_sim, depth = 2)

post_sim <- extract.samples(fit_sim)

### Fixed effects recovery -------------

(coef_recovery <- check_recovery(
  true = true_coefs,
  
  # Parameter recovery:
  post_draws = list(
    loga1 = post_sim$loga[,1],
    loga2 = post_sim$loga[,2],
    s_conv = post_sim$s_conv,
    gap_shift = post_sim$gap_shift,
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
            recovery = coef_recovery, 
            model_name = "The Structured Splitter")
save_gg("sim_contrast_density", model_id, p_sim_contrast)

## Prior predictive check -------------------------------------------------------

n_prior <- 1000
extracted_prior <- extract.prior(fit_sim, n = n_prior)

prior_pred <- map_dfr(seq_len(n_prior), function(i){
  simulate_from_priors_MDLS2v(draw_true(extracted_prior, i), shift = 1)
}, .id = "draw")

summary(prior_pred$Dv_shifted)

p_prior_pc <- prior_predictive_spaghetti(
  prior_pred, value_col = "Dv_shifted", upper_q = 0.99, model = model,
  title = "Prior predictive check", observed = dat_sim$Dv_shifted) ; p_prior_pc

save_gg("sim_prior_PC", model_id, p_prior_pc)


### WE SHOULD STOP HERE
# Don't chase the extreme mean and sd! this is a prior; data will pull it back.
# Just run SBC:


## Simulation-based calibration (SBC), via the SBC package -----------------------
# Template for the project-wide rollout (src/utils/sbc_backend_ulam.R). Every
# raw parameter gets a rank/coverage check automatically -- no extra_estimands
# list needed. b[Lo]/tr[Tr]/yr[Yr] stay out of `variables` on purpose: the
# simulator draws fresh location/tree/year offsets from sigma_loc/sigma_tr,
# not from these prior-drawn arrays, so their "true" values wouldn't match
# what actually generated the data.
source('src/utils/sbc_backend_ulam.R')
library(SBC)
future::plan(future::multisession)  # compute_SBC() then uses all available cores automatically

generate_one_MDLS2v <- function(){
  # extract.prior() does its own tiny internal sampling run and prints the
  # same divergence-style chatter as a real fit -- silence it here too,
  # same reasoning as the compute_SBC() call below.
  true_params <- suppressMessages(suppressWarnings(draw_true(extract.prior(fit_sim, n = 1, refresh = 0), 1)))[
    c("loga", "s_conv", "gap_shift", "ls0", "ls_Mg", "ls_Mo", "ls_MgMo", "sigma_loc", "sigma_tr")]
  dat <- simulate_from_priors_MDLS2v(true_params, shift = 1)
  list(variables = true_params,
       generated = as.list(dat[, c("Dv", "Mg", "Lo", "Mo", "Yr", "Tr")]))
}

# The three contrasts 5b.4's own analysis script reports (May gap, July
# gap, and the seasonal change between them) -- rdo() can't reference one
# dquant from another, so each is fully self-contained rather than reusing
# a shared subexpression.
dq_MDLS2v <- derived_quantities(
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

# Tried parallelizing this via future.chunk.size -- reverted. Its worker
# processes don't have `rethinking` (or draw_true()/fit_sim(), etc.)
# available even after listing future.globals explicitly ("could not find
# function extract.prior" even once draw_true() was fixed): every package
# the generator's call chain touches would need explicit future.packages=
# too. Not worth chasing for a step that isn't the bottleneck -- the ~100
# real refits below dominate wall-clock time by orders of magnitude more
# than serially generating 100 small datasets ever will.
datasets_MDLS2v <- generate_datasets(SBC_generator_function(generate_one_MDLS2v), n_sbc)

# refresh = 0: cmdstanr's own default per-iteration progress printing,
# multiplied across 100 refits, would flood the log far more than anything
# suppressMessages()/suppressWarnings() below are doing.
backend_MDLS2v <- SBC_backend_ulam(model_MDLS2v_16S, chains = 2, iter = n_iter,
                                    refresh = 0, control = list(adapt_delta = 0.99))

# Divergence/treedepth/E-BFMI chatter from cmdstanr comes through as
# message() conditions (yes, even though the text itself says "Warning:"
# -- confirmed empirically, suppressWarnings() alone does not catch it),
# plus whatever genuine warning()s show up elsewhere (e.g. ks.test ties).
# Safe to blanket-silence here because the same information is already
# captured numerically in $backend_diagnostics (n_divergent/
# n_max_treedepth/n_low_ebfmi) and compute_SBC()'s own end-of-run summary
# isn't part of this stream, so nothing is actually lost, only decluttered.
sbc_MDLS2v <- compute_SBC(
  datasets_MDLS2v, backend_MDLS2v, dquants = dq_MDLS2v,
  cache_mode = "results", cache_location = file.path(hiermod_out_dir, "sbc_cache_MDLS2v"),
  # Required for a custom backend under future::multisession -- worker
  # sessions don't otherwise see S3 methods defined outside the SBC package.
  globals = c("SBC_fit.SBC_backend_ulam", "SBC_fit_to_draws_matrix.ulam",
              "SBC_fit_to_diagnostics.ulam"))

# Rank histogram / ECDF-diff / coverage, one panel per parameter+contrast
p_sbc_rank  <- plot_rank_hist(sbc_MDLS2v)
p_sbc_ecdf  <- plot_ecdf_diff(sbc_MDLS2v)
p_sbc_cover <- plot_coverage(sbc_MDLS2v)

sbc_step <- paste0(model_id, "_", n_sbc, "sbc_iter")
save_gg("SBC_rank_hist", sbc_step, p_sbc_rank, width = 9, height = 7)
save_gg("SBC_ecdf_diff", sbc_step, p_sbc_ecdf, width = 9, height = 7)
save_gg("SBC_coverage",  sbc_step, p_sbc_cover, width = 9, height = 7)

# Divergences/treedepth/E-BFMI across all 100 refits
(sbc_diag_summary <- sbc_MDLS2v$backend_diagnostics %>%
   dplyr::summarise(total_divergent = sum(n_divergent), total_max_treedepth = sum(n_max_treedepth),
                     total_low_ebfmi = sum(n_low_ebfmi), n_replicates = dplyr::n()))

# _______ if need be, caclibrate more:


# # Which prior is driving the extreme tail? No stored sigma[cell] parameter
# # here (unlike MDLS2) -- derive the 4 cell sigmas from ls0/ls_Mg/ls_Mo/ls_MgMo
# # via the shared sigma_cell() helper (MDLS2v_model.R) instead.
# diagnose_extreme_tail(
#   prior_pred, value_col = "Dv_shifted",
#   candidates = list(
#     loga_max  = pmax(extracted_prior$loga[,1], extracted_prior$loga[,2]),
#     sigma_loc = extracted_prior$sigma_loc,
#     sigma_tr  = extracted_prior$sigma_tr,
#     sigma_max = pmax(
#       sigma_cell(1,1, extracted_prior$ls0, extracted_prior$ls_Mg, extracted_prior$ls_Mo, extracted_prior$ls_MgMo),
#       sigma_cell(1,2, extracted_prior$ls0, extracted_prior$ls_Mg, extracted_prior$ls_Mo, extracted_prior$ls_MgMo),
#       sigma_cell(2,1, extracted_prior$ls0, extracted_prior$ls_Mg, extracted_prior$ls_Mo, extracted_prior$ls_MgMo),
#       sigma_cell(2,2, extracted_prior$ls0, extracted_prior$ls_Mg, extracted_prior$ls_Mo, extracted_prior$ls_MgMo))
#   ))
# 
# # Biggest driver are the priors on sigma (cor is 0.9)
# # Shrink them, retry!
# 
# ## Calibration 1 ---------- 
# 
# model_ppc1 <- model
# model_ppc1$prior_ls0 <- quote(ls0        ~ dnorm(0,0.5))    # don't tighten baseline as much
# model_ppc1$prior_ls_Mg<- quote(ls_Mg     ~ dnorm(0,0.4))
# model_ppc1$prior_ls_Mo<- quote(ls_Mo     ~ dnorm(0,0.4))
# model_ppc1$prior_ls_MgMo<- quote(ls_MgMo ~ dnorm(0,0.4))
# 
# fit_sim_ppc1 <- ulam(
#   model_ppc1,
#   data = as.list(dat_sim),
#   chains = 6, cores = 6, iter = 5000,
#   control = list(adapt_delta = 0.99)
# )
# 
# precis(fit_sim_ppc1, depth = 2)
# 
# post_sim_ppc1 <- extract.samples(fit_sim_ppc1)
# 
# ### Fixed effects recovery -------------
# 
# (coef_recovery_ppc1 <- check_recovery(
#   true = true_coefs,
#   
#   # Parameter recovery:
#   post_draws = list(
#     loga1 = post_sim_ppc1$loga[,1],
#     loga2 = post_sim_ppc1$loga[,2],
#     s_conv = post_sim_ppc1$s_conv,
#     gap_shift = post_sim_ppc1$gap_shift,
#     ls0 = post_sim_ppc1$ls0, ls_Mg = post_sim_ppc1$ls_Mg,
#     ls_Mo = post_sim_ppc1$ls_Mo, ls_MgMo = post_sim_ppc1$ls_MgMo)) )
# 
# 
# extracted_prior_ppc1 <- extract.prior(fit_sim_ppc1, n = n_prior)
# 
# prior_pred_ppc1 <- map_dfr(seq_len(n_prior), function(i){
#   simulate_from_priors_MDLS2v(draw_true(extracted_prior_ppc1, i), shift = 1)
# }, .id = "draw")
# 
# p_prior_pc_ppc1 <- prior_predictive_spaghetti(
#   prior_pred_ppc1, value_col = "Dv_shifted", upper_q = 0.99, model = model_ppc1,
#   title = "Prior predictive check", observed = dat_sim$Dv_shifted) ; p_prior_pc_ppc1
# 
# save_gg("sim_prior_PC", model_id, p_prior_pc)
# 
# ## Extreme values check --------------
# diagnose_extreme_tail(
#   prior_pred_ppc1, value_col = "Dv_shifted",
#   candidates = list(
#     loga_max  = pmax(extracted_prior_ppc1$loga[,1], extracted_prior_ppc1$loga[,2]),
#     sigma_loc = extracted_prior_ppc1$sigma_loc,
#     sigma_tr  = extracted_prior_ppc1$sigma_tr,
#     sigma_max = pmax(
#       sigma_cell(1,1, extracted_prior_ppc1$ls0, extracted_prior_ppc1$ls_Mg, extracted_prior_ppc1$ls_Mo, extracted_prior_ppc1$ls_MgMo),
#       sigma_cell(1,2, extracted_prior_ppc1$ls0, extracted_prior_ppc1$ls_Mg, extracted_prior_ppc1$ls_Mo, extracted_prior_ppc1$ls_MgMo),
#       sigma_cell(2,1, extracted_prior_ppc1$ls0, extracted_prior_ppc1$ls_Mg, extracted_prior_ppc1$ls_Mo, extracted_prior_ppc1$ls_MgMo),
#       sigma_cell(2,2, extracted_prior_ppc1$ls0, extracted_prior_ppc1$ls_Mg, extracted_prior_ppc1$ls_Mo, extracted_prior_ppc1$ls_MgMo))
#   ))
# 
# 
# ## Calibration 2 --------------------------------
# 
# model_ppc2 <- model_ppc1
# model_ppc2$prior_loga <- quote(loga[Mg] ~ dnorm(4.5, 2))
# model_ppc2$pr_sigma_tr <- quote(sigma_tr ~ dexp(2.5))
# 
# fit_sim_ppc2 <- ulam(
#   model_ppc2,
#   data = as.list(dat_sim),
#   chains = 6, cores = 6, iter = 5000,
#   control = list(adapt_delta = 0.99)
# )
# 
# (precis_ppc2 <- precis(fit_sim_ppc2, depth = 2))
# 
# post_sim_ppc2 <- extract.samples(fit_sim_ppc2)
# 
# ### Fixed effects recovery -------------
# 
# (coef_recovery_ppc2 <- check_recovery(
#   true = true_coefs,
#   
#   # Parameter recovery:
#   post_draws = list(
#     loga1 = post_sim_ppc2$loga[,1],
#     loga2 = post_sim_ppc2$loga[,2],
#     s_conv = post_sim_ppc2$s_conv,
#     gap_shift = post_sim_ppc2$gap_shift,
#     ls0 = post_sim_ppc2$ls0, ls_Mg = post_sim_ppc2$ls_Mg,
#     ls_Mo = post_sim_ppc2$ls_Mo, ls_MgMo = post_sim_ppc2$ls_MgMo)) )
# 
# extracted_prior_ppc2 <- extract.prior(fit_sim_ppc2, n = n_prior)
# 
# prior_pred_ppc2 <- map_dfr(seq_len(n_prior), function(i){
#   simulate_from_priors_MDLS2v(draw_true(extracted_prior_ppc2, i), shift = 1)
# }, .id = "draw")
# 
# p_prior_pc_ppc2 <- prior_predictive_spaghetti(
#   prior_pred_ppc2, value_col = "Dv_shifted", upper_q = 0.99, model = model_ppc2,
#   title = "Prior predictive check", observed = dat_sim$Dv_shifted) ; p_prior_pc_ppc2
# 
# save_gg("sim_prior_PC", model_id, p_prior_pc)
# 
# ## Extreme values check --------------
# diagnose_extreme_tail(
#   prior_pred_ppc2, value_col = "Dv_shifted",
#   candidates = list(
#     loga_max  = pmax(extracted_prior_ppc2$loga[,1], extracted_prior_ppc2$loga[,2]),
#     sigma_loc = extracted_prior_ppc2$sigma_loc,
#     sigma_tr  = extracted_prior_ppc2$sigma_tr,
#     sigma_max = pmax(
#       sigma_cell(1,1, extracted_prior_ppc2$ls0, extracted_prior_ppc2$ls_Mg, extracted_prior_ppc2$ls_Mo, extracted_prior_ppc2$ls_MgMo),
#       sigma_cell(1,2, extracted_prior_ppc2$ls0, extracted_prior_ppc2$ls_Mg, extracted_prior_ppc2$ls_Mo, extracted_prior_ppc2$ls_MgMo),
#       sigma_cell(2,1, extracted_prior_ppc2$ls0, extracted_prior_ppc2$ls_Mg, extracted_prior_ppc2$ls_Mo, extracted_prior_ppc2$ls_MgMo),
#       sigma_cell(2,2, extracted_prior_ppc2$ls0, extracted_prior_ppc2$ls_Mg, extracted_prior_ppc2$ls_Mo, extracted_prior_ppc2$ls_MgMo))
#   ))
