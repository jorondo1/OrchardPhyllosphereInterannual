# 5b.2_MDLS2_shifted_validation.R -- MODEL 5, SHIFTED -- identical structure
# to Model 5, but the likelihood is fit to (Hill_1 - 1) instead of Hill_1,
# since Hill numbers can't go below 1. No separate "5b.1" model file --
# the alist itself never changes for the shifted variant (see means_MDLS2()/
# sim_div_MDLS2()'s own shift= argument), so this sources 5.1_MDLS2_model.R
# directly.

source('src/hiermod/ITS/0_SETUP.R')
source('src/hiermod/ITS/5.1_MDLS2_model.R') # model, means_MDLS2(), sim_div_MDLS2(), contrast_may_gap_MDLS2(), simulate_from_priors()
hiermod_out_dir <- "out/hiermod/ITS_5_lognormal_MDLS2_shifted"

## Model specification ---------------------------------------------------------

# Same as 5.1_MDLS2_model.R: likelihood stays `Dv ~ dlnorm(mu, sigma[cell])`:
# Dv is meant to be the *unshifted* (0, Inf)-support lognormal draw, and
# that's true whether or not a floor exists in the real world -- the floor
# only has to be added back (or subtracted out of real data) outside the
# likelihood. Concretely, the only two differences from plain Model 5:
#   - sim_div_MDLS2()/simulate_from_priors() called with shift = 1, so
#     dat_sim$Dv_shifted (= 1 + Dv) is also produced, purely for reporting/
#     comparison -- never fed to the likelihood.
#   - the real-data `dat$Dv` is built as `div$Hill_1 - 1` (subtracting the
#     floor to get back to (0, Inf) support), instead of Model 5's raw
#     `div$Hill_1`.
# means_MDLS2()'s `shift = 1` argument (threaded through post_full()'s ...)
# adds the floor back only when backtransforming to a raw mean/median.
# contrasts don't need it since a shared additive constant cancels
# in any group1-group2 difference.

## Parameter recovery -----------------------------------------------------------

# Same as model 5
may_conv <- 4          # hill scale
may_org <- 7           # hill scale
july_conv_shift <- 0.35 # log scale
july_org_shift <- 0.2  # log scale

true_sigma <- cv_to_sigma(c(0.6, 0.25, 0.8, 0.8))

dat_sim <- sim_div_MDLS2(
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

# Dv is what the model actually fits; Dv_shifted (= 1 + Dv) is only for
# sanity-checking that the floor looks right; its NOT passed to ulam().
hist(dat_sim$Dv_shifted, breaks = 30)

fit_sim <- ulam(
  model,
  data = as.list(dat_sim),
  chains = 6, cores = 6, iter = 10000,
  control = list(adapt_delta = 0.99)
)

precis(fit_sim, depth = 2)

### Fixed effects recovery -------------

post_sim <- extract.samples(fit_sim)

fixed_recovery <- check_recovery(
  true = c(loga1 = log(may_conv),
           loga2 = log(may_org),
           s_conv = july_conv_shift,
           gap_shift = july_org_shift),
  post_draws = list(
    loga1 = post_sim$loga[,1],
    loga2 = post_sim$loga[,2],
    s_conv = post_sim$s_conv,
    gap_shift = post_sim$gap_shift)
); fixed_recovery

### Sigma recovery -------------

sigma_recovery <- check_recovery(
  true = list(sigma1 = true_sigma[1], sigma2 = true_sigma[2],
              sigma3 = true_sigma[3], sigma4 = true_sigma[4]),
  post_draws = list(sigma1 = post_sim$sigma[,1], sigma2 = post_sim$sigma[,2],
                    sigma3 = post_sim$sigma[,3], sigma4 = post_sim$sigma[,4])
); sigma_recovery # One sigma is a bit off, nothing to panick about

### Contrast recovery ------------------------------

# shift = 1 here: raw mean/median cell values should reflect the floor.
# Contrasts wouldn't need it (cancels), but the cell-level columns
# (median_1..median_4 below) do.
pf_sim <- post_full(fit_sim, means_MDLS2, shift = 1)
names(pf_sim)

# mean_1=conv_May, mean_2=conv_July, mean_3=org_May, mean_4=org_July:
m_sim  <- pf_sim$median

pc_estimands_sim <- estimand_rows(list(
  "Median May gap (Organic - Conventional)"  = m_sim$median_3 - m_sim$median_1,
  "Median July gap (Organic - Conventional)" = m_sim$median_4 - m_sim$median_2,
  "Seasonal change in median gap (July - May)"= (m_sim$median_4 - m_sim$median_2) - (m_sim$median_3 - m_sim$median_1)
))

may_gap <- may_org-may_conv
july_gap <- may_org*exp(july_conv_shift+july_org_shift)-may_conv*exp(july_conv_shift)

true_estimands <- tribble(
  ~statistic,                                   ~value,
  "Median May gap (Organic - Conventional)",    may_gap,
  "Median July gap (Organic - Conventional)",   july_gap,
  "Seasonal change in median gap (July - May)", july_gap-may_gap,
)

p_sim_contrast <- contrast_plot_panels(
  pc_estimands_sim, quant = c(0, 0.995), scales = 'free_y',
  group_pal = Management_palette,
  true_vals = true_estimands); p_sim_contrast

save_report("sim_summary", "MDLS2_shifted", fit_sim, pc_estimands_sim, model,
            recovery = bind_rows(fixed_recovery, sigma_recovery), model_name = "The Floor Raiser")
save_gg("sim_contrast_density", "MDLS2_shifted", p_sim_contrast)

## Prior predictive check -------------------------------------------------------

n_prior <- 1000
extracted_prior <- extract.prior(fit_sim, n = n_prior)

# shift = 1 here too, so prior_pred$Dv_shifted matches what real Dv would look like.
prior_pred <- map_dfr(seq_len(n_prior), function(i){
  simulate_from_priors(draw_true(extracted_prior, i), shift = 1)
}, .id = "draw")

summary(prior_pred$Dv_shifted)

p_prior_pc <- prior_predictive_spaghetti(
  prior_pred, value_col = "Dv_shifted", upper_q = 0.99, model = model,
  title = "Prior predictive check", observed = dat_sim$Dv_shifted) ; p_prior_pc

save_gg("sim_prior_PC", "MDLS2_shifted", p_prior_pc)

## Simulation-based calibration (SBC) --------------------------------------------

# contrast_may_gap_MDLS2() (5.1_MDLS2_model.R) needs no shift-awareness
# at all -- it's a contrast (org_May - conv_May), and the shift cancels.
# simulate_fn wraps simulate_from_priors() with shift = 1 baked in, since
# run_sbc() only ever calls it as simulate_fn(true_params).

sbc_MDLS2_shifted <- run_sbc(
  model_fit   = fit_sim,
  means_fn    = means_MDLS2,
  contrast_fn = contrast_may_gap_MDLS2,
  simulate_fn = function(true_params) simulate_from_priors(true_params, shift = 1),
  n_sbc = 30, iter = 15000, n_parallel = 4, chains = 2, cores = 2,
  control = list(adapt_delta = 0.99))

sbc_out_MDLS2_shifted <- save_sbc_report(sbc_MDLS2_shifted, "MDLS2_shifted_30sbc_iter")
hist(sbc_out_MDLS2_shifted$ranks, breaks = 30)

## Model fit ----------------------------------------------------------------

dat <- list(
  Dv = div$Hill_1 - 1, # subtract the floor because Dv must be (0, Inf)-support to match the likelihood
  Mg = idx$Mg$to_index(div$Management),
  Lo = idx$Lo$to_index(div$Location),
  Mo = idx$Mo$to_index(div$Time),
  Yr = idx$Yr$to_index(div$Year),
  Tr = idx$Tr$to_index(div$Tree_id)
)
dat$cell <- (dat$Mg - 1) * 2 + dat$Mo   # 1=Conv-May, 2=Conv-July, 3=Org-May, 4=Org-July

fit_MDLS2_shifted <- ulam(
  model,
  data = dat,
  chains = 6, cores = 6, iter = 20000,
  control = list(adapt_delta = 0.99)
)
save_fit("fit", "MDLS2_shifted", fit_MDLS2_shifted)

precis(fit_MDLS2_shifted, depth = 2)

save_pdf("fit_trankplot", "MDLS2_shifted", function() trankplot(fit_MDLS2_shifted, n_cols = 6, max_rows = 10))

## Posterior predictive check --------------------------------------------------

### Overall, by Management x Season cell ----
# dat$Dv is on the same (Hill_1 - 1) scale the model was fit to, and sim()
# generates yrep on that same scale -- no +1 correction needed anywhere in
# this section, since both sides were shifted down by 1 consistently.

pp_group <- interaction(idx$Mg$to_label(dat$Mg), idx$Mo$to_label(dat$Mo), sep = " ")
(p_postpred <- plot_ppc_overlay(fit_MDLS2_shifted, dat, pp_group, xlim = c(0,150)))
save_gg("postpred_density", "MDLS2_shifted", p_postpred)

### Contrast test statistics ----

(p_ppc <- plot_ppc_season_contrast_stats(fit_MDLS2_shifted, dat))

save_gg("postpred_stat", "MDLS2_shifted", p_ppc)
