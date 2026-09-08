# 4.2_MDLS_validation.R -- MODEL 4 (Management effect by Season): parameter
# recovery, prior-predictive check, SBC, the real fit, and PPC.

source('src/hiermod/ITS/0_SETUP.R')
source('src/hiermod/ITS/4.1_MDLS_model.R') # model, means_MDLS(), sim_div_MDLS(), contrast_may_gap_MDLS(), simulate_from_priors()
hiermod_out_dir <- "out/hiermod/ITS_4_lognormal_MDLS"

## Model specification ---------------------------------------------------------
# THREE ESTIMANDS: May gap, July gap, seasonal change in gap (July - May).
# Tree random effect added (same trees sampled May+July -- pseudoreplication
# otherwise); Location-varying gap (Model 3) dropped to isolate Season; Year
# enters as a fixed effect (unbalanced sampling across years). Full why: see
# MODEL_HISTORY.md.
#
# gamma <- s + gap_shift*(Mg-1): McElreath-style two-equation interaction
# (Statistical Rethinking 7.1.2) -- gap_shift *is* the interaction estimand
# directly, not derived after the fact.
#
# mean has 4 columns here (no single group1/group2 contrast), so the
# contrast tibble/plot are built by hand in 4.3_MDLS_analysis.R, not via
# compute_contrasts().

## Parameter recovery -----------------------------------------------------------

# Let's simulate data where the baseline div is
# May Conventional:   4
# May Organic:        7
# July Conventional:  4*exp(0.35) = 5.6
# July Organic:       7*exp(0.55) = 12.1
# Therefore:
#   May gap   = 3
#   July gap  = 6.5
#   Gap shift = 3.5 (July - May)

may_conv <- 4          # hill scale
may_org <- 7           # hill sncale
july_conv_shift <- 0.35 # log scale
july_org_shift <- 0.2  # log scale

dat_sim <- sim_div_MDLS(
  N_samples = 240,
  n_loc = 4,
  # true values, in the model's own units (log-median for loga, log-shift for s/gap_shift):
  loga = log(c(may_conv, may_org)),   # passing log directly
  s_conv = july_conv_shift,           # Conventional: May->July grows  MULTIPLICATIVE
  gap_shift = july_org_shift,         # Organic's own May->July shift is s+gap_shift (increases):
  sigma = cv_to_sigma(c(0.6, 0.8)),   # Same management-specific variance within each season
  year_offset = c(0, 0.3, -0.2),  # in dnorm space
  p_dropout = 0.1
); head(dat_sim)


# View design unbalance
dat_sim %>%
  count(Lo, Yr, Mo) %>% print(n=100)

hist(dat_sim$Dv)

fit_sim <- ulam(
  model,
  data = as.list(dat_sim),
  chains = 6, cores = 6, iter = 5000 )

precis(fit_sim, depth = 2)

### Fixed effects recovery -------------

post_sim <- extract.samples(fit_sim)

fixed_recovery <- check_recovery(
  true = c(loga1 = log(may_conv),
           loga2 = log(may_org),
           s_conv = july_conv_shift,
           gap_shift = july_org_shift),
  post_draws = list(loga1 = post_sim$loga[,1], loga2 = post_sim$loga[,2],
                     s_conv = post_sim$s_conv, gap_shift = post_sim$gap_shift)
); fixed_recovery

### Contrast recovery ------------------------------

# means_MDLS() checks if params are picked up . Let's verify if the contrasts
# are picked up,  which is what we aim for.

# So here, we plot the whole posterior distributions of these estimands
# and check where the true (synthetic) data contrasts actually fall:

# list of one tibble per parameter, with one coloumn per group when applicable
pf_sim <- post_full(fit_sim, means_MDLS)
names(pf_sim) # all lists

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

save_report("sim_summary", "MDLS", fit_sim, pc_estimands_sim, model, recovery = fixed_recovery, model_name = "The Season Ticket")
save_gg("sim_contrast_density", "MDLS", p_sim_contrast)

## Prior predictive check -------------------------------------------------------

# Worth doing even though it's more setup than models 1-3: there are three
# new fixed-effect parameters (s_conv, gap_shift, yr[Yr]) that can each push
# implied Dv in ways sigma[Mg]/sigma_loc alone never had to account for.

n_prior <- 1000
extracted_prior <- extract.prior(fit_sim, n = n_prior)

# simulate_from_priors() is sourced from 4.1_MDLS_model.R above.
# Simulate from priors
prior_pred <- map_dfr(seq_len(n_prior), function(i){
  simulate_from_priors(draw_true(extracted_prior, i))
}, .id = "draw")

summary(prior_pred$Dv) # 3rd quartile around 40, mean ~330

p_prior_pc <- prior_predictive_spaghetti(
  prior_pred, upper_q = 0.99, model = model,
  title = "Prior predictive check", observed = dat_sim$Dv) ; p_prior_pc

save_gg("sim_prior_PC", "MDLS", p_prior_pc)

# all seems reasonable but which parameters are driving the tail?
# Flag draws with extreme values only
extreme_draws <- prior_pred %>%
  group_by(draw) %>%
  summarise(max_Dv = max(Dv)) %>%
  filter(max_Dv > quantile(prior_pred$Dv, 0.99)) %>%
  pull(draw) %>% as.integer()

tibble(
  draw = seq_len(n_prior),
  extreme = seq_len(n_prior) %in% extreme_draws,
  sigma_max = pmax(extracted_prior$sigma[,1], extracted_prior$sigma[,2]),
  sigma_loc = extracted_prior$sigma_loc,
  sigma_tr  = extracted_prior$sigma_tr,
  gamma_max = pmax(abs(extracted_prior$s_conv + extracted_prior$gap_shift), abs(extracted_prior$s_conv))  # unshrunk terms
) %>%
  group_by(extreme) %>%
  summarise(across(where(is.numeric), median))

# Nothing stands out (if anything, sigma is the biggest in non-extremes!)

## Simulation-based calibration (SBC) --------------------------------------------

# contrast_may_gap_MDLS() is sourced from 4.1_MDLS_model.R above --
# needed because means_MDLS()'s `mean` has 4 columns, so run_sbc()'s default
# contrast_from_means() (mean[,2]-mean[,1]) would silently compare
# conv_July-conv_May against an unrelated "true" Mg gap.

# Heads up on cost: this model has 129 Tree-level parameters, so n_sbc=100 at
# iter=10000 (MDLb's setting) may be too slow to be worth it here.

sbc_MDLS <- run_sbc(
  model_fit   = fit_sim,
  means_fn    = means_MDLS,
  contrast_fn = contrast_may_gap_MDLS,
  simulate_fn = simulate_from_priors,
  n_sbc = 100, iter = 5000, n_parallel=4, chains=2, cores=2,
  control = list(adapt_delta = 0.99))

sbc_out_MDLS <- save_sbc_report(sbc_MDLS, "MDLS_100sbc_iter")
hist(sbc_out_MDLS$ranks, breaks = 30)

# More concentrated in the middle (~60% in th emiddle 30% range).
# U-shapes mean posterior is too narrow / overconfident, missing the true
# value too often. Bells show underconfidence (posteriors are a bit too
# wide) but that's ok.

# Check aliasing, aka the collinearity counterpart for variances:
cor( post_sim$sigma_tr, post_sim$sigma)
    # 0.39 not a red flag; sigma[Mg] is within-tree noise and sigma_tr is between-tree spread
cor(post_sim$sigma_tr, post_sim$sigma_loc)
cor(post_sim$gap_shift, post_sim$sigma_tr)

# We move on!

## Model fit ----------------------------------------------------------------

dat <- list(
  Dv = div$Hill_1,
  Mg = idx$Mg$to_index(div$Management),
  Lo = idx$Lo$to_index(div$Location),
  Mo = idx$Mo$to_index(div$Time),
  Yr = idx$Yr$to_index(div$Year),
  Tr = idx$Tr$to_index(div$Tree_id)
)

fit_MDLS <- ulam(
  model,
  data = dat,
  chains = 6, cores = 6, iter = 10000,
  control = list(adapt_delta = 0.99)
)
save_fit("fit", "MDLS", fit_MDLS)

precis(fit_MDLS, depth = 2)

traceplot(fit_MDLS); trankplot(fit_MDLS)
save_pdf("fit_traceplot", "MDLS", function() traceplot(fit_MDLS))
save_pdf("fit_trankplot", "MDLS", function() trankplot(fit_MDLS))

## Posterior predictive check --------------------------------------------------

### Overall, by Management x Season cell ----
# ppc_dens_overlay_grouped() (inside plot_ppc_overlay()) isn't limited to a
# binary group, so the 4-cell interaction can be shown in one call.

pp_group <- interaction(idx$Mg$to_label(dat$Mg), idx$Mo$to_label(dat$Mo), sep = " ")
(p_postpred <- plot_ppc_overlay(fit_MDLS, dat, pp_group, xlim = c(0,150)))
save_gg("postpred_density", "MDLS", p_postpred)

### Contrast test statistics ----
# contrast_stat()/plot_ppc_contrast_stat() only handle a single binary group
# (dat$Mg); with a 4-cell Mg x Mo design, write the stat functions directly.

(p_ppc <- plot_ppc_season_contrast_stats(fit_MDLS, dat))

save_gg("postpred_stat", "MDLS", p_ppc)
