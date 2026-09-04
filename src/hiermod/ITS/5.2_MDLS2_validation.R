# 5.2_MDLS2_validation.R -- MODEL 5 (sigma[Mg] -> sigma[cell]): parameter
# recovery, prior-predictive check, SBC, the real fit, and PPC.

source('src/hiermod/ITS/0_SETUP.R')
source('src/hiermod/ITS/5.1_MDLS2_model.R') # model, means_MDLS2(), sim_div_MDLS2(), contrast_may_gap_MDLS2(), simulate_from_priors()
hiermod_out_dir <- "out/hiermod/ITS_5_lognormal_MDLS2"

## Model specification ---------------------------------------------------------

# Same three estimands as MODEL 4 (May gap, July gap, seasonal change in
# gap) and the same random-effects structure (Location, Year, Tree). The
# only change: MODEL 4's postpred density overlay showed Conventional-July
# tightly clustered near a low value while its yrep was visibly wider and
# shifted right -- sigma[Mg] was shared across both months within a
# Management group, forcing one spread to fit two seasons that likely have
# different true variance. Here, sigma varies by Management AND Month (4
# cells) via a combined index `cell = (Mg-1)*2 + Mo` -- see
# 5.1_MDLS2_model.R for why (ulam doesn't support double-indexing a
# dexp-distributed scale parameter) and for the model/means_fn/sim_div_fn
# definitions themselves.

## Parameter recovery -----------------------------------------------------------

#
# sigma now 4 values (cell order: conv_May, conv_July, org_May, org_July).
# Deliberately picked to differ a lot: conv_July tighter than conv_May,
# so here we test whether 4 distinguishable cells come back
# distinguishable, not just plausible.

# Same baseline diversity/gap values as MODEL 4:
may_conv <- 4          # hill scale
may_org <- 7           # hill scale
july_conv_shift <- 0.35 # log scale
july_org_shift <- 0.2  # log scale

true_sigma <- cv_to_sigma(c(0.6, 0.25, 0.8, 0.8)) # conv_May, conv_July, org_May, org_July

# Pinned so re-runs are comparable -- without it, every run draws a fresh
# dat_sim (new Lo x Yr grid, new loc/tree offsets), and diagnostics
# (divergences/treedepth/E-BFMI) can vary a lot run to run just from that,
# not from anything about the model itself.
set.seed(20260905)

dat_sim <- sim_div_MDLS2(
  N_samples = 240,
  n_loc = 4,
  loga = log(c(may_conv, may_org)),
  s_conv = july_conv_shift,
  gap_shift = july_org_shift,
  sigma = true_sigma,
  year_offset = c(0, 0.3, -0.2),
  p_dropout = 0.1
); head(dat_sim)

# View design unbalance
dat_sim %>%
  count(Lo, Yr, Mo, Mg) %>% print(n=100); hist(dat_sim$Dv)

fit_sim <- ulam(
  model,
  data = as.list(dat_sim),
  chains = 6, cores = 6, iter = 10000 )
# fewer iterations = more divergencres (e.g. 5000 = 5% divergences; here essentially 0%)
# 10,000 iterations works

precis(fit_sim, depth = 2)

post_sim <- extract.samples(fit_sim)

### Divergences check ------- 

# Divergences (red points)/max-treedepth/E-BFMI showed up together here,
# and got worse (not better) at higher iter -- that combination points to
# real posterior geometry (sigma_loc/sigma_tr/sigma[cell] fighting for the
# same variance, flagged as a known Model 5 weakness in MODEL_HISTORY.md),
# not an under-sampling problem. This pairs() plot is how to tell which
# parameters are actually implicated:
#  - divergent (red) points clustered in a specific corner -- e.g. small
#    sigma_loc or sigma_tr with the other parameter free to roam -- is the
#    classic funnel signature, meaning even the non-centered form isn't
#    fully escaping it here.
#  - a tight diagonal ridge between any two of these (not just clumped
#    points, an actual correlated *shape*) means those two are only weakly
#    separately identified by the data -- the sampler has to explore a
#    narrow degenerate slice, which is exactly what produces low E-BFMI.
#  - watch sigma[2] (conv_July) specifically -- it's the deliberately
#    tightest true cell (cv=0.25) with the fewest effective obs (~60), the
#    most likely one to be hard to pull apart from sigma_loc/sigma_tr.
#  - a skewed/heavy marginal density on the diagonal for any one parameter
#    can drive divergences on its own, independent of any pairwise
#    correlation with another.

# cs$draws() (not extract.samples(), which flattens chains together) keeps
# the iteration x chain x variable structure bayesplot wants -- feeding it
# a flattened data.frame is what triggered the "only one chain" warning.
# thin_for_pairs() (hiermod_core.R) caps the actual point count at 5000
# (keeping every divergent draw) instead of asking ggplot to render all
# ~45,000 -- that's what produced a 40MB PDF; points overlap heavily well
# before 5000 anyway, so nothing is visually lost. png on top of that is
# belt-and-suspenders, not the fix itself.

pairs_vars <- c("sigma_loc","sigma_tr","sigma[1]","sigma[2]","sigma[3]","sigma[4]")
cs <- attr(fit_sim, "cstanfit")
thinned <- thin_for_pairs(cs, pairs_vars)
p_pairs <- bayesplot::mcmc_pairs(thinned$draws, np = thinned$np)
save_gg("fit_pairs", "MDLS2_sim_variance", p_pairs, width = 14, height = 14, type = "png")

### Fixed effects recovery -------------

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
); fixed_recovery # s_conv doesn't 

### Sigma recovery -------------
# The actual new thing this model tests: are the 4 cells identifiable at
# this design's per-cell sample sizes, not just plausible in aggregate?

sigma_recovery <- check_recovery(
  true = list(sigma1 = true_sigma[1], sigma2 = true_sigma[2],
              sigma3 = true_sigma[3], sigma4 = true_sigma[4]),
  post_draws = list(sigma1 = post_sim$sigma[,1], sigma2 = post_sim$sigma[,2],
                    sigma3 = post_sim$sigma[,3], sigma4 = post_sim$sigma[,4])
); sigma_recovery

# Seems ok but again sigma2 is really on the margin of the 89%

### Contrast recovery ------------------------------

# means_MDLS2() checks if params are picked up. Let's verify if the
# contrasts are picked up, which is what we aim for: plot the whole
# posterior distributions of these estimands and check where the true
# (synthetic) data contrasts actually fall.

pf_sim <- post_full(fit_sim, means_MDLS2)
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

save_report("sim_summary", "MDLS2", fit_sim, pc_estimands_sim, model,
            recovery = bind_rows(fixed_recovery, sigma_recovery), model_name = "The Splitter")
save_gg("sim_contrast_density", "MDLS2", p_sim_contrast)

## Prior predictive check -------------------------------------------------------

n_prior <- 1000
extracted_prior <- extract.prior(fit_sim, n = n_prior)

# Simulate from priors
prior_pred <- map_dfr(seq_len(n_prior), function(i){
  simulate_from_priors(draw_true(extracted_prior, i))
}, .id = "draw")

summary(prior_pred$Dv)

p_prior_pc <- prior_predictive_spaghetti(
  prior_pred, upper_q = 0.99, model = model,
  title = "Prior predictive check", observed = dat_sim$Dv) ; p_prior_pc
# median is great, 3rd quartile in the same ballpark as before

save_gg("sim_prior_PC", "MDLS2", p_prior_pc)

## Simulation-based calibration (SBC) --------------------------------------------

sbc_MDLS2 <- run_sbc(
  model_fit   = fit_sim,
  means_fn    = means_MDLS2,
  contrast_fn = contrast_may_gap_MDLS2,
  simulate_fn = simulate_from_priors,
  n_sbc = 30, iter = 15000, n_parallel = 4, chains = 2, cores = 2,
  control = list(adapt_delta = 0.99))

(sbc_out_MDLS2 <- summarize_sbc(sbc_MDLS2))
save_sbc_report(sbc_out_MDLS2, "MDLS2_30sbc_iter")
hist(sbc_out_MDLS2$ranks, breaks = 30)

## Model fit ----------------------------------------------------------------

dat <- list(
  Dv = div$Hill_1,
  Mg = idx$Mg$to_index(div$Management),
  Lo = idx$Lo$to_index(div$Location),
  Mo = idx$Mo$to_index(div$Time),
  Yr = idx$Yr$to_index(div$Year),
  Tr = idx$Tr$to_index(div$Tree_id)
)
dat$cell <- (dat$Mg - 1) * 2 + dat$Mo   # 1=Conv-May, 2=Conv-July, 3=Org-May, 4=Org-July

fit_MDLS2 <- ulam(
  model,
  data = dat,
  chains = 6, cores = 6, iter = 20000,
  control = list(adapt_delta = 0.99)
)
save_fit("fit", "MDLS2", fit_MDLS2)

precis(fit_MDLS2, depth = 2)

save_pdf("fit_traceplot", "MDLS2", function() traceplot(fit_MDLS2))
save_pdf("fit_trankplot", "MDLS2", function() trankplot(fit_MDLS2))

## Posterior predictive check --------------------------------------------------

### Overall, by Management x Season cell ----

pp_group <- interaction(idx$Mg$to_label(dat$Mg), idx$Mo$to_label(dat$Mo), sep = " ")
(p_postpred <- plot_ppc_overlay(fit_MDLS2, dat, pp_group, xlim = c(0,150)))
save_gg("postpred_density", "MDLS2", p_postpred)

# Compare against fit_contrast_density_MDLS.pdf (model 4)
# Conventional-July's yrep hug the observed spike much more closely than model 4's did.

### Contrast test statistics ----

(p_ppc <- plot_ppc_season_contrast_stats(fit_MDLS2, dat))

save_gg("postpred_stat", "MDLS2", p_ppc)
