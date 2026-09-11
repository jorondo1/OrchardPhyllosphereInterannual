# 2.2_MDL_validation.R
# See how much priors need to be tightened
# informed by ITS codes

source('src/hiermod_16S/0_SETUP_16S.R')
source('src/hiermod_16S/2_MDL/2.1_MDL_16S_model.R') # model, means_MDL(), sim_div_ML(), simulate_from_priors()
hiermod_out_dir <- "out/hiermod/16S_2_lognormal_MDL"

# MODEL 2 -- Partial pooling across Location (non-centered) =================
#same rationale as ITS model 2

## Parameter recovery -----------------------------------------------------------
# Dummy values
conv <- 180
org <- 120
true_sigma <- cv_to_sigma(c(0.5, 0.8))# keep different variances
true_sigma_loc <- 0.5

dat_sim <- sim_div_ML(
  Mg = rbern(250)+1,
  # 4 locations for simulation; fewer creates a funnel (see below)
  Lo = sample(1:4, 250, replace = TRUE),
  loga = log(c(conv, org)),           # Difference of 60 in median, raw scale
  sigma = true_sigma,
  sigma_loc = true_sigma_loc) 

fit_sim <- ulam(
  model,
  data = as.list(dat_sim),
  chains = 6, cores = 6, iter = 5000 )
precis(fit_sim, depth = 2 )

### Fixed effectt recovery ---------

post_sim <- extract.samples(fit_sim)

check_recovery(
  true = c(loga1 = log(conv),
           loga2 = log(org)),
  post_draws = list(
    loga1 = post_sim$loga[,1],
    loga2 = post_sim$loga[,2])
) # all gooood

### Sigma recovery -------------

check_recovery(
  true = list(sigma1 = true_sigma[1], sigma2 = true_sigma[2]),
  post_draws = list(sigma1 = post_sim$sigma[,1], sigma2 = post_sim$sigma[,2])
) # all gooood

### Contrast recovery -------------

# means_MDL() already returns both mean and median

# True contrast: median needs no total_var (median = exp(loga), the raw
# conv/org values loga was built from); mean does (sigma[Mg]^2 +
# sigma_loc^2, matching means_MDL()'s own total_var.

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

save_report("sim_summary", "MDL", fit_sim, pc_sim, model, model_name = "The Wildcard")
save_gg("sim_contrast_density", "MDL", p_sim_contrast)

## Prior predictive checks -------------------------------------------------------

### ITS dataset priors (dnorm(2,2)) ------------------------

# Let's see what happens with the model we imported from ITS
# Likely wrong, because in general the distribution of diversity is a lot higher

n_prior <- 10000

prior <- extract.prior(fit_sim, n = n_prior)

# draw_true() is shared (sbc_helpers.R) -- generic slice of any extract.prior()
# output, no per-model rewrite needed. simulate_from_priors() is sourced from
# 2.1_MDL_model.R above.

# One simulated dataset per prior draw
prior_pred <- map_dfr(seq_len(n_prior), function(i){
  simulate_from_priors(draw_true(prior, i))
}, .id = "draw")

summary(prior_pred$Dv)  # raw per-observation draws
# Median suspiciously low
# Mean extremely high (trillions)
# sd a bit low

p_sim_spaghetti <- prior_predictive_spaghetti(
  prior_pred, upper_q = 0.99, model = model,
  title = "Prior predictive check: loga[Mg]~dnorm(2,2) (imported from ITS model))"); p_sim_spaghetti
# sd is in the billions! haha
save_gg("sim_prior_PC_default", "MDL", p_sim_spaghetti)

# So our model will produce draws with undeestimated medians and (highly) overestimated means


summary(div$Hill_1)
# Indeed, median is way to low (data is 50) and mean way too high (data is 145)
# Median is probably a loga prior problem;

# mean is probably a sigma problem, because multiplicative can randomly escalate.

div %>% group_by(Management) %>%
  summarise(mean_ = mean(Hill_1),
            sd_ = sd(Hill_1),
            median_ = median(Hill_1)) %>%
  mutate(
    logmean_ = log(mean_),
    logmedian_ = log(median_)
  )

# Which prior is driving the extreme tail?
(xlim_upper <- quantile(prior_pred$Dv, 0.99))
extreme_draws <- prior_pred %>%
  group_by(draw) %>%
  summarise(max_Dv = max(Dv)) %>%
  filter(max_Dv > xlim_upper) %>%
  pull(draw) %>% as.integer()

tibble(
  draw       = seq_len(n_prior),
  extreme    = seq_len(n_prior) %in% extreme_draws,
  sigma_loc  = prior$sigma_loc,
  sigma_max  = pmax(prior$sigma[,1], prior$sigma[,2])) %>%
  group_by(extreme) %>%
  summarise(across(c(sigma_loc, sigma_max), median))

# sigma max is quite extreme; sigma loc a little too

### Update loga prior --------------------------------------

# Sit halfway between log means , make sure means sit within 1sd of the mean
model_ppc1 <- model
model_ppc1$prior_ppc1 <- quote(loga[Mg] ~ dnorm(5,2)) #slightly more skeptical than ITS?
model_ppc1$prior_sigma <- quote(sigma[Mg] ~ dexp(2))
model_ppc1$prior_sigma_loc <- quote(sigma_loc ~ dexp(2))

fit_sim_ppc1 <- ulam(
  model_ppc1,
  data = as.list(dat_sim),
  chains = 6, cores = 6, iter = 5000 )
precis(fit_sim_ppc1, depth = 2 )

prior_ppc1 <- extract.prior(fit_sim_ppc1, n = n_prior)

# One simulated dataset per prior draw
prior_ppc1_pred <- map_dfr(seq_len(n_prior), function(i){
  simulate_from_priors(draw_true(prior_ppc1, i))
}, .id = "draw")

summary(prior_ppc1_pred$Dv) # raw per-observation draws
# Much better!
# Mean still very high (trillions), but 
# sd more reasonable and 4rd quartile also

p_sim_spaghetti_ppc1 <- prior_predictive_spaghetti(
  prior_ppc1_pred, upper_q = 0.99, model = model_ppc1, 
  title = "Prior predictive check: adapted priors)"); p_sim_spaghetti_ppc1

save_gg("sim_prior_PC_ppc1", "MDL", p_sim_spaghetti_ppc1)

# Visually, priors can easily explore the <2000 diversity 

## Parameter recovery: updated priors -----------------------------------------------------------

fit_sim_ppc1 <- ulam(
  model_ppc1,  # <- <- <- Updated model
  data = as.list(dat_sim),
  chains = 4, cores = 4, iter = 4000,
  control = list(adapt_delta = 0.99) )

precis(fit_sim_ppc1, depth = 2 ) # good r_hats

### Fixed effectt recovery ---------

post_sim_ppc1 <- extract.samples(fit_sim_ppc1)

check_recovery(
  true = c(loga1 = log(conv),
           loga2 = log(org)),
  post_draws = list(
    loga1 = post_sim_ppc1$loga[,1],
    loga2 = post_sim_ppc1$loga[,2])
) # all gooood

### Sigma recovery -------------

check_recovery(
  true = list(sigma1 = true_sigma[1], sigma2 = true_sigma[2]),
  post_draws = list(sigma1 = post_sim_ppc1$sigma[,1], sigma2 = post_sim_ppc1$sigma[,2])
) # all gooood

### Contrast recovery -------------

pf_sim_ppc1 <- post_full(fit_sim_ppc1, means_MDL)
pc_sim_ppc1 <- compute_contrasts(pf_sim_ppc1, keep = c("mean", "median"), group_levels = idx$Mg$levels)

p_sim_contrast_ppc1 <- contrast_plot_panels(
  pc_sim_ppc1, quant = c(0.00, 0.995), group_pal = Management_palette,
  true_vals = true_vals); p_sim_contrast_ppc1

save_report(
  "sim_summary", "MDL_ppc1", fit_sim_ppc1, pc_sim_ppc1, model_ppc1, model_name = "The Wildcard")
save_gg("sim_contrast_density", "MDL_ppc1", p_sim_contrast_ppc1)


# Traces look fine:
traceplot(fit_sim_ppc1)
trankplot(fit_sim_ppc1)
save_pdf("sim_traceplot", "MDLb", function() traceplot(fit_sim_ppc1))
save_pdf("sim_trankplot", "MDLb", function() trankplot(fit_sim_ppc1))

save_report("sim_summary", "MDLb", fit_sim_ppc1, pc_b_sim, model_ppc1, model_name = "The Tamed Wildcard")

pb_sim_contrast <- contrast_plot_panels(pc_b_sim, quant = c(0, 0.995), group_pal = Management_palette,
                                         true_vals = true_vals) +
  labs(x = 'Estimates for mean Hill number of order 1 and its contrast'); pb_sim_contrast
# Not incredible, but much better (depends on iteration/unstable; let's SBC!)

save_gg("sim_contrast_density", "MDLb", pb_sim_contrast)

## Simulation-based calibration (SBC) --------------------------------------------

# Does this model recover the contrast, without divergences, across many
# datasets drawn from its own priors (not just the one simulation above)?
# run_sbc() pulls "true" params via extract.prior() on fit_sim_ppc1 itself (via
# the shared draw_true()), so it always matches whatever priors model_ppc1
# actually declares.

# We reuse simulate_from_priors() from before
sbcb <- run_sbc(
  model_fit = fit_sim_ppc1,  # model_ppc1's formula + priors, both from this one fit
  means_fn = means_MDL,
  simulate_fn = simulate_from_priors,
  n_sbc = 100, iter = 10000, chains = 4,
  control = list(adapt_delta = 0.99))

sbc_out <- save_sbc_report(sbcb, "MDLb")
par(mfrow= c(1,1))
hist(sbc_out$ranks)

## Model fit ----------------------------------------------------------------

dat <- list(
  Dv = div$Hill_1,
  Mg = idx$Mg$to_index(div$Management),
  Lo = idx$Lo$to_index(div$Location)
)

fitb <- ulam(
  model_ppc1,
  data = dat,
  chains = 6, cores = 6, iter = 10000,
  control = list(adapt_delta = 0.99)
)
save_fit("fit", "MDLb", fitb)

precis(fitb, depth = 2 )

traceplot(fitb); trankplot(fitb)
save_pdf("fit_traceplot", "MDLb", function() traceplot(fitb))
save_pdf("fit_trankplot", "MDLb", function() trankplot(fitb))

## Posterior predictive check --------------------------------------------------
# Overlay the real data with data simulated from the posterior distribution.
# AFAIK, we are showing the estimate (i.e. its distribution), on which we
# overlay the real data and see if it makes sense. It doesn't have to fit
# perfectly, and the differences we see between posterior and actual data
# are a story.

### Overall ----

pb_postpred <- plot_ppc_overlay(fitb, dat, idx$Mg$to_label(dat$Mg), xlim = c(0,150)); pb_postpred
# The model seems to underestimate the center of mass for the organic group
# as well as overestimate its spread

save_gg("postpred_density", "MDLb", pb_postpred)


### By location (marginal shape) ----
postpred_long <-
  sim(fitb, dat, n = 1000) %>%
  postpred_as_long_tibble()

# n-annotated Location labels for this plot only (to_label_n(), make_index()
# in hiermod_core.R) -- postpred_long/obs_df keep plain `Lo` for joins;
# Lo_n is purely a display label carrying the per-cultivar sample size.
Loc_labels_n <- tibble(
  Lo   = idx$Lo$to_label(seq_along(idx$Lo$levels)),
  Lo_n = idx$Lo$to_label_n(seq_along(idx$Lo$levels))
)

# stratification variables for postpred_long's `obs` index
pp_joined <- tibble(
  obs = seq_len(length(dat$Mg)),
  Mg  = idx$Mg$to_label(dat$Mg),
  Lo  = idx$Lo$to_label(dat$Lo)) %>%
  left_join(postpred_long, by = "obs") %>%
  left_join(Loc_labels_n, by = "Lo")

obs_df <- tibble(
  Dv = dat$Dv,
  Mg = idx$Mg$to_label(dat$Mg),
  Lo  = idx$Lo$to_label(dat$Lo)
) %>% left_join(Loc_labels_n, by = "Lo")

# quick summaries
obs_df %>%
  group_by(Mg, Lo) %>%
  summarise(med_obs = median(Dv), MAD = mad(Dv), .groups = 'drop')

pp_joined %>%
  group_by(Mg, Lo) %>%
  summarise(med_obs = median(Dv_sim), MAD = mad(Dv_sim), .groups = 'drop')

(xlim_upper_pp <- quantile(pp_joined$Dv_sim, 0.99))
quantile(pp_joined$Dv_sim, 0.999)

p_postpred_ridges <- pp_joined %>%
  ggplot(
    aes(
      x = Dv_sim,
      y = Lo_n,
      height = after_stat(density),
      fill = Mg)) +
  ggridges::geom_density_ridges(
    stat = "density",
    alpha = 0.55,
    colour = "white",
    scale = 0.8,
    linewidth = 0.3) +

  geom_point(
    data = obs_df, inherit.aes = FALSE,
    aes(x = Dv, y = as.numeric(Lo_n) + 0.4,
        fill = Mg,
        shape = Lo),
    position = position_jitter(height = 0.1, width = 0),
    size = 2, alpha = 0.5, colour = "grey20"
  ) +

  scale_shape_manual(values = c(21:24)) +
  scale_fill_manual(values = Management_palette) +
  scale_colour_manual(values = Management_palette) +
  coord_cartesian(xlim = c(0, xlim_upper_pp)) +
  theme_light() +
  guides(colour = 'none') +
  labs(x = "Diversity", y = "Location", fill = "Management",
       caption = "Ridges = posterior predictive density; points = observed data (jittered)"); p_postpred_ridges

save_gg("postpred_ridges", "MDL", p_postpred_ridges)

## Contrast statistic (target derived quantities) ----------------------------
# The ridge/dens_overlay checks above test marginal shape (does each group's
# simulated distribution look plausible).

# For each posterior draw, simulate a full replicate dataset
# and recompute the same contrast we'd compute on real data (median gap, and
# the MAD gap that stands in for the sigma story above), then see where the
# real data's contrast falls among the replicates. Standard posterior-
# predictive test-statistic check (bayesplot::ppc_stat); yrep rows are in the
# same observation order as dat$Dv/dat$Mg, so indexing by dat$Mg inside each
# stat function is valid for every replicate row too.

p_ppc_median_contrast <- plot_ppc_contrast_stat(fitb, dat, dat$Mg, median, "median", idx$Mg$levels)
# !!!!!!! not what we'd expect?! how extreme can it get before redflagging?
p_ppc_mad_contrast <- plot_ppc_contrast_stat(fitb, dat, dat$Mg, mad, "dispersion (MAD)", idx$Mg$levels)

p_ppc_median_contrast
p_ppc_mad_contrast

save_gg("postpred_stat_median_contrast", "MDLb", p_ppc_median_contrast)
save_gg("postpred_stat_mad_contrast", "MDLb", p_ppc_mad_contrast)
