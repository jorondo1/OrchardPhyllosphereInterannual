# 2.2_MDL_validation.R -- MODEL 2 (partial pooling across Location, exploratory,
# untightened priors) -- MODEL 2B (tightened priors, actually fit to real data).
# Parameter recovery, prior-predictive check, SBC, the real fit, and PPC.

source('src/hiermod/ITS/0_SETUP.R')
source('src/hiermod/ITS/2.1_MDL_model.R') # model, means_MDL(), sim_div_ML(), simulate_from_priors()
hiermod_out_dir <- "out/hiermod/ITS_2_lognormal_MDL"

# MODEL 2 -- Partial pooling across Location (non-centered) =================

# Keeps Management-specific variance from Model 1's MDv.
# Exploratory: priors here turn out to be too loose (see Prior predictive
# check below) -- MODEL 2B is the version actually fit to real data.

## Parameter recovery -----------------------------------------------------------
# Funnel check via pairs(): if divergences cluster where sigma_loc is small,
# that's the classic non-centered-parameterization funnel. Also serves as the
# compiled `ulam` object extract.prior() needs below, in the Prior predictive
# check -- extract.prior() can't run without a compiled fit of the right
# structure, so this has to come first even though conceptually you'd want to
# sanity-check the priors before trusting any fit's recovery.

dat_sim <- sim_div_ML(
  Mg = rbern(250)+1,
  # 4 locations for simulation; fewer creates a funnel (see below)
  Lo = sample(1:4, 250, replace = TRUE),
  loga = log(c(8, 12)),           # Difference of 4 in mean, raw scale
  sigma = cv_to_sigma(c(0.5, 0.8)))  # keep different variances

fit_sim <- ulam(
  model,
  data = as.list(dat_sim),
  chains = 6, cores = 6, iter = 10000 ) # divergences!

post_counts_sim <- postcounts(fit_sim, means_fn = means_MDL)

# Mean is often apart from median (iterate to see range)

pairs_plot <- pairs(fit_sim, pars = c("sigma_loc","b[1]","b[2]"))
save_pdf("fit_pairs", "MDL", function()
  pairs(fit_sim, pars = c("sigma_loc","b[1]","b[2]")))

precis(fit_sim, depth = 2 )

p_sim_contrast <- plot_contrast_density(post_counts_sim, 0.999, group_name = 'Posterior mean'); p_sim_contrast

save_report("sim_summary", "MDL", fit_sim, post_counts_sim, model, model_name = "The Wildcard")
save_gg("sim_contrast_density", "MDL", p_sim_contrast)

## Prior predictive check -------------------------------------------------------
# How extreme can Dv get under these priors?

# E[Dv] = exp(mu + sigma^2/2) means any prior tail draw of sigma[Mg]
# or sigma_loc amplifies Dv exponentially. Some large sigma draws are
# expected under dexp(1), P(sigma_loc>3) ~ 5%. If this tail looks
# implausible for Hill_1, tighten sigma[Mg]/sigma_loc, e.g. dexp(2)
# or dexp(3)).

n_prior <- 1000
prior <- extract.prior(fit_sim, n = n_prior)

# draw_true() is shared (sbc_helpers.R) -- generic slice of any extract.prior()
# output, no per-model rewrite needed. simulate_from_priors() is sourced from
# 2.1_MDL_model.R above.

# One simulated dataset per prior draw
prior_pred <- map_dfr(seq_len(n_prior), function(i){
  simulate_from_priors(draw_true(prior, i))
}, .id = "draw")

summary(prior_pred$Dv)  # raw per-observation draws
# heavy max regardless of prior
# signal of interest: 3rd quantile

p_sim_spaghetti <- prior_predictive_spaghetti(
  prior_pred, upper_q = 0.99, model = model,
  title = "Prior predictive check -- model (dexp(1))"); p_sim_spaghetti

save_gg("sim_prior_PC", "MDL", p_sim_spaghetti)

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

# both sigma separate extreme from non-extreme sharply (~1 -> ~2.5);
# we don't check loga because it's in the dnorm space
# E[Dv] = exp(mu + σ²/2) with unbounded priors will always have some values
# fall in the far range of the tail. We can slightly tighten the sigma priors.

# Try tightening sigma priors
model_ppc1 <- model
model_ppc1$prior_sigma = sigma[Mg] ~ dexp(3)
model_ppc1$prior_sigma_loc = sigma_loc ~ dexp(2)

# Minimal, throwaway fit, not a real parameter-recovery fit
fit_ppc1 <- ulam(
  model_ppc1,
  data = as.list(dat_sim),
  chains = 1, cores = 1, iter = 100, messages = FALSE)

# Extract the priors
prior_ppc1 <- extract.prior(fit_ppc1, n = n_prior)

# simulate over each of the n_prior
prior_pred_ppc1 <- map_dfr(seq_len(n_prior), function(i){
  simulate_from_priors(draw_true(prior_ppc1, i))
}, .id = "draw")

summary(prior_pred_ppc1$Dv)


# Check distribution of parameters between extreme and non-extreme draws
xlim_upper_ppc1 <- 0.99
extreme_draws_ppc1 <- prior_pred_ppc1 %>%
  group_by(draw) %>%
  summarise(max_Dv = max(Dv)) %>%
  filter(max_Dv > xlim_upper_ppc1) %>%
  pull(draw) %>% as.integer()

tibble(
  draw       = seq_len(n_prior),
  extreme    = seq_len(n_prior) %in% extreme_draws_ppc1,
  sigma_loc  = prior_ppc1$sigma_loc,
  sigma_max  = pmax(prior_ppc1$sigma[,1], prior_ppc1$sigma[,2])
) %>%
  group_by(extreme) %>%
  summarise(across(c(sigma_loc, sigma_max), median))
# much more balanced

p_sim_spaghetti_ppc1 <- prior_predictive_spaghetti(
  prior_pred_ppc1, upper_q = xlim_upper_ppc1,
  model = model_ppc1,
  title = paste("Prior predictive check: sigma tightened to",
                deparse1(model_ppc1$prior_sigma[[3]]))) +
  labs(x = 'Mean Hill diversity'); p_sim_spaghetti_ppc1

save_gg("sim_prior_PC_ppc1", "MDL", p_sim_spaghetti_ppc1)
# More reasonable!
# sd is 4 orders-of-magnitude lower,
# 3rd quartile 25% lower
# mean is nearly 3 OOM lower

## Parameter recovery: updated priors -----------------------------------------------------------

fitb_sim <- ulam(
  model_ppc1,  # <- <- <- Updated model
  data = as.list(dat_sim),
  chains = 4, cores = 4, iter = 4000,
  control = list(adapt_delta = 0.99) ) # smaller step size

precis(fitb_sim, depth = 2 ) # good r_hats

post_countsb_sim <- postcounts(fitb_sim, means_fn = means_MDL)
# Iterate, mean sometimes way over target
mean(post_countsb_sim$contrast > 100) # should be very small (~0)

# if divergences cluster where sigma_loc is small, funnel:
pairs_plot <- pairs(fitb_sim, pars = c("sigma_loc","b[1]","b[2]"))
save_pdf("fit_pairs", "MDLb", function()
  pairs(fitb_sim, pars = c("sigma_loc","b[1]","b[2]")))
# They do, not sure if ok

# Traces look fine:
traceplot(fitb_sim)
trankplot(fitb_sim)
save_pdf("sim_traceplot", "MDLb", function() traceplot(fitb_sim))
save_pdf("sim_trankplot", "MDLb", function() trankplot(fitb_sim))

save_report("sim_summary", "MDLb", fitb_sim, post_countsb_sim, model_ppc1, model_name = "The Tamed Wildcard")

pb_sim_contrast <- plot_contrast_density(
  post_countsb_sim, 0.995, group_name = 'Posterior mean') +
  labs(x = 'Estimates for mean Hill number of order 1 and its contrast'); pb_sim_contrast
# Not incredible, but much better (depends on iteration/unstable; let's SBC!)

save_gg("sim_contrast_density", "MDLb", pb_sim_contrast)

## Simulation-based calibration (SBC) --------------------------------------------

# Does this model recover the contrast, without divergences, across many
# datasets drawn from its own priors (not just the one simulation above)?
# run_sbc() pulls "true" params via extract.prior() on fitb_sim itself (via
# the shared draw_true()), so it always matches whatever priors model_ppc1
# actually declares.

# We reuse simulate_from_priors() from before
sbcb <- run_sbc(
  model_fit = fitb_sim,  # model_ppc1's formula + priors, both from this one fit
  means_fn = means_MDL,
  simulate_fn = simulate_from_priors,
  n_sbc = 100, iter = 10000, chains = 4,
  control = list(adapt_delta = 0.99))

sbc_out <- summarize_sbc(sbcb)

save_sbc_report(sbc_out, "MDLb")
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
