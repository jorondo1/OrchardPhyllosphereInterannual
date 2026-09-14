# MODEL 2 (MDL), 16S: same structure as ITS's Model 2 (partial pooling
# across Location). 16S's raw diversity scale is much higher than ITS's, so
# this checks how much the shared naive priors need retuning here.

hiermod_marker <- "16S"
source('src/hiermod/0_SETUP.R')
source('src/hiermod/Models/MDL_model.R') # model_MDL_16S, means_MDL(), sim_div_MDL(), simulate_from_priors_MDL()
model <- model_MDL_16S

hiermod_out_dir <- "out/hiermod/16S_2_lognormal_MDL"

# MODEL 2 -- Partial pooling across Location (non-centered) =================
#same rationale as ITS model 2

## Parameter recovery -----------------------------------------------------------
# Dummy values
conv <- 180
org <- 120
true_sigma <- cv_to_sigma(c(0.5, 0.8))# keep different variances
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
  chains = 6, cores = 6, iter = 5000 )
precis(fit_sim, depth = 2 )

### Fixed effect recovery ---------

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
# output, no per-model rewrite needed. simulate_from_priors_MDL() is sourced
# from hiermod/Models/MDL_model.R above.

# One simulated dataset per prior draw
prior_pred <- map_dfr(seq_len(n_prior), function(i){
  simulate_from_priors_MDL(draw_true(prior, i))
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
diagnose_extreme_tail(
  prior_pred, value_col = "Dv",
  candidates = list(
    sigma_loc = prior$sigma_loc,
    sigma_max = pmax(prior$sigma[,1], prior$sigma[,2])))

# sigma max is quite extreme; sigma loc a little too

### Update loga prior --------------------------------------

# Sit halfway between log means , make sure means sit within 1sd of the mean
model_ppc1 <- model
model_ppc1$prior_loga <- quote(loga[Mg] ~ dnorm(6,2)) #slightly more skeptical than ITS?
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
  simulate_from_priors_MDL(draw_true(prior_ppc1, i))
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

### Fixed effect recovery ---------

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
  pc_sim_ppc1, quant = c(0.001, 0.995), group_pal = Management_palette,
  true_vals = true_vals); p_sim_contrast_ppc1

save_report(
  "sim_summary", model_id, fit_sim_ppc1, pc_sim_ppc1, model_ppc1, model_name = "The Tamed Wildcard")
save_gg("sim_contrast_density", model_id, p_sim_contrast_ppc1)


# Traces look fine:
save_pdf("sim_traceplot", model_id, function() traceplot(fit_sim_ppc1))
save_pdf("sim_trankplot", model_id, function() trankplot(fit_sim_ppc1))

## Simulation-based calibration (SBC) --------------------------------------------

# Does this model recover the contrast, without divergences, across many
# datasets drawn from its own priors (not just the one simulation above)?
# run_sbc() pulls "true" params via extract.prior() on fit_sim_ppc1 itself (via
# the shared draw_true()), so it always matches whatever priors model_ppc1
# actually declares.


ncores <- 24
nchains <- 2
n_sbc = 100
n_iter = 10000

# We reuse simulate_from_priors_MDL() from before
sbc <- run_sbc(
  model_fit = fit_sim_ppc1,  # model_ppc1's formula + priors, both from this one fit
  means_fn = means_MDL,
  simulate_fn = simulate_from_priors_MDL,
  n_sbc = 100, iter = 5000, chains = 4,
  control = list(adapt_delta = 0.99))

sbc_out <- save_sbc_report(sbc, model_id)
par(mfrow= c(1,1))
hist(sbc_out$ranks)

saveRDS(model_ppc1, file.path(hiermod_out_dir, "model_ppc1.rds"))
