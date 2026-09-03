# 2.1_MDL_model.R -- MODEL 2 (partial pooling across Location) base
# definition: model alist, means_fn, data-generating function, and the SBC
# prior-simulator. See 2.2_MDL_validation.R for the actual runs/fits/plots,
# including the in-script discovery of MODEL 2B's tightened-prior variant
# (model_ppc1 <- model; model_ppc1$prior_sigma <- ...) -- that derivation is
# part of the analysis narrative (found via that script's own prior
# predictive check), so it stays there rather than here.

## Model specification ---------------------------------------------------------
# + b[Lo]*sigma_loc

model <- alist(
  likelihood   = Dv ~ dlnorm(mu, sigma[Mg]),
  linear_model = mu <- loga[Mg] + b[Lo]*sigma_loc,
  prior_loga   = loga[Mg] ~ dnorm(2,2),
  prior_b      = b[Lo] ~ dnorm(0,1),
  prior_sigma  = sigma[Mg] ~ dexp(1),
  prior_sigma_loc = sigma_loc ~ dexp(1)
)

# Mean backtransformation wrapper
# total_var -> mean[,1:2]: sigma[Mg] (observation-level) + sigma_loc^2
# (Location mean-offset, marginalized over the Location population). Median
# needs no such adjustment (lognormal_mean() with total_var=0 is just
# exp(loga)) -- see hiermod_core.R.

means_MDL <- function(post){
  total_var1 <- post$sigma[,1]^2 + as.vector(post$sigma_loc)^2
  total_var2 <- post$sigma[,2]^2 + as.vector(post$sigma_loc)^2
  list(
    mean = cbind(
      lognormal_mean(post$loga[,1], total_var1),
      lognormal_mean(post$loga[,2], total_var2)),
    median = cbind(
      lognormal_mean(post$loga[,1], 0),
      lognormal_mean(post$loga[,2], 0))
  )
}

## Data-generating function ---------------------------------------------------
# Adds location-level log-scale offsets on top of sim_div_M()'s structure.
# Takes loga/sigma directly (the model's own units), not a raw mean_/cv_ to
# backconvert -- see sim_div_MDLS() in 4.1_MDLS_model.R for why: a
# mean_/cv_ argument here would just get converted straight back to mu/sigma2
# below, an identity round-trip. cv_to_sigma() (hiermod_core.R) is there for
# the one place that's actually useful: picking a human-readable CV by hand
# at a call site below.

sim_div_ML <- function(Mg, Lo, loga, sigma, sigma_loc = 0.5){
  N <- length(Mg)
  n_loc <- length(unique(Lo))

  loc_offset <- rnorm(n_loc, 0, sigma_loc)   # location deviations, log scale
  mu <- loga[Mg] + loc_offset[Lo]

  Dv <- rlnorm(N, meanlog = mu, sdlog = sigma[Mg])
  data.frame(Mg, Lo, Dv)
}

# Prior simulator (SBC/prior-predictive glue). sim_div_ML() already takes
# loga/sigma directly, so no mean_/cv_ conversion needed here.
simulate_from_priors <- function(true_params, N_per_group = 125, n_loc = 4){
  Mg <- rep(1:2, each = N_per_group)
  Lo <- sample(1:n_loc, length(Mg), replace = TRUE)
  sim_div_ML(Mg, Lo, loga = true_params$loga, sigma = true_params$sigma,
             sigma_loc = true_params$sigma_loc)
}
