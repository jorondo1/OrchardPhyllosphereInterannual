# MODEL 2 (MDL): adds partial pooling across Location (non-centered
# b[Lo]*sigma_loc) on top of Model 2/MDv's Management-specific variance.

## Data-generating function ---------------------------------------------------

sim_div_MDL <- function(Mg, Lo, loga, sigma, sigma_loc = 0.5){
  N <- length(Mg)
  n_loc <- length(unique(Lo))

  loc_offset <- rnorm(n_loc, 0, sigma_loc)   # location deviations, log scale
  mu <- loga[Mg] + loc_offset[Lo]

  Dv <- rlnorm(N, meanlog = mu, sdlog = sigma[Mg])
  data.frame(Mg, Lo, Dv)
}

simulate_from_priors_MDL <- function(true_params, N_per_group = 125, n_loc = 4){
  Mg <- rep(1:2, each = N_per_group)
  Lo <- sample(1:n_loc, length(Mg), replace = TRUE)
  sim_div_MDL(Mg, Lo, loga = true_params$loga, sigma = true_params$sigma,
              sigma_loc = true_params$sigma_loc)
}

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

## Model spec (first-tested, per marker) -----------------------------------

model_MDL_ITS <- alist(
  likelihood      = Dv ~ dlnorm(mu, sigma[Mg]),
  linear_model    = mu <- loga[Mg] + b[Lo]*sigma_loc,
  prior_loga      = loga[Mg] ~ dnorm(2,2),
  prior_b         = b[Lo] ~ dnorm(0,1),
  prior_sigma     = sigma[Mg] ~ dexp(1),
  prior_sigma_loc = sigma_loc ~ dexp(1)
)

# Same as ITS: no calibration for either marker yet.
model_MDL_16S <- model_MDL_ITS

# Canonical id from the real fit onward (both markers already use this) --
# calibration's own naive pre-tightening stage keeps its own local "MDL" label.
model_id <- "MDLb"
