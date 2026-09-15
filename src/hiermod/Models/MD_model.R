# MD_model.R --- MODEL 1 (MD, constant variance) and MODEL 2 (MDv,
# Mg-specific variance) definitions, shared by ITS and 16S: model alist,
# means_fn, and data-generating function for each. 
## Data-generating function ---------------------------------------------------
# Hand-rolled simulator: raw mean_/cv_ -> lognormal Dv, given a group index Mg.

sim_div_M <- function(Mg, mean_, cv_){
  N      <- length(Mg)
  Dv     <- numeric(N)

  sigma2 <- log(1 +  cv_[Mg]^2)
  mu     <- log( mean_[Mg]) - sigma2 / 2
  Dv     <- rlnorm(N, meanlog = mu, sdlog = sqrt(sigma2))
  data.frame(Mg, Dv)
}

## MODEL 1 -- Mean difference by Management, constant variance ---------------

model_MD_ITS <- alist(
  likelihood   = Dv ~ dlnorm(mu, sigma),  # mu = mean of log(Dv)
  linear_model = mu <- loga[Mg],
  prior_loga   = loga[Mg] ~ dnorm(2, 2),
  prior_sigma  = sigma ~ dexp(1)
)

# Same for 16S because no calibration so far
model_MD_16S <- model_MD_ITS

model_id_MD <- "MD"

# Direct loga/sigma parameterization (loga IS the model's own log-median,
# sigma IS its own likelihood sdlog) -- unlike sim_div_M()'s mean_/cv_
# round-trip above (kept as-is for the existing effect-size sanity check),
# this matches the convention used by every later model's own
# simulate_from_priors_X()/SBC setup, for direct comparability.
sim_div_MD <- function(Mg, loga, sigma){
  data.frame(Mg, Dv = rlnorm(length(Mg), meanlog = loga[Mg], sdlog = sigma))
}

simulate_from_priors_MD <- function(true_params, N_samples = 250){
  sim_div_MD(Mg = rbern(N_samples) + 1, loga = true_params$loga, sigma = true_params$sigma)
}

# total_var -> mean[,1:2]: a single shared scalar sigma (no per-group or pooling variance)
means_MD <- function(post){
  total_var <- as.vector(post$sigma)^2
  list(mean = cbind(
    lognormal_mean(post$loga[,1], total_var),
    lognormal_mean(post$loga[,2], total_var)
  ))
}

## MODEL 2 -- Allow Management-specific variance (heteroscedasticity)
# One variance per group: sigma[Mg] instead of a shared sigma.

model_MDv_ITS <- alist(
  likelihood   = Dv ~ dlnorm(mu, sigma[Mg]),
  linear_model = mu <- loga[Mg],
  prior_loga   = loga[Mg] ~ dnorm(2,2),
  prior_sigma  = sigma[Mg] ~ dexp(1)
)

model_MDv_16S <- model_MDv_ITS

model_id_MDv <- "MDv"

# total_var -> mean[,1:2]: sigma[Mg] differs by group, so each
# group gets its own total_var.
means_MDv <- function(post){
  list(mean = cbind(
    lognormal_mean(post$loga[,1], post$sigma[,1]^2),
    lognormal_mean(post$loga[,2], post$sigma[,2]^2)
  ))
}
