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

# 16S variant, patched to the 16S investigation's own starting point
# (dnorm(5,2)/dexp(2), not the never-validated ITS import above) -- this
# reflects where the calibration investigation STARTED, not what it
# concluded (MD needed no correction; MDv's dexp(2) here is known
# miscalibrated -- see 1.2_MDv_16S_calibration.R -- consumers needing the
# learned dhalfnorm(0,1) apply it as a local override, same as every
# calibration script already does).
model_MD_16S <- model_MD_ITS
model_MD_16S$prior_loga  <- quote(loga[Mg] ~ dnorm(5,2))
model_MD_16S$prior_sigma <- quote(sigma ~ dexp(2))

attr(model_MD_16S, "name") <- "Samwise the Steadfast"
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

# total_var -> mean[,1:2]: a single shared scalar sigma (no per-group or pooling variance).
# median added (lognormal_mean with total_var=0, same identity used by
# every later model) so compute_contrasts(keep=c("mean","median")) works
# uniformly across the whole model family.
means_MD <- function(post){
  total_var <- as.vector(post$sigma)^2
  list(
    mean = cbind(
      lognormal_mean(post$loga[,1], total_var),
      lognormal_mean(post$loga[,2], total_var)
    ),
    median = cbind(
      lognormal_mean(post$loga[,1], 0),
      lognormal_mean(post$loga[,2], 0)
    )
  )
}

# SBC estimands -- median/mean contrast, matching means_MD()'s own
# total_var convention (single shared sigma). Property of the model, not
# any one calibration script -- every MD SBC section reuses this.
dq_MD <- SBC::derived_quantities(
  median_contrast = exp(loga[2]) - exp(loga[1]),
  mean_contrast   = exp(loga[2] + sigma^2/2) - exp(loga[1] + sigma^2/2)
)

## MODEL 2 -- Allow Management-specific variance (heteroscedasticity)
# One variance per group: sigma[Mg] instead of a shared sigma.

model_MDv_ITS <- alist(
  likelihood   = Dv ~ dlnorm(mu, sigma[Mg]),
  linear_model = mu <- loga[Mg],
  prior_loga   = loga[Mg] ~ dnorm(2,2),
  prior_sigma  = sigma[Mg] ~ dexp(1)
)

model_MDv_16S <- model_MDv_ITS
model_MDv_16S$prior_loga  <- quote(loga[Mg] ~ dnorm(5,2))
model_MDv_16S$prior_sigma <- quote(sigma[Mg] ~ dexp(2))

attr(model_MDv_16S, "name") <- "Gollum the Two-Faced"
model_id_MDv <- "MDv"

# Direct loga/sigma[Mg] parameterization, same rationale as sim_div_MD() above.
sim_div_MDv <- function(Mg, loga, sigma){
  data.frame(Mg, Dv = rlnorm(length(Mg), meanlog = loga[Mg], sdlog = sigma[Mg]))
}

simulate_from_priors_MDv <- function(true_params, N_samples = 250){
  sim_div_MDv(Mg = rbern(N_samples) + 1, loga = true_params$loga, sigma = true_params$sigma)
}

# total_var -> mean[,1:2]: sigma[Mg] differs by group, so each
# group gets its own total_var. median added, same reason as means_MD().
means_MDv <- function(post){
  list(
    mean = cbind(
      lognormal_mean(post$loga[,1], post$sigma[,1]^2),
      lognormal_mean(post$loga[,2], post$sigma[,2]^2)
    ),
    median = cbind(
      lognormal_mean(post$loga[,1], 0),
      lognormal_mean(post$loga[,2], 0)
    )
  )
}

# SBC estimands -- median/mean contrast, matching means_MDv()'s own
# per-group total_var convention.
dq_MDv <- SBC::derived_quantities(
  median_contrast = exp(loga[2]) - exp(loga[1]),
  mean_contrast   = exp(loga[2] + sigma[2]^2/2) - exp(loga[1] + sigma[1]^2/2)
)
