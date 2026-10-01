# MD_model.R --- MODEL 1 (MD: Management, shared sigma) and MDv (Management-specific sigma)
# - 16S and ITS variants, means, SBC estimands, simulators

source('src/hiermod/0_INDEX.R')

## Data-generating function ---------------------------------------------------

# Raw-scale mean/CV per group -> lognormal Dv
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
attr(model_MD_ITS, "name") <- "Samwise the Steadfast"

# 16S variant: starting priors of the 16S calibration
# - MDv's dexp(2) turned out miscalibrated -> dhalfnorm(0,1), applied as a
#   local override in calibration scripts (1.2_MDv_16S_calibration.R)
model_MD_16S <- model_MD_ITS
model_MD_16S$prior_loga  <- quote(loga[Mg] ~ dnorm(5,2))
model_MD_16S$prior_sigma <- quote(sigma ~ dexp(2))

attr(model_MD_16S, "name") <- "Samwise the Steadfast"
model_id_MD <- "MD"

# Simulator on the model's own scale (loga = log-median, sigma = sdlog)
# - same convention as every later model's simulator
sim_div_MD <- function(Mg, loga, sigma){
  data.frame(Mg, Dv = rlnorm(length(Mg), meanlog = loga[Mg], sdlog = sigma))
}

# Simulate one dataset from one prior draw (draw_true() output)
simulate_from_priors_MD <- function(true_params, N_samples = 250){
  sim_div_MD(Mg = rbern(N_samples) + 1, loga = true_params$loga, sigma = true_params$sigma)
}

# Hill-scale group means (shared sigma) and medians (total_var = 0)
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

# SBC estimands: median and mean contrast (Organic - Conventional)
dq_MD <- SBC::derived_quantities(
  median_contrast = exp(loga[2]) - exp(loga[1]),
  mean_contrast   = exp(loga[2] + sigma^2/2) - exp(loga[1] + sigma^2/2)
)

## MODEL MDv -- Management-specific sigma[Mg] (heteroscedasticity)

model_MDv_ITS <- alist(
  likelihood   = Dv ~ dlnorm(mu, sigma[Mg]),
  linear_model = mu <- loga[Mg],
  prior_loga   = loga[Mg] ~ dnorm(2,2),
  prior_sigma  = sigma[Mg] ~ dexp(1)
)
attr(model_MDv_ITS, "name") <- "Gollum the Two-Faced"

model_MDv_16S <- model_MDv_ITS
model_MDv_16S$prior_loga  <- quote(loga[Mg] ~ dnorm(5,2))
model_MDv_16S$prior_sigma <- quote(sigma[Mg] ~ dexp(2))

attr(model_MDv_16S, "name") <- "Gollum the Two-Faced"
model_id_MDv <- "MDv"

# Simulator, model scale (as sim_div_MD)
sim_div_MDv <- function(Mg, loga, sigma){
  data.frame(Mg, Dv = rlnorm(length(Mg), meanlog = loga[Mg], sdlog = sigma[Mg]))
}

# Simulate one dataset from one prior draw (draw_true() output)
simulate_from_priors_MDv <- function(true_params, N_samples = 250){
  sim_div_MDv(Mg = rbern(N_samples) + 1, loga = true_params$loga, sigma = true_params$sigma)
}

# Hill-scale means (per-group sigma) and medians
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

# SBC estimands: median and mean contrast
dq_MDv <- SBC::derived_quantities(
  median_contrast = exp(loga[2]) - exp(loga[1]),
  mean_contrast   = exp(loga[2] + sigma[2]^2/2) - exp(loga[1] + sigma[1]^2/2)
)
