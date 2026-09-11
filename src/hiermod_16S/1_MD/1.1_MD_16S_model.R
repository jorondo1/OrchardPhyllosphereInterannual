# 1.1_MD_model.R --- MODEL 1 (constant variance) and MODEL 2 (Mg-specific
# variance) definitions: model alist, means_fn, and data-generating function
# for each.

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

model_MD <- alist(
  likelihood   = Dv ~ dlnorm(mu, sigma),  # mu = mean of log(Dv)
  linear_model = mu <- loga[Mg],
  prior_loga   = loga[Mg] ~ dnorm(2, 2),
  prior_sigma  = sigma ~ dexp(1)
)

# total_var -> mean[,1:2]: a single shared scalar sigma (no
# per-group or pooling variance) -- same total_var used for both groups.
means_MD <- function(post){
  total_var <- as.vector(post$sigma)^2
  list(mean = cbind(
    lognormal_mean(post$loga[,1], total_var),
    lognormal_mean(post$loga[,2], total_var)
  ))
}

## MODEL 2 -- Allow Management-specific variance (heteroscedasticity) --------
# One variance per group: sigma[Mg] instead of a shared sigma.

model_MDv <- alist(
  likelihood   = Dv ~ dlnorm(mu, sigma[Mg]),
  linear_model = mu <- loga[Mg],
  prior_loga   = loga[Mg] ~ dnorm(2,2),
  prior_sigma  = sigma[Mg] ~ dexp(1)
)

# total_var -> mean[,1:2]: sigma[Mg] differs by group, so each
# group gets its own total_var. When this model gets an SBC section, reuse
# this here too: contrast_MDv <- function(post, tp) contrast_from_means(post, tp, means_MDv)
means_MDv <- function(post){
  list(mean = cbind(
    lognormal_mean(post$loga[,1], post$sigma[,1]^2),
    lognormal_mean(post$loga[,2], post$sigma[,2]^2)
  ))
}
