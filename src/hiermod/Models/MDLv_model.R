# MODEL 3 (MDLv), ITS only: Organic-Conventional gap now varies by Location
# (g[Lo]*sigma_g), on top of Model 2's Location-varying intercept.

model_MDLv_ITS <- alist(
  likelihood   = Dv ~ dlnorm(mu, sigma[Mg]),
  linear_model = mu <- loga[Mg] + b[Lo]*sigma_loc + g[Lo]*sigma_g*(Mg-1),

  # Priors:
  prior_loga   = loga[Mg] ~ dnorm(2,2), # Estimand
  prior_b      = b[Lo]    ~ dnorm(0,1), # Location-specific parameters (centered)
  prior_g      = g[Lo]    ~ dnorm(0,1), # Location-specific gap between managements (centered)

  # Hyperpriors:
  pr_sigma     = sigma[Mg] ~ dexp(3),   # Management-specific noise
  pr_sigma_loc = sigma_loc ~ dexp(2),   # Global location spread
  pr_sigma_g   = sigma_g   ~ dexp(2)    # Spread of Organic-Conventional gap by Location
)

mdlv_labels <- c(
  median  = "Median diversity (Hill number scale)",
  mean    = "Mean diversity (Hill number scale)",
  sigma   = "Residual SD (log scale)",
  sigma_g = "Location-specific gap SD (log scale)",
  sigma_loc = "Global location SD (log scale)")
# sigma_g isn't Mg-indexed -- compute_contrasts() gives it a single
# Contrast-only row (no group1/group2 pair to difference). Per-Location
# breakdown is in g_by_loc (3.4_MDLv_analysis.R). sigma_loc isn't part of
# this panel set at all -- see precis(fitb, depth = 2) for it directly.

# Mean backtransformation: Organic carries an extra variance term (sigma_g)
# that Conventional doesn't, since g[Lo]*sigma_g only ever multiplies
# (Mg-1). Median needs no such adjustment (lognormal_mean() with
# total_var=0 is just exp(loga)) -- see hiermod_core.R. Reads
# extract.samples()'s own post$param[,i] shape; post_full() merges this
# function's mean/median output back into that same list, by name, for
# compute_contrasts() to pick up.
means_MDLv <- function(post){
  total_var1 <- post$sigma[,1]^2 + as.vector(post$sigma_loc)^2
  total_var2 <- post$sigma[,2]^2 + as.vector(post$sigma_loc)^2 + as.vector(post$sigma_g)^2
  list(
    mean = cbind(
      lognormal_mean(post$loga[,1], total_var1),
      lognormal_mean(post$loga[,2], total_var2)
    ),
    median = cbind(
      lognormal_mean(post$loga[,1], 0),
      lognormal_mean(post$loga[,2], 0)
    ))
}

## Data-generating function ---------------------------------------------------
# Same skeleton as sim_div_MDL() (MDLb), plus one more term: loc_gap[Lo], a
# per-Location deviation added *only* to the Organic linear predictor (via
# (Mg-1), same as the model above).

# Takes loga/sigma directly (the model's own units), not a raw mean_/cv_ to
# backconvert -- a mean_/cv_ argument here would just get converted straight
# back to mu/sigma2 below, an identity round-trip (see sim_div_MDLS() in
# MDLS_model.R). cv_to_sigma() (hiermod_core.R) is there for the
# one place that's actually useful: picking a human-readable CV by hand at a
# call site.
sim_div_MDLv <- function(Mg, Lo, loga, sigma, sigma_loc = 0.5, sigma_g = 0.3){
  N <- length(Mg)
  n_loc <- length(unique(Lo))

  loc_offset <- rnorm(n_loc, 0, sigma_loc)  # shared Location baseline
  loc_gap    <- rnorm(n_loc, 0, sigma_g)    # this Location's deviation in the Organic-Conventional gap

  mu <- loga[Mg] + loc_offset[Lo] + loc_gap[Lo]*(Mg-1)

  Dv <- rlnorm(N, meanlog = mu, sdlog = sigma[Mg])
  data.frame(Mg, Lo, Dv)
}

# Prior simulator (SBC/prior-predictive glue).
simulate_from_priors_MDLv <- function(true_params, N_per_group = 125, n_loc = 4){
  Mg <- rep(1:2, each = N_per_group)
  Lo <- sample(1:n_loc, length(Mg), replace = TRUE)
  sim_div_MDLv(Mg, Lo, loga = true_params$loga, sigma = true_params$sigma,
               sigma_loc = true_params$sigma_loc, sigma_g = true_params$sigma_g)
}
