# MDSY_model.R --- MODEL 4 (MDSY): MDS + Year as a fixed effect (no pooling)
source('src/hiermod/0_INDEX.R')

# Why fixed: 3 years only, too few to estimate a sigma_yr
# yr[Yr] ~ dnorm(0,1), unconstrained
# - not identified from loga[Mg]: miscalibrated -> fixed by sum-to-zero in MDSYz

model_MDSY_16S <- alist(
  likelihood = Dv ~ dlnorm(mu, sigma[Mg]),
  main_model = mu <- loga[Mg] + gamma*(Mo-1) + yr[Yr],
  gamma_def  = gamma <- s_conv + gap_shift*(Mg-1), # interactive Season effect

  prior_loga = loga[Mg]  ~ dnorm(5,2),
  prior_s    = s_conv    ~ dnorm(0,1),
  prior_gs   = gap_shift ~ dnorm(0,1),
  pr_sigma   = sigma[Mg] ~ dhalfnorm(0,1),

  prior_yr   = yr[Yr]    ~ dnorm(0,1)
)

attr(model_MDSY_16S, "name") <- "Elrond the Ageless"
model_id_MDSY <- "MDSY"

## Backtransform wrapper ------------------------------------------------------
# Same as means_MDS(): Year is location-only, reported at the average year

means_MDSY <- function(post, shift = 0){
  total_var_conv <- post$sigma[,1]^2
  total_var_org  <- post$sigma[,2]^2

  s_conv    <- as.vector(post$s_conv)
  gap_shift <- as.vector(post$gap_shift)

  mu_conv_May  <- post$loga[,1]
  mu_conv_July <- mu_conv_May + s_conv
  mu_org_May   <- post$loga[,2]
  mu_org_July  <- mu_org_May + s_conv + gap_shift

  list(
    mean = cbind(
      lognormal_mean(mu_conv_May,  total_var_conv, shift = shift),
      lognormal_mean(mu_conv_July, total_var_conv, shift = shift),
      lognormal_mean(mu_org_May,   total_var_org,  shift = shift),
      lognormal_mean(mu_org_July,  total_var_org,  shift = shift)
    ),
    median = cbind(
      lognormal_mean(mu_conv_May,  0, shift = shift),
      lognormal_mean(mu_conv_July, 0, shift = shift),
      lognormal_mean(mu_org_May,   0, shift = shift),
      lognormal_mean(mu_org_July,  0, shift = shift)
    )
  )
}

# SBC estimands: same as dq_MDS
dq_MDSY <- SBC::derived_quantities(
  may_gap =
    exp(loga[2] + sigma[2]^2 / 2) -
    exp(loga[1] + sigma[1]^2 / 2),
  july_gap =
    exp(loga[2] + s_conv + gap_shift + sigma[2]^2 / 2) -
    exp(loga[1] + s_conv +             sigma[1]^2 / 2),
  seasonal_change =
    (exp(loga[2] + s_conv + gap_shift + sigma[2]^2 / 2) -
       exp(loga[1] + s_conv +             sigma[1]^2 / 2)) -
    (exp(loga[2] + sigma[2]^2 / 2) -
       exp(loga[1] + sigma[1]^2 / 2))
)

## Data-generating function ---------------------------------------------------
# As sim_div_MDS() + a random year per unit

sim_div_MDSY <- function(N_samples, loga, s_conv, gap_shift, sigma, yr, shift = NULL){
  n_unit <- N_samples %/% 2 # 2 rows/unit (May + July)
  n_yr   <- length(yr)

  units <- tibble(Un = seq_len(n_unit), Mg = rbern(n_unit) + 1, Yr = sample(n_yr, n_unit, replace = TRUE))
  dat <- units %>% crossing(Mo = 1:2) %>% arrange(Un)

  gamma <- s_conv + gap_shift*(dat$Mg - 1)
  mu <- loga[dat$Mg] + gamma*(dat$Mo - 1) + yr[dat$Yr]

  dat$Dv <- rlnorm(nrow(dat), meanlog = mu, sdlog = sigma[dat$Mg])
  if (!is.null(shift)) dat$Dv_shifted <- shift + dat$Dv
  dat
}

# Simulate one dataset from one prior draw (draw_true() output)
simulate_from_priors_MDSY <- function(true_params, N_samples = 250, shift = NULL){
  sim_div_MDSY(
    N_samples = N_samples,
    loga = true_params$loga,
    s_conv = true_params$s_conv,
    gap_shift = true_params$gap_shift,
    sigma = true_params$sigma,
    yr = true_params$yr,
    shift = shift
  )
}
