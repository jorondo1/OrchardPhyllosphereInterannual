# MDS_model.R --- MODEL 2 (MDS): Management x Season interaction
# - no random effects; validated alone before adding any

source('src/hiermod/0_INDEX.R')

# Priors
# - loga, sigma: as validated in MDv
# - new: s_conv (conventional May -> July shift), gap_shift (organic extra shift) ~ dnorm(0,1)

model_MDS_16S <- alist(
  likelihood = Dv ~ dlnorm(mu, sigma[Mg]),
  main_model = mu <- loga[Mg] + gamma*(Mo-1),
  gamma_def  = gamma <- s_conv + gap_shift*(Mg-1), # interactive Season effect

  prior_loga = loga[Mg]  ~ dnorm(5,2),
  prior_s    = s_conv    ~ dnorm(0,1),
  prior_gs   = gap_shift ~ dnorm(0,1),
  pr_sigma   = sigma[Mg] ~ dhalfnorm(0,1)
)

attr(model_MDS_16S, "name") <- "Strider the Unrooted"
model_id_MDS <- "MDS"

## ITS variant -----------------------------------------------------------------
# Only the intercept prior changes: fungal Hill numbers are much lower
# - other priors are log-scale/relative quantities, kept from 16S
model_MDS_ITS <- model_MDS_16S
model_MDS_ITS$prior_loga <- quote(loga[Mg] ~ dnorm(2,2))
attr(model_MDS_ITS, "name") <- "Strider the Unrooted"

## Backtransform wrapper ------------------------------------------------------
# Hill-scale means/medians for the 4 Mg x Mo cells

means_MDS <- function(post, shift = 0){
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

# SBC estimands: May gap, July gap, seasonal change in gap (mean-based)
dq_MDS <- SBC::derived_quantities(
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
# Balanced Mg x Mo design, one May + one July row per unit
# - simulator matches the model's own likelihood (needed for SBC)

sim_div_MDS <- function(N_samples, loga, s_conv, gap_shift, sigma, shift = NULL){
  n_unit <- N_samples %/% 2 # 2 rows/unit (May + July)

  units <- tibble(Un = seq_len(n_unit), Mg = rbern(n_unit) + 1)
  dat <- units %>% crossing(Mo = 1:2) %>% arrange(Un)

  gamma <- s_conv + gap_shift*(dat$Mg - 1)
  mu <- loga[dat$Mg] + gamma*(dat$Mo - 1)

  dat$Dv <- rlnorm(nrow(dat), meanlog = mu, sdlog = sigma[dat$Mg])
  if (!is.null(shift)) dat$Dv_shifted <- shift + dat$Dv
  dat
}

# Simulate one dataset from one prior draw (draw_true() output)
simulate_from_priors_MDS <- function(true_params, N_samples = 250, shift = NULL){
  sim_div_MDS(
    N_samples = N_samples,
    loga = true_params$loga,
    s_conv = true_params$s_conv,
    gap_shift = true_params$gap_shift,
    sigma = true_params$sigma,
    shift = shift
  )
}
