# MDS_model.R --- MODEL 2 (MDS, Management x Season interaction), 16S: the
# from-scratch rebuild's next step after MD/MDv. Adapted directly from the
# archived MDS2's parameterization (src/hiermod/16S/archive/3_MDS2/,
# src/hiermod/Models/MDS2_model.R) with every Tree/Year term removed --
# just the Mg x Mo interaction estimand, no random effects at all, so it
# can be freshly SBC-validated on its own before any random effect is
# reintroduced.
#

source('src/hiermod/0_INDEX.R')

# Priors here hardcode what Model 1 (MD/MDv) actually learned, not a fresh
# starting guess: loga[Mg] ~ dnorm(5,2) and sigma[Mg] ~ dhalfnorm(0,1) are
# MDv's own validated answer (see 1.2_MDv_16S_calibration.R), carried
# forward as this model's starting point. s_conv/gap_shift are new
# parameters this family hasn't calibrated before -- dnorm(0,1) is MDS2's
# own (untested by our own SBC) choice, kept as the starting hypothesis to
# validate in 2.2_MDS_16S_calibration.R.

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
# Built by referencing the 16S object directly and overriding only what's
# genuinely scale-dependent: loga[Mg] ~ dnorm(2,2) is ITS's own already-
# established Model 1 starting point (MD_ITS/MDv_ITS, MD_model.R), since
# ITS's own Hill_1 values sit on a much smaller scale than 16S's (fungal
# vs bacterial diversity). s_conv/gap_shift/sigma[Mg] carry over unchanged
# from model_MDS_16S -- these are additive log-scale/CV-like quantities,
# not tied to the raw Hill_1 baseline, so 16S's own validated dnorm(0,1)/
# dhalfnorm(0,1) choices are a reasonable starting hypothesis here too.
# A guess to start from, not a conclusion -- expect this to get refined as
# real ITS SBC results come in, same as every 16S model in this family was.
model_MDS_ITS <- model_MDS_16S
model_MDS_ITS$prior_loga <- quote(loga[Mg] ~ dnorm(2,2))
attr(model_MDS_ITS, "name") <- "Strider the Unrooted"

## Backtransform wrapper ------------------------------------------------------
# Same 4-cell shape as means_MDS2(), minus the sigma_tr term (no Tree here).

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

# SBC estimands -- may_gap/july_gap/seasonal_change, matching means_MDS()'s
# own total_var convention (sigma[Mg]^2 only, no Tree term).
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
# Balanced Mg x Mo design, one row per (unit, Mo) pair -- no Tree/Year
# offsets, so every row is drawn independently given (Mg, Mo) alone,
# matching this model's own likelihood exactly (a valid SBC self-consistency
# check needs the simulator to match the model, not the real design).

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
