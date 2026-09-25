# MDSY_model.R --- MODEL 4 (MDSY), 16S: MDS plus Year as a FIXED effect
# (yr[Yr] ~ dnorm(0,1), no sigma_yr, no pooling) -- deliberately not random.
source('src/hiermod/0_INDEX.R')

# Only 3 years exist. A hierarchical/pooled treatment would need sigma_yr
# estimated from just 3 group means -- barely more informative than fitting
# 3 independent offsets directly, while adding another entangled scale
# parameter of exactly the kind that just caused MDST's funnel. Fixed costs
# nothing here and sidesteps that risk entirely. (This also matches the
# archived MDS2's own original treatment of Year, before it was later
# "promoted" to pooled much further down that lineage.)
#
# loga[Mg]/s_conv/gap_shift/sigma[Mg] are MDS's own validated answer,
# hardcoded as this model's starting point. yr[Yr] ~ dnorm(0,1) is the one
# new assumption -- unlike Tree's tr[Tr] (cardinality scales with N_samples,
# so excluded from SBC's `keep`), Year has a small, fixed cardinality (3)
# that doesn't depend on sample size, so it can be fully tracked and
# calibrated via SBC just like any other fixed-effect parameter.
#
# Note: yr[Yr] isn't sum-to-zero constrained, so it's not perfectly
# identified separately from loga[Mg]'s own overall level -- same
# non-identifiability MDS2 already lived with. Not a problem for the
# combined mu each row actually needs; means_MDSY() below reports the
# Mg x Mo estimand at the (implicit) average-across-years level, matching
# MDS2's own convention of leaving yr[]/tr[] out of the reported means.

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
# Identical to means_MDS() -- yr[Yr] is a location-only fixed effect, no
# variance term to add, and the reported cells are implicitly averaged
# across years via loga[Mg]'s own role (same convention as MDS2's means_fn,
# which also left yr[]/tr[] out of the reported means).

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

# SBC estimands -- identical formulas to dq_MDS (Year doesn't enter the
# Mg x Mo estimand or its variance).
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
# Same balanced Mg x Mo design as sim_div_MDS() (no Tree here, so no
# tree-pairing bookkeeping needed), plus an independent Yr draw per unit and
# its additive yr[Yr] offset.

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
