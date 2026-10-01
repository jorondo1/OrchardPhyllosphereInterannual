# MDSYC_model.R --- MODEL 5 (MDSYC): MDSYz + 3 standardised covariates (additive slopes)
# - deg_h_z: degree-hours before sampling (centred within season)
# - precip_72h_z: precipitation, 72 h before sampling
# - seq_depth_z: log read count (technical: deeper sequencing -> more taxa)
# - computed in 0.3.2_Metadata_phyloseq.R

source('src/hiermod/0_INDEX.R')

# - no Tree yet (separate branch, merged in model 7)
# - no Location: confounded with Management (A/C single-management)
# Priors: as MDSYz; new b_deg, b_precip, b_seq ~ dnorm(0,1)

model_MDSYC_16S <- alist(
  likelihood = Dv ~ dlnorm(mu, sigma[Mg]),
  main_model = mu <- loga[Mg] + gamma*(Mo-1) + yr_eff +
    b_deg*deg_h_z + b_precip*precip_72h_z + b_seq*seq_depth_z,
  gamma_def  = gamma <- s_conv + gap_shift*(Mg-1), # interactive Season effect

  prior_loga = loga[Mg]  ~ dnorm(5,2),
  prior_s    = s_conv    ~ dnorm(0,1),
  prior_gs   = gap_shift ~ dnorm(0,1),
  pr_sigma   = sigma[Mg] ~ dhalfnorm(0,1),

  # Year, sum-to-zero: yr1/yr2 free, yr3 = -(yr1+yr2)
  yr_eff_def = yr_eff <- yr1*(Yr==1) + yr2*(Yr==2) - (yr1+yr2)*(Yr==3),
  prior_yr1  = yr1 ~ dnorm(0,1),
  prior_yr2  = yr2 ~ dnorm(0,1),

  # Covariates (standardized)
  prior_deg    = b_deg    ~ dnorm(0,1),
  prior_precip = b_precip ~ dnorm(0,1),
  prior_seq    = b_seq    ~ dnorm(0,1)
)

attr(model_MDSYC_16S, "name") <- "Radagast the Grower"
model_id_MDSYC <- "MDSYC"

## ITS variant -----------------------------------------------------------------
# Intercept prior only (as MDS_ITS)
model_MDSYC_ITS <- model_MDSYC_16S
model_MDSYC_ITS$prior_loga <- quote(loga[Mg] ~ dnorm(2,2))
attr(model_MDSYC_ITS, "name") <- "Radagast the Grower"

## Backtransform wrapper ------------------------------------------------------
# 4-cell means/medians at reference covariates (z = 0 by default)

means_MDSYC <- function(post, shift = 0, deg_h_z = 0, precip_72h_z = 0, seq_depth_z = 0){
  total_var_conv <- post$sigma[,1]^2
  total_var_org  <- post$sigma[,2]^2

  s_conv    <- as.vector(post$s_conv)
  gap_shift <- as.vector(post$gap_shift)
  covariate_offset <- as.vector(post$b_deg)*deg_h_z + as.vector(post$b_precip)*precip_72h_z +
    as.vector(post$b_seq)*seq_depth_z

  mu_conv_May  <- post$loga[,1] + covariate_offset
  mu_conv_July <- mu_conv_May + s_conv
  mu_org_May   <- post$loga[,2] + covariate_offset
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

# SBC estimands: same as dq_MDSYz (covariates at z = 0 cancel)
dq_MDSYC <- SBC::derived_quantities(
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

## Variance partition / Bayesian R2 -------------------------------------------
# Sequential (type I) partition, build order: Mg x Season -> Year -> Covariates
# - each share = variance added by that group; residual = size-weighted sigma[Mg]^2
# - order-dependent; increments can be slightly negative

variance_partition_MDSYC <- function(post, dat){
  loga_obs   <- post$loga[, dat$Mg]                                # n_draws x N
  gamma      <- as.vector(post$s_conv) + outer(as.vector(post$gap_shift), dat$Mg - 1)
  gamma_term <- sweep(gamma, 2, dat$Mo - 1, "*")
  yr3        <- -(as.vector(post$yr1) + as.vector(post$yr2))
  yr_obs     <- cbind(post$yr1, post$yr2, yr3)[, dat$Yr]           # n_draws x N
  covariates <- outer(as.vector(post$b_deg), dat$deg_h_z) +
    outer(as.vector(post$b_precip), dat$precip_72h_z) +
    outer(as.vector(post$b_seq), dat$seq_depth_z)

  fixed_MgMo  <- loga_obs + gamma_term
  fixed_MgMoY <- fixed_MgMo + yr_obs
  fixed_full  <- fixed_MgMoY + covariates

  var_MgMo  <- apply(fixed_MgMo,  1, var)
  var_Y     <- apply(fixed_MgMoY, 1, var) - var_MgMo
  var_Cov   <- apply(fixed_full,  1, var) - apply(fixed_MgMoY, 1, var)
  explained <- apply(fixed_full,  1, var) # length n_draws, == var_MgMo+var_Y+var_Cov

  mg_n <- as.integer(table(factor(dat$Mg, levels = 1:2)))
  residual_var <- as.vector((post$sigma^2) %*% (mg_n / sum(mg_n)))

  total <- explained + residual_var

  bind_rows(
    tibble(statistic = "Variance partition", group = "Management x Season", value = var_MgMo / total),
    tibble(statistic = "Variance partition", group = "Year",                 value = var_Y / total),
    tibble(statistic = "Variance partition", group = "Covariates",           value = var_Cov / total),
    tibble(statistic = "Variance partition", group = "Residual",             value = residual_var / total)
  )
}

## Data-generating function ---------------------------------------------------
# As sim_div_MDSYz() + 3 independent N(0,1) covariates per row

sim_div_MDSYC <- function(N_samples, loga, s_conv, gap_shift, sigma, yr1, yr2,
                           b_deg, b_precip, b_seq, shift = NULL){
  n_unit <- N_samples %/% 2 # 2 rows/unit (May + July)
  yr_vec <- c(yr1, yr2, -(yr1 + yr2))

  units <- tibble(Un = seq_len(n_unit), Mg = rbern(n_unit) + 1, Yr = sample(3, n_unit, replace = TRUE))
  dat <- units %>% crossing(Mo = 1:2) %>% arrange(Un)

  gamma <- s_conv + gap_shift*(dat$Mg - 1)
  mu_structural <- loga[dat$Mg] + gamma*(dat$Mo - 1) + yr_vec[dat$Yr]

  dat$deg_h_z      <- rnorm(nrow(dat))
  dat$precip_72h_z <- rnorm(nrow(dat))
  dat$seq_depth_z  <- rnorm(nrow(dat))

  mu <- mu_structural + b_deg*dat$deg_h_z + b_precip*dat$precip_72h_z + b_seq*dat$seq_depth_z

  dat$Dv <- rlnorm(nrow(dat), meanlog = mu, sdlog = sigma[dat$Mg])
  if (!is.null(shift)) dat$Dv_shifted <- shift + dat$Dv
  dat
}

# Simulate one dataset from one prior draw (draw_true() output)
simulate_from_priors_MDSYC <- function(true_params, N_samples = 250, shift = NULL){
  sim_div_MDSYC(
    N_samples = N_samples,
    loga = true_params$loga,
    s_conv = true_params$s_conv,
    gap_shift = true_params$gap_shift,
    sigma = true_params$sigma,
    yr1 = true_params$yr1,
    yr2 = true_params$yr2,
    b_deg = true_params$b_deg,
    b_precip = true_params$b_precip,
    b_seq = true_params$b_seq,
    shift = shift
  )
}
