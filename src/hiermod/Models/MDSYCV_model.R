# MDSYCV_model.R --- MODEL 6 (MDSYCV, "Bombadil the Eldest"), 16S: MDSYC
# plus Cultivar as a FIXED, sum-to-zero effect (5 levels: Cortland, Liberty,
# Paulared, Honeycrisp, Spartan).

source('src/hiermod/0_INDEX.R')

# Motivation: MDSYC found sigma[Mg] meaningfully higher for Organic than
# Conventional (89% contrast [0.16, 0.48], excludes 0) even after season,
# year, and weather/sequencing covariates. Management and Cultivar aren't
# perfectly balanced in the real data (Liberty: 35 Conventional vs 18
# Organic; Honeycrisp: 20 vs 24) -- not a hard confound like Location (every
# cultivar has substantial presence in both groups), but real correlation
# that could plausibly explain some of that residual asymmetry if cultivars
# themselves differ in their own diversity variability.
#
# Sum-to-zero from the start this time (4 free scalars + 1 derived),
# applying what MDSY's own severe loga leak taught us, rather than the old
# lineage's plain unconstrained cv[Cv] ~ dnorm(0,1) (MDLSY_model.R) -- no
# reason to risk rediscovering that bug.
#
# Cv index order (idx$Cv$levels): 1=Cortland, 2=Liberty, 3=Paulared,
# 4=Honeycrisp, 5=Spartan. Liberty has the most combined observations across
# both Management groups (53, vs 44-51 for the other four), so it's the
# derived slot -- same "put the data-richest level where the construction's
# extra prior variance matters least" logic used for Year (see
# MDSYz_model.R). Free parameters are named cv_1/cv_3/cv_4/cv_5 (matching
# their own Cv index directly, cv_2/Liberty deliberately skipped) so the
# mapping stays unambiguous everywhere this model is used.
#
# loga[Mg]/s_conv/gap_shift/sigma[Mg]/yr1/yr2/b_deg/b_precip/b_seq are
# MDSYC's own validated answer, hardcoded as this model's starting point.
# cv_1/cv_3/cv_4/cv_5 ~ dnorm(0,1) is the one new assumption to validate.

model_MDSYCV_16S <- alist(
  likelihood = Dv ~ dlnorm(mu, sigma[Mg]),
  main_model = mu <- loga[Mg] + gamma*(Mo-1) + yr_eff + cv_eff +
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
  prior_seq    = b_seq    ~ dnorm(0,1),

  # Cultivar, sum-to-zero: cv_1/cv_3/cv_4/cv_5 free, cv_2 (Liberty) derived
  cv_eff_def = cv_eff <- cv_1*(Cv==1) + cv_3*(Cv==3) + cv_4*(Cv==4) + cv_5*(Cv==5) - (cv_1+cv_3+cv_4+cv_5)*(Cv==2),
  prior_cv1  = cv_1 ~ dnorm(0,1),
  prior_cv3  = cv_3 ~ dnorm(0,1),
  prior_cv4  = cv_4 ~ dnorm(0,1),
  prior_cv5  = cv_5 ~ dnorm(0,1)
)

attr(model_MDSYCV_16S, "name") <- "Bombadil the Eldest"
model_id_MDSYCV <- "MDSYCV"

## ITS variant -----------------------------------------------------------------
# Same rationale as MDS_ITS. cv_1/cv_3/cv_4/cv_5 ~ dnorm(0,1) carry over
# unchanged -- same Cultivar factor (idx$Cv, shared across Kingdoms), same
# additive log-scale-offset reasoning as Year.
model_MDSYCV_ITS <- model_MDSYCV_16S
model_MDSYCV_ITS$prior_loga <- quote(loga[Mg] ~ dnorm(2,2))
attr(model_MDSYCV_ITS, "name") <- "Bombadil the Eldest"

## means_MDSYCV()/dq_MDSYCV ----------------------------------------------------
# Same formulas as means_MDSYC()/dq_MDSYC -- Cultivar, like Year, doesn't
# enter the reported Mg x Mo estimand or its variance (assigned
# independently of Mg x Mo in the simulator, reported at the default/average
# level on real data).

means_MDSYCV <- function(post, shift = 0, deg_h_z = 0, precip_72h_z = 0, seq_depth_z = 0){
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

dq_MDSYCV <- SBC::derived_quantities(
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

## variance_partition_MDSYCV() --------------------------------------------------
# Same as variance_partition_MDSYC() plus Cultivar as its own group,
# sequentially last (Model 6's own addition, after Covariates). cv_derived
# (Liberty) computed the same way the model itself does; cv_mat's columns
# are built directly in Cv's own 1..5 index order so dat$Cv can index it
# with no remapping. See variance_partition_MDSYC()'s own comment for the
# sequential-decomposition rationale/telescoping property.

variance_partition_MDSYCV <- function(post, dat){
  loga_obs   <- post$loga[, dat$Mg]
  gamma      <- as.vector(post$s_conv) + outer(as.vector(post$gap_shift), dat$Mg - 1)
  gamma_term <- sweep(gamma, 2, dat$Mo - 1, "*")
  yr3        <- -(as.vector(post$yr1) + as.vector(post$yr2))
  yr_obs     <- cbind(post$yr1, post$yr2, yr3)[, dat$Yr]

  cv_derived <- -(as.vector(post$cv_1) + as.vector(post$cv_3) + as.vector(post$cv_4) + as.vector(post$cv_5))
  cv_mat     <- cbind(post$cv_1, cv_derived, post$cv_3, post$cv_4, post$cv_5) # Cv index order: 1..5
  cv_obs     <- cv_mat[, dat$Cv]

  covariates <- outer(as.vector(post$b_deg), dat$deg_h_z) +
    outer(as.vector(post$b_precip), dat$precip_72h_z) +
    outer(as.vector(post$b_seq), dat$seq_depth_z)

  fixed_MgMo    <- loga_obs + gamma_term
  fixed_MgMoY   <- fixed_MgMo + yr_obs
  fixed_MgMoYCo <- fixed_MgMoY + covariates
  fixed_full    <- fixed_MgMoYCo + cv_obs

  var_MgMo  <- apply(fixed_MgMo,    1, var)
  var_Y     <- apply(fixed_MgMoY,   1, var) - var_MgMo
  var_Cov   <- apply(fixed_MgMoYCo, 1, var) - apply(fixed_MgMoY, 1, var)
  var_Cv    <- apply(fixed_full,    1, var) - apply(fixed_MgMoYCo, 1, var)
  explained <- apply(fixed_full,    1, var)

  mg_n <- as.integer(table(factor(dat$Mg, levels = 1:2)))
  residual_var <- as.vector((post$sigma^2) %*% (mg_n / sum(mg_n)))

  total <- explained + residual_var

  bind_rows(
    tibble(statistic = "Variance partition", group = "Management x Season", value = var_MgMo / total),
    tibble(statistic = "Variance partition", group = "Year", value = var_Y / total),
    tibble(statistic = "Variance partition", group = "Covariates", value = var_Cov / total),
    tibble(statistic = "Variance partition", group = "Cultivar", value = var_Cv / total),
    tibble(statistic = "Variance partition", group = "Residual", value = residual_var / total)
  )
}

## Data-generating function ---------------------------------------------------
# Same balanced Mg x Mo x Year design as sim_div_MDSYC(), plus an
# independent Cv draw per unit (5 categories) and its cv_eff offset.
# Cultivar assigned independently of Mg here (unlike the real data's mild
# imbalance) -- tests whether cv_1..cv_4 are recoverable in principle, not
# how identifiable they are under the real design's modest correlation with
# Management.
#
# rho_deg_season/rho_precip_season/rho_seq_mg (all default 0, i.e. unchanged
# original behaviour -- independent draws) optionally correlate deg_h_z/
# precip_72h_z with Season and seq_depth_z with Management instead, at
# approximately the given correlation -- a stress test for identifiability
# under realistic collinearity, not just independent-covariate simulation.
# Standard target-correlation construction: z = rho*scale(x) + sqrt(1-rho^2)*noise
# (same template as the ITS lineage's own MDLSYC_model.R). Real data has all
# three: cor(deg_h_z, Mo)=0.73, cor(precip_72h_z, Mo)=0.55,
# cor(seq_depth_z, Mg)=0.31 -- see 6.5_MDSYCV_16S_collinearity_check.R.

sim_div_MDSYCV <- function(N_samples, loga, s_conv, gap_shift, sigma, yr1, yr2,
                            cv_1, cv_3, cv_4, cv_5, b_deg, b_precip, b_seq, shift = NULL,
                            rho_deg_season = 0, rho_precip_season = 0, rho_seq_mg = 0){
  n_unit <- N_samples %/% 2 # 2 rows/unit (May + July)
  yr_vec <- c(yr1, yr2, -(yr1 + yr2))
  cv_vec <- c(cv_1, -(cv_1 + cv_3 + cv_4 + cv_5), cv_3, cv_4, cv_5) # Cv index order: 1..5

  units <- tibble(
    Un = seq_len(n_unit), Mg = rbern(n_unit) + 1,
    Yr = sample(3, n_unit, replace = TRUE), Cv = sample(5, n_unit, replace = TRUE))
  dat <- units %>% crossing(Mo = 1:2) %>% arrange(Un)

  gamma <- s_conv + gap_shift*(dat$Mg - 1)
  mu_structural <- loga[dat$Mg] + gamma*(dat$Mo - 1) + yr_vec[dat$Yr] + cv_vec[dat$Cv]

  season_z <- as.vector(scale(dat$Mo - 1))
  dat$deg_h_z <- if (rho_deg_season == 0) rnorm(nrow(dat)) else
    rho_deg_season * season_z + sqrt(1 - rho_deg_season^2) * rnorm(nrow(dat))
  dat$precip_72h_z <- if (rho_precip_season == 0) rnorm(nrow(dat)) else
    rho_precip_season * season_z + sqrt(1 - rho_precip_season^2) * rnorm(nrow(dat))

  mg_z <- as.vector(scale(dat$Mg - 1))
  dat$seq_depth_z <- if (rho_seq_mg == 0) rnorm(nrow(dat)) else
    rho_seq_mg * mg_z + sqrt(1 - rho_seq_mg^2) * rnorm(nrow(dat))

  mu <- mu_structural + b_deg*dat$deg_h_z + b_precip*dat$precip_72h_z + b_seq*dat$seq_depth_z

  dat$Dv <- rlnorm(nrow(dat), meanlog = mu, sdlog = sigma[dat$Mg])
  if (!is.null(shift)) dat$Dv_shifted <- shift + dat$Dv
  dat
}

simulate_from_priors_MDSYCV <- function(true_params, N_samples = 250, shift = NULL,
                                         rho_deg_season = 0, rho_precip_season = 0, rho_seq_mg = 0){
  sim_div_MDSYCV(
    N_samples = N_samples,
    loga = true_params$loga,
    s_conv = true_params$s_conv,
    gap_shift = true_params$gap_shift,
    sigma = true_params$sigma,
    yr1 = true_params$yr1,
    yr2 = true_params$yr2,
    cv_1 = true_params$cv_1,
    cv_3 = true_params$cv_3,
    cv_4 = true_params$cv_4,
    cv_5 = true_params$cv_5,
    b_deg = true_params$b_deg,
    b_precip = true_params$b_precip,
    b_seq = true_params$b_seq,
    shift = shift,
    rho_deg_season = rho_deg_season,
    rho_precip_season = rho_precip_season,
    rho_seq_mg = rho_seq_mg
  )
}
