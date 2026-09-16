# MDSYC_model.R --- MODEL 5 (MDSYC), 16S: MDSYz plus three standardized
# control covariates (deg_h_z, precip_72h_z, seq_depth_z) as additive fixed
# slopes -- degree-hours and 72h precipitation before sampling (recent
# growing conditions), and log sequencing depth (technical: deeper
# sequencing detects more taxa, inflating diversity metrics independent of
# any real biological effect). All three already computed in 0_SETUP.R.
#
# Continues the MDS -> MDSYz branch specifically (Year, sum-to-zero fixed),
# not yet merged with the separate MDST (Tree) branch -- one thing at a
# time, same discipline as the rest of this rebuild. Location stays out
# entirely for now, deferred to its own later sensitivity check (see
# MDS2_model.R's header on why: Location A/C are single-management,
# confounded with Management itself, not a random-effect identifiability
# problem like Year's was).
#
# Direct precedent: the ITS lineage's own Model 7 (MDLSYC,
# src/hiermod/Models/MDLSYC_model.R) added the identical three covariates
# the identical way (b_deg/b_precip/b_seq ~ dnorm(0,1), additive in mu).
# That model's own "K stays at 4" comment applies here too, even more
# simply: these are additive mu-level fixed effects, not new summed
# variance terms, so sigma[Mg]'s own already-validated dhalfnorm(0,1) prior
# needs no rescaling at all.
#
# loga[Mg]/s_conv/gap_shift/sigma[Mg]/yr1/yr2 are MDSYz's own validated
# answer, hardcoded as this model's starting point. b_deg/b_precip/b_seq
# ~ dnorm(0,1) is the one new assumption -- MDLSYC's own choice, carried
# forward as the hypothesis to validate here.

source('src/hiermod/Models/MDSYz_model.R') # model_MDSYz_16S, means_MDSYz(), dq_MDSYz

model_MDSYC_16S <- model_MDSYz_16S
model_MDSYC_16S$main_model <- quote(
  mu <- loga[Mg] + gamma*(Mo-1) + yr_eff +
    b_deg*deg_h_z + b_precip*precip_72h_z + b_seq*seq_depth_z
)
model_MDSYC_16S$prior_deg    <- quote(b_deg    ~ dnorm(0,1))
model_MDSYC_16S$prior_precip <- quote(b_precip ~ dnorm(0,1))
model_MDSYC_16S$prior_seq    <- quote(b_seq    ~ dnorm(0,1))

attr(model_MDSYC_16S, "name") <- "Radagast the Grower"
model_id_MDSYC <- "MDSYC"

## ITS variant -----------------------------------------------------------------
# Same rationale as MDS_ITS. b_deg/b_precip/b_seq ~ dnorm(0,1) carry over
# unchanged -- these are additive log-scale slopes on standardized (z-score)
# covariates, not tied to Hill_1's own baseline scale, and the old ITS
# lineage's own MDLSYC_model.R already used the identical dnorm(0,1) choice.
model_MDSYC_ITS <- model_MDSYC_16S
model_MDSYC_ITS$prior_loga <- quote(loga[Mg] ~ dnorm(2,2))
attr(model_MDSYC_ITS, "name") <- "Radagast the Grower"

## Backtransform wrapper ------------------------------------------------------
# Same 4-cell shape as means_MDSYz(), plus a covariate_offset term
# (deg_h_z/precip_72h_z/seq_depth_z default to 0, i.e. this sample's own
# average weather and sequencing depth -- same convention as MDLSYC's own
# means_fn).

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

# SBC estimands -- identical to dq_MDSYz. b_deg/b_precip/b_seq cancel in the
# Mg x Mo contrast at the default (z=0) reference level, same reasoning as
# MDLSYC's own dq_MDLSYC.
dq_MDSYC <- dq_MDSYz

## Variance partition / Bayesian R2 -------------------------------------------
# Adapted from the ITS lineage's own variance_partition_MDLSYC()
# (MDLSYC_model.R). Residual is sigma[Mg]^2 -- genuinely heteroscedastic by
# Management, unlike the old lineage's single pooled/cell-level sigma --
# weighted by each group's share of the sample (mirrors the old
# cell_n/sum(cell_n) weighting).
#
# Fixed effects split by GROUP (Management x Season, Year, Covariates), not
# lumped into one "Explained" bucket -- via a sequential (Type I)
# decomposition: each group's share is the variance ADDED by including it,
# in this model family's own build order (Mg x Season, Model 2 -> Year,
# Model 4 -> Covariates, Model 5). This telescopes exactly to the same total
# as the old lumped version (Var(A) + [Var(A+B)-Var(A)] + [Var(A+B+C)-
# Var(A+B)] = Var(A+B+C)), so it's a genuine refinement, not a different
# number -- order-dependent in principle if predictors are correlated
# (close to orthogonal here by design), and an individual increment can
# come out slightly negative if a later term happens to reduce a given
# draw's cumulative variance.

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
# Same balanced Mg x Mo x Year design as sim_div_MDSYz(), plus three
# independent standard-normal covariate draws per row (tests whether
# b_deg/b_precip/b_seq are recoverable in principle, not how identifiable
# they are under any real-world collinearity with season/structural mu --
# MDLSYC_model.R's own rho_deg_season/rho_seq_mu stress test is the
# template if that's worth adding here later).

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
