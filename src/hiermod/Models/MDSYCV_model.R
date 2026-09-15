# MDSYCV_model.R --- MODEL 6 (MDSYCV, "Bombadil the Eldest"), 16S: MDSYC
# plus Cultivar as a FIXED, sum-to-zero effect (5 levels: Cortland, Liberty,
# Paulared, Honeycrisp, Spartan).
#
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

source('src/hiermod/Models/MDSYC_model.R') # model_MDSYC_16S, means_MDSYC(), dq_MDSYC

model_MDSYCV_16S <- model_MDSYC_16S
model_MDSYCV_16S$main_model <- quote(
  mu <- loga[Mg] + gamma*(Mo-1) + yr_eff + cv_eff +
    b_deg*deg_h_z + b_precip*precip_72h_z + b_seq*seq_depth_z
)
model_MDSYCV_16S$cv_eff_def <- quote(
  cv_eff <- cv_1*(Cv==1) + cv_3*(Cv==3) + cv_4*(Cv==4) + cv_5*(Cv==5) - (cv_1+cv_3+cv_4+cv_5)*(Cv==2)
)
model_MDSYCV_16S$prior_cv1 <- quote(cv_1 ~ dnorm(0,1))
model_MDSYCV_16S$prior_cv3 <- quote(cv_3 ~ dnorm(0,1))
model_MDSYCV_16S$prior_cv4 <- quote(cv_4 ~ dnorm(0,1))
model_MDSYCV_16S$prior_cv5 <- quote(cv_5 ~ dnorm(0,1))

model_id_MDSYCV <- "MDSYCV"

## means_MDSYCV()/dq_MDSYCV ----------------------------------------------------
# Identical to means_MDSYC()/dq_MDSYC -- Cultivar, like Year, doesn't enter
# the reported Mg x Mo estimand or its variance (assigned independently of
# Mg x Mo in the simulator, reported at the default/average level on real
# data).
means_MDSYCV <- means_MDSYC
dq_MDSYCV    <- dq_MDSYC

## variance_partition_MDSYCV() --------------------------------------------------
# Same as variance_partition_MDSYC() plus cv_eff in the fixed-effects sum.
# cv_derived (Liberty) computed the same way the model itself does; cv_mat's
# columns are built directly in Cv's own 1..5 index order so dat$Cv can
# index it with no remapping.

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

  fixed_mu  <- loga_obs + gamma_term + yr_obs + cv_obs + covariates
  explained <- apply(fixed_mu, 1, var)

  mg_n <- as.integer(table(factor(dat$Mg, levels = 1:2)))
  residual_var <- as.vector((post$sigma^2) %*% (mg_n / sum(mg_n)))

  total <- explained + residual_var

  bind_rows(
    tibble(statistic = "Variance partition", group = "Explained (fixed effects)", value = explained / total),
    tibble(statistic = "Variance partition", group = "Residual",                  value = residual_var / total)
  )
}

## Data-generating function ---------------------------------------------------
# Same balanced Mg x Mo x Year design as sim_div_MDSYC(), plus an
# independent Cv draw per unit (5 categories) and its cv_eff offset.
# Cultivar assigned independently of Mg here (unlike the real data's mild
# imbalance) -- tests whether cv_1..cv_4 are recoverable in principle, not
# how identifiable they are under the real design's modest correlation with
# Management.

sim_div_MDSYCV <- function(N_samples, loga, s_conv, gap_shift, sigma, yr1, yr2,
                            cv_1, cv_3, cv_4, cv_5, b_deg, b_precip, b_seq, shift = NULL){
  n_unit <- N_samples %/% 2 # 2 rows/unit (May + July)
  yr_vec <- c(yr1, yr2, -(yr1 + yr2))
  cv_vec <- c(cv_1, -(cv_1 + cv_3 + cv_4 + cv_5), cv_3, cv_4, cv_5) # Cv index order: 1..5

  units <- tibble(
    Un = seq_len(n_unit), Mg = rbern(n_unit) + 1,
    Yr = sample(3, n_unit, replace = TRUE), Cv = sample(5, n_unit, replace = TRUE))
  dat <- units %>% crossing(Mo = 1:2) %>% arrange(Un)

  gamma <- s_conv + gap_shift*(dat$Mg - 1)
  mu_structural <- loga[dat$Mg] + gamma*(dat$Mo - 1) + yr_vec[dat$Yr] + cv_vec[dat$Cv]

  dat$deg_h_z      <- rnorm(nrow(dat))
  dat$precip_72h_z <- rnorm(nrow(dat))
  dat$seq_depth_z  <- rnorm(nrow(dat))

  mu <- mu_structural + b_deg*dat$deg_h_z + b_precip*dat$precip_72h_z + b_seq*dat$seq_depth_z

  dat$Dv <- rlnorm(nrow(dat), meanlog = mu, sdlog = sigma[dat$Mg])
  if (!is.null(shift)) dat$Dv_shifted <- shift + dat$Dv
  dat
}

simulate_from_priors_MDSYCV <- function(true_params, N_samples = 250, shift = NULL){
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
    shift = shift
  )
}
