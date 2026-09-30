# MDSTYCL_model.R --- MODEL 9 (MDSTYCL, "Faramir the Judicious"), 16S:
# MDSTYCV with Cultivar's fixed effect (cv[Cv]) REPLACED by Location's
# fixed effect (lo[Lo]), same sum-to-zero recipe.

source('src/hiermod/0_INDEX.R')

# Why this exists: Model 8 (MDSTYCVr, "Gimli the Greedy") tried making
# Cultivar a RANDOM effect instead and failed badly -- real SBC
# (n_sbc=500): 297 divergences, 154/500 (30.8%) fits Rhat>1.01, and
# loga[1]/loga[2] themselves severely MISCALIBRATED (z=-11.96/-14.22),
# not just sigma_cv (z=-5.76). Nesting sigma_cv with sigma_tr over the
# same 129 trees didn't just fail to identify sigma_cv cleanly -- it
# leaked into the headline estimand. Back to a fixed-effect design instead
# of trying to fix that nesting.
#
# Location is a natural next fixed effect to test (never in this rebuild's
# family before -- Tree is deterministically nested in Location, same
# 129/129 mapping as Cultivar, so structurally this is the same "does
# Tree's random effect coexist with a coarser fixed grouping over the same
# trees" question MDSTYCV already answered cleanly for Cultivar).
#
# Cultivar can NOT coexist with Location here, though -- confirmed by
# direct query of the real data, not just assumed: restricted to the two
# Locations with both Management levels present (B, D -- see this
# model's own calibration/fit scripts for why), Location D has ONLY
# Honeycrisp and Spartan (0 Cortland/Liberty/Paulared) -- i.e. within that
# subset, 3 of Cultivar's 5 levels are perfectly aliased with "Location B",
# not just correlated with it. Cultivar is dropped entirely for this
# model, not added alongside Location.
#
# Real-data motivation for the subset (see 9.2/9.3's own headers for the
# full numbers): Location x Management in the full data is A=Conventional-
# only (50), C=Organic-only (51), B and D have both -- fitting Location
# on the full 4-level factor would confound Location with Management for
# A/C. Model 9's real fit (9.3, not yet built) will restrict to the B/D
# subset (141/242 rows) for exactly this reason. This calibration-stage
# model file itself doesn't know about that restriction -- Lo is just a
# 2-level index, assigned independently of Mg in the simulator, same
# "recoverable in principle" scope as every other calibration script in
# this family.
#
# loga[Mg]/s_conv/gap_shift/sigma[Mg]/yr1/yr2/tr[Tr]*sigma_tr/b_deg/
# b_precip/b_seq are MDSTYCV's own validated answer, hardcoded as this
# model's starting point. lo1 ~ dnorm(0,1) is the one new assumption to
# validate -- same sum-to-zero recipe as yr1/yr2 and cv_1/cv_3/cv_4/cv_5,
# just 2 levels (1 free parameter, 1 derived as its negative) instead of
# 3 or 5.
#
# ITS variant deliberately NOT built yet -- same discipline as Model 8:
# validate on 16S first.

model_MDSTYCL_16S <- alist(
  likelihood = Dv ~ dlnorm(mu, sigma[Mg]),
  main_model = mu <- loga[Mg] + gamma*(Mo-1) + yr_eff + lo_eff + tr[Tr]*sigma_tr +
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

  # Tree, non-centered random effect
  prior_tr    = tr[Tr]   ~ dnorm(0,1),
  pr_sigma_tr = sigma_tr ~ dhalfnorm(0,1),

  # Location (2 levels within the real B/D subset), sum-to-zero: lo1 free,
  # lo2 = -lo1 -- same construction as Year/Cultivar, N-1=1 free here
  lo_eff_def = lo_eff <- lo1*(Lo==1) - lo1*(Lo==2),
  prior_lo1  = lo1 ~ dnorm(0,1)
)

attr(model_MDSTYCL_16S, "name") <- "Faramir the Judicious"
model_id_MDSTYCL <- "MDSTYCL"

## means_MDSTYCL()/dq_MDSTYCL ----------------------------------------------------
# Same formulas as means_MDSTYCV()/dq_MDSTYCV -- Location, like Year, is a
# fixed effect held at its own observed-level average and doesn't enter
# the reported Mg x Mo estimand or its variance (same reasoning as
# means_MDSYCV() dropping Cultivar).

means_MDSTYCL <- function(post, shift = 0, deg_h_z = 0, precip_72h_z = 0, seq_depth_z = 0){
  total_var_conv <- post$sigma[,1]^2 + as.vector(post$sigma_tr)^2
  total_var_org  <- post$sigma[,2]^2 + as.vector(post$sigma_tr)^2
  
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

dq_MDSTYCL <- SBC::derived_quantities(
  may_gap =
    exp(loga[2] + (sigma[2]^2 + sigma_tr^2) / 2) -
    exp(loga[1] + (sigma[1]^2 + sigma_tr^2) / 2),
  july_gap =
    exp(loga[2] + s_conv + gap_shift + (sigma[2]^2 + sigma_tr^2) / 2) -
    exp(loga[1] + s_conv +             (sigma[1]^2 + sigma_tr^2) / 2),
  seasonal_change =
    (exp(loga[2] + s_conv + gap_shift + (sigma[2]^2 + sigma_tr^2) / 2) -
       exp(loga[1] + s_conv +             (sigma[1]^2 + sigma_tr^2) / 2)) -
    (exp(loga[2] + (sigma[2]^2 + sigma_tr^2) / 2) -
       exp(loga[1] + (sigma[1]^2 + sigma_tr^2) / 2))
)

## variance_partition_MDSTYCL() --------------------------------------------------
# Same as variance_partition_MDSTYCV() (see its own comment: Shapley/LMG
# default, effect-coded Management / Season / interaction split) -- Location
# swapped in for Cultivar, in the same slot (coarser grouping before the
# Tree it nests, same 129/129 deterministic mapping as Cultivar).

variance_partition_MDSTYCL <- function(post, dat, method = c("lmg", "margin"), split_mgmo = TRUE){
  method <- match.arg(method)
  yr3    <- -(as.vector(post$yr1) + as.vector(post$yr2))
  yr_obs <- cbind(post$yr1, post$yr2, yr3)[, dat$Yr]

  lo2    <- -(as.vector(post$lo1))
  # NB: was `[dat$Lo]` (no comma) -- linear/column-major indexing on a
  # matrix, not column selection by dat$Lo like every other realized-level
  # term here (yr_obs, cv_obs). That silently returned near-constant
  # garbage instead of each observation's actual Location coefficient --
  # the likely cause of Location's R2 share collapsing to ~0 (and sometimes
  # negative).
  lo_obs <- cbind(post$lo1, lo2)[, dat$Lo]

  tree_obs <- sweep(post$tr[, dat$Tr], 1, as.vector(post$sigma_tr), "*")

  mg_n <- as.integer(table(factor(dat$Mg, levels = 1:2)))
  residual_var <- as.vector((post$sigma^2) %*% (mg_n / sum(mg_n)))

  loga_obs   <- post$loga[, dat$Mg]
  gamma      <- as.vector(post$s_conv) + outer(as.vector(post$gap_shift), dat$Mg - 1)
  gamma_term <- sweep(gamma, 2, dat$Mo - 1, "*")

  terms <- list(
    "Reads count"          = outer(as.vector(post$b_seq), dat$seq_depth_z),
    "Degree-hours"         = outer(as.vector(post$b_deg), dat$deg_h_z),
    "Precipitation"        = outer(as.vector(post$b_precip), dat$precip_72h_z),
    "Year"                 = yr_obs,
    "Location"             = lo_obs,
    "Tree"                 = tree_obs
  )
  if (split_mgmo) terms <- c(terms, mgmo_effect_terms(loga_obs + gamma_term, dat$Mg, dat$Mo))
  else terms[["Management x Season"]] <- loga_obs + gamma_term

  if (method == "lmg") variance_partition_lmg(terms, residual_var)
  else variance_partition_panels(terms, residual_var)
}

## Data-generating function ---------------------------------------------------
# Same Tree-as-study-unit design as sim_div_MDSTYCV(), Lo assigned per
# tree (2 levels) exactly like Cv was -- independently of Mg here (tests
# recoverability in principle; the real B/D subset's actual Location x
# Cultivar aliasing is a real-fit-stage concern, not a calibration one,
# see this file's own header).

sim_div_MDSTYCL <- function(N_samples, loga, s_conv, gap_shift, sigma, yr1, yr2,
                             lo1, sigma_tr, b_deg, b_precip, b_seq, shift = NULL){
  n_tree <- N_samples %/% 2 # 2 rows/tree (May + July)
  yr_vec <- c(yr1, yr2, -(yr1 + yr2))
  lo_vec <- c(lo1, -lo1) # Lo index order: 1, 2

  trees <- tibble(
    Tr = seq_len(n_tree), Mg = rbern(n_tree) + 1,
    Yr = sample(3, n_tree, replace = TRUE), Lo = sample(2, n_tree, replace = TRUE))
  dat <- trees %>% crossing(Mo = 1:2) %>% arrange(Tr)

  tree_offset <- rnorm(n_tree, 0, sigma_tr)

  gamma <- s_conv + gap_shift*(dat$Mg - 1)
  mu_structural <- loga[dat$Mg] + gamma*(dat$Mo - 1) + yr_vec[dat$Yr] +
    lo_vec[dat$Lo] + tree_offset[dat$Tr]

  dat$deg_h_z      <- rnorm(nrow(dat))
  dat$precip_72h_z <- rnorm(nrow(dat))
  dat$seq_depth_z  <- rnorm(nrow(dat))

  mu <- mu_structural + b_deg*dat$deg_h_z + b_precip*dat$precip_72h_z + b_seq*dat$seq_depth_z

  dat$Dv <- rlnorm(nrow(dat), meanlog = mu, sdlog = sigma[dat$Mg])
  if (!is.null(shift)) dat$Dv_shifted <- shift + dat$Dv
  dat
}

simulate_from_priors_MDSTYCL <- function(true_params, N_samples = 250, shift = NULL){
  sim_div_MDSTYCL(
    N_samples = N_samples,
    loga = true_params$loga, s_conv = true_params$s_conv, gap_shift = true_params$gap_shift,
    sigma = true_params$sigma, yr1 = true_params$yr1, yr2 = true_params$yr2,
    lo1 = true_params$lo1, sigma_tr = true_params$sigma_tr,
    b_deg = true_params$b_deg, b_precip = true_params$b_precip, b_seq = true_params$b_seq,
    shift = shift
  )
}
