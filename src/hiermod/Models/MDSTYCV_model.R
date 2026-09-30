# MDSTYCV_model.R --- MODEL 7 (MDSTYCV), 16S: MDSYCV plus Tree, merging the
# two branches this rebuild has kept separate since Model 3 (MDST, Tree)
# and Model 4 (MDSYz, Year, later +Covariates +Cultivar).

source('src/hiermod/0_INDEX.R')

# Nested Tree-Location-Cultivar structure, testing because sigma_tree
# was identified as fragile (MDST's own SBC:
# divergences + a U-shaped rank histogram tied to sparse per-tree N, see
# MDST_model.R and 3.2_MDST_16S_calibration.R). Calibration is made
# to test whether that nestedness carries the same fragility.
#
# loga[Mg]/s_conv/gap_shift/sigma[Mg]/yr1/yr2/cv_1/cv_3/cv_4/cv_5/b_deg/
# b_precip/b_seq are MDSYCV's own validated answer, hardcoded as this
# model's starting point. 

model_MDSTYCV_16S <- alist(
  likelihood = Dv ~ dlnorm(mu, sigma[Mg]),
  main_model = mu <- loga[Mg] + gamma*(Mo-1) + yr_eff + cv_eff + tr[Tr]*sigma_tr +
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
  prior_cv5  = cv_5 ~ dnorm(0,1),

  # Tree, non-centered random effect
  prior_tr    = tr[Tr]   ~ dnorm(0,1),
  pr_sigma_tr = sigma_tr ~ dhalfnorm(0,1)
)

attr(model_MDSTYCV_16S, "name") <- "Saruman the Fool"
model_id_MDSTYCV <- "MDSTYCV"

## ITS variant -----------------------------------------------------------------
# Same rationale as MDS_ITS/MDST_ITS. tr[Tr]/sigma_tr carry over unchanged
# -- same Tree factor (idx$Tr), same non-centered form, same CV-like
# sigma_tr scale reasoning.
model_MDSTYCV_ITS <- model_MDSTYCV_16S
model_MDSTYCV_ITS$prior_loga <- quote(loga[Mg] ~ dnorm(2,2))
attr(model_MDSTYCV_ITS, "name") <- "Saruman the Fool"

## means_MDSTYCV() -------------------------------------------------------------
# Same as means_MDSYCV(), plus sigma_tr^2 folded into total_var -- matching
# MDST's own means_MDST() convention (sigma[Mg]^2 + sigma_tr^2).

means_MDSTYCV <- function(post, shift = 0, deg_h_z = 0, precip_72h_z = 0, seq_depth_z = 0){
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

# SBC estimands -- same formula as MDST's own dq_MDST (sigma[Mg]^2 +
# sigma_tr^2; yr/cv/covariates all cancel in the Mg x Mo contrast, same
# reasoning as dq_MDSYCV).
dq_MDSTYCV <- SBC::derived_quantities(
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

## variance_partition_MDSTYCV() --------------------------------------------------
# Bayesian R2 partition (see doc/R2_methods.txt). Default: Shapley/LMG shares
# (variance_partition_lmg(), non-negative, sum to the explained fraction) with
# Management x Season split into effect-coded Management / Season /
# interaction terms (mgmo_effect_terms(), postcontrast_helpers.R).
#
# Why effect coding: the model's own Management (loga[Mg]) and interaction
# (gap_shift*(Mg-1)*(Mo-1)) terms are 0/1-dummy coded, so they overlap
# heavily -- splitting them as-is gave a strongly negative "by margin"
# interaction share (median ~ -0.33 in the real 16S fit), and even Shapley
# over-credits the interaction under dummy coding. Recentring each piece
# around the grand mean makes the three orthogonal, so the split no longer
# depends on the reference level.
#
# method = "margin" gives the older no-refit shortcut (variance_partition_panels(),
# can go negative); split_mgmo = FALSE keeps Management x Season as one term.
#
# Tree is a per-observation term (realized tr[Tr] draws x sigma_tr), like
# Year/Cultivar's own realized-level construction.

variance_partition_MDSTYCV <- function(post, dat, method = c("lmg", "margin"), split_mgmo = TRUE){
  method <- match.arg(method)
  yr3    <- -(as.vector(post$yr1) + as.vector(post$yr2))
  yr_obs <- cbind(post$yr1, post$yr2, yr3)[, dat$Yr]

  cv_derived <- -(as.vector(post$cv_1) + as.vector(post$cv_3) + as.vector(post$cv_4) + as.vector(post$cv_5))
  cv_mat     <- cbind(post$cv_1, cv_derived, post$cv_3, post$cv_4, post$cv_5) # Cv index order: 1..5
  cv_obs     <- cv_mat[, dat$Cv]

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
    "Cultivar"             = cv_obs,
    "Tree"                 = tree_obs
  )
  if (split_mgmo) terms <- c(terms, mgmo_effect_terms(loga_obs + gamma_term, dat$Mg, dat$Mo))
  else terms[["Management x Season"]] <- loga_obs + gamma_term

  if (method == "lmg") variance_partition_lmg(terms, residual_var)
  else variance_partition_panels(terms, residual_var)
}

## Data-generating function ---------------------------------------------------
# Tree is the study unit (as in sim_div_MDST()): each tree gets one
# Management, one Cultivar, one Year, and exactly one May + one July row.
# Cv and Yr are assigned PER TREE, not per row (real, confirmed nesting for
# Cv; already the de facto behaviour for Yr in sim_div_MDSYCV()'s own
# per-unit assignment, just renamed now that "unit" is explicitly Tree).
# This is the one genuine design change from a mechanical merge: MDSYCV's
# own simulator drew Cv independently per unit because Tree wasn't in that
# model at all, so there was no nesting to represent.

sim_div_MDSTYCV <- function(N_samples, loga, s_conv, gap_shift, sigma, yr1, yr2,
                            cv_1, cv_3, cv_4, cv_5, sigma_tr,
                            b_deg, b_precip, b_seq, shift = NULL){
  n_tree <- N_samples %/% 2 # 2 rows/tree (May + July)
  yr_vec <- c(yr1, yr2, -(yr1 + yr2))
  cv_vec <- c(cv_1, -(cv_1 + cv_3 + cv_4 + cv_5), cv_3, cv_4, cv_5) # Cv index order: 1..5
  
  trees <- tibble(
    Tr = seq_len(n_tree), Mg = rbern(n_tree) + 1,
    Yr = sample(3, n_tree, replace = TRUE), Cv = sample(5, n_tree, replace = TRUE))
  dat <- trees %>% crossing(Mo = 1:2) %>% arrange(Tr)
  
  tree_offset <- rnorm(n_tree, 0, sigma_tr)
  
  gamma <- s_conv + gap_shift*(dat$Mg - 1)
  mu_structural <- loga[dat$Mg] + gamma*(dat$Mo - 1) + yr_vec[dat$Yr] + cv_vec[dat$Cv] + tree_offset[dat$Tr]
  
  dat$deg_h_z      <- rnorm(nrow(dat))
  dat$precip_72h_z <- rnorm(nrow(dat))
  dat$seq_depth_z  <- rnorm(nrow(dat))
  
  mu <- mu_structural + b_deg*dat$deg_h_z + b_precip*dat$precip_72h_z + b_seq*dat$seq_depth_z
  
  dat$Dv <- rlnorm(nrow(dat), meanlog = mu, sdlog = sigma[dat$Mg])
  if (!is.null(shift)) dat$Dv_shifted <- shift + dat$Dv
  dat
}

## Prior simulation wrapper ---------------
simulate_from_priors_MDSTYCV <- function(true_params, N_samples = 250, shift = NULL){
  sim_div_MDSTYCV(
    N_samples = N_samples,
    loga = true_params$loga, s_conv = true_params$s_conv, gap_shift = true_params$gap_shift,
    sigma = true_params$sigma, yr1 = true_params$yr1, yr2 = true_params$yr2,
    cv_1 = true_params$cv_1, cv_3 = true_params$cv_3, cv_4 = true_params$cv_4, cv_5 = true_params$cv_5,
    sigma_tr = true_params$sigma_tr,
    b_deg = true_params$b_deg, b_precip = true_params$b_precip, b_seq = true_params$b_seq,
    shift = shift
  )
}
