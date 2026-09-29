# MDSTYCV_model.R --- MODEL 7 (MDSTYCV), 16S: MDSYCV plus Tree, merging the
# two branches this rebuild has kept separate since Model 3 (MDST, Tree)
# and Model 4 (MDSYz, Year, later +Covariates +Cultivar).

source('src/hiermod/0_INDEX.R')

# This is not a purely mechanical merge. Direct query of the real data
# confirms Tree is DETERMINISTICALLY NESTED in Cultivar (129/129 trees map
# to exactly one cultivar) and in Location (129/129, though Location still
# isn't in any model here). Cultivar's fixed effect (cv_1..cv_4, ~26 trees
# each) and Tree's own random effect (tr[Tr]*sigma_tr, 2 obs/tree) now
# share the exact same 129 trees for the first time -- worth testing
# directly rather than assuming it's fine, especially since sigma_tr
# already has documented, real fragility of its own (MDST's own SBC:
# divergences + a U-shaped rank histogram tied to sparse per-tree N, see
# MDST_model.R's header and 3.2_MDST_16S_calibration.R). The calibration
# script (7.2) is built around exactly this question: does sigma_tr's
# known fragility get better, worse, or stay the same once Cultivar/Year/
# covariates are also in the model, and is it identifiable jointly with
# cv_1..cv_4 (never tested together before)?
#
# loga[Mg]/s_conv/gap_shift/sigma[Mg]/yr1/yr2/cv_1/cv_3/cv_4/cv_5/b_deg/
# b_precip/b_seq are MDSYCV's own validated answer, hardcoded as this
# model's starting point. sigma_tr ~ dhalfnorm(0,1) is MDST's own validated
# answer for that parameter (the tightened dlnorm(log(0.15),0.5) attempt
# didn't fully resolve MDST's own fragility either, so no reason to prefer
# it over the simpler, equally-imperfect dhalfnorm(0,1) as a starting
# point here) -- tr[Tr] ~ dnorm(0,1) is the standard non-centered form.

source('src/hiermod/Models/MDSYCV_model.R') # model_MDSYCV_16S, means_MDSYCV(), dq_MDSYCV

model_MDSTYCV_16S <- model_MDSYCV_16S
model_MDSTYCV_16S$main_model <- quote(
  mu <- loga[Mg] + gamma*(Mo-1) + yr_eff + cv_eff + tr[Tr]*sigma_tr +
    b_deg*deg_h_z + b_precip*precip_72h_z + b_seq*seq_depth_z
)
model_MDSTYCV_16S$prior_tr    <- quote(tr[Tr]   ~ dnorm(0,1))
model_MDSTYCV_16S$pr_sigma_tr <- quote(sigma_tr ~ dhalfnorm(0,1))

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
# Reports a marginal ("by margin", order-free) decomposition via the shared
# variance_partition_panels() engine (postcontrast_helpers.R).
#
# Management x Season is kept as ONE combined term (loga[Mg] + gamma_term),
# not split into Management/Season/interaction -- splitting them was tried
# and produced a strongly negative "by margin" share for the interaction
# (median ~ -0.33 in the real 16S fit). That's not confounding in the usual
# ANOVA sense -- it's that this function's "marginal" share is a no-refit
# shortcut (var(full) - var(full minus term), holding every OTHER
# coefficient fixed at its joint full-model value), which lacks the
# projection-based non-negativity guarantee a real refit-based ANOVA/
# PERMANOVA Type III SS has. Management (loga[Mg]) and the interaction
# (gap_shift*(Mg-1)*(Mo-1)) are raw 0/1-indicator coded (unlike Year/
# Cultivar's sum-to-zero coding), so they aren't orthogonal by construction
# -- in this fit, gap_shift and the loga[Mg] contrast have opposite signs
# (Organic starts lower in May but climbs more steeply by July), which
# drives a large negative cross term. Combining them back into one term
# sidesteps the issue entirely (a term's marginal share vs. itself is just
# its own variance). method = "lmg" gives the Shapley/LMG version instead
# (variance_partition_lmg(), non-negative shares).
#
# Tree is a proper per-observation term (realized tr[Tr] draws x sigma_tr,
# like Year/Cultivar's own realized-level construction) instead of the old
# bulk sigma_tr^2/total bolt-on that was commented out of this model's total
# entirely -- Tree is included in the explained/total variance here for the
# first time, matching variance_partition_MDSTYCL()'s own (already-included)
# treatment for the first time.

variance_partition_MDSTYCV <- function(post, dat, method = c("margin", "lmg")){
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
    "Tree"                 = tree_obs,
    "Management x Season"  = loga_obs + gamma_term
  )

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
