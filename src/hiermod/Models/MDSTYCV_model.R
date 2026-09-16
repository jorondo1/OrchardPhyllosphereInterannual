# MDSTYCV_model.R --- MODEL 7 (MDSTYCV), 16S: MDSYCV plus Tree, merging the
# two branches this rebuild has kept separate since Model 3 (MDST, Tree)
# and Model 4 (MDSYz, Year, later +Covariates +Cultivar).
#
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

model_id_MDSTYCV <- "MDSTYCV"

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
