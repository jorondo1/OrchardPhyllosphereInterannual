# MDST_model.R --- MODEL 3 (MDST, "Treebeard the Skeptic"): MDS + Tree random effect
# - non-centered: tr[Tr]*sigma_tr
# Why: each tree gives a May and a July sample, treated as independent in MDS
# - overstates effective N; tree-to-tree variation ends up in sigma[Mg]
# Known fragility: small sigma_tr + 1-2 obs/tree -> funnel (divergences, Rhat)

source('src/hiermod/0_INDEX.R')

# Priors: as MDS; new sigma_tr ~ dhalfnorm(0,1) (same family as sigma[Mg])

model_MDST_16S <- alist(
  likelihood = Dv ~ dlnorm(mu, sigma[Mg]),
  main_model = mu <- loga[Mg] + gamma*(Mo-1) + tr[Tr]*sigma_tr,
  gamma_def  = gamma <- s_conv + gap_shift*(Mg-1), # interactive Season effect

  prior_loga  = loga[Mg]  ~ dnorm(5,2),
  prior_s     = s_conv    ~ dnorm(0,1),
  prior_gs    = gap_shift ~ dnorm(0,1),
  pr_sigma    = sigma[Mg] ~ dhalfnorm(0,1),

  prior_tr    = tr[Tr]    ~ dnorm(0,1),
  pr_sigma_tr = sigma_tr  ~ dhalfnorm(0,1)
)

attr(model_MDST_16S, "name") <- "Treebeard the Skeptic"
model_id_MDST <- "MDST"

## ITS variant -----------------------------------------------------------------
# Intercept prior only (as MDS_ITS)
model_MDST_ITS <- model_MDST_16S
model_MDST_ITS$prior_loga <- quote(loga[Mg] ~ dnorm(2,2))
attr(model_MDST_ITS, "name") <- "Treebeard the Skeptic"

## Backtransform wrapper ------------------------------------------------------
# 4-cell means/medians; mean variance = sigma[Mg]^2 + sigma_tr^2

means_MDST <- function(post, shift = 0){
  total_var_conv <- post$sigma[,1]^2 + as.vector(post$sigma_tr)^2
  total_var_org  <- post$sigma[,2]^2 + as.vector(post$sigma_tr)^2

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

# SBC estimands: May gap, July gap, seasonal change (variance incl. sigma_tr^2)
dq_MDST <- SBC::derived_quantities(
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
# Tree = study unit: one Management, one May + one July row, shared tree offset

sim_div_MDST <- function(N_samples, loga, s_conv, gap_shift, sigma, sigma_tr, shift = NULL){
  n_tree <- N_samples %/% 2 # 2 rows/tree (May + July)

  trees <- tibble(Tr = seq_len(n_tree), Mg = rbern(n_tree) + 1)
  dat <- trees %>% crossing(Mo = 1:2) %>% arrange(Tr)

  tree_offset <- rnorm(n_tree, 0, sigma_tr)

  gamma <- s_conv + gap_shift*(dat$Mg - 1)
  mu <- loga[dat$Mg] + gamma*(dat$Mo - 1) + tree_offset[dat$Tr]

  dat$Dv <- rlnorm(nrow(dat), meanlog = mu, sdlog = sigma[dat$Mg])
  if (!is.null(shift)) dat$Dv_shifted <- shift + dat$Dv
  dat
}

# Simulate one dataset from one prior draw (draw_true() output)
simulate_from_priors_MDST <- function(true_params, N_samples = 250, shift = NULL){
  sim_div_MDST(
    N_samples = N_samples,
    loga = true_params$loga,
    s_conv = true_params$s_conv,
    gap_shift = true_params$gap_shift,
    sigma = true_params$sigma,
    sigma_tr = true_params$sigma_tr,
    shift = shift
  )
}
