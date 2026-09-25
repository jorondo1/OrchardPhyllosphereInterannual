# MDST_model.R --- MODEL 3 (MDST, "Treebeard the Skeptic"), 16S: MDS plus a
# Tree random effect, non-centered (tr[Tr]*sigma_tr), mirroring the archived
# MDLS_model.R's proven form for isolating Tree cleanly (historically only
# mild underconfidence there, not the severe bias Location caused).
#
# Motivation: every tree contributes one May row and one July row that MDS
# treats as independent draws from dlnorm(mu, sigma[Mg]).
# Not modelling that overstates effective N and risks absorbing real
# tree-to-tree variation into sigma[Mg] as pure noise. 

source('src/hiermod/0_INDEX.R')

# Priors: loga[Mg]/s_conv/gap_shift/sigma[Mg] are hardcoded at MDS's own
# validated answer (this model's starting point, per the rolling
# "each model bakes in what the previous one learned" convention) --
# sigma_tr ~ dhalfnorm(0,1) is new territory, but uses the same validated
# scale-parameter family (light-tailed, mode off zero, extract.prior()-safe)
# rather than reaching back for dexp() (MDS2's own original, untested,
# choice for sigma_tr was dexp(2)).

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
# Same rationale as MDS_ITS (MDS_model.R): only loga[Mg] is scale-dependent
# and gets ITS's own dnorm(2,2) starting point; sigma_tr ~ dhalfnorm(0,1)
# carries over unchanged (a CV-like quantity, not baseline-dependent).
model_MDST_ITS <- model_MDST_16S
model_MDST_ITS$prior_loga <- quote(loga[Mg] ~ dnorm(2,2))
attr(model_MDST_ITS, "name") <- "Treebeard the Skeptic"

## Backtransform wrapper ------------------------------------------------------
# Same 4-cell shape as means_MDS(), total_var now includes sigma_tr^2
# (matching MDS2's own total_var convention: sigma[Mg]^2 + sigma_tr^2).

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

# SBC estimands -- may_gap/july_gap/seasonal_change, matching means_MDST()'s
# own total_var convention (sigma[Mg]^2 + sigma_tr^2).
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
# Tree is the study unit: each gets one Management + exactly one May and one
# July row (the actual repeated-measures design), plus its own offset
# tree_offset[Tr] shared by both its rows -- the non-independence MDS's own
# simulator deliberately didn't have.

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
