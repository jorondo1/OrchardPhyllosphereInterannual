# MODEL 6 (MDLSY), ITS only: pools yr[Yr] into a non-centered random effect
# (was a fixed dnorm(0,1) effect) and adds Cultivar (cv[Cv]) as a fixed
# effect, on top of Model 5.

source('src/hiermod/Models/MDLS2_model.R')

model_MDLSY_ITS <- model_MDLS2_ITS
model_MDLSY_ITS$main_model <- quote(
  mu <- loga[Mg] + gamma*(Mo-1) + b[Lo]*sigma_loc + yr[Yr]*sigma_yr +
    tr[Tr]*sigma_tr + cv[Cv]
)
model_MDLSY_ITS$prior_cv    <- quote(cv[Cv]   ~ dnorm(0,1)) # fixed/unpooled, same as loga[Mg]
model_MDLSY_ITS$pr_sigma_yr <- quote(sigma_yr ~ dexp(2))    # same rate as sigma_loc/sigma_tr

model_id <- "MDLSY"

# Backtransforming function --------------------------------------------------

# means_MDLSY(): total_var gains sigma_yr^2 -- Year is now a proper
# marginalized-over population (same logic as Location/Tree, see
# MODEL_HISTORY.md). cv[Cv] deliberately left OUT of mu here -- no sigma_cv
# to marginalize over (unpooled fixed effect), same precedent yr[Yr] set
# pre-Model-6 -- reported instead via its own posterior panel in
# 6.4_MDLSY_ITS_analysis.R.
means_MDLSY <- function(post, shift = 0){
  total_var <- post$sigma^2 + as.vector(post$sigma_loc)^2 +
    as.vector(post$sigma_tr)^2 + as.vector(post$sigma_yr)^2

  s_conv    <- as.vector(post$s_conv)
  gap_shift <- as.vector(post$gap_shift)

  mu_conv_May  <- post$loga[,1]
  mu_conv_July <- mu_conv_May + s_conv
  mu_org_May   <- post$loga[,2]
  mu_org_July  <- mu_org_May + s_conv + gap_shift

  list(
    mean = cbind(
      lognormal_mean(mu_conv_May,  total_var[,1], shift = shift),
      lognormal_mean(mu_conv_July, total_var[,2], shift = shift),
      lognormal_mean(mu_org_May,   total_var[,3], shift = shift),
      lognormal_mean(mu_org_July,  total_var[,4], shift = shift)
    ),
    median = cbind(
      lognormal_mean(mu_conv_May,  0, shift = shift),
      lognormal_mean(mu_conv_July, 0, shift = shift),
      lognormal_mean(mu_org_May,   0, shift = shift),
      lognormal_mean(mu_org_July,  0, shift = shift)
    )
  )
}

# Variance partition / Bayesian R2 -----------------------------------------

# For each posterior draw: var() across the *actual observations* of the
# fixed-effect part of mu (Management+Season+Cultivar, at each row's real
# covariate values) vs. each random effect's population variance vs. the
# (cell-membership-weighted) residual variance -- five components, all
# sharing one total, so each is a genuine share of it (sums to 1 per draw).
# "Explained" here is exactly Bayesian R2 (var(fixed)/var(total)); the other
# four are the variance-partition coefficients (VPC) for what's left.
variance_partition_MDLSY <- function(post, dat){
  N <- length(dat$Mg)

  loga_obs <- post$loga[, dat$Mg]                                   # n_draws x N
  gamma    <- as.vector(post$s_conv) + outer(as.vector(post$gap_shift), dat$Mg - 1)
  gamma_term <- sweep(gamma, 2, dat$Mo - 1, "*")
  cv_obs   <- post$cv[, dat$Cv]                                     # n_draws x N

  fixed_mu  <- loga_obs + gamma_term + cv_obs
  explained <- apply(fixed_mu, 1, var)                              # length n_draws

  cell_n <- as.integer(table(factor(dat$cell, levels = 1:4)))
  residual_var <- as.vector((post$sigma^2) %*% (cell_n / sum(cell_n)))

  sigma_loc_sq <- as.vector(post$sigma_loc)^2
  sigma_tr_sq  <- as.vector(post$sigma_tr)^2
  sigma_yr_sq  <- as.vector(post$sigma_yr)^2

  total <- explained + sigma_loc_sq + sigma_tr_sq + sigma_yr_sq + residual_var

  bind_rows(
    tibble(statistic = "Variance partition", group = "Explained (fixed effects)", value = explained / total),
    tibble(statistic = "Variance partition", group = "Location",                  value = sigma_loc_sq / total),
    tibble(statistic = "Variance partition", group = "Tree",                      value = sigma_tr_sq / total),
    tibble(statistic = "Variance partition", group = "Year",                      value = sigma_yr_sq / total),
    tibble(statistic = "Variance partition", group = "Residual",                  value = residual_var / total)
  )
}

## Data-generating function ---------------------------------------------------
# Same skeleton as sim_div_MDLS2() (MDLS2_model.R), plus: year_offset is
# now drawn internally from sigma_yr (like loc_offset/tree_offset), not a
# caller-supplied fixed vector -- required for SBC to test what it's
# actually supposed to (recovering sigma_yr, not an arbitrary fixed vector).
# Cultivar (Cv) assigned per Tree, cv[] used as a direct true-value lookup
# (unpooled, like loga).
sim_div_MDLSY <- function(
    N_samples, loga, s_conv, gap_shift, sigma, cv,
    sigma_loc = 0.5, sigma_tr = 0.3, sigma_yr = 0.3,
    n_loc, p_dropout = 0, shift = NULL){
  n_tree <- N_samples %/% 2
  n_yr   <- 3   # fixed at the real design's Year count, not a free arg --
                # pooling is only meaningfully tested against the real number
                # of Year levels (weak regularizing power by design).

  loc_yr_grid <- expand.grid(Lo = seq_len(n_loc), Yr = seq_len(n_yr)) %>%
    filter(runif(n()) > p_dropout)
  stopifnot(nrow(loc_yr_grid) > 0)

  trees <- tibble(
    Tr = seq_len(n_tree),
    cell = sample(nrow(loc_yr_grid), n_tree, replace = TRUE),
    Lo   = loc_yr_grid$Lo[cell],
    Yr   = loc_yr_grid$Yr[cell],
    Mg   = rbern(n_tree) + 1,
    Cv   = sample(seq_along(cv), n_tree, replace = TRUE)
  ) %>% dplyr::select(-cell)

  dat <- trees %>% crossing(Mo = 1:2) %>% arrange(Tr)
  dat$cell <- (dat$Mg - 1) * 2 + dat$Mo

  loc_offset  <- rnorm(n_loc,  0, sigma_loc)
  tree_offset <- rnorm(n_tree, 0, sigma_tr)
  year_offset <- rnorm(n_yr,   0, sigma_yr)

  gamma <- s_conv + gap_shift*(dat$Mg - 1)
  mu <- loga[dat$Mg] + gamma*(dat$Mo - 1) +
    loc_offset[dat$Lo] + year_offset[dat$Yr] + tree_offset[dat$Tr] + cv[dat$Cv]

  dat$Dv <- rlnorm(nrow(dat), meanlog = mu, sdlog = sigma[dat$cell])
  if (!is.null(shift)) dat$Dv_shifted <- shift + dat$Dv
  dat
}

# SBC contrast_fn -- same shape as contrast_may_gap_MDLS2() (MDLS2_model.R).
# cv doesn't need to enter here: Cv is assigned independently of Mg in the
# simulator, so it cancels in the org-conv contrast specifically.
contrast_may_gap_MDLSY <- function(post, true_params, means_fn){
  m <- means_fn(post)$mean
  post_contrast <- m[,3] - m[,1]

  total_var_conv_May <- true_params$sigma[1]^2 + true_params$sigma_loc^2 +
    true_params$sigma_tr^2 + true_params$sigma_yr^2
  total_var_org_May  <- true_params$sigma[3]^2 + true_params$sigma_loc^2 +
    true_params$sigma_tr^2 + true_params$sigma_yr^2
  true_conv_May <- lognormal_mean(true_params$loga[1], total_var_conv_May)
  true_org_May  <- lognormal_mean(true_params$loga[2], total_var_org_May)

  list(post_contrast = post_contrast, true_contrast = true_org_May - true_conv_May)
}

# Prior simulator (SBC/prior-predictive glue).
simulate_from_priors_MDLSY <- function(true_params, N_samples = 250,
                                 n_loc = 4, p_dropout = 0.1, shift = NULL){
  sim_div_MDLSY(
    N_samples = N_samples,
    loga = true_params$loga,
    s_conv = true_params$s_conv,
    gap_shift = true_params$gap_shift,
    sigma = true_params$sigma,
    cv = true_params$cv,
    sigma_loc = true_params$sigma_loc,
    sigma_tr = true_params$sigma_tr,
    sigma_yr = true_params$sigma_yr,
    n_loc = n_loc,
    p_dropout = p_dropout,
    shift = shift
  )
}
