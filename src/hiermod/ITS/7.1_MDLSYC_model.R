# 7.1_MDLSYC_model.R

# Adds three control variables on top of Model 6: 
# deg_h_z, precip_72h_z and seq_depth_z

## Model definition ---------------------------------------------------

source('src/hiermod/ITS/6.1_MDLSY_model.R')

model$main_model <- quote(
  mu <- loga[Mg] + gamma*(Mo-1) + b[Lo]*sigma_loc + yr[Yr]*sigma_yr +
    tr[Tr]*sigma_tr + cv[Cv] + 
    # Here:
    b_deg*deg_h_z + b_precip*precip_72h_z + b_seq*seq_depth_z
)
# #all with standardized-scale slopes:
model$prior_deg    <- quote(b_deg    ~ dnorm(0,1))
model$prior_precip <- quote(b_precip ~ dnorm(0,1))
model$prior_seq    <- quote(b_seq    ~ dnorm(0,1))

# Sigma priors: Model 6's variance-budget-calibrated rates (scale_dexp_rate(3|2, 3, 4)
# in 6.2_MDLSY_validation.R), hardcoded here since K stays at 4 (no new summed
# variance term) and the calibrated model there only exists as a local
# `model_vbc` variable, not something to source from a validation script.
#     sigma  sigma_loc   sigma_tr   sigma_yr
#      3.46       2.31       2.31       2.31
model$pr_sigma     <- quote(sigma[cell] ~ dexp(3.46))
model$pr_sigma_loc <- quote(sigma_loc   ~ dexp(2.31))
model$pr_sigma_tr  <- quote(sigma_tr    ~ dexp(2.31))
model$pr_sigma_yr  <- quote(sigma_yr    ~ dexp(2.31))

# Backtransforming function --------------------------------------------------

# means_MDLSYC(): add covariate_offset term
# (deg_h_z/precip_72h_z/seq_depth_z all default to 0, i.e. this sample's own
# average weather and sequencing depth).
means_MDLSYC <- function(post, shift = 0, deg_h_z = 0, precip_72h_z = 0, seq_depth_z = 0){
  total_var <- post$sigma^2 + as.vector(post$sigma_loc)^2 +
    as.vector(post$sigma_tr)^2 + as.vector(post$sigma_yr)^2

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

# Add three covariates into "Explained" alongside Management/Season/Cultivar.
# kept as one lumped amount rather than split out per covariate, same
# reasoning as 6.1_MDLSY_model.R 

variance_partition_MDLSYC <- function(post, dat){
  N <- length(dat$Mg)

  loga_obs <- post$loga[, dat$Mg] # n_draws x N
  gamma    <- as.vector(post$s_conv) + outer(as.vector(post$gap_shift), dat$Mg - 1)
  gamma_term <- sweep(gamma, 2, dat$Mo - 1, "*")
  cv_obs   <- post$cv[, dat$Cv] # n_draws x N
  covariates <- outer(as.vector(post$b_deg), dat$deg_h_z) +
    outer(as.vector(post$b_precip), dat$precip_72h_z) +
    outer(as.vector(post$b_seq), dat$seq_depth_z)

  fixed_mu  <- loga_obs + gamma_term + cv_obs + covariates
  explained <- apply(fixed_mu, 1, var) # length n_draws

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

## Data-generating function --------------------------------------------

# Add deg_h_z/precip_72h_z/seq_depth_z, each  independent rnorm(0,1) draws.
#  NOT reproducing the real deg_h/precip_72h ~ Season/Year
# correlation (see 7.2_MDLSYC_validation.R); this tests whether b_deg/
# b_precip/b_seq are recoverable in principle, not how identifiable they are
# under the real design's collinearity.
sim_div_MDLSYC <- function(
    N_samples, loga, s_conv, gap_shift, sigma, b_deg, b_precip, b_seq, cv,
    sigma_loc = 0.5, sigma_tr = 0.3, sigma_yr = 0.3,
    n_loc, p_dropout = 0, shift = NULL){
  n_tree <- N_samples %/% 2
  n_yr   <- 3

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

  dat$deg_h_z      <- rnorm(nrow(dat))
  dat$precip_72h_z <- rnorm(nrow(dat))
  dat$seq_depth_z  <- rnorm(nrow(dat))

  loc_offset  <- rnorm(n_loc,  0, sigma_loc)
  tree_offset <- rnorm(n_tree, 0, sigma_tr)
  year_offset <- rnorm(n_yr,   0, sigma_yr)

  gamma <- s_conv + gap_shift*(dat$Mg - 1)
  mu <- loga[dat$Mg] + gamma*(dat$Mo - 1) +
    loc_offset[dat$Lo] + year_offset[dat$Yr] + tree_offset[dat$Tr] +
    cv[dat$Cv] + b_deg*dat$deg_h_z + b_precip*dat$precip_72h_z + b_seq*dat$seq_depth_z

  dat$Dv <- rlnorm(nrow(dat), meanlog = mu, sdlog = sigma[dat$cell])
  if (!is.null(shift)) dat$Dv_shifted <- shift + dat$Dv
  dat
}

# SBC contrast_fn -- same shape as contrast_may_gap_MDLSY() (6.1_MDLSY_model.R).
# b_deg/b_precip/b_seq/cv don't need to enter here: covariate_offset in
# means_MDLSYC() is identical for mean_1/mean_3 (both use the default
# deg_h_z=precip_72h_z=seq_depth_z=0 reference), and Cv is assigned
# independently of Mg in the simulator, so all cancel in the org-conv
# contrast specifically.
contrast_may_gap_MDLSYC <- function(post, true_params, means_fn){
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
simulate_from_priors <- function(true_params, N_samples = 250,
                                 n_loc = 4, p_dropout = 0.1, shift = NULL){
  sim_div_MDLSYC(
    N_samples = N_samples,
    loga = true_params$loga,
    s_conv = true_params$s_conv,
    gap_shift = true_params$gap_shift,
    sigma = true_params$sigma,
    b_deg = true_params$b_deg,
    b_precip = true_params$b_precip,
    b_seq = true_params$b_seq,
    cv = true_params$cv,
    sigma_loc = true_params$sigma_loc,
    sigma_tr = true_params$sigma_tr,
    sigma_yr = true_params$sigma_yr,
    n_loc = n_loc,
    p_dropout = p_dropout,
    shift = shift
  )
}
