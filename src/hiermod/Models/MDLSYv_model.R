# MODEL 6v (MDLSYv), 16S only: pools yr[Yr] into a non-centered random
# effect and adds Cultivar (cv[Cv]), same jump ITS made at its own Model 6
# (MDLSY) -- but built on MDLS2v's structured log-linear residual variance
# instead of MDLS2's independent sigma[cell]. Covariates (ITS Model 7)
# are a separate, later step, not this one.

source('src/hiermod/Models/MDLS2v_model.R') # model_MDLS2v_16S, sigma_cell(), true_sigma_from_ls()
source('src/hiermod/Models/MDLSY_model.R')  # sim_div_MDLSY(), contrast_may_gap_MDLSY() -- marker-agnostic, reused unchanged

model_MDLSYv_16S <- model_MDLS2v_16S
model_MDLSYv_16S$main_model <- quote(
  mu <- loga[Mg] + gamma*(Mo-1) + b[Lo]*sigma_loc + yr[Yr]*sigma_yr +
    tr[Tr]*sigma_tr + cv[Cv]
)
model_MDLSYv_16S$prior_cv    <- quote(cv[Cv]   ~ dnorm(0,1))  # fixed/unpooled, same as loga[Mg]
model_MDLSYv_16S$pr_sigma_yr <- quote(sigma_yr ~ dexp(2))     # starting guess -- VBC-calibrated in calibration script

model_id <- "MDLSYv"

# Backtransforming function --------------------------------------------------

# Same math as means_MDLS2v(), plus means_MDLSY()'s sigma_yr^2 term in
# total_var. cv[Cv] stays out of mu here, same precedent as means_MDLSY().
means_MDLSYv <- function(post, shift = 0){
  ls0     <- as.vector(post$ls0)
  ls_Mg   <- as.vector(post$ls_Mg)
  ls_Mo   <- as.vector(post$ls_Mo)
  ls_MgMo <- as.vector(post$ls_MgMo)

  sigma_sq <- cbind(
    sigma_cell(1,1, ls0, ls_Mg, ls_Mo, ls_MgMo),
    sigma_cell(1,2, ls0, ls_Mg, ls_Mo, ls_MgMo),
    sigma_cell(2,1, ls0, ls_Mg, ls_Mo, ls_MgMo),
    sigma_cell(2,2, ls0, ls_Mg, ls_Mo, ls_MgMo)
  )^2

  total_var <- sigma_sq + as.vector(post$sigma_loc)^2 + as.vector(post$sigma_tr)^2 +
    as.vector(post$sigma_yr)^2

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
# Same structure as variance_partition_MDLSY(), but residual_var comes from
# sigma_cell()-evaluated per-cell sigma instead of a stored post$sigma
# matrix -- same adaptation means_MDLSYv() makes relative to means_MDLSY().
# dat$cell isn't assumed to exist (MDLSYv's own fit script doesn't build
# one, since the likelihood indexes by Mg/Mo directly) -- derived locally.
variance_partition_MDLSYv <- function(post, dat){
  cell <- (dat$Mg - 1) * 2 + dat$Mo

  loga_obs   <- post$loga[, dat$Mg]
  gamma      <- as.vector(post$s_conv) + outer(as.vector(post$gap_shift), dat$Mg - 1)
  gamma_term <- sweep(gamma, 2, dat$Mo - 1, "*")
  cv_obs     <- post$cv[, dat$Cv]

  fixed_mu  <- loga_obs + gamma_term + cv_obs
  explained <- apply(fixed_mu, 1, var)

  ls0     <- as.vector(post$ls0)
  ls_Mg   <- as.vector(post$ls_Mg)
  ls_Mo   <- as.vector(post$ls_Mo)
  ls_MgMo <- as.vector(post$ls_MgMo)
  sigma_sq <- cbind(
    sigma_cell(1,1, ls0, ls_Mg, ls_Mo, ls_MgMo),
    sigma_cell(1,2, ls0, ls_Mg, ls_Mo, ls_MgMo),
    sigma_cell(2,1, ls0, ls_Mg, ls_Mo, ls_MgMo),
    sigma_cell(2,2, ls0, ls_Mg, ls_Mo, ls_MgMo)
  )^2

  cell_n <- as.integer(table(factor(cell, levels = 1:4)))
  residual_var <- as.vector(sigma_sq %*% (cell_n / sum(cell_n)))

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

# Derives the 4 true cell-level sigmas from ls0/ls_Mg/ls_Mo/ls_MgMo, then
# delegates to sim_div_MDLSY() unchanged -- it only needs a length-4 sigma
# vector, doesn't care where it came from.
simulate_from_priors_MDLSYv <- function(true_params, N_samples = 250,
                                         n_loc = 4, p_dropout = 0.1, shift = NULL){
  true_sigma <- true_sigma_from_ls(true_params)

  sim_div_MDLSY(
    N_samples = N_samples,
    loga = true_params$loga,
    s_conv = true_params$s_conv,
    gap_shift = true_params$gap_shift,
    sigma = true_sigma,
    cv = true_params$cv,
    sigma_loc = true_params$sigma_loc,
    sigma_tr = true_params$sigma_tr,
    sigma_yr = true_params$sigma_yr,
    n_loc = n_loc,
    p_dropout = p_dropout,
    shift = shift
  )
}

# SBC contrast_fn -- derive the sigma vector, delegate to
# contrast_may_gap_MDLSY() (MDLSY_model.R) unchanged.
contrast_may_gap_MDLSYv <- function(post, true_params, means_fn){
  true_params$sigma <- true_sigma_from_ls(true_params)
  contrast_may_gap_MDLSY(post, true_params, means_fn)
}

# SBC estimands -- may_gap/july_gap/seasonal_change, same three contrasts
# as MDLS2v's own SBC, with sigma_yr^2 added into every total_var term
# (matching means_MDLSYv()'s own total_var). Cell sigmas inlined via the
# log-linear ls0/ls_Mg/ls_Mo/ls_MgMo formula (derived_quantities() formulas
# must be self-contained -- no calling sigma_cell()/true_sigma_from_ls()
# from inside).
dq_MDLSYv <- SBC::derived_quantities(
  may_gap =
    exp(loga[2] + (exp(ls0 + ls_Mg)^2 + sigma_loc^2 + sigma_tr^2 + sigma_yr^2) / 2) -
    exp(loga[1] + (exp(ls0)^2         + sigma_loc^2 + sigma_tr^2 + sigma_yr^2) / 2),
  july_gap =
    exp(loga[2] + s_conv + gap_shift + (exp(ls0 + ls_Mg + ls_Mo + ls_MgMo)^2 + sigma_loc^2 + sigma_tr^2 + sigma_yr^2) / 2) -
    exp(loga[1] + s_conv +             (exp(ls0 + ls_Mo)^2                 + sigma_loc^2 + sigma_tr^2 + sigma_yr^2) / 2),
  seasonal_change =
    (exp(loga[2] + s_conv + gap_shift + (exp(ls0 + ls_Mg + ls_Mo + ls_MgMo)^2 + sigma_loc^2 + sigma_tr^2 + sigma_yr^2) / 2) -
       exp(loga[1] + s_conv +             (exp(ls0 + ls_Mo)^2                 + sigma_loc^2 + sigma_tr^2 + sigma_yr^2) / 2)) -
    (exp(loga[2] + (exp(ls0 + ls_Mg)^2 + sigma_loc^2 + sigma_tr^2 + sigma_yr^2) / 2) -
       exp(loga[1] + (exp(ls0)^2         + sigma_loc^2 + sigma_tr^2 + sigma_yr^2) / 2))
)
