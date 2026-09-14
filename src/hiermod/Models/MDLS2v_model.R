# MODEL 5v (MDLS2v), 16S only: alternative to Model 5 (MDLS2) that
# reparametrizes sigma[cell] (4 independent dexp-distributed values) as a
# log-linear function of Management/Season main effects + interaction,
# exactly mirroring how the mean model already avoids one free parameter
# per cell (loga[Mg] + gamma*(Mo-1)). Same 4 degrees of freedom as today's
# sigma[cell], but on the unconstrained log scale with normal priors --
# no boundary-at-zero mode the way dexp() has, which is what let a handful
# of posterior draws collapse a cell's sigma toward zero and blow up the
# likelihood (see Model 5's divergence/PSIS investigation).

source('src/hiermod/Models/MDLS2_model.R') # sim_div_MDLS2(), contrast_may_gap_MDLS2()

model_MDLS2v_16S <- alist(
  likelihood   = Dv ~ dlnorm(mu, sigma),   # no [cell] index -- sigma is per-observation
  main_model   = mu <- loga[Mg] + gamma*(Mo-1) + b[Lo]*sigma_loc + yr[Yr] + tr[Tr]*sigma_tr,
  gamma_def    = gamma <- s_conv + gap_shift*(Mg-1),
  # Sigma now varies by Mg-Mo combination:
  sigma_def    = sigma <- exp(ls0 + ls_Mg*(Mg-1) + ls_Mo*(Mo-1) + ls_MgMo*(Mg-1)*(Mo-1)),

  prior_loga = loga[Mg]   ~ dnorm(5,2),    # unchanged from model_MDLS2_16S
  prior_b    = b[Lo]      ~ dnorm(0,1),
  prior_s    = s_conv     ~ dnorm(0,1),
  prior_gs   = gap_shift  ~ dnorm(0,1),
  prior_yr   = yr[Yr]     ~ dnorm(0,1),
  prior_tr   = tr[Tr]     ~ dnorm(0,1),

  prior_ls0     = ls0      ~ dnorm(0,0.5),    # starting guesses; calibration tunes these
  prior_ls_Mg   = ls_Mg    ~ dnorm(0,0.5),
  prior_ls_Mo   = ls_Mo    ~ dnorm(0,0.5),
  prior_ls_MgMo = ls_MgMo  ~ dnorm(0,0.5),

  pr_sigma_loc = sigma_loc ~ dexp(2),     # tighten
  pr_sigma_tr  = sigma_tr  ~ dexp(2)        # tighten more 
)

model_id <- "MDLS2v_shifted"

# Backtransforming function --------------------------------------------------

# Shared by means_MDLS2v()/true_sigma_from_ls(): the log-linear cell-sigma
# formula, taking the 4 coefficients as explicit arguments rather than
# reaching into the caller's scope by name -- that's what makes it safe to
# reuse for both scalar true values and vectorized posterior draws, since
# R's arithmetic vectorizes automatically either way.
sigma_cell <- function(Mg, Mo, ls0, ls_Mg, ls_Mo, ls_MgMo){
  exp(ls0 + ls_Mg*(Mg-1) + ls_Mo*(Mo-1) + ls_MgMo*(Mg-1)*(Mo-1))
}

# Same math as means_MDLS2(), but the 4 cell-level sigmas are first derived
# from ls0/ls_Mg/ls_Mo/ls_MgMo (evaluated at each of the 4 (Mg,Mo)
# combinations) rather than read off a post$sigma[,1:4] group-indexed param.
# Cell order matches sim_div_MDLS2()'s own convention: 1=conv_May,
# 2=conv_July, 3=org_May, 4=org_July.
means_MDLS2v <- function(post, shift = 0){
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

  total_var <- sigma_sq + as.vector(post$sigma_loc)^2 + as.vector(post$sigma_tr)^2

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

# Shared by simulate_from_priors_MDLS2v()/contrast_may_gap_MDLS2v(): derives
# the length-4 true cell-level sigma vector (cell order matching
# sim_div_MDLS2()'s convention: 1=conv_May, 2=conv_July, 3=org_May,
# 4=org_July) from a true_params list holding ls0/ls_Mg/ls_Mo/ls_MgMo.
true_sigma_from_ls <- function(true_params){
  with(true_params, c(
    sigma_cell(1,1, ls0, ls_Mg, ls_Mo, ls_MgMo),
    sigma_cell(1,2, ls0, ls_Mg, ls_Mo, ls_MgMo),
    sigma_cell(2,1, ls0, ls_Mg, ls_Mo, ls_MgMo),
    sigma_cell(2,2, ls0, ls_Mg, ls_Mo, ls_MgMo)
  ))
}

# Derives the 4 true cell-level sigmas from ls0/ls_Mg/ls_Mo/ls_MgMo, then
# delegates to sim_div_MDLS2() -- it only needs a length-4 sigma vector and
# doesn't care whether those 4 values came from independent parameters or
# from evaluating this log-linear formula.
simulate_from_priors_MDLS2v <- function(true_params, N_samples = 250,
                                         year_offset = c(0, 0.3, -0.2),
                                         n_loc = 4, p_dropout = 0.1, shift = NULL){
  true_sigma <- true_sigma_from_ls(true_params)

  sim_div_MDLS2(
    N_samples = N_samples,
    loga = true_params$loga,
    s_conv = true_params$s_conv,
    gap_shift = true_params$gap_shift,
    sigma = true_sigma,
    year_offset = year_offset,
    sigma_loc = true_params$sigma_loc,
    sigma_tr = true_params$sigma_tr,
    n_loc = n_loc,
    p_dropout = p_dropout,
    shift = shift
  )
}

# SBC contrast_fn -- run_sbc() draws true_params straight from the prior,
# which for this model means ls0/ls_Mg/ls_Mo/ls_MgMo, not a sigma[1:4]
# vector. Derive that vector first, then delegate to
# contrast_may_gap_MDLS2() (MDLS2_model.R) unchanged.
contrast_may_gap_MDLS2v <- function(post, true_params, means_fn){
  true_params$sigma <- true_sigma_from_ls(true_params)
  contrast_may_gap_MDLS2(post, true_params, means_fn)
}
