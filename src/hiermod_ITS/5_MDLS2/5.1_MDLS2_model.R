# 5.1_MDLS2_model.R -- MODEL 5: sigma[Mg] -> sigma[cell], residual SD
# varies by Management AND Season (4 cells). Why: see MODEL_HISTORY.md.
# cell = (Mg-1)*2 + Mo (1-based, 1..4), precomputed in every data list --
# order (1=Conv-May, 2=Conv-July, 3=Org-May, 4=Org-July) already matches
# means_MDLS()'s own mean/median column order, no downstream remapping.

## Model definition -------------------------------------------------

source('src/hiermod_ITS/4_MDLS/4.1_MDLS_model.R')

# Only the likelihood's scale index and its prior change; sourcing model 4's
# alist and patching these two named elements directly avoids re-declaring
# everything else (main_model/gamma_def/priors are all unchanged). Elements
# must exist and be named in the base alist for this to actually replace
# rather than silently append -- true here (4.1_MDLS_model.R names every
# element), but worth remembering if either model's alist changes shape.

model$likelihood <- quote(Dv ~ dlnorm(mu, sigma[cell]))   # was sigma[Mg]
model$pr_sigma   <- quote(sigma[cell] ~ dexp(3))           # was sigma[Mg] ~ dexp(3);
# TODO: revisit via prior-predictive check like model 4

# Backtransforming function --------------------------------------------------

# means_MDLS2(): post$sigma now comes back as an n_draws x 4 matrix (cell
# order = conv_May, conv_July, org_May, org_July), so total_var collapses
# from two vectors (means_MDLS()'s total_var_conv/total_var_org) to one
# matrix, one column per cell.

# Also we add the shift parameter so this is defined for model 6 (MDLSY) already
means_MDLS2 <- function(post, shift = 0){
  total_var <- post$sigma^2 + as.vector(post$sigma_loc)^2 + as.vector(post$sigma_tr)^2
  s_conv    <- as.vector(post$s_conv)
  gap_shift <- as.vector(post$gap_shift)

  mu_conv_May  <- post$loga[,1]
  mu_conv_July <- mu_conv_May + s_conv
  mu_org_May   <- post$loga[,2]
  mu_org_July  <- mu_org_May + s_conv + gap_shift

  list(
    mean = cbind( # HERE >>>>
      lognormal_mean(mu_conv_May,  total_var[,1], shift = shift),
      lognormal_mean(mu_conv_July, total_var[,2], shift = shift),
      lognormal_mean(mu_org_May,   total_var[,3], shift = shift),
      lognormal_mean(mu_org_July,  total_var[,4], shift = shift)
    ),
    median = cbind(
      lognormal_mean(mu_conv_May, 0, shift = shift),
      lognormal_mean(mu_conv_July, 0, shift = shift),
      lognormal_mean(mu_org_May, 0, shift = shift),
      lognormal_mean(mu_org_July, 0, shift = shift)
    )
  )
}


## Data-generating function ---------------------------------------------------
# Same design-generation logic as sim_div_MDLS() (4.1_MDLS_model.R) -- Tree
# owns a fixed Location/Management and gets exactly one May + one July row.
# Two changes from model 4: a `cell` column (matching the model's own
# sigma[cell] indexing), and `sigma` taken as a length-4 vector (cell order)
# instead of length-2.

sim_div_MDLS2 <- function(
    N_samples, loga, s_conv, gap_shift, sigma, year_offset,
    sigma_loc = 0.5, sigma_tr = 0.3,
    n_loc, p_dropout = 0,
    shift = NULL) # Shifting parameter, also used by model 6 (MDLSY)
  {
  n_tree <- N_samples %/% 2
  n_yr   <- length(year_offset)

  loc_yr_grid <- expand.grid(Lo = seq_len(n_loc), Yr = seq_len(n_yr)) %>%
    filter(runif(n()) > p_dropout)
  stopifnot(nrow(loc_yr_grid) > 0)

  trees <- tibble(
    Tr = seq_len(n_tree),
    cell = sample(nrow(loc_yr_grid), n_tree, replace = TRUE),
    Lo   = loc_yr_grid$Lo[cell],
    Yr   = loc_yr_grid$Yr[cell],
    Mg   = rbern(n_tree) + 1
  ) %>% dplyr::select(-cell)

  dat <- trees %>% crossing(Mo = 1:2) %>% arrange(Tr)
  dat$cell <- (dat$Mg - 1) * 2 + dat$Mo   # 1=Conv-May, 2=Conv-July, 3=Org-May, 4=Org-July

  loc_offset  <- rnorm(n_loc,  0, sigma_loc)
  tree_offset <- rnorm(n_tree, 0, sigma_tr)

  gamma <- s_conv + gap_shift*(dat$Mg - 1)
  mu <- loga[dat$Mg] + gamma*(dat$Mo - 1) +
    loc_offset[dat$Lo] + year_offset[dat$Yr] + tree_offset[dat$Tr]

  dat$Dv <- rlnorm(nrow(dat), meanlog = mu, sdlog = sigma[dat$cell])

  if(!is.null(shift))  # >>> SHIFTED VALUES on hill scale[1,Inf]
  {dat$Dv_shifted <- shift + dat$Dv}

  dat
}

# SBC contrast_fn -- same shape as contrast_may_gap_MDLS() (4.1_MDLS_model.R),
# testing May gap specifically. true_params$sigma is now length 4 (cell
# order), so the true side pulls indices [1] (conv_May) and [3] (org_May)
# instead of model 4's [1]/[2].
contrast_may_gap_MDLS2 <- function(post, true_params, means_fn){
  m <- means_fn(post)$mean
  post_contrast <- m[,3] - m[,1]

  total_var_conv_May <- true_params$sigma[1]^2 + true_params$sigma_loc^2 + true_params$sigma_tr^2
  total_var_org_May  <- true_params$sigma[3]^2 + true_params$sigma_loc^2 + true_params$sigma_tr^2
  true_conv_May <- lognormal_mean(true_params$loga[1], total_var_conv_May)
  true_org_May  <- lognormal_mean(true_params$loga[2], total_var_org_May)

  list(post_contrast = post_contrast, true_contrast = true_org_May - true_conv_May)
}

# Prior simulator (SBC/prior-predictive glue) -- identical to model 4's
# except the inner call target (sim_div_MDLS2() instead of sim_div_MDLS()).
# shift = NULL (default): plain model 5, no Dv_shifted column, matches
# sim_div_MDLS2()'s own default -- pass shift = 1 (e.g. at the run_sbc()/
# prior-predictive call site) to use this for the shifted variant instead;
# no separate model file needed for that, since the model alist itself
# never changes (see means_MDLS2()/sim_div_MDLS2()'s own shift= argument
# and the "Model specification" note in 5b.2_MDLS2_shifted_validation.R).
simulate_from_priors <- function(true_params, N_samples = 250,
                                 year_offset = c(0, 0.3, -0.2),
                                 n_loc = 4, p_dropout = 0.1, shift = NULL){
  sim_div_MDLS2(
    N_samples = N_samples,
    loga = true_params$loga,
    s_conv = true_params$s_conv,
    gap_shift = true_params$gap_shift,
    sigma = true_params$sigma,
    year_offset = year_offset,
    sigma_loc = true_params$sigma_loc,
    sigma_tr = true_params$sigma_tr,
    n_loc = n_loc,
    p_dropout = p_dropout,
    shift = shift
  )
}
