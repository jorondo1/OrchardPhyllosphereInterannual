# 4.1_MDLS_model.R -- MODEL 4 (Management effect by Season) definition:
# model alist, means_fn, data-generating function, the custom SBC contrast_fn
# (needed because means_MDLS()'s 4-column `mean` doesn't fit run_sbc()'s
# default 2-column contrast_from_means()), and the prior-simulator. See
# 4.2_MDLS_validation.R for the actual runs/fits/plots and the
# design-rationale comments (three estimands, the McElreath two-equation
# interaction, Tree random effects) -- that narrative stays there.

model <- alist(
  likelihood   = Dv ~ dlnorm(mu, sigma[Mg]), # Link function

  # Model:
  main_model = mu <- loga[Mg] + gamma*(Mo-1) + b[Lo]*sigma_loc + yr[Yr] + tr[Tr]*sigma_tr,
  gamma_def  = gamma <- s_conv + gap_shift*(Mg-1), # interactive gap
                # we use Mo-1 and Mg-1 because they are index variables [1,2] but
                # here are just building a linear combination of parameters
                # in order to get an intercept for each Mo-Mg combination

  # Old Priors:
  prior_loga = loga[Mg]   ~ dnorm(2,2),  # May-baseline median diversity per Management
  prior_b    = b[Lo]      ~ dnorm(0,1),  # Location-specific parameters (centered)

  # New priors:
  prior_s    = s_conv     ~ dnorm(0,1),  # Season main effect, at Conventional (Mg=1)
  prior_gs   = gap_shift  ~ dnorm(0,1),  # Interaction estimand (change in Mg gap from May to July)
  prior_yr   = yr[Yr]     ~ dnorm(0,1),  # Year-specific parameters (centered)
  prior_tr   = tr[Tr]     ~ dnorm(0,1),  # Tree-specific parameters

  # Hyperpriors:
  pr_sigma     = sigma[Mg] ~ dexp(3), # Management-specific noise
  pr_sigma_loc = sigma_loc ~ dexp(2), # Global location spread
  pr_sigma_tr  = sigma_tr  ~ dexp(2)  # Global tree spread
)

# Backtransform wrapper
# Specific to this model, needs to produce the target estimands
# i-e compute means for each managament-month pair
#
means_MDLS <- function(post){
  # total_var is shared across both months within a Management group; only
  # sigma[Mg] differs by month; sigma_loc/sigma_tr apply identically to every
  # cell. yr[Yr] deliberately NOT included -- no sigma_yr exists to marginalize
  # it properly (dnorm(0,1) is a fixed, unpooled prior, not a population
  # scale); "mean" here is for an average Location/Tree, reference Year.
  total_var_conv <- post$sigma[,1]^2 + as.vector(post$sigma_loc)^2 + as.vector(post$sigma_tr)^2
  total_var_org  <- post$sigma[,2]^2 + as.vector(post$sigma_loc)^2 + as.vector(post$sigma_tr)^2

  # Build the means for each MAnagement-Month pair:
  s_conv    <- as.vector(post$s_conv)
  gap_shift <- as.vector(post$gap_shift)

  mu_conv_May  <- post$loga[,1]
  mu_conv_July <- mu_conv_May + s_conv            # Mean seasonal shift at conventional
  mu_org_May   <- post$loga[,2]
  mu_org_July  <- mu_org_May + s_conv + gap_shift # Mean seasonal shift at organic

  # order: conv_May, conv_July, org_May, org_July -- fixed, since post_full()
  # overwrites these colnames with mean_1..mean_4 anyway (see below)
  list(
    mean = cbind(
      lognormal_mean(mu_conv_May,  total_var_conv),
      lognormal_mean(mu_conv_July, total_var_conv),
      lognormal_mean(mu_org_May,   total_var_org),
      lognormal_mean(mu_org_July,  total_var_org)
    ),
    median = cbind(
      lognormal_mean(mu_conv_May,  0),
      lognormal_mean(mu_conv_July, 0),
      lognormal_mean(mu_org_May, 0),
      lognormal_mean(mu_org_July,  0)
    )
  )
}

## Data-generating function ---------------------------------------------------

# Mirrors the model's own mu line directly: every fixed-effect argument
# below is given in the *same units the model would report it in* (loga is
# log-median, sigma is the model's own sigma[Mg], not a raw mean_/cv_ to
# backconvert like before; that round-trip is an identity that just cancels out, see
# sim_div_ML()/MLvary() in 2.1_MDL_model.R and 3.1_MDLv_model.R for the
# same fix). cv_to_sigma() (hiermod_core.R) is there for picking a
# human-readable CV by hand below.
#
# Only provide N and desired parameters, the rest gets built from within to
# ensure a good reflection of the actual design.
# Tree is the study unit: each gets one Location x Management combination,
# and exactly one May row + one July row.
#
# p_dropout randomly knocks out some Location x Year combinations before
# assigning trees, to unbalance the design.

sim_div_MDLS <- function(N_samples, loga, s_conv, gap_shift, sigma, year_offset,
                          sigma_loc = 0.5, sigma_tr = 0.3,
                          n_loc, p_dropout = 0){
  n_tree <- N_samples %/% 2       # 2 rows/tree (May + July)
  n_yr   <- length(year_offset)

  # Complete Location X year grid
  loc_yr_grid <- expand.grid(Lo = seq_len(n_loc), Yr = seq_len(n_yr)) %>%
    filter(runif(n()) > p_dropout)                    # drop some Lo x Yr cells
  stopifnot(nrow(loc_yr_grid) > 0)                     # unlucky draw wiped every cell

  # Design tibble;
  trees <- tibble(                                     # one row per Tree...
    Tr = seq_len(n_tree),
    cell = sample(nrow(loc_yr_grid), n_tree, replace = TRUE),
    Lo   = loc_yr_grid$Lo[cell],
    Yr   = loc_yr_grid$Yr[cell],
    Mg   = rbern(n_tree) + 1
  ) %>% dplyr::select(-cell)

  # expand to May+July pair (double rows):
  dat <- trees %>% crossing(Mo = 1:2) %>%
    arrange(Tr)

  # draw offsets:
  loc_offset  <- rnorm(n_loc,  0, sigma_loc)
  tree_offset <- rnorm(n_tree, 0, sigma_tr)

  # Rewrite the model equation to determine mu
  gamma <- s_conv + gap_shift*(dat$Mg - 1)
  mu <- loga[dat$Mg] + gamma*(dat$Mo - 1) +
    loc_offset[dat$Lo] + year_offset[dat$Yr] + tree_offset[dat$Tr]

  dat$Dv <- rlnorm(nrow(dat), meanlog = mu, sdlog = sigma[dat$Mg])
  dat
}

# SBC contrast_fn -- run_sbc()'s default contrast_fn (contrast_from_means())
# assumes `mean` is a 2-column matrix and tests column2-column1. means_MDLS()'s
# `mean` has 4 columns, so the default would silently compare
# conv_July-conv_May against an unrelated "true" Mg gap. Testing May gap
# here: structurally the same single Mg contrast models 1-3 already
# validate, and gap_shift's own recovery is separately checked by the
# Parameter recovery section's fixed_recovery table in the numbered script.
contrast_may_gap_MDLS <- function(post, true_params, means_fn){
  m <- means_fn(post)$mean               # mean_1=conv_May, mean_3=org_May
  post_contrast <- m[,3] - m[,1]

  # true side: same total_var convention as means_MDLS() itself (+ sigma_loc^2
  # + sigma_tr^2) -- otherwise the two sides aren't on matched footing even
  # once the right columns are compared.
  total_var_conv <- true_params$sigma[1]^2 + true_params$sigma_loc^2 + true_params$sigma_tr^2
  total_var_org  <- true_params$sigma[2]^2 + true_params$sigma_loc^2 + true_params$sigma_tr^2
  true_conv_May <- lognormal_mean(true_params$loga[1], total_var_conv)
  true_org_May  <- lognormal_mean(true_params$loga[2], total_var_org)

  list(post_contrast = post_contrast, true_contrast = true_org_May - true_conv_May)
}

# yr[Yr] is a fixed dnorm(0,1) prior, not an estimated population scale
# worth testing calibration of here.
simulate_from_priors <- function(true_params, N_samples = 250,
                                  year_offset = c(0, 0.3, -0.2),
                                  n_loc = 4, p_dropout = 0.1){
  sim_div_MDLS(
    N_samples = N_samples,
    loga = true_params$loga,
    s_conv = true_params$s_conv,
    gap_shift = true_params$gap_shift,
    sigma = true_params$sigma,
    year_offset = year_offset,
    sigma_loc = true_params$sigma_loc,
    sigma_tr = true_params$sigma_tr,
    n_loc = n_loc,
    p_dropout = p_dropout
  )
}
