# MDSTYCVr_model.R --- MODEL 8 (MDSTYCVr, "Gimli the Greedy"): MDSTYCV with Cultivar
# as a random effect (cv[Cv]*sigma_cv, non-centered, like Tree)

source('src/hiermod/0_INDEX.R')

# Why: the 5 cultivars are a choice among many -> arguably exchangeable
# Risk: Tree nested in Cultivar -> sigma_cv and sigma_tr compete for the same variance
# - can sigma_cv be identified from 5 levels? (pairs check in 8.2: cor(sigma_tr, sigma_cv))
# Priors: as MDSTYCV; cv[Cv] ~ dnorm(0,1), sigma_cv ~ dhalfnorm(0,1) (as Tree)
# Outcome: failed calibration (divergences, miscalibrated intercepts) -> not used

model_MDSTYCVr_16S <- alist(
  likelihood = Dv ~ dlnorm(mu, sigma[Mg]),
  main_model = mu <- loga[Mg] + gamma*(Mo-1) + yr_eff + cv[Cv]*sigma_cv + tr[Tr]*sigma_tr +
    b_deg*deg_h_z + b_precip*precip_72h_z + b_seq*seq_depth_z,
  gamma_def  = gamma <- s_conv + gap_shift*(Mg-1), # interactive Season effect

  prior_loga = loga[Mg]  ~ dnorm(5,2),
  prior_s    = s_conv    ~ dnorm(0,1),
  prior_gs   = gap_shift ~ dnorm(0,1),
  pr_sigma   = sigma[Mg] ~ dhalfnorm(0,1),

  # Year, sum-to-zero: yr1/yr2 free, yr3 = -(yr1+yr2)
  yr_eff_def = yr_eff <- yr1*(Yr==1) + yr2*(Yr==2) - (yr1+yr2)*(Yr==3),
  prior_yr1  = yr1 ~ dnorm(0,1),
  prior_yr2  = yr2 ~ dnorm(0,1),

  # Covariates (standardized)
  prior_deg    = b_deg    ~ dnorm(0,1),
  prior_precip = b_precip ~ dnorm(0,1),
  prior_seq    = b_seq    ~ dnorm(0,1),

  # Tree, non-centered random effect
  prior_tr    = tr[Tr]   ~ dnorm(0,1),
  pr_sigma_tr = sigma_tr ~ dhalfnorm(0,1),

  # Cultivar, non-centered random effect
  prior_cv    = cv[Cv]   ~ dnorm(0,1),
  pr_sigma_cv = sigma_cv ~ dhalfnorm(0,1)
)

attr(model_MDSTYCVr_16S, "name") <- "Gimli the Greedy"
model_id_MDSTYCVr <- "MDSTYCVr"

## ITS variant -----------------------------------------------------------------
# Intercept prior only (as MDS_ITS)
model_MDSTYCVr_ITS <- model_MDSTYCVr_16S
model_MDSTYCVr_ITS$prior_loga <- quote(loga[Mg] ~ dnorm(2,2))
attr(model_MDSTYCVr_ITS, "name") <- "Gimli the Greedy"

## means_MDSTYCVr() -------------------------------------------------------------
# As means_MDSTYCV(); mean variance also includes sigma_cv^2 (cultivar now random)

means_MDSTYCVr <- function(post, shift = 0, deg_h_z = 0, precip_72h_z = 0, seq_depth_z = 0){
  total_var_conv <- post$sigma[,1]^2 + as.vector(post$sigma_tr)^2 + as.vector(post$sigma_cv)^2
  total_var_org  <- post$sigma[,2]^2 + as.vector(post$sigma_tr)^2 + as.vector(post$sigma_cv)^2

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

# SBC estimands: as dq_MDSTYCV, variance incl. sigma_cv^2
dq_MDSTYCVr <- SBC::derived_quantities(
  may_gap =
    exp(loga[2] + (sigma[2]^2 + sigma_tr^2 + sigma_cv^2) / 2) -
    exp(loga[1] + (sigma[1]^2 + sigma_tr^2 + sigma_cv^2) / 2),
  july_gap =
    exp(loga[2] + s_conv + gap_shift + (sigma[2]^2 + sigma_tr^2 + sigma_cv^2) / 2) -
    exp(loga[1] + s_conv +             (sigma[1]^2 + sigma_tr^2 + sigma_cv^2) / 2),
  seasonal_change =
    (exp(loga[2] + s_conv + gap_shift + (sigma[2]^2 + sigma_tr^2 + sigma_cv^2) / 2) -
       exp(loga[1] + s_conv +             (sigma[1]^2 + sigma_tr^2 + sigma_cv^2) / 2)) -
    (exp(loga[2] + (sigma[2]^2 + sigma_tr^2 + sigma_cv^2) / 2) -
       exp(loga[1] + (sigma[1]^2 + sigma_tr^2 + sigma_cv^2) / 2))
)

## variance_partition_MDSTYCVr() -------------------------------------------------
# Sequential partition as MDSYCV; Cultivar and Tree as variance shares
# (sigma_cv^2, sigma_tr^2 over total)

variance_partition_MDSTYCVr <- function(post, dat){
  loga_obs   <- post$loga[, dat$Mg]
  gamma      <- as.vector(post$s_conv) + outer(as.vector(post$gap_shift), dat$Mg - 1)
  gamma_term <- sweep(gamma, 2, dat$Mo - 1, "*")
  yr3        <- -(as.vector(post$yr1) + as.vector(post$yr2))
  yr_obs     <- cbind(post$yr1, post$yr2, yr3)[, dat$Yr]

  covariates <- outer(as.vector(post$b_deg), dat$deg_h_z) +
    outer(as.vector(post$b_precip), dat$precip_72h_z) +
    outer(as.vector(post$b_seq), dat$seq_depth_z)

  fixed_MgMo    <- loga_obs + gamma_term
  fixed_MgMoY   <- fixed_MgMo + yr_obs
  fixed_MgMoYCo <- fixed_MgMoY + covariates

  var_MgMo  <- apply(fixed_MgMo,    1, var)
  var_Y     <- apply(fixed_MgMoY,   1, var) - var_MgMo
  var_Cov   <- apply(fixed_MgMoYCo, 1, var) - apply(fixed_MgMoY, 1, var)
  explained <- apply(fixed_MgMoYCo, 1, var)

  sigma_tr_sq <- as.vector(post$sigma_tr)^2
  sigma_cv_sq <- as.vector(post$sigma_cv)^2

  mg_n <- as.integer(table(factor(dat$Mg, levels = 1:2)))
  residual_var <- as.vector((post$sigma^2) %*% (mg_n / sum(mg_n)))

  total <- explained + sigma_tr_sq + sigma_cv_sq + residual_var

  bind_rows(
    tibble(statistic = "Variance partition", group = "Management x Season", value = var_MgMo / total),
    tibble(statistic = "Variance partition", group = "Year",                 value = var_Y / total),
    tibble(statistic = "Variance partition", group = "Covariates",           value = var_Cov / total),
    tibble(statistic = "Variance partition", group = "Cultivar",             value = sigma_cv_sq / total),
    tibble(statistic = "Variance partition", group = "Tree",                 value = sigma_tr_sq / total),
    tibble(statistic = "Variance partition", group = "Residual",             value = residual_var / total)
  )
}

## Data-generating function ---------------------------------------------------
# As sim_div_MDSTYCV(), cultivar offsets drawn from N(0, sigma_cv)

sim_div_MDSTYCVr <- function(N_samples, loga, s_conv, gap_shift, sigma, yr1, yr2,
                              sigma_cv, sigma_tr, b_deg, b_precip, b_seq, shift = NULL){
  n_tree <- N_samples %/% 2 # 2 rows/tree (May + July)
  yr_vec <- c(yr1, yr2, -(yr1 + yr2))

  trees <- tibble(
    Tr = seq_len(n_tree), Mg = rbern(n_tree) + 1,
    Yr = sample(3, n_tree, replace = TRUE), Cv = sample(5, n_tree, replace = TRUE))
  dat <- trees %>% crossing(Mo = 1:2) %>% arrange(Tr)

  tree_offset <- rnorm(n_tree, 0, sigma_tr)
  cv_offset   <- rnorm(5, 0, sigma_cv) # 5 Cultivar levels, random effect

  gamma <- s_conv + gap_shift*(dat$Mg - 1)
  mu_structural <- loga[dat$Mg] + gamma*(dat$Mo - 1) + yr_vec[dat$Yr] +
    cv_offset[dat$Cv] + tree_offset[dat$Tr]

  dat$deg_h_z      <- rnorm(nrow(dat))
  dat$precip_72h_z <- rnorm(nrow(dat))
  dat$seq_depth_z  <- rnorm(nrow(dat))

  mu <- mu_structural + b_deg*dat$deg_h_z + b_precip*dat$precip_72h_z + b_seq*dat$seq_depth_z

  dat$Dv <- rlnorm(nrow(dat), meanlog = mu, sdlog = sigma[dat$Mg])
  if (!is.null(shift)) dat$Dv_shifted <- shift + dat$Dv
  dat
}

# Simulate one dataset from one prior draw (draw_true() output)
simulate_from_priors_MDSTYCVr <- function(true_params, N_samples = 250, shift = NULL){
  sim_div_MDSTYCVr(
    N_samples = N_samples,
    loga = true_params$loga, s_conv = true_params$s_conv, gap_shift = true_params$gap_shift,
    sigma = true_params$sigma, yr1 = true_params$yr1, yr2 = true_params$yr2,
    sigma_cv = true_params$sigma_cv, sigma_tr = true_params$sigma_tr,
    b_deg = true_params$b_deg, b_precip = true_params$b_precip, b_seq = true_params$b_seq,
    shift = shift
  )
}
