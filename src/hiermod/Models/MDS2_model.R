# MODEL 3 (MDS2), 16S only: Management x Season (4-cell estimand: conv_May,
# conv_July, org_May, org_July), with Tree random effects (repeated
# measures) and Year as a fixed effect -- on top of Model 2 (MDL), but with
# Location dropped entirely rather than patched.
#
# Location A/C in the real data are single-management (100% Conventional /
# 100% Organic respectively; only B/D have both), so no amount of
# hierarchical adjustment on b[Lo] can cleanly separate a "Management
# effect" from "being at that particular Location" -- SBC on MDLS2/MDLS2vz
# confirmed this shows up as systematic loga miscalibration that a
# sum-to-zero constraint only partially addresses. Rather than continuing
# to engineer around Location's low cardinality, 16S's model family drops
# Location from here on: fit the full dataset without it (a transparent,
# explicitly-caveated-as-confounded estimate), and separately replicate on
# just Location B/D (the only locations with within-location Management
# variation) as the causally clean check -- restrict-and-compare instead of
# statistical adjustment. This is a divergence from ITS's own model lineage
# (which keeps Location throughout), not a naming coincidence: 16S's Model
# 3 (MDS2) is NOT the same lineage as ITS's Model 3 (MDLv).
#
# "2" keeps MDLS2's naming (4-cell mean structure: gamma <- s_conv +
# gap_shift*(Mg-1), i.e. an interactive Season effect, not just an
# additive one) without yet adopting MDLS2's OTHER defining feature --
# sigma varying by cell (Mg x Mo). sigma stays Management-only (sigma[Mg])
# for now, matching Model 2's simplicity; revisit the cell-level/log-linear
# variance structure only after this simpler version calibrates cleanly.
#
# Real data uses Dv = Hill_1 - 1 (dlnorm needs (0,Inf) support; Hill_1's
# natural floor is 1) -- same fix Model 5 (now archived) already needed.
# `shift = 1` in the simulator/means function converts back to the Hill_1
# scale for interpretable plots/contrasts, matching that same precedent.

model_MDS2_16S <- alist(
  likelihood = Dv ~ dlnorm(mu, sigma[Mg]),
  main_model = mu <- loga[Mg] + gamma*(Mo-1) + yr[Yr] + tr[Tr]*sigma_tr,
  gamma_def  = gamma <- s_conv + gap_shift*(Mg-1), # interactive gap

  prior_loga = loga[Mg]  ~ dnorm(5,2),
  prior_s    = s_conv    ~ dnorm(0,1),
  prior_gs   = gap_shift ~ dnorm(0,1),
  prior_yr   = yr[Yr]    ~ dnorm(0,1),
  prior_tr   = tr[Tr]    ~ dnorm(0,1),

  pr_sigma    = sigma[Mg] ~ dexp(2),
  pr_sigma_tr = sigma_tr  ~ dexp(2)
)

model_id <- "MDS2"

## Backtransform wrapper ------------------------------------------------------
# Same shape as means_MDLS() (MDLS_model.R) minus sigma_loc; shift arg
# added directly since 16S needs it from the start here (unlike MDLS(),
# which only grew a shift-aware sibling later at MDLS2()).

means_MDS2 <- function(post, shift = 0){
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

# SBC estimands -- may_gap/july_gap/seasonal_change, matching means_MDS2()'s
# own total_var convention (sigma[Mg]^2 + sigma_tr^2, no Location term).
dq_MDS2 <- SBC::derived_quantities(
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
# Tree is the study unit: each gets one Management + Year, and exactly one
# May + one July row. No Location, so no Lo x Yr dropout grid needed
# (that was purely to unbalance the Location design) -- Year is sampled
# directly per tree.

sim_div_MDS2 <- function(N_samples, loga, s_conv, gap_shift, sigma, year_offset,
                          sigma_tr = 0.3, shift = NULL){
  n_tree <- N_samples %/% 2 # 2 rows/tree (May + July)
  n_yr   <- length(year_offset)

  trees <- tibble(
    Tr = seq_len(n_tree),
    Yr = sample(n_yr, n_tree, replace = TRUE),
    Mg = rbern(n_tree) + 1
  )

  dat <- trees %>% crossing(Mo = 1:2) %>% arrange(Tr)

  tree_offset <- rnorm(n_tree, 0, sigma_tr)

  gamma <- s_conv + gap_shift*(dat$Mg - 1)
  mu <- loga[dat$Mg] + gamma*(dat$Mo - 1) + year_offset[dat$Yr] + tree_offset[dat$Tr]

  dat$Dv <- rlnorm(nrow(dat), meanlog = mu, sdlog = sigma[dat$Mg])
  if (!is.null(shift)) dat$Dv_shifted <- shift + dat$Dv
  dat
}

simulate_from_priors_MDS2 <- function(true_params, N_samples = 250,
                                       year_offset = c(0, 0.3, -0.2), shift = NULL){
  sim_div_MDS2(
    N_samples = N_samples,
    loga = true_params$loga,
    s_conv = true_params$s_conv,
    gap_shift = true_params$gap_shift,
    sigma = true_params$sigma,
    year_offset = year_offset,
    sigma_tr = true_params$sigma_tr,
    shift = shift
  )
}
