# MDSYz_model.R --- MODEL 4 (MDSYz), 16S: MDSY + a sum-to-zero constraint on
# yr[Yr]. SBC on MDSY found loga[1]/loga[2] severely miscalibrated (mean
# rank-fraction ~0.22-0.24, biased HIGH) while yr[1]/yr[2]/yr[3] were all
# miscalibrated in the opposite direction (~0.75, biased LOW) -- a
# one-directional leak between the two, exactly the same mechanism (and
# same fix) as Location's own historical problem: with no sum-to-zero
# constraint, yr[Yr] isn't separately identified from loga[Mg]'s overall
# level, and that ambiguity gets resolved inconsistently across replicates.

source('src/hiermod/0_INDEX.R')

# Built as its own sibling model rather than editing MDSY in place, so the
# two stay independently comparable (same convention as MDLS2v/MDLS2vz).
#
# Stan's native sum_to_zero_vector isn't supported by this rethinking::ulam()
# version. What does work: N-1 free scalar parameters,
# with the Nth level's effect computed per-row as the negative sum of the
# others via plain arithmetic. guarantees sum=0 exactly, no vector/array
# construct involved. For Year (3 levels, one fewer than Location's 4):
# yr1/yr2 free, yr3 implied as -(yr1+yr2).
#
# Known, accepted prior asymmetry (same caveat as MDLS2vz): yr1/yr2 are iid
# dnorm(0,1), so the derived yr3 has ~2x their prior variance a priori --
# not perfectly exchangeable before the likelihood, though MDLS2vz found in
# practice the shared constraint + likelihood regularized all levels to
# comparable posterior SDs anyway. Worth re-checking here if SBC still
# shows any residual asymmetry specifically on yr3.

model_MDSYz_16S <- alist(
  likelihood = Dv ~ dlnorm(mu, sigma[Mg]),
  main_model = mu <- loga[Mg] + gamma*(Mo-1) + yr_eff,
  gamma_def  = gamma <- s_conv + gap_shift*(Mg-1), # interactive Season effect

  prior_loga = loga[Mg]  ~ dnorm(5,2),
  prior_s    = s_conv    ~ dnorm(0,1),
  prior_gs   = gap_shift ~ dnorm(0,1),
  pr_sigma   = sigma[Mg] ~ dhalfnorm(0,1),

  # Year, sum-to-zero: yr1/yr2 free, yr3 = -(yr1+yr2)
  yr_eff_def = yr_eff <- yr1*(Yr==1) + yr2*(Yr==2) - (yr1+yr2)*(Yr==3),
  prior_yr1  = yr1 ~ dnorm(0,1),
  prior_yr2  = yr2 ~ dnorm(0,1)
)

attr(model_MDSYz_16S, "name") <- "Elrond the Ageless"
model_id_MDSYz <- "MDSYz"

## ITS variant -----------------------------------------------------------------
# Referencing model_MDSYz_16S directly (not MDSY_16S) picks up the FULL,
# already-fixed structure in one step -- sum-to-zero yr_eff_def, yr1/yr2
# priors, everything -- so the naive/miscalibrated MDSY intermediate stage
# (4.2_MDSY_16S_calibration.R's own floor test) doesn't need re-running for
# ITS at all; that bug is already understood and fixed structurally, not a
# per-Kingdom finding. Same loga-only override as MDS_ITS/MDST_ITS.
model_MDSYz_ITS <- model_MDSYz_16S
model_MDSYz_ITS$prior_loga <- quote(loga[Mg] ~ dnorm(2,2))
attr(model_MDSYz_ITS, "name") <- "Elrond the Ageless"

## means_MDSYz()/dq_MDSYz -----------------------------------------------------
# Same formulas as means_MDSY()/dq_MDSY -- Year still doesn't enter the
# reported Mg x Mo estimand or its variance, sum-to-zero or not.

means_MDSYz <- function(post, shift = 0){
  total_var_conv <- post$sigma[,1]^2
  total_var_org  <- post$sigma[,2]^2

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

dq_MDSYz <- SBC::derived_quantities(
  may_gap =
    exp(loga[2] + sigma[2]^2 / 2) -
    exp(loga[1] + sigma[1]^2 / 2),
  july_gap =
    exp(loga[2] + s_conv + gap_shift + sigma[2]^2 / 2) -
    exp(loga[1] + s_conv +             sigma[1]^2 / 2),
  seasonal_change =
    (exp(loga[2] + s_conv + gap_shift + sigma[2]^2 / 2) -
       exp(loga[1] + s_conv +             sigma[1]^2 / 2)) -
    (exp(loga[2] + sigma[2]^2 / 2) -
       exp(loga[1] + sigma[1]^2 / 2))
)

## Data-generating function ---------------------------------------------------
# Same balanced Mg x Mo design as sim_div_MDSY(), but takes yr1/yr2 directly
# (matching the model's own free parameters) and derives yr3 the same way
# the model does, so the simulator and the likelihood agree exactly.

sim_div_MDSYz <- function(N_samples, loga, s_conv, gap_shift, sigma, yr1, yr2, shift = NULL){
  n_unit <- N_samples %/% 2 # 2 rows/unit (May + July)
  yr_vec <- c(yr1, yr2, -(yr1 + yr2))

  units <- tibble(Un = seq_len(n_unit), Mg = rbern(n_unit) + 1, Yr = sample(3, n_unit, replace = TRUE))
  dat <- units %>% crossing(Mo = 1:2) %>% arrange(Un)

  gamma <- s_conv + gap_shift*(dat$Mg - 1)
  mu <- loga[dat$Mg] + gamma*(dat$Mo - 1) + yr_vec[dat$Yr]

  dat$Dv <- rlnorm(nrow(dat), meanlog = mu, sdlog = sigma[dat$Mg])
  if (!is.null(shift)) dat$Dv_shifted <- shift + dat$Dv
  dat
}

simulate_from_priors_MDSYz <- function(true_params, N_samples = 250, shift = NULL){
  sim_div_MDSYz(
    N_samples = N_samples,
    loga = true_params$loga,
    s_conv = true_params$s_conv,
    gap_shift = true_params$gap_shift,
    sigma = true_params$sigma,
    yr1 = true_params$yr1,
    yr2 = true_params$yr2,
    shift = shift
  )
}
