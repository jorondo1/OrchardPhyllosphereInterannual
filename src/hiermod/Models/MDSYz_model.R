# MDSYz_model.R --- MODEL 4 (MDSYz), 16S: MDSY + a sum-to-zero constraint on
# yr[Yr]. SBC on MDSY found loga[1]/loga[2] severely miscalibrated (mean
# rank-fraction ~0.22-0.24, biased HIGH) while yr[1]/yr[2]/yr[3] were all
# miscalibrated in the opposite direction (~0.75, biased LOW) -- a
# one-directional leak between the two, exactly the same mechanism (and
# same fix) as Location's own historical problem: with no sum-to-zero
# constraint, yr[Yr] isn't separately identified from loga[Mg]'s overall
# level, and that ambiguity gets resolved inconsistently across replicates.
#
# Built as its own sibling model rather than editing MDSY in place, so the
# two stay independently comparable (same convention as MDLS2v/MDLS2vz).
#
# Stan's native sum_to_zero_vector isn't supported by this rethinking::ulam()
# version (confirmed when Location hit this same problem -- see
# src/hiermod/Models/archive/MDLS2vz_model.R's own header for the two
# approaches that didn't work). What does work: N-1 free scalar parameters,
# with the Nth level's effect computed per-row as the negative sum of the
# others via plain arithmetic -- guarantees sum=0 exactly, no vector/array
# construct involved. For Year (3 levels, one fewer than Location's 4):
# yr1/yr2 free, yr3 implied as -(yr1+yr2).
#
# Known, accepted prior asymmetry (same caveat as MDLS2vz): yr1/yr2 are iid
# dnorm(0,1), so the derived yr3 has ~2x their prior variance a priori --
# not perfectly exchangeable before the likelihood, though MDLS2vz found in
# practice the shared constraint + likelihood regularized all levels to
# comparable posterior SDs anyway. Worth re-checking here if SBC still
# shows any residual asymmetry specifically on yr3.

source('src/hiermod/Models/MDSY_model.R') # model_MDSY_16S, means_MDSY(), dq_MDSY

model_MDSYz_16S <- model_MDSY_16S
model_MDSYz_16S$main_model <- quote(
  mu <- loga[Mg] + gamma*(Mo-1) + yr_eff
)
model_MDSYz_16S$yr_eff_def <- quote(
  yr_eff <- yr1*(Yr==1) + yr2*(Yr==2) - (yr1+yr2)*(Yr==3)
)
model_MDSYz_16S$prior_yr  <- NULL # drop MDSY's old yr[Yr] ~ dnorm(0,1)
model_MDSYz_16S$prior_yr1 <- quote(yr1 ~ dnorm(0,1))
model_MDSYz_16S$prior_yr2 <- quote(yr2 ~ dnorm(0,1))

model_id_MDSYz <- "MDSYz"

## means_MDSYz()/dq_MDSYz -----------------------------------------------------
# Identical to means_MDSY()/dq_MDSY -- Year still doesn't enter the reported
# Mg x Mo estimand or its variance, sum-to-zero or not.
means_MDSYz <- means_MDSY
dq_MDSYz    <- dq_MDSY

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
