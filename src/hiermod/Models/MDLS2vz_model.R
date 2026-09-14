# MODEL 5c (MDLS2vz), 16S only: MDLS2v + a sum-to-zero constraint on the
# Location random effect b[Lo]. SBC on MDLS2v found loga[1]/loga[2]
# severely miscalibrated (mean rank-fraction ~0.16, not ~0.5), correlating
# strongly with sigma_loc (r=0.75) and sigma_tr (r=0.54) but not with ls0
# (r~0). With only n_loc=4 locations, the 4 b[Lo] draws don't reliably
# average to zero on their own, and that chance leftover gets partly
# absorbed into loga instead of staying in the random effect. This is
# worse for Location than Tree specifically because Location has far
# fewer levels (4 vs ~120+). A sum-to-zero constraint removes that
# arithmetic ambiguity, but it does NOT fix the separate, real confound that
# 2 of the 4 locations are single-management (Location B/D-only sensitivity
# check planned for Model 7).
#
# Built as its own sibling model rather than editing MDLS2v in place, so
# the two stay independently comparable.

source('src/hiermod/Models/MDLS2v_model.R') # model_MDLS2v_16S, means_MDLS2v(), sigma_cell(), true_sigma_from_ls(), contrast_may_gap_MDLS2v()

# Stan's native sum_to_zero_vector type (would give all 4 locations equal
# marginal variance via an isometric-log-ratio transform) is NOT supported
# by this version of rethinking::ulam() -- confirmed directly: its
# declare_templates list only has real/vector/row_vector/matrix/int/
# int_array/corr_matrix/cholesky_factor_corr/ordered/simplex. The next
# option (3 free values + a transpars> vector[4] <- append_row(...) whole-
# vector assignment) also failed: ulam()'s compose_assignment() always
# wraps a "vector"-typed transformed-parameter assignment in a per-index
# loop (b[i] = RHS) regardless of whether RHS is itself already a whole
# vector -- confirmed via its own generated Stan code, which tried
# `b[i] = append_row(b_raw, -sum(b_raw))` (a type error: vector assigned to
# a scalar slot). Only "matrix" types are exempted from that loop.
#
# What DOES work, using only constructs already proven inside this exact
# model (main_model/gamma_def are themselves per-row "local" formulas):
# 3 independent scalar parameters b1/b2/b3, with the 4th location's effect
# computed per-row as -(b1+b2+b3) via plain arithmetic (Lo==k comparisons
# are just 0/1-valued). This guarantees sum(b)=0 exactly (verified: max
# abs deviation across posterior draws = 0), using nothing but scalar
# priors and a per-row arithmetic formula -- no vector/array declaration
# involved at all, so none of the above landmines apply.
#
# Known, accepted PRIOR asymmetry: b1/b2/b3 are i.i.d. dnorm(0,1), so the
# derived 4th value (-(b1+b2+b3)) has 3x their prior variance -- unlike
# Stan's native transform, this naive construction isn't fully exchangeable
# a priori. In practice this washed out in 5c.2's own calibration check:
# with real data/likelihood feeding in, the 4 posterior SDs came back
# comparable (~0.33-0.41, not a 3x split) -- the shared sum-to-zero
# constraint and the likelihood itself regularize all 4 together. Still
# worth re-checking whenever true sigma_loc is large, since that's exactly
# the regime the prior asymmetry would matter most.
model_MDLS2vz_16S <- model_MDLS2v_16S
model_MDLS2vz_16S$main_model <- quote(
  mu <- loga[Mg] + gamma*(Mo-1) + loc_eff*sigma_loc + yr[Yr] + tr[Tr]*sigma_tr
)
model_MDLS2vz_16S$loc_eff_def <- quote(
  loc_eff <- b1*(Lo==1) + b2*(Lo==2) + b3*(Lo==3) - (b1+b2+b3)*(Lo==4)
)
model_MDLS2vz_16S$prior_b  <- NULL  # drop MDLS2v's old b[Lo] ~ dnorm(0,1)
model_MDLS2vz_16S$prior_b1 <- quote(b1 ~ dnorm(0,1))
model_MDLS2vz_16S$prior_b2 <- quote(b2 ~ dnorm(0,1))
model_MDLS2vz_16S$prior_b3 <- quote(b3 ~ dnorm(0,1))
# Tried reporting all 4 location effects as gq> generated quantities for
# analysis-script convenience -- ulam()'s gq> handling fails even for the
# most trivial scalar-from-scalar case here ("Unable to determine type and
# dimensions", confirmed in isolation, not an artifact of the arithmetic).
# Not needed for the actual fix: b_loc4 is trivially reconstructed in R
# post-hoc from extract.samples() as -(post$b1 + post$b2 + post$b3).

model_id <- "MDLS2vz"

# Data-generating function ---------------------------------------------
# Copy of sim_div_MDLS2() (MDLS2_model.R) with exactly one change:
# loc_offset is demeaned so it sums to zero, matching the fitting model's
# constrained b[Lo] -- required for SBC to test the model against
# self-consistent data. sim_div_MDLS2() itself is left untouched (still
# used by MDLS2/MDLS2v/MDLSY/MDLSYv's own calibration) rather than edited
# in place, matching this model's own "sibling, not in-place" scope.
sim_div_MDLS2vz <- function(
    N_samples, loga, s_conv, gap_shift, sigma, year_offset,
    sigma_loc = 0.5, sigma_tr = 0.3,
    n_loc, p_dropout = 0,
    shift = NULL){
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
  dat$cell <- (dat$Mg - 1) * 2 + dat$Mo

  loc_offset_raw <- rnorm(n_loc, 0, sigma_loc)
  loc_offset  <- loc_offset_raw - mean(loc_offset_raw)  # <<< the fix: exact sum-to-zero
  tree_offset <- rnorm(n_tree, 0, sigma_tr)

  gamma <- s_conv + gap_shift*(dat$Mg - 1)
  mu <- loga[dat$Mg] + gamma*(dat$Mo - 1) +
    loc_offset[dat$Lo] + year_offset[dat$Yr] + tree_offset[dat$Tr]

  dat$Dv <- rlnorm(nrow(dat), meanlog = mu, sdlog = sigma[dat$cell])
  if (!is.null(shift)) dat$Dv_shifted <- shift + dat$Dv
  dat
}

simulate_from_priors_MDLS2vz <- function(true_params, N_samples = 250,
                                          year_offset = c(0, 0.3, -0.2),
                                          n_loc = 4, p_dropout = 0.1, shift = NULL){
  true_sigma <- true_sigma_from_ls(true_params)

  sim_div_MDLS2vz(
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

# means_MDLS2v() and contrast_may_gap_MDLS2v() (MDLS2v_model.R) are reused
# unchanged -- neither references b[Lo] directly, only sigma_loc (the
# scale hyperparameter), which this change doesn't touch.
