# MDSTYCL_model.R --- MODEL 9 (MDSTYCL, "Faramir the Judicious"), 16S:
# MDSTYCV with Cultivar's fixed effect (cv[Cv]) REPLACED by Location's
# fixed effect (lo[Lo]), same sum-to-zero recipe.

source('src/hiermod/0_INDEX.R')

# Why this exists: Model 8 (MDSTYCVr, "Gimli the Greedy") tried making
# Cultivar a RANDOM effect instead and failed badly -- real SBC
# (n_sbc=500): 297 divergences, 154/500 (30.8%) fits Rhat>1.01, and
# loga[1]/loga[2] themselves severely MISCALIBRATED (z=-11.96/-14.22),
# not just sigma_cv (z=-5.76). Nesting sigma_cv with sigma_tr over the
# same 129 trees didn't just fail to identify sigma_cv cleanly -- it
# leaked into the headline estimand. Back to a fixed-effect design instead
# of trying to fix that nesting.
#
# Location is a natural next fixed effect to test (never in this rebuild's
# family before -- Tree is deterministically nested in Location, same
# 129/129 mapping as Cultivar, so structurally this is the same "does
# Tree's random effect coexist with a coarser fixed grouping over the same
# trees" question MDSTYCV already answered cleanly for Cultivar).
#
# Cultivar can NOT coexist with Location here, though -- confirmed by
# direct query of the real data, not just assumed: restricted to the two
# Locations with both Management levels present (B, D -- see this
# model's own calibration/fit scripts for why), Location D has ONLY
# Honeycrisp and Spartan (0 Cortland/Liberty/Paulared) -- i.e. within that
# subset, 3 of Cultivar's 5 levels are perfectly aliased with "Location B",
# not just correlated with it. Cultivar is dropped entirely for this
# model, not added alongside Location.
#
# Real-data motivation for the subset (see 9.2/9.3's own headers for the
# full numbers): Location x Management in the full data is A=Conventional-
# only (50), C=Organic-only (51), B and D have both -- fitting Location
# on the full 4-level factor would confound Location with Management for
# A/C. Model 9's real fit (9.3, not yet built) will restrict to the B/D
# subset (141/242 rows) for exactly this reason. This calibration-stage
# model file itself doesn't know about that restriction -- Lo is just a
# 2-level index, assigned independently of Mg in the simulator, same
# "recoverable in principle" scope as every other calibration script in
# this family.
#
# loga[Mg]/s_conv/gap_shift/sigma[Mg]/yr1/yr2/tr[Tr]*sigma_tr/b_deg/
# b_precip/b_seq are MDSTYCV's own validated answer, hardcoded as this
# model's starting point. lo1 ~ dnorm(0,1) is the one new assumption to
# validate -- same sum-to-zero recipe as yr1/yr2 and cv_1/cv_3/cv_4/cv_5,
# just 2 levels (1 free parameter, 1 derived as its negative) instead of
# 3 or 5.
#
# ITS variant deliberately NOT built yet -- same discipline as Model 8:
# validate on 16S first.

source('src/hiermod/Models/MDSTYCV_model.R') # model_MDSTYCV_16S, means_MDSTYCV(), dq_MDSTYCV

model_MDSTYCL_16S <- model_MDSTYCV_16S
model_MDSTYCL_16S$main_model <- quote(
  mu <- loga[Mg] + gamma*(Mo-1) + yr_eff + lo_eff + tr[Tr]*sigma_tr +
    b_deg*deg_h_z + b_precip*precip_72h_z + b_seq*seq_depth_z
)
# Drop MDSTYCV's fixed Cultivar construction -- can't coexist with
# Location in the real B/D subset (see header).
model_MDSTYCL_16S$cv_eff_def <- NULL
model_MDSTYCL_16S$prior_cv1  <- NULL
model_MDSTYCL_16S$prior_cv3  <- NULL
model_MDSTYCL_16S$prior_cv4  <- NULL
model_MDSTYCL_16S$prior_cv5  <- NULL

# Location: 2 levels within the real B/D subset -- 1 free scalar, 1
# derived as its negative, same construction as every other sum-to-zero
# fixed effect in this family (Year: N-1=2 free; Cultivar: N-1=4 free;
# here N-1=1 free).
model_MDSTYCL_16S$lo_eff_def <- quote(
  lo_eff <- lo1*(Lo==1) - lo1*(Lo==2)
)
model_MDSTYCL_16S$prior_lo1 <- quote(lo1 ~ dnorm(0,1))

attr(model_MDSTYCL_16S, "name") <- "Faramir the Judicious"
model_id_MDSTYCL <- "MDSTYCL"

## means_MDSTYCL()/dq_MDSTYCL ----------------------------------------------------
# Identical to means_MDSTYCV()/dq_MDSTYCV -- Location, like Year, is a
# fixed effect held at its own observed-level average and doesn't enter
# the reported Mg x Mo estimand or its variance (same reasoning as
# means_MDSYCV() dropping Cultivar).
means_MDSTYCL <- means_MDSTYCV
dq_MDSTYCL    <- dq_MDSTYCV

## variance_partition_MDSTYCL() --------------------------------------------------
# Same redesign as variance_partition_MDSTYCV() (see its own comment for the
# by-margin rationale and why Management x Season stays one combined term) --
# Location swapped in for Cultivar, in the same slot (coarser grouping
# before the Tree it nests, same 129/129 deterministic mapping as Cultivar).
# Tree was already folded into this model's own total (unlike MDSTYCV's,
# until its own redesign) via the bulk sigma_tr^2 bolt-on; now a proper
# per-observation term like every other entry here, via
# variance_partition_panels()'s shared engine.

variance_partition_MDSTYCL <- function(post, dat){
  yr3    <- -(as.vector(post$yr1) + as.vector(post$yr2))
  yr_obs <- cbind(post$yr1, post$yr2, yr3)[, dat$Yr]

  lo2    <- -(as.vector(post$lo1))
  # NB: was `[dat$Lo]` (no comma) -- linear/column-major indexing on a
  # matrix, not column selection by dat$Lo like every other realized-level
  # term here (yr_obs, cv_obs). That silently returned near-constant
  # garbage instead of each observation's actual Location coefficient --
  # the likely cause of Location's R2 share collapsing to ~0 (and sometimes
  # negative).
  lo_obs <- cbind(post$lo1, lo2)[, dat$Lo]

  tree_obs <- sweep(post$tr[, dat$Tr], 1, as.vector(post$sigma_tr), "*")

  mg_n <- as.integer(table(factor(dat$Mg, levels = 1:2)))
  residual_var <- as.vector((post$sigma^2) %*% (mg_n / sum(mg_n)))

  loga_obs   <- post$loga[, dat$Mg]
  gamma      <- as.vector(post$s_conv) + outer(as.vector(post$gap_shift), dat$Mg - 1)
  gamma_term <- sweep(gamma, 2, dat$Mo - 1, "*")

  terms <- list(
    "Reads count"          = outer(as.vector(post$b_seq), dat$seq_depth_z),
    "Degree-hours"         = outer(as.vector(post$b_deg), dat$deg_h_z),
    "Precipitation"        = outer(as.vector(post$b_precip), dat$precip_72h_z),
    "Year"                 = yr_obs,
    "Location"             = lo_obs,
    "Tree"                 = tree_obs,
    "Management x Season"  = loga_obs + gamma_term
  )

  variance_partition_panels(terms, residual_var)
}

## Data-generating function ---------------------------------------------------
# Same Tree-as-study-unit design as sim_div_MDSTYCV(), Lo assigned per
# tree (2 levels) exactly like Cv was -- independently of Mg here (tests
# recoverability in principle; the real B/D subset's actual Location x
# Cultivar aliasing is a real-fit-stage concern, not a calibration one,
# see this file's own header).

sim_div_MDSTYCL <- function(N_samples, loga, s_conv, gap_shift, sigma, yr1, yr2,
                             lo1, sigma_tr, b_deg, b_precip, b_seq, shift = NULL){
  n_tree <- N_samples %/% 2 # 2 rows/tree (May + July)
  yr_vec <- c(yr1, yr2, -(yr1 + yr2))
  lo_vec <- c(lo1, -lo1) # Lo index order: 1, 2

  trees <- tibble(
    Tr = seq_len(n_tree), Mg = rbern(n_tree) + 1,
    Yr = sample(3, n_tree, replace = TRUE), Lo = sample(2, n_tree, replace = TRUE))
  dat <- trees %>% crossing(Mo = 1:2) %>% arrange(Tr)

  tree_offset <- rnorm(n_tree, 0, sigma_tr)

  gamma <- s_conv + gap_shift*(dat$Mg - 1)
  mu_structural <- loga[dat$Mg] + gamma*(dat$Mo - 1) + yr_vec[dat$Yr] +
    lo_vec[dat$Lo] + tree_offset[dat$Tr]

  dat$deg_h_z      <- rnorm(nrow(dat))
  dat$precip_72h_z <- rnorm(nrow(dat))
  dat$seq_depth_z  <- rnorm(nrow(dat))

  mu <- mu_structural + b_deg*dat$deg_h_z + b_precip*dat$precip_72h_z + b_seq*dat$seq_depth_z

  dat$Dv <- rlnorm(nrow(dat), meanlog = mu, sdlog = sigma[dat$Mg])
  if (!is.null(shift)) dat$Dv_shifted <- shift + dat$Dv
  dat
}

simulate_from_priors_MDSTYCL <- function(true_params, N_samples = 250, shift = NULL){
  sim_div_MDSTYCL(
    N_samples = N_samples,
    loga = true_params$loga, s_conv = true_params$s_conv, gap_shift = true_params$gap_shift,
    sigma = true_params$sigma, yr1 = true_params$yr1, yr2 = true_params$yr2,
    lo1 = true_params$lo1, sigma_tr = true_params$sigma_tr,
    b_deg = true_params$b_deg, b_precip = true_params$b_precip, b_seq = true_params$b_seq,
    shift = shift
  )
}
