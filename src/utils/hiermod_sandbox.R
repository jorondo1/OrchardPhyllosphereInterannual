# ==============================================================================
# hiermod_sandbox.R -- archived / not sourced by default
# ==============================================================================
# Retired helpers from hiermod_utils.R, kept here for reference. Not sourced
# by hiermod_setup_ITS.R or any model script -- source this file explicitly
# if you actually want one of these back.
# ==============================================================================

# ---- postpred_density_by_Mg() -- archived [was in hiermod_utils.R] ----------
# Retired in favour of hardcoding each model's posterior predictive check
# inline (see model_MDLb's PPC section in hiermod_lognormal_ITS_3.R for a
# worked example -- observed vs. posterior-predictive medians by Location,
# not just by Management) or, for a quick generic density overlay, reaching
# for bayesplot (already a project dependency) instead of this:
#
#   bayesplot::ppc_dens_overlay(y, yrep)          -- pooled overlay
#   bayesplot::ppc_dens_overlay_grouped(y, yrep, group)  -- faceted, e.g. by Mg
#
# where yrep is a [S draws] x [N obs] matrix, e.g.
#   yrep <- sim(fit_MDC, dat_MDC, n = 500)
#   bayesplot::ppc_dens_overlay_grouped(dat_MDC$Dv, yrep, group = idx$Mg$to_label(dat_MDC$Mg))
#
# One gotcha to carry over if you use it: bayesplot's default per-curve
# bw = "nrd0" has the exact same outlier-blowup failure mode diagnosed for
# the prior-predictive spaghetti plots (a single extreme sim() draw inflates
# that curve's bandwidth and smears it flat). Pass an explicit shared `bw`
# (e.g. bw.nrd0() computed on the pooled, 99th-percentile-clipped yrep) if a
# fit's posterior still puts real mass on extreme values. Unlike the prior-
# predictive case, yrep here is a genuine matrix (S x N, same N real
# observations replicated every draw) so global winsorizing (pmin(yrep, cap))
# keeps it rectangular -- that's what makes bayesplot a good fit for PPC
# specifically, and an awkward one for prior predictive checks (independent,
# variable-length replicate datasets with no shared observation identity;
# see prior_predictive_spaghetti() in hiermod_utils.R for that case instead).

postpred_density_by_Mg <- function(fit, div, management_idx, data_list, shift = 0){
  Mg <- data_list$Mg
  Dv_sim <- shift + sim(fit, data = data_list)  # draws x N matrix

  plot_df <- map_dfr(seq_along(management_idx$levels), function(m){
    bind_rows(
      tibble(value = div$Hill_1[Mg == m],         source = "Observed",  Management = management_idx$levels[m]),
      tibble(value = as.vector(Dv_sim[, Mg == m]), source = "Simulated", Management = management_idx$levels[m])
    )
  })

  # clip the view to where the simulated mass actually is; the density
  # itself is still computed on the full (unclipped) data above
  xlim_upper <- quantile(plot_df$value[plot_df$source == "Simulated"], 0.99)

  ggplot(plot_df, aes(x = value, colour = source, fill = source)) +
    geom_density(alpha = 0.25, linewidth = 1) +
    coord_cartesian(xlim = c(0, xlim_upper)) +
    facet_wrap(~ Management) +
    theme_light() +
    labs(x = "Diversity (Hill_1)", y = "Density",
         title = "Posterior predictive check")
}
