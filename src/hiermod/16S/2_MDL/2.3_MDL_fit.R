# MODEL 2B (MDLb), 16S: Model 2's priors retuned for 16S's much higher raw
# diversity scale (loga dnorm(2,2) -> dnorm(6,2); sigma dexp(1) -> dexp(2)/
# dexp(2)) per 2.2_MDL_calibration.R's prior predictive check. Real fit + PPC.

hiermod_marker <- "16S"
source('src/hiermod/0_SETUP.R')
source('src/hiermod/Models/MDL_model.R') # means_MDL(), sim_div_MDL()

hiermod_out_dir <- "out/hiermod/16S_2_lognormal_MDL"
model_ppc1 <- readRDS(file.path(hiermod_out_dir, "model_ppc1.rds"))

## Model fit ----------------------------------------------------------------

dat <- list(
  Dv = div$Hill_1,
  Mg = idx$Mg$to_index(div$Management),
  Lo = idx$Lo$to_index(div$Location)
)

fitb <- ulam(
  model_ppc1,
  data = dat,
  chains = 6, cores = 6, iter = 10000,
  control = list(adapt_delta = 0.99)
)
save_fit("fit", "MDLb", fitb)

precis(fitb, depth = 2 )

save_pdf("fit_traceplot", "MDLb", function() traceplot(fitb))
save_pdf("fit_trankplot", "MDLb", function() trankplot(fitb))

## Posterior predictive check --------------------------------------------------
# Overlay the real data with data simulated from the posterior distribution.
# AFAIK, we are showing the estimate (i.e. its distribution), on which we
# overlay the real data and see if it makes sense. It doesn't have to fit
# perfectly, and the differences we see between posterior and actual data
# are a story.

### Overall ----

pb_postpred <- plot_ppc_overlay(fitb, dat, idx$Mg$to_label(dat$Mg),
                                xlim = c(0,2000)); pb_postpred
# Not bad!
# Conventional has a hump around 800-900 which the model doesn't see.

save_gg("postpred_density", "MDLb", pb_postpred)


### By location (marginal shape) ----
postpred_long <-
  sim(fitb, dat, n = 1000) %>%
  postpred_as_long_tibble()

# n-annotated Location labels for this plot only (to_label_n(), make_index()
# in hiermod_core.R) -- postpred_long/obs_df keep plain `Lo` for joins;
# Lo_n is purely a display label carrying the per-cultivar sample size.
Loc_labels_n <- tibble(
  Lo   = idx$Lo$to_label(seq_along(idx$Lo$levels)),
  Lo_n = idx$Lo$to_label_n(seq_along(idx$Lo$levels))
)

# stratification variables for postpred_long's `obs` index
pp_joined <- tibble(
  obs = seq_len(length(dat$Mg)),
  Mg  = idx$Mg$to_label(dat$Mg),
  Lo  = idx$Lo$to_label(dat$Lo)) %>%
  left_join(postpred_long, by = "obs") %>%
  left_join(Loc_labels_n, by = "Lo")

obs_df <- tibble(
  Dv = dat$Dv,
  Mg = idx$Mg$to_label(dat$Mg),
  Lo  = idx$Lo$to_label(dat$Lo)
) %>% left_join(Loc_labels_n, by = "Lo")

# quick summaries
obs_df %>%
  group_by(Mg, Lo) %>%
  summarise(med_obs = median(Dv), MAD = mad(Dv), .groups = 'drop')

pp_joined %>%
  group_by(Mg, Lo) %>%
  summarise(med_obs = median(Dv_sim), MAD = mad(Dv_sim), .groups = 'drop')

(xlim_upper_pp <- quantile(pp_joined$Dv_sim, 0.99))
quantile(pp_joined$Dv_sim, 0.999)

p_postpred_ridges <- pp_joined %>%
  ggplot(
    aes(
      x = Dv_sim,
      y = Lo_n,
      height = after_stat(density),
      fill = Mg)) +
  ggridges::geom_density_ridges(
    stat = "density",
    alpha = 0.55,
    colour = "white",
    scale = 0.8,
    linewidth = 0.3) +

  geom_point(
    data = obs_df, inherit.aes = FALSE,
    aes(x = Dv, y = as.numeric(Lo_n) + 0.4,
        fill = Mg,
        shape = Lo),
    position = position_jitter(height = 0.1, width = 0),
    size = 2, alpha = 0.5, colour = "grey20"
  ) +

  scale_shape_manual(values = c(21:24)) +
  scale_fill_manual(values = Management_palette) +
  scale_colour_manual(values = Management_palette) +
  coord_cartesian(xlim = c(0, xlim_upper_pp)) +
  theme_light() +
  guides(colour = 'none') +
  labs(x = "Diversity", y = "Location", fill = "Management",
       caption = "Ridges = posterior predictive density; points = observed data (jittered)"); p_postpred_ridges

save_gg("postpred_ridges", "MDL", p_postpred_ridges)

### Contrast test statistics ----
# Model 2 is a clean 2-group (Mg) split, so plot_ppc_contrast_stat() (not
# plot_ppc_season_contrast_stats(), which is specifically for the Mg x Mo
# design from Model 4 onward) is the right generic helper here. Median and
# mean contrast checked against the real data's own gap (paralleling the
# mean+median pair already reported in Contrast recovery above); MAD
# checks the heteroscedasticity assumption itself (sigma[Mg] differing by
# group is the whole point of Model 2 over Model 1), not redundant with
# the other two. All sit comfortably within their predictive spread (see
# discussion) -- not a red flag.

p_ppc_median_contrast <- plot_ppc_contrast_stat(fitb, dat, dat$Mg, median, "median", idx$Mg$levels)
p_ppc_mean_contrast   <- plot_ppc_contrast_stat(fitb, dat, dat$Mg, mean, "mean", idx$Mg$levels)
p_ppc_mad_contrast    <- plot_ppc_contrast_stat(fitb, dat, dat$Mg, mad, "dispersion (MAD)", idx$Mg$levels)

# Median is spot on at the peak; others are ok but slightly tailed, has to do with variance
(p_ppc <- p_ppc_median_contrast / p_ppc_mean_contrast / p_ppc_mad_contrast +
    plot_layout(guides = 'collect'))
save_gg("postpred_stat", "MDLb", p_ppc)
