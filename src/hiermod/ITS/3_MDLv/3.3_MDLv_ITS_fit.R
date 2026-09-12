# MODEL 3 (MDLv): Organic-Conventional gap varies by Location. Real fit and PPC.

hiermod_marker <- "ITS"
source('src/hiermod/0_SETUP.R')
source('src/hiermod/Models/MDLv_model.R') # model_MDLv_ITS, means_MDLv(), mdlv_labels, sim_div_MDLv()
model <- model_MDLv_ITS

hiermod_out_dir <- "out/hiermod/ITS_3_lognormal_MDLv"

## Model fit ----------------------------------------------------------------

dat <- list(
  Dv = div$Hill_1,
  Mg = idx$Mg$to_index(div$Management),
  Lo = idx$Lo$to_index(div$Location)
)

fitb <- ulam(
  model,
  data = dat,
  chains = 6, cores = 6, iter = 10000,
  control = list(adapt_delta = 0.99)
)
save_fit("fit", "MDLv", fitb)

precis(fitb, depth = 2)

traceplot(fitb); trankplot(fitb)
save_pdf("fit_traceplot", "MDLv", function() traceplot(fitb))
save_pdf("fit_trankplot", "MDLv", function() trankplot(fitb))

## Posterior predictive check --------------------------------------------------

pb_postpred <- plot_ppc_overlay(fitb, dat, idx$Mg$to_label(dat$Mg), xlim = c(0,150)); pb_postpred

save_gg("postpred_density", "MDLv", pb_postpred)

### By Location ---------------------------

postpred_long <- sim(fitb, dat, n = 1000) %>% postpred_as_long_tibble()

Loc_labels_n <- tibble(
  Lo   = idx$Lo$to_label(seq_along(idx$Lo$levels)),
  Lo_n = idx$Lo$to_label_n(seq_along(idx$Lo$levels))
)

pp_joined <- tibble(
  obs = seq_len(length(dat$Mg)),
  Mg  = idx$Mg$to_label(dat$Mg),
  Lo  = idx$Lo$to_label(dat$Lo)) %>%
  left_join(postpred_long, by = "obs") %>%
  left_join(Loc_labels_n, by = "Lo")

obs_df <- tibble(
  Dv = dat$Dv,
  Mg = idx$Mg$to_label(dat$Mg),
  Lo = idx$Lo$to_label(dat$Lo)
) %>% left_join(Loc_labels_n, by = "Lo")

xlim_upper_pp <- quantile(pp_joined$Dv_sim, 0.98)

p_postpred_ridges <- pp_joined %>%
  ggplot(aes(x = Dv_sim, y = Lo_n, height = after_stat(density), fill = Mg)) +
  ggridges::geom_density_ridges(
    stat = "density", alpha = 0.55, colour = "white", scale = 0.8, linewidth = 0.3) +
  geom_point(
    data = obs_df, inherit.aes = FALSE,
    shape = 21, stroke = 0.2, size = 2, alpha = 0.6,
    aes(x = Dv, y = as.numeric(Lo_n) - 0.15, fill = Mg),
    position = position_jitter(height = 0.1, width = 0)) +
  scale_fill_manual(values = Management_palette) +
  scale_colour_manual(values = Management_palette) +
  coord_cartesian(xlim = c(0, xlim_upper_pp)) +
  theme_light() +
  guides(colour = 'none') +
  labs(x = "Diversity", y = "Location", fill = "Management",
       caption = "Ridges = posterior predictive density; points = observed data (jittered)"); p_postpred_ridges

save_gg("postpred_ridges", "MDLv", p_postpred_ridges)

## Contrast statistic ------------------------------
# Same idea as MDLb's: does the model reproduce the specific gap we're
# reporting, not just plausible marginal shapes per group?

p_ppc_median_contrast <- plot_ppc_contrast_stat(fitb, dat, dat$Mg, median, "median", idx$Mg$levels)
p_ppc_mean_contrast   <- plot_ppc_contrast_stat(fitb, dat, dat$Mg, mean, "mean", idx$Mg$levels)
p_ppc_mad_contrast    <- plot_ppc_contrast_stat(fitb, dat, dat$Mg, mad, "dispersion (MAD)", idx$Mg$levels)

(p_ppc <- p_ppc_median_contrast / p_ppc_mean_contrast / p_ppc_mad_contrast)
save_gg("postpred_stat", "MDLv", p_ppc)

# Next step: Tree ID + Season/Year. Deliberately not bundled in here -- see
# write-up. All Dv rows currently get treated as independent even though many
# trees contribute two rows (May + July) within a year; that's the next
# thing to model properly, on its own, not tacked onto this one.
