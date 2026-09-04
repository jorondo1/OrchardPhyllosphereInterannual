# 3.3_MDLv_analysis.R -- MODEL 3 (MDLv): posterior contrast, run against the
# fit saved by 3.2_MDLv_validation.R -- no refit needed.

source('src/hiermod/ITS/0_SETUP.R')
source('src/hiermod/ITS/3.1_MDLv_model.R') # model, means_MDLv(), mdlv_labels
hiermod_out_dir <- "out/hiermod/ITS_3_lognormal_MDLv"

fitb <- readRDS(file.path(hiermod_out_dir, "fit_MDLv.rds"))

dat <- list(
  Dv = div$Hill_1,
  Mg = idx$Mg$to_index(div$Management),
  Lo = idx$Lo$to_index(div$Location)
)

## Posterior contrast ---------------------------------------------------------

pf <- post_full(fitb, means_MDLv)      # raw draws + mean/median; pf$g feeds g_by_loc below

pc_full <- compute_contrasts(
  pf, labels = mdlv_labels,
  group_levels = idx$Mg$levels)

p_contrasts_panel <- contrast_plot_panels(
  pc_full,
  quant = c(0.001,0.995),
  group_pal = Management_palette); p_contrasts_panel

# proportion below zero?
pc_full %>%
  filter(group == "Contrast") %>%
  group_by(statistic) %>%
  summarise(n_neg = sum(value<=0), n = n()) %>%
  mutate(prop_neg = n_neg/n) # ~25%

save_report("fit_summary", "MDLv", fitb, pc_full, model, model_name = "The Copycat")
save_gg("fit_contrasts_panel", "MDLv", p_contrasts_panel)

# Does the gap actually vary by Location, or was MDLb's shared-slope
# assumption fine all along? sigma_g sits comfortably above 0, therefore
# Location-varying Management is a real pattern here, not just added flexibility.

# g[Lo] per Location: how far that Location's gap sits from the population
# average. Saint-Benoît/Windsor are Conventional-only/Organic-only -- see
# the model definition's caveat above before reading too much into those two.
g_mat <- as.matrix(pf$g)
g_by_loc <- tibble(
  Lo       = idx$Lo$levels,
  median_g = apply(g_mat, 2, median),
  lo89     = apply(g_mat, 2, function(x) PI(x)[1]),
  hi89     = apply(g_mat, 2, function(x) PI(x)[2])
); g_by_loc


# g_i -> Location name; row_number() within each group re-creates the draw
# index that post_full()'s >2-category branch doesn't carry explicitly
# (rows are still in original draw order, just not labelled as such)
g_wide <- pc_full %>%
  filter(statistic == "g") %>%
  mutate(Lo = idx$Lo$levels[as.integer(str_remove(group, "g_"))]) %>%
  group_by(Lo) %>%
  mutate(draw = row_number()) %>%
  ungroup() %>%
  select(draw, Lo, value) %>%
  pivot_wider(names_from = Lo, values_from = value)   # one column per Location, one row per draw

# every pairwise contrast (col2 - col1), long format -- a Location-pairwise
# breakdown, not the same thing as the (dropped) Management forest plot.
loc_pairs <- combn(idx$Lo$levels, 2, simplify = FALSE)

g_contrasts <- map_dfr(loc_pairs, function(pair){
  tibble(
    contrast = paste(pair[2], "-", pair[1]),
    value    = g_wide[[pair[2]]] - g_wide[[pair[1]]]
  )
})

# summarize + simple point-range plot -- 10 pairs is a lot for stacked
# densities, a point-range per pair reads faster
p_loc_pairs <- g_contrasts %>%
  group_by(contrast) %>%
  summarise(median = median(value), lo89 = PI(value)[1], hi89 = PI(value)[2]) %>%
  ggplot(aes(x = median, y = contrast)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_pointrange(aes(xmin = lo89, xmax = hi89)) +
  labs(x = "Location-gap contrast (g[Lo] difference)", y = NULL,
       title = "Pairwise contrasts of the Organic-Conventional gap by Location")

save_gg("fit_posterior_locations", "MDLv", p_loc_pairs)
