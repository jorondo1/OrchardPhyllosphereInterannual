# 5.3_MDLS2_analysis.R -- MODEL 5: posterior contrast and variance
# components, run against the fit saved by 5.2_MDLS2_validation.R -- no
# refit needed.

source('src/hiermod/ITS/0_SETUP.R')
source('src/hiermod/ITS/5.1_MDLS2_model.R') # model, means_MDLS2()
hiermod_out_dir <- "out/hiermod/ITS_5_lognormal_MDLS2"

fit_MDLS2 <- readRDS(file.path(hiermod_out_dir, "fit_MDLS2.rds"))

pf <- post_full(fit_MDLS2, means_MDLS2)
m  <- pf$mean
md <- pf$median

## Management x Season contrasts -----------------------------------------------

pc_estimands_means <- estimand_panels(
  pairs = list(`May mean` = list(m$mean_1, m$mean_3),
               `July mean` = list(m$mean_2, m$mean_4)),
  extra = list("Change in management contrast, from May to July" =
                 (m$mean_4 - m$mean_2) - (m$mean_3 - m$mean_1)),
  group_levels = c("Conventional", "Organic")
)

# fit_MDLS2@formula (not the freshly-sourced `model`) -- the real fit was
# built from model_vbc (5.2_MDLS2_validation.R's variance-budget-calibrated
# priors), which this script never sees; the fit's own compiled formula is
# the only copy guaranteed to match what was actually fit.
save_report("fit_summary", "MDLS2", fit_MDLS2, pc_estimands_means, fit_MDLS2@formula, model_name = "The Splitter")

p_contrast_mean <- contrast_plot_panels(
  pc_estimands_means, quant = c(0.005, 0.995), scales = 'free_y',
  group_pal = Management_palette,
  legend_title = "Posteriors (population means)"); p_contrast_mean

save_gg("fit_contrast_mean", "MDLS2", p_contrast_mean)

pc_estimands_medians <- estimand_panels(
  pairs = list(`May median` = list(md$median_1, md$median_3),
               `July median`= list(md$median_2, md$median_4)),
  extra = list("Change in management contrast, from May to July" =
                 (md$median_4 - md$median_2) - (md$median_3 - md$median_1)),
  group_levels = c("Conventional", "Organic")
)

p_contrast_median <- contrast_plot_panels(
  pc_estimands_medians, quant = c(0.005,.995), scales = 'free_y',
  group_pal = Management_palette,
  legend_title = "Posteriors (population medians)"); p_contrast_median
save_gg("fit_contrast_median", "MDLS2", p_contrast_median, width = 10, height = 6)

## Variance components -----------------------------------------------------------

# estimand_panels() only builds paired 2-group panels (May vs July, etc.) --
# Location (4) and Year (3) aren't a 2-group contrast, so this uses
# compute_contrasts()'s already-existing >2-column branch instead: one
# density per column, group = that column's own name, no fabricated
# Contrast row. b[Lo] is the raw non-centered z-score, not interpretable on
# its own -- b[Lo]*sigma_loc (the actual log-scale Location offset, additive
# to mu identically across every Mg/Mo/Yr cell) is what's plotted: each
# density is that Location's typical deviation from the population-average
# Location, in log-diversity units, i.e. "how much impact does being at
# this orchard have, relative to a typical one."
#
# variance_component_panels() (postcontrast_helpers.R) gives each statistic
# its own patchwork panel with its own right-side legend, instead of one
# shared bottom legend -- with 4 Locations + 3 Years + sigma_loc/sigma_tr,
# a single combined legend was getting unreadably long.

stat_levels <- c("Location effect (log scale)",
                  "Location vs Tree effect SD (population)",
                  "Year effect (log scale)")

pf_re <- list(
  `Location effect (log scale)` = as_tibble(as.matrix(pf$b) * pf$sigma_loc[[1]]) %>%
    setNames(as.character(idx$Lo$to_label_n(seq_along(idx$Lo$levels)))),
  `Year effect (log scale)` = as_tibble(pf$yr) %>%
    setNames(as.character(idx$Yr$to_label_n(seq_along(idx$Yr$levels))))
)

# sigma_loc vs sigma_tr, overlaid in one panel -- built by hand rather than
# through compute_contrasts(), since a 2-column entry there hits the paired
# group1/group2/Contrast branch, which would fabricate a "sigma_tr -
# sigma_loc" difference that isn't a meaningful quantity (these are two
# separate random effects' population SDs, not two levels of one factor).
sigma_re <- bind_rows(
  tibble(statistic = stat_levels[2], group = "sigma_loc", value = pf$sigma_loc$sigma_loc),
  tibble(statistic = stat_levels[2], group = "sigma_tr",  value = pf$sigma_tr$sigma_tr)
)

# One combined palette: idx$Lo/idx$Yr's own colour schemes (palette_n(),
# hiermod_core.R) + sigma_loc/sigma_tr for the population-SD panel.
# sigma_tr kept neutral/grey -- it's here for scale context (how big is
# Tree variance relative to Location variance), not as a headline estimate.
re_pal <- c(
  idx$Lo$palette_n(),
  idx$Yr$palette_n(),
  sigma_loc = Management_palette[["Population"]],
  sigma_tr  = "grey40"
)

pc_random_effects <- bind_rows(compute_contrasts(pf_re, keep = names(pf_re)), sigma_re) %>%
  mutate(statistic = factor(statistic, levels = stat_levels))

p_random_effects <- variance_component_panels(
  pc_random_effects, quant = c(0.005, 0.995), palette = re_pal); p_random_effects

save_gg("fit_location_year_effects", "MDLS2", p_random_effects, width = 8, height = 10)
