# Shared preamble for every hiermod script (ITS and 16S). Set
# hiermod_marker <- "ITS" or "16S" before sourcing this file.

pacman::p_load(rethinking, tidyverse, bayesm, bayesplot, ggridges, magrittr, patchwork, rlang, scales, posterior, kableExtra)
source('src/utils/hiermod_core.R')
source('src/utils/saver_functions.R')
source('src/utils/postcontrast_helpers.R')
source('src/utils/predictive_checks.R')
source('src/utils/sbc_workflow.R') # SBC-package-based workflow (draw_true/SBC_backend_ulam/make_sbc_generator/run_sbc_pipeline/plot_sbc_diagnostics/save_sbc_health_report)

div_all <- read_rds('data/diversity_data.rds')
div <- if (hiermod_marker == "ITS") div_all$Fungi$alpha else div_all$Bacteria$alpha

message("Make sure to define hiermod_out_dir <- out/hiermod/<...>")

# Category <-> index codebooks (make_index(), hiermod_core.R) -- shared here so the index<->category mapping can't drift between scripts.

fill_loc  <- c("#5DB63B", "#2C9EE3", "gold", "#F8A11C")
fill_cult <- c("#7DB16B", "#45818E", "#AB4F84", "#AE9FCB", "#99CFE1")

idx <- list(
  Mg = make_index(div$Management, levels = c("Conventional", "Organic"),
                   palette = c(Conventional = "#F28E2B", Organic = "#499894")),
  Lo = make_index(div$Location, palette = fill_loc),
  Tr = make_index(div$Tree_id),
  Cv = make_index(div$Cultivar, levels = c("Cortland", "Liberty", "Paulared", "Honeycrisp", "Spartan"),
                   palette = fill_cult),
  Mo = make_index(div$Time, levels = c("May", "July")),
  Yr = make_index(div$Year, levels = c("2022", "2023", "2024"))
)

# Built from idx$Mg$palette; Contrast/Population are generic roles layered on top, not Management levels.
Management_palette <- c(idx$Mg$palette, Contrast = "#98494d", Population = "#895a92")

# Covariate effect panel labels/colours (Model 5+).
cov_labels <- c("Degree-hours", "Precipitation (72h)", "Seq. depth")
cov_pal    <- setNames(scales::hue_pal()(3), cov_labels)

# Variance-partition panel colours (Model 5+, variance_partition_*()'s own sequential decomposition).
Variance_partition_palette <- c(
  "Management x Season" = "#4C72B0",
  "Year"                = "#55A868",
  "Covariates"          = "#8172B2",
  "Cultivar"            = "#CCB974",
  "Tree"                = "#DD8452",
  "Residual"            = "grey50"
)

# Fold-change panel colours (Model 2+, contrast_plot_panels()'s ratio_stats=). Reference these names dynamically
# (names(Fold_change_palette)) in each script's own extra= list rather than retyping them -- a mismatch here silently
# breaks the combined ratio panel (it did once already).
Fold_change_palette <- c(
  "May fold difference"  = "#CC8FBB",
  "July fold difference" = "#7AAB32"
)

# deg_h_z is centered WITHIN Season (not globally) so the season term (s_conv/gap_shift) keeps the real, reliably
# monotonic (July always warmer) weather-driven seasonal signal instead of it leaking into the covariate slope.
# precip_72h_z is centered globally -- unlike temperature, precip has no reliable May-vs-July direction (e.g. 2022
# reverses it), so within-season centering would just bake a given year's idiosyncratic rain pattern into the season
# term instead of a real seasonal identity. seq_depth_z stays globally centered too -- its Management correlation is
# a confound we want removed, not preserved.
deg_h_season_mean <- tapply(div$deg_h, div$Time, mean)

div$deg_h_z <- (div$deg_h - deg_h_season_mean[div$Time])
div$deg_h_z <- div$deg_h_z / sd(div$deg_h_z)

div$precip_72h_z <- (div$precip_72h - mean(div$precip_72h)) / sd(div$precip_72h)

# Seq_depth logged first (right-skewed, ~20x range), then standardized.
log_seq_depth   <- log(div$Seq_depth)
div$seq_depth_z <- (log_seq_depth - mean(log_seq_depth)) / sd(log_seq_depth)

# themes
ggplot2::theme_set(
  ggplot2::theme_light() +
    theme(
      strip.text = element_text(colour = 'black'),
      strip.background = element_rect(
        fill = 'grey90',  linewidth = 0.2),
      panel.border = element_rect(colour = 'grey90'),
      panel.grid = element_blank(),
      panel.spacing = unit(0, "lines")
    ))
