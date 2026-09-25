# Shared preamble for every hiermod script (ITS and 16S). Set

pacman::p_load(rethinking, tidyverse, bayesm, bayesplot, ggridges, magrittr, patchwork, rlang, scales, posterior, kableExtra)
source('src/utils/hiermod_core.R')
source('src/utils/saver_functions.R')
source('src/utils/postcontrast_helpers.R')
source('src/utils/predictive_checks.R')
source('src/utils/sbc_workflow.R') # SBC-package-based workflow (draw_true/SBC_backend_ulam/make_sbc_generator/run_sbc_pipeline/plot_sbc_diagnostics/save_sbc_health_report)

# Level -> colour palettes (make_index(), hiermod_core.R). Named literals, no
# `div` dependency -- lets this file stay usable in marker-agnostic contexts
# (e.g. 1.X_Fig_hiermod.R) before a single hiermod_marker is picked. The
# actual div/idx codebooks (which DO need a marker) live in 0_INDEX.R,
# sourced by each model file. Named (not positional) so any script can
# safely subset a palette by label -- e.g. Model 9's Location-only-B/D
# subset uses fill_loc[c("B", "D")] rather than relying on level order.
fill_mg   <- c(Conventional = "#F28E2B", Organic = "#499894")
fill_loc  <- c(A = "#5DB63B", B = "#2C9EE3", C = "gold", D = "#F8A11C")
fill_cult <- c(Cortland = "#7DB16B", Liberty = "#45818E", Paulared = "#AB4F84",
               Honeycrisp = "#AE9FCB", Spartan = "#99CFE1")

# Contrast/Population are generic roles layered on top, not Management levels.
Management_palette <- c(fill_mg, Contrast = "#98494d", Population = "#895a92")

# Covariate effect panel labels/colours (Model 5+).
cov_labels <- c("Degree-hours", "Precipitation (72h)", "Seq. depth")
cov_pal    <- setNames(scales::hue_pal()(3), cov_labels)

# Variance-partition panel colours (Model 5+, variance_partition_*()'s own marginal decomposition,
# variance_partition_panels() in postcontrast_helpers.R). Management x Season stays one combined
# term -- splitting it into main effects + interaction produced a strongly negative "by margin"
# share for the interaction (see MDSTYCV_model.R's own comment).
Variance_partition_palette <- c(
  "Management x Season" = "#4C72B0",
  "Year"                = "#A95571",
  "Precipitation"       = "#7D702E",
  "Reads count"         = "#008170",
  "Degree-hours"        = "#785186",
  "Cultivar"            = "#B17259",
  "Location"            = "#B17259",
  "Tree"                = "#3D6936",
  "Residual"            = "grey50"
)

# Fold-change panel colours (Model 2+, contrast_plot_panels()'s ratio_stats=). Reference these names dynamically
# (names(Fold_change_palette)) in each script's own extra= list rather than retyping them -- a mismatch here silently
# breaks the combined ratio panel (it did once already).
Fold_change_palette <- c(
  "May fold difference"  = "#CC8FBB",
  "July fold difference" = "#7AAB32"
)

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
