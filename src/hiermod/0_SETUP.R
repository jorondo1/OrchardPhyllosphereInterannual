# Shared config for every hiermod script: packages, helpers, palettes

pacman::p_load(update = FALSE,
  rethinking, tidyverse, bayesm, bayesplot, ggridges, magrittr, patchwork, rlang, scales, posterior, kableExtra)
source('src/utils/hiermod_core.R')
source('src/utils/saver_functions.R')
source('src/utils/postcontrast_helpers.R')
source('src/utils/predictive_checks.R')
source('src/0.0_Config.R')

# Contrast panels
Management_palette <- c(fill_mg[c("Conventional", "Organic")], Contrast = "#98494d", Population = "#895a92")

# Covariate effect panels (models 5+); colours from env_var_colors (0.0_Config.R)
cov_labels <- c("Degree-hours", "Precipitation (72h)", "Seq. depth")
cov_pal    <- setNames(unname(env_var_colors[c("deg_h", "precip_72h", "seq_depth")]), cov_labels)

# Variance-partition panels
Variance_partition_palette <- c(
  "Management"          = "#2E4A7D",
  "Season"              = "#C9A227",
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

# Fold-change panels
Fold_change_palette <- c(
  "May fold difference"  = "#CC8FBB",
  "July fold difference" = "#7AAB32"
)