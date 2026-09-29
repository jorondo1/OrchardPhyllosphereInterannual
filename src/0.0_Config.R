# Config: colour/shape/label palettes 

source('src/0.0_ggplot_themes.R')


## Cross-package function-name conflicts ---------------------------------------
# Resolves dplyr/base functions masked by other loaded packages (e.g.
# MASS::select, stats::filter) .quiet = TRUE since this is deliberate
pacman::p_load(conflicted, update = FALSE)
conflicted::conflicts_prefer(base::intersect,  .quiet = TRUE)
conflicted::conflicts_prefer(dplyr::filter,    .quiet = TRUE)
conflicted::conflicts_prefer(dplyr::select,    .quiet = TRUE)
conflicted::conflicts_prefer(dplyr::rename,    .quiet = TRUE)
conflicted::conflicts_prefer(dplyr::slice,     .quiet = TRUE)
conflicted::conflicts_prefer(dplyr::combine,   .quiet = TRUE)
conflicted::conflicts_prefer(dplyr::desc,      .quiet = TRUE)
conflicted::conflicts_prefer(dplyr::count,     .quiet = TRUE)
conflicted::conflicts_prefer(dplyr::first,     .quiet = TRUE)
conflicted::conflicts_prefer(dplyr::mutate,    .quiet = TRUE)
conflicted::conflicts_prefer(dplyr::arrange,   .quiet = TRUE)

## Year -----------------------------------------------------------------------
color_year <- c(`2022` = "black", `2023` = "black", `2024` = "black")
shape_year <- c(`2022` = 21, `2023` = 21, `2024` = 21)
fill_year  <- c(`2022` = "#3B839E", `2023` = "#BB4D3C", `2024` = "#9AB88A") # first two colors adopted from Sophie

## Time (May/July) --------------------------------------------------------------
color_time  <- c(May = "black", July = "black", other = "grey80")
shape_time  <- c(May = 21, July = 21, other = 21)
shape_time2 <- c(May = 21, July = 22) # 2-level variant (no "other"), for plots covering only May/July
fill_time   <- c(May = "#CC8FBB", July = "#7AAB32", other = "grey90") # colors adopted from Sophie

## Management ---------------------------------------------------------------------
# fill_mg matches src/hiermod/0_SETUP.R's own values; "other" fallback kept
# for consumers with a 3rd level (hiermod's own 2-level Mg ignores it).

color_mg <- c(Conventional = "black", Organic = "black", other = "grey80")
shape_mg <- c(Conventional = 21, Organic = 21, other = 21)
fill_mg  <- c(Conventional = "#F28E2B", Organic = "#499894", other = "grey90")

## Time x Management -----------------------------------------------------------------

color_timman <- c(
  `May, conventional`  = "black", 
  `May, organic`       = "black",
  `July, conventional` = "black", 
  `July, organic`      = "black", 
  other = "grey80")

fill_timman  <- c(
  `May, conventional` = "#D89C60", 
  `May, organic`      = "#8F83AA",
  `July, conventional`= "#AFAA1B", 
  `July, organic`      = "#669165", 
  other = "grey80")

# https://colorkit.co/color-mixer/?mix=7aab32-e4a804&ratios=1-1&space=rgb

## Site -------------------------------------------------------------------------------

color_site <- c(A = "black", B1 = "black", B2 = "black", C = "black", D1 = "black", D2 = "black")
shape_site <- c(A = 22, B1 = 21, B2 = 21, C = 22, D1 = 23, D2 = 23)
fill_site  <- c(A = "#D17913", B1 = "#F3A44A", B2 = "#096EA4",
                C = "#1C9EE4", D1 = "#FFC787", D2 = "#89CEF3") # colors adopted from Sophie

## Location -----------------------------------------------------------------------------
color_loc <- c(A = "black", B = "black", C = "black", D = "black", other = "black")
# shape_loc: distinct shapes (unlike most other shape_* palettes, which are
# uniform placeholders) -- needed wherever Location has to be told apart by
# shape alone, e.g. src/2.5_betadiv_envfit_fig.R's Time-coloured panel.
shape_loc <- c(A = 21, B = 22, C = 23, D = 24, other = 25)
fill_loc  <- c(A = "#5DB63B", B = "#2C9EE3", C = "gold", D = "#F8A11C", other = "grey90")

## Cultivar -----------------------------------------------------------------------------
# Canonical level order: Cortland, Liberty, Paulared, Honeycrisp, Spartan --
# matches src/hiermod/0_SETUP.R's idx$Cv order (Cv index 1..5).
color_cult <- c(Cortland = "black", Liberty = "black", Paulared = "black",
                Honeycrisp = "black", Spartan = "black", other = "black")
shape_cult <- c(Cortland = 21, Liberty = 21, Paulared = 21,
                Honeycrisp = 21, Spartan = 21, other = 21)
fill_cult  <- c(Cortland = "#7DB16B", Liberty = "#45818E", Paulared = "#AB4F84",
                Honeycrisp = "#AE9FCB", Spartan = "#99CFE1", other = "grey90") # first three colors adopted from Sophie

## Gradient (continuous colour scales) -------------------------------------------------
grad_low <- "red"
grad_mid <- "white"
grad_hi  <- "blue"
grad_na  <- "white"


## Environmental/covariate variables ------------------------------------------------------
# Shared by src/1.4_Fig_BetaDiv_Envfit.R (raw mean_temp/precip_72h/deg_h vs.
# PCoA axes) and src/hiermod/0_SETUP.R's cov_pal (same deg_h/precip_72h plus
# seq_depth instead of mean_temp, as z-scored regression coefficients).
# Keyed by raw column name; each consumer re-keys to its own display labels
# (see cov_pal). Colours: ColorBrewer Set1.
env_var_colors <- c(
  mean_temp  = "#e41a1c",
  precip_72h = "#377eb8",
  deg_h      = "#4daf4a",
  seq_depth  = "#984ea3"
)

## ps_objects_full.rds display labels ----------------------------------------------------
# The 6-item ps.ls (src/0.3.2.Metadata_phyloseq.R) uses $-safe identifiers as
# its keys (Bacteria/Fungi/Bacteria_3y/Bacteria_2y/Fungi_3y/Fungi_2y) so they
# never need backtick-quoting -- this maps them to a human-readable label for
# kable tables / plot facet strips.
ps_dataset_labels <- c(
  Bacteria = "Bacteria", Fungi = "Fungi",
  Bacteria_3y = "Bacteria (3-year subset)", Bacteria_2y = "Bacteria (2-year subset)",
  Fungi_3y = "Fungi (3-year subset)", Fungi_2y = "Fungi (2-year subset)"
)

# Same idea for the Dataset column itself (diversity_data.rds / ps sample
# data: "2-year"/"3-year") when a longer facet-strip label is wanted instead
# of the raw filtering value.
dataset_facet_labels <- c(`2-year` = "2-year dataset", `3-year` = "3-year dataset")
