# Config file
# colour/shape/labels mappings
# package conflict resolution
# Species cluster names overrides,

source('src/0.0_ggplot_themes.R')

## Cross-package function-name conflicts ---------------------------------------
# Resolves dplyr/base functions masked by other loaded packages 
pacman::p_load(conflicted, update = FALSE)
conflicted::conflicts_prefer(
  base::intersect,
  base::match,
  bayesplot::rhat,
  stats::sd,
  purrr::map,
  stats::var,
  dplyr::filter,  
  dplyr::select,  
  dplyr::rename,  
  dplyr::slice,   
  dplyr::combine, 
  dplyr::desc,    
  dplyr::count,   
  dplyr::first,   
  dplyr::mutate,  
  dplyr::arrange, .quiet = TRUE)

# Fungi Species clusters overrides, determined from BLAST
fungi_label_overrides <- c( 
  "NA_sp_clust_4"                              = "Cladosporium_4*",
  "Ascomycota_sp_clust_5"                      = "Didymellaceae_5*",
  "Pleosporales_gen_Incertae_sedis_sp_clust_7" = "Alternaria_7*",
  "Ascomycota_sp_clust_14"                     = "Melanommataceae_14*",
  "NA_sp_clust_15"                             = "Filobasidium_15*",
  "Helotiales_sp_clust_17"                     = "Lemonniera_17*"
)

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
# fill_mg: also used by hiermod; "other" for consumers with a 3rd level

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
fill_site  <- c(A = "#5DB63B", B1 = "#2C9EE3", B2 = "#096EA4",
                C = "gold", D1 = "#FCB452", D2 = "#F8A11C") # colors adopted from Sophie

## Location -----------------------------------------------------------------------------
color_loc <- c(A = "black", B = "black", C = "black", D = "black", other = "black")
# shape_loc: distinct shapes, for location told apart by shape alone
shape_loc <- c(A = 21, B = 22, C = 23, D = 24, other = 25)
fill_loc  <- c(A = "#5DB63B", B = "#1B86C4", C = "gold", D = "#FAAB37", other = "grey90")

## Cultivar -----------------------------------------------------------------------------
# Level order = hiermod's idx$Cv (Cv index 1..5)
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
# Keyed by raw column name; consumers re-key to their labels (e.g. hiermod's cov_pal)
# Colours: ColorBrewer Set1
env_var_colors <- c(
  mean_temp  = "#e41a1c",
  precip_72h = "#377eb8",
  deg_h      = "#4daf4a",
  seq_depth  = "#984ea3"
)

## ps_objects_full.rds display labels ----------------------------------------------------
# ps.ls keys ($-safe) -> readable labels for tables / facet strips
ps_dataset_labels <- c(
  Bacteria = "Bacteria", Fungi = "Fungi",
  Bacteria_3y = "Bacteria (3-year subset)", Bacteria_2y = "Bacteria (2-year subset)",
  Fungi_3y = "Fungi (3-year subset)", Fungi_2y = "Fungi (2-year subset)"
)

# Dataset column values -> longer facet-strip labels
dataset_facet_labels <- c(`2-year` = "2-year dataset", `3-year` = "3-year dataset")
