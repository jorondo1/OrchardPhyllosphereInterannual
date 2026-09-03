# Shared preamble for every hiermod ITS script

pacman::p_load(rethinking, tidyverse, bayesm, bayesplot, ggridges, magrittr, patchwork, rlang, scales, posterior)
source('~/Repos/orchardPhyllosphere2/src/hiermod/utils/hiermod_core.R')
source('~/Repos/orchardPhyllosphere2/src/hiermod/utils/saver_functions.R')
source('~/Repos/orchardPhyllosphere2/src/hiermod/utils/postcontrast_helpers.R')
source('~/Repos/orchardPhyllosphere2/src/hiermod/utils/predictive_checks.R')
source('~/Repos/orchardPhyllosphere2/src/hiermod/utils/sbc_helpers.R')

div.ITS <- read_rds('~/Repos/orchardPhyllosphere2/data/diversity_data.rds')
div <- div.ITS$Fungi$alpha

message("Make sure to define hiermod_out_dir <- out/hiermod/<...>")

# ---- Category <-> index codebooks (see make_index() in hiermod_core.R) ----
# Defined once, here. Every Mg/Lo/Tr/Cv conversion in every model script --
# both building a model's data list (to_index) and recovering labels for
# postpred checks (to_label) -- goes through these, so the index<->category
# mapping can't drift between model scripts (it previously did once: one
# model's data list had Conventional/Organic flipped relative to the rest).
# palette= carries each index's own colour scheme -- see $palette/$palette_n()
# in make_index() -- so scripts stop hand-building scales::hue_pal() palettes
# inline (as every variance-component panel used to).

fill_loc  <- c("#5DB63B", "#2C9EE3", "gold", "#F8A11C")
fill_cult <- c("#7DB16B", "#45818E", "#AB4F84", "#AE9FCB", "#99CFE1")

idx <- list(
  Mg = make_index(div$Management, levels = c("Conventional", "Organic"),
                   palette = c(Conventional = "#F28E2B", Organic = "#499894")),
  Lo = make_index(div$Location, levels = c("Saint-Benoît", "Compton", "Windsor", "Milton"),
                   palette = fill_loc),
  Tr = make_index(div$Tree_id),
  Cv = make_index(div$Cultivar, levels = c("Cortland", "Liberty", "Paulared", "Honeycrisp", "Spartan"),
                   palette = fill_cult),
  Mo = make_index(div$Time, levels = c("May", "July")),
  Yr = make_index(div$Year, levels = c("2022", "2023", "2024"))
)

# Colours -- Management_palette built FROM idx$Mg$palette (same hex values as
# before) so every existing `group_pal = Management_palette` call site keeps
# working unchanged; Contrast/Population are generic roles, not Management
# levels, so they're layered on top rather than folded into idx$Mg itself.
Management_palette <- c(idx$Mg$palette, Contrast = "#98494d", Population = "#895a92")

# ---- Standardized weather covariates (Model 6 (MDLSY) onward) ----
# Keep the raw mean/sd as named scalars (not just baked into `div`) so any
# script that needs to report/predict on the raw deg_h/precip_72h scale can
# invert the z-score consistently, rather than re-deriving mean/sd locally.
deg_h_mean      <- mean(div$deg_h);      deg_h_sd      <- sd(div$deg_h)
precip_72h_mean <- mean(div$precip_72h); precip_72h_sd <- sd(div$precip_72h)

div$deg_h_z      <- (div$deg_h      - deg_h_mean)      / deg_h_sd
div$precip_72h_z <- (div$precip_72h - precip_72h_mean) / precip_72h_sd

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
