# Shared preamble for every hiermod script (ITS and 16S). Set
# hiermod_marker <- "ITS" or "16S" before sourcing this file.

pacman::p_load(rethinking, tidyverse, bayesm, bayesplot, ggridges, magrittr, patchwork, rlang, scales, posterior)
source('src/utils/hiermod_core.R')
source('src/utils/saver_functions.R')
source('src/utils/postcontrast_helpers.R')
source('src/utils/predictive_checks.R')
source('src/utils/sbc_workflow.R') # SBC-package-based workflow (draw_true/SBC_backend_ulam/make_sbc_generator/run_sbc_pipeline/plot_sbc_diagnostics/save_sbc_health_report)

div_all <- read_rds('data/diversity_data.rds')
div <- if (hiermod_marker == "ITS") div_all$Fungi$alpha else div_all$Bacteria$alpha

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
  Lo = make_index(div$Location, palette = fill_loc),
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

# Covariate effect panel labels/colours (Model 5 (MDSYC) onward) -- same
# 3 labels/colours every fit script built inline for its own
# fit_covariate_effects panel.
cov_labels <- c("Degree-hours", "Precipitation (72h)", "Seq. depth")
cov_pal    <- setNames(scales::hue_pal()(3), cov_labels)

# Variance-partition panel palette (Model 5 (MDSYC) onward) -- one colour
# per fixed-effect group (variance_partition_*()'s own sequential
# decomposition, see MDSYC_model.R), plus Tree (random effect, Model 7
# only) and Residual. Shared across every model's own variance-partition
# plot so the same group always gets the same colour; a model that doesn't
# have a given group (e.g. Models 5/6 have no Tree) just never uses that
# name.
Variance_partition_palette <- c(
  "Management x Season" = "#4C72B0",
  "Year"                = "#55A868",
  "Covariates"          = "#8172B2",
  "Cultivar"            = "#CCB974",
  "Tree"                = "#DD8452",
  "Residual"            = "grey50"
)

# ---- Standardized control covariates (Model 7 (MDLSYC) onward) ----
# deg_h_z/precip_72h_z are centered WITHIN each Season (mean-subtracted per
# Time level), not globally -- Season (Mo) is already a predictor, and
# degree-hours/precipitation are strongly collinear with it in the real
# data (r=0.735/0.547, see TODO.md). Global centering would let
# gamma/s_conv/gap_shift (the May->July contrast terms) report the
# seasonal shift net of the portion explained by b_deg/b_precip -- "the
# seasonal change if temperature had been constant," which never happens
# and isn't the estimand of interest. Within-Season centering makes each
# covariate orthogonal to Mo by construction (so it still absorbs
# day-to-day sampling-date weather jitter within a season), while the real
# between-season temperature/precipitation difference flows entirely into
# the season term, where it belongs. seq_depth_z below stays globally
# centered on purpose: its correlation with Management is a technical
# sequencing-depth artifact we want removed from the Management contrast,
# not a substantive quantity to preserve.
# *_season_mean are named vectors (keyed by Time level, not a single
# scalar) so any script needing the raw scale back can invert the z-score
# per season; *_sd is the pooled within-season residual SD.
deg_h_season_mean      <- tapply(div$deg_h,      div$Time, mean)
precip_72h_season_mean <- tapply(div$precip_72h, div$Time, mean)

div$deg_h_z      <- (div$deg_h      - deg_h_season_mean[div$Time])
div$precip_72h_z <- (div$precip_72h - precip_72h_season_mean[div$Time])

deg_h_sd      <- sd(div$deg_h_z);      div$deg_h_z      <- div$deg_h_z      / deg_h_sd
precip_72h_sd <- sd(div$precip_72h_z); div$precip_72h_z <- div$precip_72h_z / precip_72h_sd

# Seq_depth (raw pre-rarefaction read count) is right-skewed and spans
# ~20x (see MODEL_HISTORY.md Model 7) -- logged first so a linear slope
# matches a saturating detection-effort effect, then standardized like the
# weather covariates above.
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
