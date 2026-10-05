# - needs hiermod_marker ("16S"/"ITS") set beforehand
# - sourced by every Models/*.R; sources 0_SETUP.R itself (safe to re-source)
source('src/hiermod/0_SETUP.R')

div_all <- read_rds('data/diversity_data.rds')
div <- if (hiermod_marker == "ITS") div_all$Fungi$alpha else div_all$Bacteria$alpha

# Codebooks (make_index(), hiermod_core.R), one shared mapping for all scripts
idx <- list(
  Mg = make_index(div$Management, levels = names(fill_mg),   palette = fill_mg),
  Lo = make_index(div$Location,   levels = names(fill_loc),  palette = fill_loc),
  Tr = make_index(div$Tree_id),
  Cv = make_index(div$Cultivar,   levels = names(fill_cult), palette = fill_cult),
  Mo = make_index(div$Time, levels = c("May", "July")),
  Yr = make_index(div$Year, levels = c("2022", "2023", "2024"))
)
