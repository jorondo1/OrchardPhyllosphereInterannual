# 0_INDEX.R -- marker-specific data + category<->index codebooks. Sourced by
# every model file (Models/*.R), not by 0_SETUP.R itself -- 0_SETUP.R stays
# usable in marker-agnostic contexts (e.g. 1.X_Fig_hiermod.R) that haven't
# picked a single hiermod_marker yet. Self-sources 0_SETUP.R so a model file
# is fully runnable on its own regardless of what the calling script already
# sourced (re-sourcing 0_SETUP.R is idempotent -- cheap package/theme/palette
# rebinds, no side effects).
source('src/hiermod/0_SETUP.R')

div_all <- read_rds('data/diversity_data.rds')
div <- if (hiermod_marker == "ITS") div_all$Fungi$alpha else div_all$Bacteria$alpha

# Category <-> index codebooks (make_index(), hiermod_core.R) -- shared here so the index<->category mapping can't drift between scripts.
idx <- list(
  Mg = make_index(div$Management, levels = names(fill_mg),   palette = fill_mg),
  Lo = make_index(div$Location,   levels = names(fill_loc),  palette = fill_loc),
  Tr = make_index(div$Tree_id),
  Cv = make_index(div$Cultivar,   levels = names(fill_cult), palette = fill_cult),
  Mo = make_index(div$Time, levels = c("May", "July")),
  Yr = make_index(div$Year, levels = c("2022", "2023", "2024"))
)
