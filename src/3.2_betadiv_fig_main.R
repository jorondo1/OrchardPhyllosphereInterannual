# PCoA (weighted UniFrac) figures, both barcodes, 2- and 3-year subsets.
# Main figure: points/ellipses coloured by Time x Management.
# Supplementary figures, same ordinations: coloured by Year; and Cultivar /
# Location highlighted one month at a time (panel A = May, B = July).

pacman::p_load(tidyverse, vegan, patchwork, update = FALSE)

source('src/0.0_Config.R')

supp_dir <- "out/manuscript/supp"
dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(230726)

betadiv <- readRDS('data/diversity_subsets.rds')

# Helper: scores as df, merge metadata
score_tibble <- function(pcoa, meta) {
  scores(pcoa, display = "sites") %>%
    as.data.frame() %>%
    rownames_to_column("Sample") %>%
    merge(meta, by = "Sample") %>%
    mutate(
      t_m = factor(
        paste0(Time, ", ", str_to_lower(Management)),
        levels = names(fill_timman)))
}

# Helper: PCoA of one subset's weighted UniFrac, as scores + % variance per axis
run_pcoa <- function(subset, subtitle) {
  pcoa <- capscale(subset$Dist$wuf ~ 1, distance = "wunifrac")
  list(
    scores   = score_tibble(pcoa, subset$Meta),
    var      = round(100 * pcoa$CA$eig / sum(pcoa$CA$eig), 2),
    subtitle = subtitle)
}

# Ordinations -------------------------------------------------------------------
# Computed once, shared by every figure below.

ordinations <- list(
  Bacteria_3y = run_pcoa(betadiv$Bacteria$y3, "Bacteria, 3-year samples"),
  Fungi_3y    = run_pcoa(betadiv$Fungi$y3,    "Fungi, 3-year samples"),
  Bacteria_2y = run_pcoa(betadiv$Bacteria$y2, "Bacteria, 2-year samples"),
  Fungi_2y    = run_pcoa(betadiv$Fungi$y2,    "Fungi, 2-year samples")
)

# Plot --------------------------------------------------------------------------

# Helper: one PCoA panel coloured by `colour_by`. Optional: legend breaks (to
# hide levels like "other"), a legend title, no ellipses.
pcoa_panel <- function(ord, colour_by, color_pal, fill_pal, limits,
                       breaks = limits, legend_title = NULL, ellipses = TRUE) {
  p <- ggplot(ord$scores, aes(x = MDS1, y = MDS2,
                              color = .data[[colour_by]], fill = .data[[colour_by]]))

  # show.legend = TRUE keeps keys for levels absent from this panel (e.g. no
  # 2022 in the 2-year subsets), so all panels' legends are identical and merge
  p <- p + geom_point(shape = 21, stroke = 0.2, size = 3, show.legend = TRUE)

  if (ellipses) p <- p + stat_ellipse(level = 0.95, geom = "polygon", alpha = .1,
                                      linewidth = 0.2, show.legend = TRUE)
  p +
    labs(
      x = paste0("PCo1 ", ord$var[1], " %"),
      y = paste0("PCo2 ", ord$var[2], " %"),
      subtitle = ord$subtitle) +

    scale_color_manual(values = color_pal, limits = limits, breaks = breaks, name = legend_title) +
    scale_fill_manual(values = fill_pal, limits = limits, breaks = breaks, name = legend_title)
}

# Helper: assembled 4-panel figure (Panel A-D: 3-year on top, 2-year below).
# Legend limits = levels present in any panel, so all four panels share one
# legend (e.g. the 2-year subsets have no 2022 samples).
pcoa_figure <- function(colour_by, color_pal, fill_pal) {
  present <- unique(unlist(map(ordinations, \(o) as.character(o$scores[[colour_by]]))))
  limits  <- intersect(names(fill_pal), present)

  panels <- map(ordinations, pcoa_panel, colour_by, color_pal, fill_pal, limits)

  (panels$Bacteria_3y + panels$Fungi_3y) / (panels$Bacteria_2y + panels$Fungi_2y) +
    plot_layout(guides = "collect") +
    plot_annotation(tag_levels = "A") &
    theme_pcoa &
    theme(legend.title = element_blank(),
          legend.position = 'bottom')
}

# Month-highlight palettes: config colours, with the non-focal month's samples
# ("other") de-emphasised in grey.
highlight_pal <- function(color_pal, fill_pal) {
  list(color = c(color_pal[names(color_pal) != "other"], other = "grey70"),
       fill  = fill_pal)
}
month_pals <- list(
  Cultivar = highlight_pal(color_cult, fill_cult),
  Location = highlight_pal(color_loc,  fill_loc))

# Helper: one ordination, `var` shown for one month's samples only; the other
# month's samples stay in grey ("other", drawn underneath, not in the legend).
# No title: months/subsets are identified in the figure caption.
month_panel <- function(ord, var, month, limits) {
  pals <- month_pals[[var]]
  ord$scores <- ord$scores %>%
    mutate(highlight = factor(if_else(Time == month, as.character(.data[[var]]), "other"),
                              levels = names(pals$fill))) %>%
    arrange(highlight != "other")
  ord$subtitle <- NULL

  pcoa_panel(ord, "highlight", pals$color, pals$fill, limits,
             breaks = setdiff(limits, "other"), legend_title = var, ellipses = FALSE)
}

# Helper: all 8 month-highlight ordinations of one barcode, May left / July right:
# A/B Cultivar 3Y, C/D Cultivar 2Y, E/F Location 3Y, G/H Location 2Y.
# Each row is its own patchwork with its legend collected to the right, so each
# legend sits beside (and lists only the levels of) its own plot pair.
month_figure <- function(barcode) {
  rows <- tribble(
    ~var,       ~subset,
    "Cultivar", "3y",
    "Cultivar", "2y",
    "Location", "3y",
    "Location", "2y")

  row_plots <- pmap(rows, \(var, subset) {
    ord    <- ordinations[[paste0(barcode, "_", subset)]]
    limits <- intersect(names(month_pals[[var]]$fill),
                        c(as.character(ord$scores[[var]]), "other"))
    month_panel(ord, var, "May", limits) + month_panel(ord, var, "July", limits) +
      plot_layout(guides = "collect")
  })

  wrap_plots(row_plots, ncol = 1) +
    plot_annotation(tag_levels = "A") &
    theme_pcoa &
    theme(legend.position = "right",
          legend.justification = "left",
          # tag top-aligned with the panel: with no titles, the top of the
          # "plot" region (which excludes the outer margin) is the panel top
          plot.tag.location = "plot",
          plot.tag.position = "topleft",
          plot.tag = element_text(vjust = 1, hjust = 0))
}

# Figures -----------------------------------------------------------------------

# Main figure: Time x Management
fig_main <- pcoa_figure("t_m", color_timman, fill_timman)

# Supplementary figure: Year
fig_year <- pcoa_figure("Year", color_year, fill_year)

# Supplementary figures: Cultivar / Location by month, one per barcode
fig_month <- map(c(Bacteria = "Bacteria", Fungi = "Fungi"), month_figure)

# Save --------------------------------------------------------------------------

ggsave("out/manuscript/2_PCoA_wUF.pdf", fig_main,
       bg = 'white', width = 2800, height = 2900, dpi = 300, units = "px")
ggsave(file.path(supp_dir, "PCoA_wUF_year.pdf"), fig_year,
       bg = 'white', width = 2800, height = 2900, dpi = 300, units = "px")

iwalk(fig_month, \(p, barcode)
  ggsave(file.path(supp_dir, sprintf("PCoA_wUF_cultivar_location_%s.pdf", barcode)), p,
         bg = 'white', width = 2800, height = 4500, dpi = 300, units = "px"))
