source('src/1.0_Community_barplots_config.R')

# 1. One tibble: Kingdom x Dataset x Time x SiteYear x aggTaxo ------------------

family_summary <- family_class %>%
  map(\(k) map(k, bind_rows) %>% 
        bind_rows(.id = "Dataset")) %>% # Time is already a column
  bind_rows(.id = "Kingdom") %>%
  group_by(Kingdom, Dataset, Time, Sample, Code, Year, aggTaxo) %>%
  summarise(relAb = sum(relAb, na.rm = TRUE), .groups = "drop") %>%
  group_by(Kingdom, Dataset, Time, Code, Year, aggTaxo) %>%
  summarise(
    mean_rel_abund = mean(relAb, na.rm = TRUE),
    sd_rel_abund   = sd(relAb, na.rm = TRUE),
    n              = sum(!is.na(relAb)),
    se_rel_abund   = sd_rel_abund / sqrt(n),
    .groups = "drop"
  ) %>%
  mutate(
    Dataset = factor(Dataset, c("3-year", "2-year")),
    Time = factor(Time, c("May", "July")))

# 2. One 2x2 faceted plot per Kingdom (rows = Dataset, cols = Time) -------------

legend_titles <- c(Bacteria = "Bacterial families", Fungi = "Fungal families")

# Stacked family barplot for one barcode, families in fixed order
plot_family_bar <- function(df) {
  
  kingdom <- df$Kingdom[1]
  fam_order <- family_order[[kingdom]]
  
  df %>%
    mutate(aggTaxo = factor(aggTaxo, levels = fam_order),
           Dataset = fct_recode(
             Dataset, 
             !!!c(`3-year samples` = '3-year', `2-year samples` = '2-year'))) %>%
    ggplot(aes(x = Year, y = mean_rel_abund, fill = aggTaxo)) +
    geom_col(
      width = 0.8, color = bar_outline_color, 
      linewidth = bar_outline_width,
      show.legend = TRUE) + # draw keys for families absent from this panel
    ggh4x::facet_nested(
      Dataset ~ Time + Code,
      scales = "free_x", space = "free_x") +
    scale_fill_manual(
      values = fill_fam, 
      breaks = fam_order,
      limits = fam_order, # identical legends in 3-year and 2-year plots -> collected into one
      labels = italicize_label(fam_order), drop = FALSE) +
    scale_y_continuous(
      breaks = seq(0, 100, by = 20),
      labels = function(x) paste0(x, "%")) +
    coord_cartesian(ylim = c(0, 100), expand = FALSE) +
    labs(x = NULL, y = "Mean relative abundance", fill = legend_titles[[kingdom]]) +
    guides(fill = guide_legend(ncol = 1)) +
    theme(
      axis.ticks = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.title = element_text(face = "bold"),
      legend.text = ggtext::element_markdown(size = 12),
      legend.key.height = grid::unit(0.6, "cm"),
      legend.key.width = grid::unit(0.6, "cm"),
      legend.key.spacing.y = grid::unit(legend_entry_spacing, "cm"),
      strip.text = element_text(colour = 'grey20', size = month_title_size),
      strip.background = element_rect(
        fill = 'grey95',  linewidth = bar_outline_width, color = bar_outline_color),
      panel.border = element_rect(linewidth = 0),
      panel.spacing = unit(0.2, "lines"),
      panel.grid = element_blank()
      
    )
}

# One plot per Kingdom x Dataset, so each only has its own sites on the x axis
lst <- split(family_summary, ~ Kingdom + Dataset, sep = "_")

plots <- map2(
  lst,
  seq_along(lst),
  ~ plot_family_bar(.x) + labs(subtitle = LETTERS[.y])
)
# 3. Combined plot ----------------------------------------------------------------
# Per Kingdom: 3-year over 2-year, identical legends collected into one.
# Fungi y axis dropped, as in the original.

kingdom_panel <- function(kingdom, ...) {
  plots[paste0(kingdom, c("_3-year", "_2-year"))] %>%
    map(\(p) p + theme(...)) %>%
    wrap_plots(ncol = 1, guides = "collect")
}

family_panel_All <- kingdom_panel("Bacteria") |
  kingdom_panel("Fungi", axis.title.y = element_blank(), axis.text.y = element_blank())

family_panel_All

ggsave('out/manuscript/1_Comm_Barplot_SiteYear.pdf', family_panel_All,
       bg = 'white', width = 1700, height = 1200, dpi = 120, units = 'px')
