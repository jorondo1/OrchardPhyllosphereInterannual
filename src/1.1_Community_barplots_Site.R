source('src/1.0_Community_barplots_config.R')

# 1. Create May and July plot data tables --------------------------

# function to Summarize abundance by Site

summarize_family_plot <- function(df) {
  
  df %>%
    group_by(Sample, Code, aggTaxo) %>%
    summarise(relAb = sum(relAb, na.rm = TRUE), .groups = "drop") %>%
    group_by(Code, aggTaxo) %>%
    summarise(
      mean_rel_abund = mean(relAb, na.rm = TRUE),
      sd_rel_abund = sd(relAb, na.rm = TRUE),
      n = sum(!is.na(relAb)),
      se_rel_abund = sd_rel_abund / sqrt(n),
      .groups = "drop"
    )
}

family_summary <- map_depth(family_class, 3, summarize_family_plot)

# 7. Plot function -----------------------------------------

plot_family_bar <- function(df, family_order, month) {
  
  legend_labels <- italicize_label(family_order)
  names(legend_labels) <- family_order
  legend_title <- "Family"
  
  df %>%
    mutate(aggTaxo = factor(aggTaxo, levels = family_order)) %>% 
    ggplot(aes(x = Code, y = mean_rel_abund, fill = aggTaxo)) +
    
    geom_col(
      width = 0.8,
      color = bar_outline_color,
      linewidth = bar_outline_width
    ) +
    
    scale_fill_manual(
      values = fill_fam,
      breaks = family_order,
      labels = legend_labels,
      drop = FALSE
    ) +
    
    scale_y_continuous(
      breaks = seq(0, 100, by = 20),
      labels = function(x) paste0(x, "%"),
      expand = expansion(mult = c(0, 0))
    ) +
    
    coord_cartesian(ylim = c(0, 100), expand = FALSE) +
    
    labs(x = "Site", 
         y = "Mean relative abundance",
         fill = legend_title,
         title = month) +
    
    guides(fill = guide_legend(ncol = 1, byrow = FALSE))  +
    
    theme(
      legend.text = ggtext::element_markdown(),
      legend.key.height = grid::unit(0.45, "cm"),
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, size = month_title_size),
      panel.border = element_rect(colour = 'black'),
      axis.ticks = element_blank()
    )
}

## 7.1. Create plots for every Kingdom > Duration > Month combination -----

family_bar <- imap(family_summary, function(duration_list, kingdom) {
  order <- family_order[[kingdom]]
  imap(duration_list, function(month_list, duration) {
    imap(month_list, function(df, month) {
      plot_family_bar(df, family_order = order, month = month)
    })
  })
})

# 8. Legends ----------------

make_shared_family_legend <- function(family_order, legend_title) {
  
  legend_plot <- ggplot(
    tibble(
      aggTaxo = factor(family_order, levels = family_order),
      x = seq_along(family_order),
      y = 1
    ),
    aes(x = x, y = y, fill = aggTaxo)
  ) +
    
    geom_col(
      show.legend = TRUE,
      color = bar_outline_color,
      linewidth = bar_outline_width
    ) +
    
    scale_fill_manual(
      values = fill_fam,
      breaks = family_order,
      labels = italicize_label(family_order),
      drop = FALSE
    ) +
    
    labs(fill = legend_title) +
    
    guides(fill = guide_legend(ncol = 1, byrow = FALSE)) +
    
    theme_void() +
    
    theme(
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      legend.text = ggtext::element_markdown(size = 12),
      legend.key.height = grid::unit(0.6, "cm"),
      legend.key.width  = grid::unit(0.6, "cm"),
      legend.key.spacing.y = grid::unit(legend_entry_spacing, "cm")
    )
  
  shared_legend <- ggpubr::get_legend(legend_plot)
  
  return(shared_legend)
}


## 8.1. Create shared legend per Kingdom -------
legend_titles <- list(Bacteria = "Bacterial families", Fungi = "Fungal families")

shared_legends <- imap(family_order, ~ make_shared_family_legend(.x, legend_titles[[.y]]))

# Measure each legend's NATURAL width
legend_widths <- map(shared_legends, ~ grid::unit(1, "grobwidth", .x) + grid::unit(0.3, "cm"))


# 9. Formatting ---------------------------------

# Year row labels
make_year_label <- function(label) {
  
  ggplot() +
    annotate("text", x = 0.5, y = 0.4, label = label,  size = year_label_size) +
    xlim(0, 1) + ylim(0, 1) + theme_void()
}

label_3year <- make_year_label("3-year")
label_2year <- make_year_label("2-year")

## 9.1. Bacteria --------------------

# Per-position theme overrides (row = Duration, col = Month). Manual because
# asymmetric depending on plot section
theme_overrides <- list(
  Bacteria = list(
    `3-year` = list(
      May  = theme(
        axis.title.x = element_blank()),
      July = theme(
        axis.title.x = element_blank(),
        axis.title.y = element_blank(), 
        axis.text.y = element_blank())
    ),
    `2-year` = list(
      May  = theme(),
      July = theme(
        axis.title.y = element_blank(), 
        axis.text.y = element_blank())
    )
  ),
  Fungi = list(
    `3-year` = list(
      May  = theme(
        axis.title = element_blank(),
        axis.text.y = element_blank()),
      July = theme(
        axis.title.x = element_blank(),
        axis.title.y = element_blank(), 
        axis.text.y = element_blank())
    ),
    `2-year` = list(
      May  = theme(
        axis.title.y = element_blank(), 
        axis.text.y = element_blank()),
      July = theme(
        axis.title.y = element_blank(), 
        axis.text.y = element_blank())
    )
  )
)

family_bar_themed <- imap(family_bar, function(duration_list, kingdom) {
  imap(duration_list, function(month_list, duration) {
    imap(month_list, function(plot, month) {
      plot + theme_overrides[[kingdom]][[duration]][[month]] +
        theme(legend.position = "none")
    })
  })
})

# Assemble each Kingdom's 2x2 (Duration x Month) panel
kingdom_graphs <- map(family_bar_themed, function(d) {
  (
    label_3year /
      (d$`3-year`$May | d$`3-year`$July) /
      label_2year /
      (d$`2-year`$May | d$`2-year`$July)
  ) +
    plot_layout(heights = c(0.1, 1.1, 0.1, 1.1))
})

# Add each Kingdom's shared legend
family_panel_list <- imap(kingdom_graphs, function(g, kingdom) {
  (g | wrap_elements(full = shared_legends[[kingdom]])) +
    plot_layout(widths = grid::unit.c(grid::unit(1, "null"), legend_widths[[kingdom]]))
})
family_panel_Bacteria <- family_panel_list$Bacteria
family_panel_Fungi    <- family_panel_list$Fungi

# 10. Combined plot -----------------

family_panel_All <- (
  kingdom_graphs$Bacteria |
    wrap_elements(full = shared_legends$Bacteria) |
    kingdom_graphs$Fungi |
    wrap_elements(full = shared_legends$Fungi)
) +
  plot_layout(
    widths = grid::unit.c(
      grid::unit(1, "null"),
      legend_widths$Bacteria,
      grid::unit(1, "null"),
      legend_widths$Fungi))

family_panel_All

ggsave('out/manuscript/1_Comm_Barplot.pdf',
       bg = 'white', width = 1700, height = 1200, dpi = 120, units = 'px')

# 11. Export abundance table -----------------------------
# Add Year / Month / Domain identifiers to each summary table

add_ids <- function(df, year_group, month, domain) {
  df %>%
    mutate(
      Year_Group = year_group,
      Month = month,
      Domain = domain,
      .before = 1
    )
}

# build table

abundance_table_all <- imap_dfr(family_summary, function(duration_list, kingdom) {
  imap_dfr(duration_list, function(month_list, duration) {
    imap_dfr(month_list, function(df, month) {
      add_ids(df, str_to_title(gsub("-", " ", duration)), month, kingdom)
    })
  })
}) %>%
  
  rename(
    Family = aggTaxo,
    Mean_Percent_Abundance = mean_rel_abund
  ) %>%
  
  select(
    Domain, Year_Group, Month, Code, Family, Mean_Percent_Abundance
  ) %>%
  
  arrange(Domain, Year_Group, Month, Code, desc(Mean_Percent_Abundance))

abundance_table_all <- abundance_table_all %>%
  mutate(
    Mean_Percent_Abundance = round(Mean_Percent_Abundance, 2)
  )

abundance_table_all

# Save as CSV
write.csv(
  abundance_table_all,
  "family_abundance_summary.csv",
  row.names = FALSE
)

