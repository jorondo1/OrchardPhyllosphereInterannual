# ==========================================================
# FAMILY RELATIVE ABUNDANCE
# ==========================================================

# ==========================================================
# 1. Load packages
# ==========================================================

library(pacman)

p_load(
  tidyverse,
  RColorBrewer,
  phyloseq,
  patchwork,
  magrittr,
  ggpubr,
  ggtext  
)

# ==========================================================
# 2. Load datasets
# ==========================================================

Years3 <- readRDS("Years3.rds")
Years2 <- readRDS("Years2.rds")

# ==========================================================
# 3. Relative-abundance cutoff
# ==========================================================

family_cutoff <- 1

# ==========================================================
# 4. Prepare Family relative-abundance data
# ==========================================================

prepare_family_data <- function(ps) {
  
  ps_rel <- transform_sample_counts(
    ps,
    function(x) {
      (x / sum(x)) * 100
    }
  )
  
  df <- psmelt(ps_rel)
  
  df <- df %>%
    mutate(
      Family_clean = case_when(
        is.na(Family) ~ "Unclassified",
        Family == "" ~ "Unclassified",
        Family == "NA" ~ "Unclassified",
        grepl(
          "Incertae|Unclassified|uncultured|unknown",
          Family,
          ignore.case = TRUE
        ) ~ "Unclassified",
        TRUE ~ as.character(Family)
      )
    )
  
  df <- df %>%
    group_by(
      Sample,
      Code,
      Time,
      Year,
      Family_clean
    ) %>%
    summarise(
      relAb = sum(Abundance, na.rm = TRUE),
      .groups = "drop"
    )
  
  return(df)
}

# ==========================================================
# 5. Prepare all datasets
# ==========================================================

family_data_Years3_Bacteria <- prepare_family_data(Years3$Bacteria)
family_data_Years3_Fungi    <- prepare_family_data(Years3$Fungi)
family_data_Years2_Bacteria <- prepare_family_data(Years2$Bacteria)
family_data_Years2_Fungi    <- prepare_family_data(Years2$Fungi)


# ==========================================================
# 6. Find Families with average relative abundance >= 1%
# ==========================================================

get_month_families <- function(df, time_value, cutoff = 1) {
  
  abundant <- df %>%
    filter(Time == time_value) %>%
    filter(Family_clean != "Unclassified") %>%
    group_by(Sample, Family_clean) %>%
    summarise(relAb = sum(relAb, na.rm = FALSE), .groups = "drop") %>%
    group_by(Family_clean) %>%
    summarise(mean_rel_abund = mean(relAb, na.rm = FALSE), .groups = "drop") %>%
    filter(mean_rel_abund >= cutoff) %>%
    pull(Family_clean)
  
  return(abundant)
}


# ==========================================================
# 7. Make a SEPARATE >= 1% Family list for every graph
# ==========================================================

families_Years3_May_Bacteria  <- get_month_families(family_data_Years3_Bacteria, "May",  cutoff = family_cutoff)
families_Years3_July_Bacteria <- get_month_families(family_data_Years3_Bacteria, "July", cutoff = family_cutoff)

families_Years3_May_Fungi  <- get_month_families(family_data_Years3_Fungi, "May",  cutoff = family_cutoff)
families_Years3_July_Fungi <- get_month_families(family_data_Years3_Fungi, "July", cutoff = family_cutoff)

families_Years2_May_Bacteria  <- get_month_families(family_data_Years2_Bacteria, "May",  cutoff = family_cutoff)
families_Years2_July_Bacteria <- get_month_families(family_data_Years2_Bacteria, "July", cutoff = family_cutoff)

families_Years2_May_Fungi  <- get_month_families(family_data_Years2_Fungi, "May",  cutoff = family_cutoff)
families_Years2_July_Fungi <- get_month_families(family_data_Years2_Fungi, "July", cutoff = family_cutoff)


# ==========================================================
# 9. Classify Families
# ==========================================================

classify_month_families <- function(df, time_value, abundant_families) {
  
  df_out <- df %>%
    filter(Time == time_value) %>%
    mutate(
      aggTaxo = case_when(
        Family_clean == "Unclassified" ~ "Unclassified",
        Family_clean %in% abundant_families ~ Family_clean,
        TRUE ~ "Others"
      )
    )
  
  return(df_out)
}


# ==========================================================
# 10. Apply the CORRECT Family list to each graph
# ==========================================================

family_class_Years3_May_Bacteria  <- classify_month_families(family_data_Years3_Bacteria, "May",  families_Years3_May_Bacteria)
family_class_Years3_July_Bacteria <- classify_month_families(family_data_Years3_Bacteria, "July", families_Years3_July_Bacteria)

family_class_Years3_May_Fungi  <- classify_month_families(family_data_Years3_Fungi, "May",  families_Years3_May_Fungi)
family_class_Years3_July_Fungi <- classify_month_families(family_data_Years3_Fungi, "July", families_Years3_July_Fungi)

family_class_Years2_May_Bacteria  <- classify_month_families(family_data_Years2_Bacteria, "May",  families_Years2_May_Bacteria)
family_class_Years2_July_Bacteria <- classify_month_families(family_data_Years2_Bacteria, "July", families_Years2_July_Bacteria)

family_class_Years2_May_Fungi  <- classify_month_families(family_data_Years2_Fungi, "May",  families_Years2_May_Fungi)
family_class_Years2_July_Fungi <- classify_month_families(family_data_Years2_Fungi, "July", families_Years2_July_Fungi)

# ==========================================================
# 11. Calculate month-specific legend ordering stats
# ==========================================================

get_family_legend_stats <- function(df) {
  
  sample_family <- df %>%
    group_by(Sample, aggTaxo) %>%
    summarise(relAb = sum(relAb, na.rm = TRUE), .groups = "drop")
  
  stats <- sample_family %>%
    group_by(aggTaxo) %>%
    summarise(
      overall_mean = mean(relAb, na.rm = TRUE),
      overall_sd = sd(relAb, na.rm = TRUE),
      n = sum(!is.na(relAb)),
      .groups = "drop"
    )
  
  return(stats)
}

# ==========================================================
# 12. Calculate May and July legend stats separately
# ==========================================================

legend_Years3_May_Bacteria  <- get_family_legend_stats(family_class_Years3_May_Bacteria)
legend_Years3_July_Bacteria <- get_family_legend_stats(family_class_Years3_July_Bacteria)

legend_Years3_May_Fungi  <- get_family_legend_stats(family_class_Years3_May_Fungi)
legend_Years3_July_Fungi <- get_family_legend_stats(family_class_Years3_July_Fungi)

legend_Years2_May_Bacteria  <- get_family_legend_stats(family_class_Years2_May_Bacteria)
legend_Years2_July_Bacteria <- get_family_legend_stats(family_class_Years2_July_Bacteria)

legend_Years2_May_Fungi  <- get_family_legend_stats(family_class_Years2_May_Fungi)
legend_Years2_July_Fungi <- get_family_legend_stats(family_class_Years2_July_Fungi)

# ==========================================================
# 13. Summarize abundance by Site
# ==========================================================

summarize_family_plot <- function(df) {
  
  plot_df <- df %>%
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
  
  return(plot_df)
}

# ==========================================================
# 14. Create May and July plotting tables
# ==========================================================

family_Years3_May_Bacteria  <- summarize_family_plot(family_class_Years3_May_Bacteria)
family_Years3_July_Bacteria <- summarize_family_plot(family_class_Years3_July_Bacteria)

family_Years3_May_Fungi  <- summarize_family_plot(family_class_Years3_May_Fungi)
family_Years3_July_Fungi <- summarize_family_plot(family_class_Years3_July_Fungi)

family_Years2_May_Bacteria  <- summarize_family_plot(family_class_Years2_May_Bacteria)
family_Years2_July_Bacteria <- summarize_family_plot(family_class_Years2_July_Bacteria)

family_Years2_May_Fungi  <- summarize_family_plot(family_class_Years2_May_Fungi)
family_Years2_July_Fungi <- summarize_family_plot(family_class_Years2_July_Fungi)

# ==========================================================
# 15. CREATE ONE FIXED FAMILY COLOR MAP FOR ALL GRAPHS
# ==========================================================

# ----------------------------------------------------------
# Plot text sizes
# ----------------------------------------------------------

family_legend_title_size <- 16
family_legend_text_size  <- 14

year_label_size <- 7

month_title_size <- 16

# ----------------------------------------------------------
# Bar outline settings
# ----------------------------------------------------------

bar_outline_color <- "black"
bar_outline_width <- 0.2

# ----------------------------------------------------------
# Panel border settings
# ----------------------------------------------------------

panel_border_color <- "black"
panel_border_width  <- 1.5

# ----------------------------------------------------------
# Legend entry spacing
# ----------------------------------------------------------

legend_entry_spacing <- 0.15

# ----------------------------------------------------------
# Family names and their assigned colors
# ----------------------------------------------------------

labels_fam <- c(
  "Acetobacteraceae","Bacillaceae","Beijerinckiaceae","Comamonadaceae",
  "Deinococcaceae","Enterobacteriaceae","Erwiniaceae","Geodermatophilaceae",
  "Hymenobacteraceae","Lactobacillaceae","Microbacteriaceae","Micrococcaceae",
  "Nocardioidaceae","Oxalobacteraceae","Peptostreptococcaceae","Pseudomonadaceae",
  "Pseudonocardiaceae","Roseiflexaceae","Sphingomonadaceae","Spirosomataceae",
  "Weeksellaceae","Botryosphaeriaceae","Buckleyzymaceae","Bulleraceae",
  "Bulleribasidiaceae","Cladosporiaceae","Cryptococcaceae","Didymellaceae",
  "Didymosphaeriaceae","Erysiphaceae","Erythrobasidiaceae","Filobasidiaceae",
  "Intrasporangiaceae","Kineosporiaceae","Mycosphaerellaceae","Nectriaceae",
  "Phaeosphaeriaceae","Polyporaceae","Pseudeurotiaceae","Saccharimonadaceae",
  "Saccotheciaceae","Sclerotiniaceae","Sporidiobolaceae","Sporocadaceae",
  "Taphrinaceae","Venturiaceae","Others","Unclassified"
)


# ----------------------------------------------------------
# Fixed colors
# ----------------------------------------------------------

fill_fam <- c(
  "#F1F42C","#F4D801","#F8B11B","#E17919","#E13E11","#711939","#C85A78",
  "#DA1C91","#AA9ADD","#AA5ADD","#A16591","#A1C991","#A3B78C","#A3E18C",
  "#37A346","#43BCCF","#118998FF","#1F677A","#1F549A","#4A559A","#244154",
  "#F6F199","#F7E42C","#E8D101","#E8B11B","#F9B36D","#D17919","#D13E11",
  "#611939","#B85A78","#E35959","#CA1C91","#D5AEE6","#A99ADD","#A95ADD",
  "#A19591","#A1B991","#A3A78C","#A3D18C","#379346","#43ACCF","#117998FF",
  "#1E677A","#1E549A","#3A559A","#144154","grey40","grey90"
)


# ----------------------------------------------------------
#  Create named color vector
# ----------------------------------------------------------

family_colors <- fill_fam
names(family_colors) <- labels_fam

# ----------------------------------------------------------
# italicize Family names in legends
# ----------------------------------------------------------

italicize_label <- function(x) {
  ifelse(
    x %in% c("Others", "Unclassified"),
    x,
    paste0("*", x, "*")
  )
}

# ==========================================================
# Create shared Family orders for Bacteria and Fungi
# ==========================================================

bacteria_family_order <- sort(
  unique(
    c(
      legend_Years3_May_Bacteria$aggTaxo,
      legend_Years3_July_Bacteria$aggTaxo,
      legend_Years2_May_Bacteria$aggTaxo,
      legend_Years2_July_Bacteria$aggTaxo
    )
  )
)

bacteria_family_order <- c(
  intersect(c("Others", "Unclassified"), bacteria_family_order),
  setdiff(bacteria_family_order, c("Others", "Unclassified"))
)

fungi_family_order <- sort(
  unique(
    c(
      legend_Years3_May_Fungi$aggTaxo,
      legend_Years3_July_Fungi$aggTaxo,
      legend_Years2_May_Fungi$aggTaxo,
      legend_Years2_July_Fungi$aggTaxo
    )
  )
)

fungi_family_order <- c(
  intersect(c("Others", "Unclassified"), fungi_family_order),
  setdiff(fungi_family_order, c("Others", "Unclassified"))
)



# ==========================================================
# Plot function
# ==========================================================

plot_family_bar <- function(df, legend_stats, shared_family_order = NULL) {
  
  df <- df %>%
    left_join(
      legend_stats %>% select(aggTaxo, overall_mean),
      by = "aggTaxo"
    )
  
  if (!is.null(shared_family_order)) {
    
    family_order <- shared_family_order
    
  } else {
    
    family_order <- legend_stats %>%
      arrange(overall_mean) %>%
      pull(aggTaxo)
    
    family_order <- c(
      intersect(c("Others", "Unclassified"), family_order),
      setdiff(family_order, c("Others", "Unclassified"))
    )
  }
  
  df <- df %>%
    mutate(aggTaxo = factor(aggTaxo, levels = family_order))
  
  legend_labels <- italicize_label(family_order)
  names(legend_labels) <- family_order
  legend_title <- "Family"
  
  ggplot(
    df,
    aes(x = Code, y = mean_rel_abund, fill = aggTaxo)
  ) +
    
    geom_col(
      width = 0.7,
      color = bar_outline_color,
      linewidth = bar_outline_width
    ) +
    
    scale_fill_manual(
      values = family_colors,
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
    
    labs(x = "Site", y = "Mean Relative Abundance", fill = legend_title) +
    
    guides(fill = guide_legend(ncol = 1, byrow = FALSE)) +
    
    theme_classic() +
    
    theme(
      axis.title.x = element_text(face = "bold", size = 16),
      axis.text.x = element_text(size = 14),
      axis.title.y = element_text(face = "bold", size = 16),
      axis.text.y = element_text(size = 14),
      legend.title = element_text(size = 16),
      legend.text = ggtext::element_markdown(size = 14),   
      legend.key.height = grid::unit(0.45, "cm"),
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, face = "bold", size = month_title_size),
      panel.border = element_rect(
        colour = panel_border_color,
        fill = NA,
        linewidth = panel_border_width
      )
    )
}

# ==========================================================
# 17. Create BACTERIA plots
# ==========================================================

family_bar_Years3_May_Bacteria <- plot_family_bar(
  family_Years3_May_Bacteria, legend_Years3_May_Bacteria,
  shared_family_order = bacteria_family_order
) + labs(title = "3 Year - May")

family_bar_Years3_July_Bacteria <- plot_family_bar(
  family_Years3_July_Bacteria, legend_Years3_July_Bacteria,
  shared_family_order = bacteria_family_order
) + labs(title = "3 Year - July")

family_bar_Years2_May_Bacteria <- plot_family_bar(
  family_Years2_May_Bacteria, legend_Years2_May_Bacteria,
  shared_family_order = bacteria_family_order
) + labs(title = "2 Year - May")

family_bar_Years2_July_Bacteria <- plot_family_bar(
  family_Years2_July_Bacteria, legend_Years2_July_Bacteria,
  shared_family_order = bacteria_family_order
) + labs(title = "2 Year - July")



# ==========================================================
# 18. Create FUNGI plots
# ==========================================================

family_bar_Years3_May_Fungi <- plot_family_bar(
  family_Years3_May_Fungi, legend_Years3_May_Fungi,
  shared_family_order = fungi_family_order
) + labs(title = "3 Year - May")

family_bar_Years3_July_Fungi <- plot_family_bar(
  family_Years3_July_Fungi, legend_Years3_July_Fungi,
  shared_family_order = fungi_family_order
) + labs(title = "3 Year - July")

family_bar_Years2_May_Fungi <- plot_family_bar(
  family_Years2_May_Fungi, legend_Years2_May_Fungi,
  shared_family_order = fungi_family_order
) + labs(title = "2 Year - May")

family_bar_Years2_July_Fungi <- plot_family_bar(
  family_Years2_July_Fungi, legend_Years2_July_Fungi,
  shared_family_order = fungi_family_order
) + labs(title = "2 Year - July")



# ==========================================================
# 19. CREATE COMPLETE SHARED LEGENDS
# ==========================================================

make_shared_family_legend <- function(family_order, legend_title) {
  
  legend_df <- tibble(
    aggTaxo = factor(family_order, levels = family_order),
    x = seq_along(family_order),
    y = 1
  )
  
  legend_plot <- ggplot(
    legend_df,
    aes(x = x, y = y, fill = aggTaxo)
  ) +
    
    geom_col(
      show.legend = TRUE,
      color = bar_outline_color,
      linewidth = bar_outline_width
    ) +
    
    scale_fill_manual(
      values = family_colors,
      breaks = family_order,
      labels = italicize_label(family_order),
      drop = FALSE
    ) +
    
    labs(fill = legend_title) +
    
    guides(fill = guide_legend(ncol = 1, byrow = FALSE)) +
    
    theme_void() +
    
    theme(
      legend.position = "right",
      legend.title = element_text(face = "bold", size = family_legend_title_size),
      legend.text = ggtext::element_markdown(size = family_legend_text_size),  # CHANGED
      legend.key.height = grid::unit(0.6, "cm"),
      legend.key.width  = grid::unit(0.6, "cm"),
      legend.key.spacing.y = grid::unit(legend_entry_spacing, "cm")
    )
  
  shared_legend <- ggpubr::get_legend(legend_plot)
  
  return(shared_legend)
}

# ==========================================================
# Create BACTERIA shared legend
# ==========================================================

bacteria_shared_legend <- make_shared_family_legend(
  bacteria_family_order,
  legend_title = "Bacterial Families"
)

# ==========================================================
# Create FUNGI shared legend
# ==========================================================

fungi_shared_legend <- make_shared_family_legend(
  fungi_family_order,
  legend_title = "Fungal Families"
)

# ==========================================================
# Measure each legend's NATURAL width
# ==========================================================

bacteria_legend_width <- grid::unit(1, "grobwidth", bacteria_shared_legend) +
  grid::unit(0.3, "cm")

fungi_legend_width <- grid::unit(1, "grobwidth", fungi_shared_legend) +
  grid::unit(0.3, "cm")

# ==========================================================
# 20. CHANGE INDIVIDUAL GRAPH TITLES
# ==========================================================

# --------------------------
# BACTERIA
# --------------------------

bac_3may <- family_bar_Years3_May_Bacteria +
  labs(title = "May") +
  theme(legend.position = "none", axis.title.x = element_blank())

bac_3july <- family_bar_Years3_July_Bacteria +
  labs(title = "July") +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

bac_2may <- family_bar_Years2_May_Bacteria +
  labs(title = "May") +
  theme(legend.position = "none")

bac_2july <- family_bar_Years2_July_Bacteria +
  labs(title = "July") +
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

# --------------------------
# FUNGI
# --------------------------

fun_3may <- family_bar_Years3_May_Fungi +
  labs(title = "May") +
  theme(legend.position = "none", axis.title.x = element_blank())

fun_3july <- family_bar_Years3_July_Fungi +
  labs(title = "July") +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

fun_2may <- family_bar_Years2_May_Fungi +
  labs(title = "May") +
  theme(legend.position = "none")

fun_2july <- family_bar_Years2_July_Fungi +
  labs(title = "July") +
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

# ==========================================================
# 21. FUNCTION FOR YEAR ROW LABEL
# ==========================================================

make_year_label <- function(label) {
  
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = label, fontface = "bold", size = year_label_size) +
    xlim(0, 1) +
    ylim(0, 1) +
    theme_void()
}

# ==========================================================
# 22. FUNCTION FOR BLACK DIVIDER LINE
# ==========================================================

make_divider <- function() {
  
  ggplot() +
    geom_hline(yintercept = 0.5, linewidth = 1, color = "black") +
    xlim(0, 1) +
    ylim(0, 1) +
    theme_void()
}

# ==========================================================
# 23. CREATE ROW LABELS
# ==========================================================

label_3year <- make_year_label("3 Year Data")
label_2year <- make_year_label("2 Year Data")
divider_line <- make_divider()

# ==========================================================
# 24. CREATE BACTERIA GRAPH BLOCK
# ==========================================================

bacteria_graphs <- (
  label_3year /
    (bac_3may | bac_3july) /
    divider_line /
    label_2year /
    (bac_2may | bac_2july)
) +
  plot_layout(heights = c(0.12, 1, 0.035, 0.12, 1))

# ==========================================================
# 25. CREATE FINAL BACTERIA PANEL
# ==========================================================

family_panel_Bacteria <- (
  bacteria_graphs | wrap_elements(full = bacteria_shared_legend)
) +
  plot_layout(
    widths = grid::unit.c(
      grid::unit(1, "null"),
      bacteria_legend_width
    )
  )

# ==========================================================
# 26. FUNGI GRAPH BLOCK
# ==========================================================

fungi_graphs <- (
  label_3year /
    (fun_3may | fun_3july) /
    divider_line /
    label_2year /
    (fun_2may | fun_2july)
) +
  plot_layout(heights = c(0.12, 1, 0.035, 0.12, 1))

# ==========================================================
# 27. CREATE FINAL FUNGI PANEL
# ==========================================================

family_panel_Fungi <- (
  fungi_graphs | wrap_elements(full = fungi_shared_legend)
) +
  plot_layout(
    widths = grid::unit.c(
      grid::unit(1, "null"),
      fungi_legend_width
    )
  )

# ==========================================================
# 28. COMBINE BACTERIA + FUNGI
# ==========================================================

family_panel_All <- (
  bacteria_graphs |
    wrap_elements(full = bacteria_shared_legend) |
    fungi_graphs |
    wrap_elements(full = fungi_shared_legend)
) +
  plot_layout(
    widths = grid::unit.c(
      grid::unit(1, "null"),
      bacteria_legend_width,
      grid::unit(1, "null"),
      fungi_legend_width
    )
  )


family_panel_All

# ==========================================================
# 32. SAVE TABLE OF % ABUNDANCE 
# ==========================================================

# ----------------------------------------------------------
# Add Year / Month / Domain identifiers to each summary table
# ----------------------------------------------------------

add_ids <- function(df, year_group, month, domain) {
  df %>%
    mutate(
      Year_Group = year_group,
      Month = month,
      Domain = domain,
      .before = 1
    )
}

abundance_table_all <- bind_rows(
  
  add_ids(family_Years3_May_Bacteria,  "3 Year", "May",  "Bacteria"),
  add_ids(family_Years3_July_Bacteria, "3 Year", "July", "Bacteria"),
  add_ids(family_Years2_May_Bacteria,  "2 Year", "May",  "Bacteria"),
  add_ids(family_Years2_July_Bacteria, "2 Year", "July", "Bacteria"),
  
  add_ids(family_Years3_May_Fungi,  "3 Year", "May",  "Fungi"),
  add_ids(family_Years3_July_Fungi, "3 Year", "July", "Fungi"),
  add_ids(family_Years2_May_Fungi,  "2 Year", "May",  "Fungi"),
  add_ids(family_Years2_July_Fungi, "2 Year", "July", "Fungi")
  
) %>%
  
  rename(
    Family = aggTaxo,
    Mean_Percent_Abundance = mean_rel_abund
  ) %>%
  
  select(
    Domain,
    Year_Group,
    Month,
    Code,
    Family,
    Mean_Percent_Abundance
  ) %>%
  
  arrange(
    Domain,
    Year_Group,
    Month,
    Code,
    desc(Mean_Percent_Abundance)
  )

abundance_table_all <- abundance_table_all %>%
  mutate(
    Mean_Percent_Abundance = round(Mean_Percent_Abundance, 2)
  )

abundance_table_all

# ----------------------------------------------------------
# Save as CSV
# ----------------------------------------------------------

write.csv(
  abundance_table_all,
  "family_abundance_summary.csv",
  row.names = FALSE
)



# ==========================================================
# 33. TOTAL UNIQUE TAXA PER TAXONOMIC LEVEL
# ==========================================================

clean_taxa_vector <- function(x) {
  x <- as.character(x)
  case_when(
    is.na(x) ~ "Unclassified",
    x == "" ~ "Unclassified",
    x == "NA" ~ "Unclassified",
    grepl("Incertae|Unclassified|uncultured|unknown", x, ignore.case = TRUE) ~ "Unclassified",
    TRUE ~ x
  )
}

get_unique_taxa_list <- function(ps) {
  
  tax_df <- as.data.frame(as(tax_table(ps), "matrix"), stringsAsFactors = FALSE)
  
  rank_list <- purrr::map(
    names(tax_df),
    function(rank) {
      cleaned <- clean_taxa_vector(tax_df[[rank]])
      cleaned <- cleaned[cleaned != "Unclassified"]
      unique(cleaned)
    }
  )
  
  names(rank_list) <- names(tax_df)
  
  # ADDED: ASV-level count
  rank_list[["ASV"]] <- taxa_names(ps)
  
  return(rank_list)
}

# ----------------------------------------------------------
# Get the unique-taxa lists for each of the 4 datasets
# ----------------------------------------------------------

taxa_Years3_Bacteria <- get_unique_taxa_list(Years3$Bacteria)
taxa_Years3_Fungi    <- get_unique_taxa_list(Years3$Fungi)
taxa_Years2_Bacteria <- get_unique_taxa_list(Years2$Bacteria)
taxa_Years2_Fungi    <- get_unique_taxa_list(Years2$Fungi)

count_by_rank <- function(taxa_list, year_group, domain) {
  tibble(
    Domain = domain,
    Year_Group = year_group,
    Rank = names(taxa_list),
    Unique_Count = purrr::map_int(taxa_list, length)
  )
}

unique_taxa_by_dataset <- bind_rows(
  count_by_rank(taxa_Years3_Bacteria, "3 Year", "Bacteria"),
  count_by_rank(taxa_Years3_Fungi,    "3 Year", "Fungi"),
  count_by_rank(taxa_Years2_Bacteria, "2 Year", "Bacteria"),
  count_by_rank(taxa_Years2_Fungi,    "2 Year", "Fungi")
)

union_by_rank <- function(list1, list2) {
  ranks <- union(names(list1), names(list2))
  out <- purrr::map(ranks, function(r) {
    v1 <- if (r %in% names(list1)) list1[[r]] else character(0)
    v2 <- if (r %in% names(list2)) list2[[r]] else character(0)
    union(v1, v2)
  })
  names(out) <- ranks
  return(out)
}

# ----------------------------------------------------------
# Domain totals — Bacteria (Year2+Year3), Fungi (Year2+Year3)
# ----------------------------------------------------------

taxa_Bacteria_AllYears <- union_by_rank(taxa_Years3_Bacteria, taxa_Years2_Bacteria)
taxa_Fungi_AllYears    <- union_by_rank(taxa_Years3_Fungi,    taxa_Years2_Fungi)

domain_totals <- bind_rows(
  tibble(
    Domain = "Bacteria", Year_Group = "All Years",
    Rank = names(taxa_Bacteria_AllYears),
    Unique_Count = purrr::map_int(taxa_Bacteria_AllYears, length)
  ),
  tibble(
    Domain = "Fungi", Year_Group = "All Years",
    Rank = names(taxa_Fungi_AllYears),
    Unique_Count = purrr::map_int(taxa_Fungi_AllYears, length)
  )
)

# ----------------------------------------------------------
# Combine into one summary table
# ----------------------------------------------------------

unique_taxa_long <- bind_rows(
  unique_taxa_by_dataset %>% mutate(Group_Type = "By Dataset", .before = 1),
  domain_totals          %>% mutate(Group_Type = "By Domain (All Years)", .before = 1)
)

standard_rank_order <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "ASV")
rank_levels <- c(
  intersect(standard_rank_order, unique(unique_taxa_long$Rank)),
  setdiff(unique(unique_taxa_long$Rank), standard_rank_order)
)

unique_taxa_summary <- unique_taxa_long %>%
  mutate(Rank = factor(Rank, levels = rank_levels)) %>%
  pivot_wider(names_from = Rank, values_from = Unique_Count) %>%
  select(Group_Type, Domain, Year_Group, all_of(rank_levels)) %>%
  arrange(Group_Type, Domain, Year_Group)

unique_taxa_summary

# ----------------------------------------------------------
# Save as CSV
# ----------------------------------------------------------

write.csv(
  unique_taxa_summary,
  "unique_taxa_summary.csv",
  row.names = FALSE
)

