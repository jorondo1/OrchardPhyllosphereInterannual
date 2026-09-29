# FAMILY RELATIVE ABUNDANCE

# Author: Amy Heim
# Revisions-consolidation by Jonathan Rondeau-Leclaire
# Optimized by Claude Code 

# 1. Load packages -------------------

pacman::p_load(
  tidyverse,  phyloseq,  mgx.tools,  patchwork,  ggpubr,  ggtext,  update = FALSE
)

# 2. Load datasets --------------------

# Nested Kingdom > Duration list of phyloseq objects 
ps.ls <- readRDS('data/ps_objects_full.rds')
ps_nested <- list(
  Bacteria = list(`3-year` = ps.ls$Bacteria_3y, `2-year` = ps.ls$Bacteria_2y),
  Fungi    = list(`3-year` = ps.ls$Fungi_3y,    `2-year` = ps.ls$Fungi_2y)
)

# 3. Prepare Family relative-abundance data --------------

prepare_family_data <- function(ps) {
  
  transform_sample_counts(
    ps,  function(x) {(x / sum(x)) * 100}
  ) %>% 
    psflashmelt() %>%
    mutate(
      Family_clean = case_when(
        is.na(Family) ~ "Unclassified",
        Family == "" ~ "Unclassified",
        Family == "NA" ~ "Unclassified",
        grepl(
          "Incertae|Unclassified|uncultured|unknown",
          Family,
          ignore.case = TRUE ) ~ "Unclassified",
        TRUE ~ as.character(Family)
      )
    ) %>%
    group_by(Sample,Code, Time, Year, Family_clean) %>%
    summarise(
      relAb = sum(Abundance, na.rm = TRUE),
      .groups = "drop"
    )
}

# Apply to every Kingdom > Duration combination:
family_data <- map_depth(ps_nested, 2, prepare_family_data)

# 4. Classify Families  ------------------------------

# function for family classification by month
classify_month_families <- function(df, time_value, cutoff = 1) { 
  
  # find Families with average relative abundance >= 1%
  abundant_families <- df %>%
    filter(Time == time_value, Family_clean != "Unclassified") %>%
    group_by(Sample, Family_clean) %>%
    summarise(relAb = sum(relAb, na.rm = FALSE), .groups = "drop") %>%
    group_by(Family_clean) %>%
    summarise(mean_rel_abund = mean(relAb, na.rm = FALSE), .groups = "drop") %>%
    filter(mean_rel_abund >= cutoff) %>%
    pull(Family_clean)
  
  # create aggregate taxonimy
  df %>%
    filter(Time == time_value) %>%
    mutate(
      aggTaxo = case_when(
        Family_clean == "Unclassified" ~ "Unclassified",
        Family_clean %in% abundant_families ~ Family_clean,
        TRUE ~ "Others"
      ), .keep = 'unused'
    )
}


# Apply the Family classification to each plot This stage adds the Month
# level, so map_depth (fixed depth) doesn't apply directly
months <- c(May = "May", July = "July")

family_class <- map(family_data, function(duration_list) {
  map(duration_list, function(df) {
    map(months, ~ classify_month_families(df, .x))
  })
})


# 6. Plot parameters --------------------------

family_legend_title_size <- 16
family_legend_text_size  <- 14
year_label_size <- 5
month_title_size <- 12
bar_outline_color <- "black"
bar_outline_width <- 0.1
legend_entry_spacing <- 0.15

# Family names and their assigned colors
fill_fam <- c(
  # Bacteria:
  Acetobacteraceae = "#F1F42C", Bacillaceae = "#F4D801", 
  Beijerinckiaceae = "#F8B11B", Comamonadaceae = "#E17919",
  Deinococcaceae = "#E13E11", Enterobacteriaceae = "#711939",
  Erwiniaceae = "#C85A78", Geodermatophilaceae = "#DA1C91",
  Hymenobacteraceae = "#AA9AFF", Kineosporiaceae = "#A99ADD",
  Lactobacillaceae = "#AA5ADD", 
  Microbacteriaceae = "#A16591", Micrococcaceae = "#A1C991",
  Nocardioidaceae = "#A3B78C", Oxalobacteraceae = "#A3E18C", 
  Peptostreptococcaceae = "#37A346", Pseudomonadaceae = "#43BCCF",
  Pseudonocardiaceae = "#118998FF", Roseiflexaceae = "#1F677A",
  Sphingomonadaceae = "#1F549A", Spirosomataceae = "#4A559A", 
  Weeksellaceae = "#244154",
  #Fungi
  Botryosphaeriaceae = "#F6F199", Buckleyzymaceae = "#F7E42C", 
  Bulleraceae = "#E8D101", Bulleribasidiaceae = "#E8B11B",
  Cladosporiaceae = "#F9B36D", Cryptococcaceae = "#D17919", 
  Didymellaceae = "#E13E11", Didymosphaeriaceae = "#611939",
  Erysiphaceae = "#B85A78", Erythrobasidiaceae = "#E35959", 
  Filobasidiaceae = "#CA1C91", Intrasporangiaceae = "#D5AEE6",
  Mycosphaerellaceae = "#A95ADD",
  Nectriaceae = "#A19591", Phaeosphaeriaceae = "#A1B991",
  Polyporaceae = "#A3A78C", Pseudeurotiaceae = "#A3D18C", 
  Saccharimonadaceae = "#379346", Saccotheciaceae = "#43ACCF",
  Sclerotiniaceae = "#117998FF", Sporidiobolaceae = "#1E677A", 
  Sporocadaceae = "#1E549A", Taphrinaceae = "#3A559A", 
  Venturiaceae = "#144154",
  Others = "grey50", Unclassified = "grey90"
)

color_fam <- setNames(rep("black", length(fill_fam)), names(fill_fam))

extra_cats <-  c("Others", "Unclassified")

# italicize Family names in legends
italicize_label <- function(x) {
  ifelse(x %in% extra_cats, x, paste0("*", x, "*"))
}

# Create shared Family order per Kingdom (union of aggTaxo across all
# Duration x Month combinations for that Kingdom).
family_order <- map(family_class, function(duration_list) {
  duration_list %>%
    map_depth(2, "aggTaxo") %>%
    unlist(use.names = FALSE) %>%
    unique() %>% sort() %>%
    setdiff(extra_cats) %>%
    c(extra_cats, .)
})


