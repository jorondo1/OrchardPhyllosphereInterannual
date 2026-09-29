## ANCOMBC
# Process output

pacman::p_load(
  mgx.tools,
  tidyverse,
  kableExtra,
  patchwork,
  phyloseq,
  update = FALSE
)

# ---- 0. Load ------------------------------------------------

ANCOMResults    <- readRDS("data/ancom/ANCOMResults.rds")
res   <- ANCOMResults$ancom_results
gr_var <- ANCOMResults$group_variables
alpha <- ANCOMResults$alpha # p=0.01

# 1. RESULTS TABLE — filtered to taxa passing both
#    q-value < alpha AND the pseudo-count sensitivity test

extract_full_ancom_table <- function(ancom_result, ps, year_dataset, dataset_name) {
  
  res <- as.data.frame(ancom_result$res, stringsAsFactors = FALSE)
  if (!"taxon" %in% colnames(res)) res$taxon <- rownames(res)
  tax <- as.data.frame(tax_table(ps), stringsAsFactors = FALSE)
  tax$taxon <- rownames(tax)
  
  # Discover terms actually present (excludes Intercept)
  lfc_cols <- grep("^lfc_", colnames(res), value = TRUE)
  terms <- sub("^lfc_", "", lfc_cols)
  terms <- terms[terms != "(Intercept)"]
  
  purrr::map_dfr(terms, function(term) {
    get_col <- function(prefix) {
      col <- paste0(prefix, "_", term)
      if (col %in% colnames(res)) res[[col]] else NA
    }
    data.frame(
      taxon        = res$taxon,
      Year_Dataset = year_dataset,
      Dataset      = dataset_name,
      Term         = term,
      LFC          = get_col("lfc"),
      SE           = get_col("se"),
      W_stat       = get_col("W"),
      p_val        = get_col("p"),
      q_val        = get_col("q"),
      diff_abn     = get_col("diff"),
      passed_ss    = get_col("passed_ss"),
      stringsAsFactors = FALSE
    )
  }) %>%
    left_join(tax, by = "taxon") %>%
    filter(
      !is.na(Species_cluster),
      trimws(Species_cluster) != "",
      tolower(trimws(Species_cluster)) != "overall",
      !is.na(q_val), q_val < alpha,
      !is.na(passed_ss), passed_ss == TRUE
    )
}

all_ancom_table <- purrr::map_dfr(names(res), function(yd) {
  purrr::map_dfr(names(res[[yd]]), function(dn) {
    extract_full_ancom_table(
      ancom_result = res[[yd]][[dn]],
      ps           = ANCOMResults$ps_final[[yd]][[dn]],
      year_dataset = yd,
      dataset_name = dn
    )
  })
})

write_csv(all_ancom_table, "data/ancom/all_ancom_results.csv")


# 2. EXTRACTOR — pulls significant hits + full stats per taxon

get_heatmap_data_sd <- function(ancom_result, ps, gr_var,
                                year_dataset, dataset_name) {
  
  if (is.null(ancom_result) || is.null(ps)) return(NULL)
  res <- as.data.frame(ancom_result$res, stringsAsFactors = FALSE)
  if (!"taxon" %in% colnames(res)) res$taxon <- rownames(res)
  tax <- as.data.frame(tax_table(ps), stringsAsFactors = FALSE)
  tax$taxon <- rownames(tax)
  
  # Relative abundance (%) per sample
  ps_rel <- transform_sample_counts(ps, function(x) {
    if (sum(x) == 0) return(rep(0, length(x)))
    x / sum(x) * 100
  })
  otu <- as(otu_table(ps_rel), "matrix")
  if (!taxa_are_rows(ps_rel)) otu <- t(otu)
  
  sample_info <- data.frame(sample_data(ps_rel), stringsAsFactors = FALSE)
  
  summarize_period <- function(period_samples) {
    period_samples <- base::intersect(colnames(otu), period_samples)
    if (length(period_samples) == 0) {
      out <- data.frame(mean = rep(NA_real_, nrow(otu)), sd = rep(NA_real_, nrow(otu)))
      rownames(out) <- rownames(otu)
      return(out)
    }
    sub <- otu[, period_samples, drop = FALSE]
    out <- data.frame(
      mean = rowMeans(sub, na.rm = TRUE),
      sd   = if (ncol(sub) > 1) apply(sub, 1, sd, na.rm = TRUE) else rep(NA_real_, nrow(sub))
    )
    rownames(out) <- rownames(otu)
    out
  }
  
  may_samples  <- rownames(sample_info)[as.character(sample_info$Time) == "May"]
  july_samples <- rownames(sample_info)[as.character(sample_info$Time) == "July"]
  
  may_stats  <- summarize_period(may_samples)
  july_stats <- summarize_period(july_samples)
  
  abundance_table <- data.frame(
    taxon = rownames(otu),
    Mean_Relative_Abundance_Total = as.numeric(rowMeans(otu, na.rm = TRUE)),
    Mean_Relative_Abundance_May   = as.numeric(may_stats[rownames(otu), "mean"]),
    SD_Relative_Abundance_May     = as.numeric(may_stats[rownames(otu), "sd"]),
    Mean_Relative_Abundance_July  = as.numeric(july_stats[rownames(otu), "mean"]),
    SD_Relative_Abundance_July    = as.numeric(july_stats[rownames(otu), "sd"]),
    stringsAsFactors = FALSE
  )
  
  by_group <- lapply(gr_var, function(group_var) {
    
    diff_cols <- grep(paste0("^diff_", group_var), colnames(res), value = TRUE)
    if (length(diff_cols) == 0) return(NULL)
    
    inner <- lapply(diff_cols, function(diff_col) {
      
      term    <- sub("^diff_", "", diff_col)
      lfc_col <- paste0("lfc_", term)
      se_col  <- paste0("se_", term)
      q_col   <- paste0("q_", term)
      ss_col  <- paste0("passed_ss_", term)
      
      if (!lfc_col %in% colnames(res)) return(NULL)
      
      keep <- res[[diff_col]] == TRUE
      keep[is.na(keep)] <- FALSE
      if (!any(keep)) return(NULL)
      
      temp <- data.frame(
        taxon      = res$taxon[keep],
        Comparison = term,
        LFC        = as.numeric(res[[lfc_col]][keep]),
        SE         = if (se_col %in% colnames(res)) as.numeric(res[[se_col]][keep]) else NA_real_,
        q_val      = if (q_col  %in% colnames(res)) as.numeric(res[[q_col]][keep])  else NA_real_,
        passed_ss  = if (ss_col %in% colnames(res)) res[[ss_col]][keep]             else NA,
        stringsAsFactors = FALSE
      )
      
      n_before_join <- nrow(temp)
      
      temp <- temp %>%
        left_join(tax, by = "taxon") %>%
        left_join(abundance_table, by = "taxon") %>%
        mutate(Species_cluster = as.character(Species_cluster)) %>%
        filter(
          !is.na(Species_cluster),
          trimws(Species_cluster) != "",
          tolower(trimws(Species_cluster)) != "overall",
          !is.na(LFC),
          is.finite(LFC)
        )
      
      cat(
        "  ", year_dataset, dataset_name, term,
        "- sig hits:", n_before_join,
        "| kept after taxonomy join/filter:", nrow(temp), "\n"
      )
      
      if (nrow(temp) == 0) return(NULL)
      
      temp$Group <- group_var
      temp
    })
    
    inner <- inner[!vapply(inner, is.null, logical(1))]
    if (length(inner) == 0) return(NULL)
    bind_rows(inner)
  })
  
  by_group <- by_group[!vapply(by_group, is.null, logical(1))]
  if (length(by_group) == 0) return(NULL)
  
  bind_rows(by_group) %>%
    mutate(
      Year_Dataset    = year_dataset,
      Dataset         = dataset_name,
      Species_cluster = as.character(Species_cluster),
      Comparison      = as.character(Comparison)
    ) %>%
    group_by(Species_cluster, Comparison) %>%
    dplyr::slice(1) %>%
    ungroup()
}


# 3. RUN EXTRACTOR ACROSS ALL DATASETS

heatmap_list <- list()
counter <- 1

for (year_dataset in names(res)) {
  for (dataset_name in names(res[[year_dataset]])) {
    heat_data <- get_heatmap_data_sd(
      ancom_result    = res[[year_dataset]][[dataset_name]],
      ps              = ANCOMResults$ps_final[[year_dataset]][[dataset_name]],
      gr_var = gr_var,
      year_dataset    = year_dataset,
      dataset_name    = dataset_name
    )
    if (!is.null(heat_data)) {
      heatmap_list[[counter]] <- heat_data
      counter <- counter + 1
    }
  }
}

heatmap_all <- if (length(heatmap_list) > 0) bind_rows(heatmap_list) else NULL

if (is.null(heatmap_all) || nrow(heatmap_all) == 0) {
  stop("No significant taxa found (diff_ == TRUE) — nothing to filter or plot.")
}


# 4. "DEFINITE HIT"  q-value < alpha (0.01) AND passed the
#    pseudo-count sensitivity test. 


n_start <- n_distinct(heatmap_all$Species_cluster)

heatmap_all <- heatmap_all %>%
  filter(!is.na(q_val), q_val < alpha)
n_after_q <- n_distinct(heatmap_all$Species_cluster)

heatmap_all <- heatmap_all %>%
  filter(!is.na(passed_ss), passed_ss == TRUE)
n_final <- n_distinct(heatmap_all$Species_cluster)

cat(
  "\nFilter cascade (unique Species_clusters):\n",
  "  starting (diff_==TRUE hits):      ", n_start, "\n",
  "  + q-value <", alpha, ":                  ", n_after_q, "\n",
  "  + passed sensitivity test:        ", n_final, "\n"
)

if (nrow(heatmap_all) == 0) {
  stop("No taxa survive the q-value + sensitivity filter.")
}


# 4A. HEATMAP ABUNDANCE FILTER
#     Keep Species_clusters  >= 1% relative abundence in May or July


heatmap_abundance_cutoff <- 1

abundant_taxa <- heatmap_all %>%
  group_by(Species_cluster) %>%
  summarise(
    max_May_abundance  = max(Mean_Relative_Abundance_May, na.rm = TRUE),
    max_July_abundance = max(Mean_Relative_Abundance_July, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(
    max_May_abundance >= heatmap_abundance_cutoff |
      max_July_abundance >= heatmap_abundance_cutoff
  ) %>%
  pull(Species_cluster)

cat(
  "\nHeatmap abundance filter:\n",
  "  Cutoff: >=", heatmap_abundance_cutoff, "%\n",
  "  Species_clusters before abundance filter:",
  n_distinct(heatmap_all$Species_cluster), "\n",
  "  Species_clusters passing abundance filter:",
  length(abundant_taxa), "\n"
)

# Apply the filter only to the heatmap data.
heatmap_all <- heatmap_all %>%
  filter(Species_cluster %in% abundant_taxa)

if (nrow(heatmap_all) == 0) {
  stop("No Species_clusters meet the >= 1% May/July relative abundance cutoff.")
}


# 5. HEATMAP COLUMN LABELS 


year_suffix <- c(Years3 = "(3)", Years2 = "(2)")

heatmap_all <- heatmap_all %>%
  mutate(
    Level = str_remove(Comparison, paste0("^", Group)),
    Heatmap_Column = paste0(Level, " ", year_suffix[Year_Dataset])
  )


# 6. FINAL DIFFERENTIAL ABUNDANCE TABLE
#    (includes SE, q-value, sensitivity flag for reporting)


taxonomy_cols <- intersect(c("Species_cluster"), colnames(heatmap_all))

final_table <- heatmap_all %>%
  transmute(
    Year_Dataset,
    Dataset,
    Group,
    Heatmap_Column,
    across(all_of(taxonomy_cols)),
    LFC         = round(LFC, 3),
    SE          = round(SE, 3),
    `q-value`   = signif(q_val, 3),
    `Passed sensitivity test` = passed_ss,
    `May %  (mean +/- SD)`  = sprintf("%.2f +/- %.2f", Mean_Relative_Abundance_May, SD_Relative_Abundance_May),
    `July % (mean +/- SD)` = sprintf("%.2f +/- %.2f", Mean_Relative_Abundance_July, SD_Relative_Abundance_July)
  ) %>%
  arrange(Year_Dataset, Dataset, Group, desc(abs(LFC)))

tables_by_dataset <- split(final_table, final_table$Dataset)

for (d in names(tables_by_dataset)) {
  cat("\n==================== TABLE:", d, "====================\n")
  
  print(
    tables_by_dataset[[d]] %>%
      kbl(
        caption = paste(
          "Differentially abundant Species_clusters -", d,"(q <", alpha, ", sensitivity-passed)"
        )
      ) %>%
      kable_styling(
        bootstrap_options = c("striped", "hover", "condensed"),
        full_width = FALSE
      )
  )
}

write_csv(final_table, "data/ancom/differential_abundance_table.csv")


# 7. TAXON LABEL FOR HEATMAP Determined From Blast


fungi_label_overrides <- c(
  "NA_sp_clust_4"                              = "Cladosporium_4(B)",
  "Ascomycota_sp_clust_5"                      = "Didymellaceae_5(B)",
  "Pleosporales_gen_Incertae_sedis_sp_clust_7" = "Alternaria_7(B)",
  "Ascomycota_sp_clust_14"                     = "Melanommataceae_14(B)",
  "NA_sp_clust_15"                             = "Filobasidium_15(B)"
  )

heatmap_all <- heatmap_all %>%
  mutate(
    Taxon_Label = case_when(
      Dataset == "Fungi" & Species_cluster %in% names(fungi_label_overrides) ~
        fungi_label_overrides[Species_cluster],
      TRUE ~ str_replace(Species_cluster, "_sp_clust_", "_")
    )
  )


# 8. BUILD + PLOT HEATMAP


plot_heatmap <- function(df, dataset_name) {
  
  taxon_order <- df %>%
    group_by(Taxon_Label) %>%
    summarise(max_abs_lfc = max(abs(LFC), na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(max_abs_lfc)) %>%
    pull(Taxon_Label)
  
  # Column order: all Years3 columns before Years2 columns; within each
  # year, Time (month) columns before Management columns; within each
  # Group, alphabetical (e.g. May before July).
  column_order <- df %>%
    distinct(Year_Dataset, Group, Heatmap_Column) %>%
    mutate(
      Year_Dataset = factor(Year_Dataset, levels = c("Years3", "Years2")),
      Group        = factor(Group, levels = c("Time", "Management"))
    ) %>%
    arrange(Group, Year_Dataset, Heatmap_Column) %>%
    pull(Heatmap_Column)
  
  complete_grid <- expand_grid(
    Taxon_Label = unique(df$Taxon_Label),
    Heatmap_Column = column_order
  )
  
  df_complete <- complete_grid %>%
    left_join(
      df %>%
        group_by(Taxon_Label, Heatmap_Column) %>%
        summarise(LFC = ifelse(all(is.na(LFC)), NA_real_, mean(LFC, na.rm = TRUE)), .groups = "drop"),
      by = c("Taxon_Label", "Heatmap_Column")
    )
  
  df_complete$Taxon_Label   <- factor(df_complete$Taxon_Label, levels = rev(taxon_order))
  df_complete$Heatmap_Column <- factor(df_complete$Heatmap_Column, levels = column_order)
  
  max_abs_lfc <- max(abs(df_complete$LFC), na.rm = TRUE)
  if (!is.finite(max_abs_lfc) || max_abs_lfc == 0) max_abs_lfc <- 1
  
  ggplot(df_complete, aes(x = Heatmap_Column, y = Taxon_Label, fill = LFC)) +
    geom_tile(color = "black", linewidth = 0.5) +
    geom_text(aes(label = ifelse(is.na(LFC), "", sprintf("%.2f", LFC))),  color = "black") +
    scale_fill_gradient2(
      low = "tomato", mid = "white", high = "royalblue2", midpoint = 0,
      limits = c(-max_abs_lfc, max_abs_lfc), name = "LFC", na.value = "white"
    ) +
    labs(title = dataset_name, x = NULL, y = NULL) +
    theme_minimal() +
    theme(
      plot.title   = element_text(hjust = 0.5),
      axis.text.x  = element_text(angle = 45, hjust = 1),
      panel.grid   = element_blank(),
      panel.border = element_blank(),
      axis.ticks   = element_blank()
    )
}

datasets <- unique(heatmap_all$Dataset)
plots_by_dataset <- list()

for (d in datasets) {
  df_sub <- heatmap_all %>% filter(Dataset == d)
  p <- plot_heatmap(df_sub, d)
  plots_by_dataset[[d]] <- p
  print(p)
}

plots_by_dataset$Bacteria + plots_by_dataset$Fungi



ggsave('out/manuscript/DA_amys.pdf', bg = 'white', 
       width = 1700, height = 1700, units = 'px', dpi = 200)
