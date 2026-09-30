## ANCOMBC
# Process output (compact rewrite of 3.2_Ancom_process.R; same tables + figure)

pacman::p_load(
  mgx.tools,
  tidyverse,
  patchwork,
  phyloseq,
  update = FALSE
)

source('src/0.0_Config.R')

ANCOMResults <- readRDS("data/ancom/ANCOMResults.rds")

# 1. Format ANCOM -----------------------

# Mean/SD relative abundance (%) per taxon, overall and by sampling month
abund_summary <- function(ps) {
  otu <- as(otu_table(ps), "matrix")
  if (!taxa_are_rows(ps)) otu <- t(otu)
  rel <- sweep(otu, 2, colSums(otu), "/") * 100
  rel[, colSums(otu) == 0] <- 0
  time <- as.character(data.frame(sample_data(ps))[colnames(rel), "Time"])
  msd <- function(m) if (ncol(m) > 1) apply(m, 1, sd, na.rm = TRUE) else NA_real_
  may  <- rel[, time == "May",  drop = FALSE]
  july <- rel[, time == "July", drop = FALSE]
  tibble(
    taxon = rownames(rel),
    mean_relab = rowMeans(rel, na.rm = TRUE),
    mean_relab_may   = if (ncol(may))  rowMeans(may,  na.rm = TRUE) else NA_real_,
    sd_relab_may     = msd(may),
    mean_relab_july  = if (ncol(july)) rowMeans(july, na.rm = TRUE) else NA_real_,
    sd_relab_july    = msd(july)
  )
}

# Parse ANCOM output into long tibble
tidy_fit <- function(fit, ps, year_dataset, dataset_name) {
  
  res <- as.data.frame(fit$res)
  terms <- setdiff(
    sub("^lfc_", "", grep("^lfc_", names(res), value = TRUE)), 
    "(Intercept)")
  
  res %>%
    select(taxon, matches("^(lfc|se|W|p|q|diff|passed_ss)_")) %>%
    pivot_longer(-taxon, names_to = c("stat", "Term"),
                 names_pattern = "^(lfc|se|W|p|q|diff|passed_ss)_(.*)$") %>%
    filter(Term %in% terms) %>% # drops e.g. diff_robust_* parsed as term "robust_*"
    pivot_wider(names_from = stat, values_from = value) %>%
    filter(as.logical(passed_ss)) %>% 
    arrange(match(Term, terms)) %>%
    transmute(
      taxon,  Term,
      Year_Dataset = year_dataset, 
      Dataset = dataset_name, 
      LFC = lfc, SE = se, q_val = q) %>%
    # add abundance for filtering
    left_join(abund_summary(ps), 'taxon') %>%
    filter(
      !is.na(taxon),
      trimws(taxon) != "",
      tolower(trimws(taxon)) != "overall",
      q_val < ANCOMResults$alpha) 
}

# Execute
ancom_long <- imap_dfr(
  ANCOMResults$ancom_results, \(fits, yd)
  imap_dfr(
    fits, \(fit, dn) tidy_fit(fit, ANCOMResults$ps_final[[yd]][[dn]], yd, dn))) %>% 
  filter(Term %in% c('TimeMay', 'ManagementOrganic'))

# 2. Write full table ----------------------------------

all_ancom_table <- ancom_long %>% select(-starts_with(c("mean_relab", "sd_relab")))
write_csv(all_ancom_table, "out/manuscript/all_ancom_results.csv")

# 3. Plot data -------------------------

heatmap_all <- ancom_long %>%
  mutate(
    X_label = case_when(
      Term == 'TimeMay' & Year_Dataset == 'Years3' ~ 'July\n(3Y)',
      Term == 'TimeMay' & Year_Dataset == 'Years2' ~ 'July\n(2Y)',
      Term == 'ManagementOrganic' & Year_Dataset == 'Years3' ~ 'NOT_PLOTTED',
      Term == 'ManagementOrganic' & Year_Dataset == 'Years2' ~ 'Organic\n(2Y)'),
    LFC = case_when(Term=='TimeMay' ~ -LFC, TRUE~LFC),
    taxLabel = if_else(
      Dataset == "Fungi" & taxon %in% names(fungi_label_overrides),
      fungi_label_overrides[taxon],
      str_replace(taxon, "_sp_clust_", "_"))) %>% 
  # Remove very low abundance ones;  >= 1% mean abundance in May or July 
  #    (abundance max taken per taxon across all datasets/year subsets)
  filter(max(mean_relab_may, mean_relab_july) >= 1,
         .by = taxon) 

cat("Heatmap Species_clusters per dataset:\n")
print(heatmap_all %>% distinct(Dataset, taxon) %>% count(Dataset))

# 4. HEATMAP -----------------------------------------------------------------------

plot_heatmap <- function(df, title) {
  
  lim <- max(abs(df$LFC), na.rm = TRUE)
  
  
  df %>%
    # Force fully square data:
    complete(taxLabel, X_label) %>% 
    mutate(
      taxLabel = factor(taxLabel, levels = sort(unique(taxLabel), decreasing = TRUE )),
      X_label = factor(X_label, levels = c("July\n(3Y)", "July\n(2Y)", 'Organic\n(2Y)' ) )) %>% 
    ggplot(
      aes(x = X_label, y = taxLabel, fill = LFC)) +
    geom_tile(color = "black", linewidth = 0.2) +
    geom_text(
      aes(label = ifelse(is.na(LFC), "", sprintf("%.2f", LFC))), 
      color = "black", size = 3) +
    scale_fill_gradient2(
      low = "tomato", mid = "white", high = "royalblue2", midpoint = 0,
      limits = c(-lim, lim), name = "LFC", na.value = "white") +
    #    scale_x_discrete(labels = c(Years3 = "3Y", Years2 = "2Y")) +
    labs(x = NULL, y = NULL) +
    theme(
      axis.text.x = element_text(face = 'bold'),
      panel.grid   = element_blank(),
      panel.border = element_blank(),
      axis.ticks   = element_blank()) +
    scale_x_discrete(position = "top")
  
}

plots_by_dataset <- heatmap_all %>%
  filter(!X_label == 'NOT_PLOTTED') %>% 
  split(.$Dataset) %>%
  imap(plot_heatmap)

wrap_plots(plots_by_dataset) +
  plot_layout(guides = 'collect') +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag.position = c(0, 0.98))


ggsave('out/manuscript/4_DA.pdf', bg = 'white', 
       width = 2500, height = 2000, units = 'px', dpi = 300)



# 5. REPORTING TABLE -------------------------------------------------------------

final_table <- heatmap_all %>%
  filter(!X_label == 'NOT_PLOTTED') %>%
  # Wald 95% CI, LFC +/- 1.96 SE, from unrounded values (LFC already sign-flipped
  # for TimeMay above; SE is unaffected by the flip)
  mutate(ci_low  = LFC - qnorm(0.975) * SE,
         ci_high = LFC + qnorm(0.975) * SE) %>%
  transmute(
    X_label, Dataset, taxLabel,
    LFC = round(LFC, 2),
    #SE  = round(SE, 3),
    `95% CI` = sprintf("[%.2f, %.2f]", ci_low, ci_high),
    #`q-value` = signif(q_val, 3),
    `May %  (mean +/- SD)` = sprintf("%.2f +/- %.2f", mean_relab_may, sd_relab_may),
    `July % (mean +/- SD)` = sprintf("%.2f +/- %.2f", mean_relab_july, sd_relab_july)
  ) %>%
  arrange(Dataset, X_label, taxLabel, desc(abs(LFC)))

writexl::write_xlsx(final_table, "out/manuscript/supp/differential_abundance_table.xlsx")
