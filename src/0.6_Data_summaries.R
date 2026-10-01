pacman::p_load(tidyverse, here, phyloseq, mgx.tools, kableExtra, magrittr, update = FALSE)
source('src/0.0_Config.R') # ps_dataset_labels

ps.ls <- read_rds('data/ps_objects_full.rds')

# Sample count viz --------------

dat <- rbind(
  samdat_as_tibble(ps.ls$Fungi) %>% mutate(Barcode = 'Fungi'),
  samdat_as_tibble(ps.ls$Bacteria) %>% mutate(Barcode = 'Bacteria')
)

count_dat <- dat %>%
  count(Barcode, Dataset, Year, Time, Cultivar, Code, Management) %>%
  rename(N_samples = n)

count_dat %>% 
  ggplot(aes(x = Time, y = N_samples, fill = Code)) +
  geom_col(position = "dodge") +
  ggh4x::facet_nested(
    Barcode+Year ~ Dataset + Management + Cultivar) +  # Facet by 2 variables
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 5),
    legend.position = 'bottom',
    panel.grid = element_blank()) +
  guides(fill = guide_legend(nrow = 1)) +
  scale_fill_manual(values = fill_site) +
  labs(fill = 'Orchard')

ggsave('out/manuscript/supp/sample_count_by_metadata.pdf',
       bg = 'white', width = 2200, height = 1700, 
       units = 'px', dpi = 220)

# Classification rates -------------------

classrates <- compute_classification_rates(ps.ls[c('Bacteria','Fungi') ])

plot_class_rates(classrates)

ggsave('out/summaries/classification_rates.pdf',
       bg = 'white', width = 2000, height = 2000, 
       units = 'px', dpi = 300)

# Summarise
classrates %>% 
  filter(taxRank == "Genus") %>% 
  group_by(Dataset, proportion_type) %>% 
  mutate(non_classified = 1-Classification_rate) %>% 
  summarise(
    mean_rate = mean(non_classified),
    sd_rate = sd(non_classified),
    .groups = 'drop')


# ASV stats table  -----------------------------

# Rounding functions
fmt_int  <- function(x) format(round(x), big.mark = ",")
# Number formatting: 1 and 2 decimals
fmt_dec <-  function(x) format(round(x, 1), nsmall = 1)
fmt_prop <- function(x) format(round(x, 2), nsmall = 2)

# Stats per object (ps.ls already has all 6 entries -- see load above)
ps.stats <- imap(ps.ls, function(ps, name) {
  asv <- otu_table(ps)
  seq_per_sam <- rowSums(asv)
  asv_per_sam <- rowSums(asv > 0)
  asv_prevalence <- colSums(asv > 0)
  
  tibble(
    Dataset = ps_dataset_labels[[name]],
    Seq  = fmt_int(sum(asv)),
    ASVs = fmt_int(ncol(asv)),
    N    = fmt_int(nrow(asv)),
    
    `Seq per sample`   = paste0(fmt_int(mean(seq_per_sam)), " ± ", fmt_int(sd(seq_per_sam))),
    `Seq range`        = paste0("[", fmt_int(min(seq_per_sam)), " - ", fmt_int(max(seq_per_sam)), "]"),
    
    `ASV per sample`   = paste0(fmt_dec(mean(asv_per_sam)), " ± ", fmt_dec(sd(asv_per_sam))),
    `ASV range`        = paste0("[", fmt_int(min(asv_per_sam)), " - ", fmt_int(max(asv_per_sam)), "]"),
    
    `Prevalence`       = paste0(fmt_dec(mean(asv_prevalence)), " ± ", fmt_dec(sd(asv_prevalence))),
    `Prevalence range` = paste0("[", fmt_int(min(asv_prevalence)), " - ", fmt_int(max(asv_prevalence)), "]")
  )
}) %>% list_rbind()

ps.stats %>% 
  arrange(Dataset) %>% 
  kable("html", align = "l") %>%
  kable_styling(full_width = FALSE) %>%
  add_header_above(c(
    "Dataset" = 1,
    "Reads\n(1,000)" = 1,
    "ASVs" = 1,
    "Samples" = 1,
    "Mean ± SD" = 1,
    "[Min, Max]" = 1,
    "Mean ± SD" = 1,
    "[Min, Max]" = 1,
    "Mean ± SD" = 1,
    "[Min, Max]" = 1
  )) %>%
  add_header_above(c(
    " " = 4,
    "Reads per sample" = 2,
    "ASVs per sample" = 2,
    "ASV prevalence" = 2
  )) %>%
  row_spec(0, extra_css = "display: none;")  %T>%
  
  # NOTE : ITS filtered excludes all non-AMF fungi!
  save_kable(file = "out/manuscript/supp/asv_summary.html")

# Unique taxa per taxonomic rank -----------------------------------------------
# Same taxa-cleaning convention as amy_code/Family_Site.R.

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

# Unique classified taxa per taxonomic rank
get_unique_taxa_list <- function(ps) {
  tax_df <- as.data.frame(as(tax_table(ps), "matrix"), stringsAsFactors = FALSE)
  rank_list <- purrr::map(names(tax_df), function(rank) {
    cleaned <- clean_taxa_vector(tax_df[[rank]])
    unique(cleaned[cleaned != "Unclassified"])
  })
  names(rank_list) <- names(tax_df)
  rank_list[["ASV"]] <- taxa_names(ps) # ASV pseudo-rank
  rank_list
}

standard_rank_order <- c("Phylum", "Class", "Order", "Family", "Genus", "Species_cluster", "ASV")

unique_taxa_summary <- imap(ps.ls, function(ps, name) {
  taxa_list <- get_unique_taxa_list(ps)
  taxa_list <- taxa_list[standard_rank_order]
  tibble(
    Dataset = ps_dataset_labels[[name]],
    Rank = factor(names(taxa_list), levels = intersect(standard_rank_order, names(taxa_list))),
    Unique_Count = purrr::map_int(taxa_list, length)
  )
}) %>%
  list_rbind() %>%
  pivot_wider(names_from = Rank, values_from = Unique_Count)

unique_taxa_summary %>% 
  arrange(Dataset) %>% 
  kable("html", align = "lccccccc") %>%
  kable_styling(full_width = FALSE) %T>%
  save_kable(file = "out/manuscript/supp/unique_taxa_summary.html")

