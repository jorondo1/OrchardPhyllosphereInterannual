pacman::p_load(tidyverse,  here, phyloseq, mgx.tools, kableExtra, magrittr)
ps.ls <- read_rds('data/ps_objects_full.rds')


# Sample count viz --------------

dat <- rbind(
  samdat_as_tibble(ps.ls$Fungi) %>% mutate(Barcode = 'Fungi'),
  samdat_as_tibble(ps.ls$Bacteria) %>% mutate(Barcode = 'Bacteria')
) 


count_dat <- dat %>% 
  mutate(Dataset = case_when(
    Cultivar %in% c('Honeycrisp', 'Spartan') ~ '2-year dataset',
    TRUE ~ '3-year dataset'
  )) %>% 
  count(Barcode, Dataset, Year, Time, Cultivar, Code, Management) %>%
  rename(N_samples = n)

count_dat %>% 
  ggplot(aes(x = Time, y = N_samples, fill = Code)) +
  geom_col(position = "dodge") +
  ggh4x::facet_nested(Barcode+Year ~ Dataset + Management + Cultivar) +  # Facet by 2 variables
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_light()  +
  theme(
    legend.position = 'bottom',
    panel.grid = element_blank()) +
  guides(fill = guide_legend(nrow = 1)) +
  labs(fill = 'Orchard')

ggsave('out/summaries/sample_count_by_metadata.pdf',
       bg = 'white', width = 2200, height = 2000, 
       units = 'px', dpi = 220)

# Classification rates -------------------

classrates <- compute_classification_rates(ps.ls[c('Bacteria','Fungi') ])
plot_class_rates(classrates)
ggsave('out/summaries/classification_rates.pdf',
       bg = 'white', width = 2000, height = 2000, 
       units = 'px', dpi = 220)


# ASV stats table  -----------------------------

# Rounding functions
fmt_int  <- function(x) format(round(x), big.mark = ",")
fmt_dec <-  function(x) format(round(x, 1), nsmall = 1)
fmt_prop <- function(x) format(round(x, 2), nsmall = 2)

# per-dataset 
ps.ls$`Bacteria (3-year subset)` <- ps.ls$Bacteria %>% 
  phyloseq::subset_samples(Dataset == "3-year") %>% 
  prune_taxa(taxa_sums(.)>0, .)

ps.ls$`Bacteria (2-year subset)` <- ps.ls$Bacteria %>% 
  phyloseq::subset_samples(Dataset == "2-year") %>% 
  prune_taxa(taxa_sums(.)>0, .)

ps.ls$`Fungi (3-year subset)` <- ps.ls$Fungi %>% 
  phyloseq::subset_samples(Dataset == "3-year") %>% 
  prune_taxa(taxa_sums(.)>0, .)

ps.ls$`Fungi (2-year subset)` <- ps.ls$Fungi %>% 
  phyloseq::subset_samples(Dataset == "2-year") %>% 
  prune_taxa(taxa_sums(.)>0, .)

# Stats per object
ps.stats <- imap(ps.ls, function(ps, barcode) {
  asv <- otu_table(ps)
  seq_per_sam <- rowSums(asv)
  asv_per_sam <- rowSums(asv > 0)
  asv_prevalence <- colSums(asv > 0)

  tibble(
    Dataset = barcode,
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

kable(ps.stats, "html", align = "l") %>%
  kable_styling(full_width = FALSE) %>%
  add_header_above(c(
    "Dataset" = 1,
    "Sequences" = 1,
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
    "Sequences per sample" = 2,
    "ASVs per sample" = 2,
    "ASV prevalence" = 2
    )) %>%
  row_spec(0, extra_css = "display: none;")  %T>% 
  
  # NOTE : ITS filtered excludes all non-AMF fungi!
  save_kable(file = "out/manuscript/supp/asv_summary.html")

