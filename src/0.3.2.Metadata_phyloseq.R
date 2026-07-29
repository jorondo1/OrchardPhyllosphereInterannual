# Parse metadata and add it to the phyloseq objects that have trees

# Author: Jonathan Rondeau-Leclaire

pacman::p_load(readxl, tidyverse, phyloseq)
meta_raw <- read_xlsx("data/Meta_interannual.xlsx", sheet = "Combined")

# fix metadata naming and typos ------------------------------------------
meta_formatted <- meta_raw %>% 
  dplyr::mutate(
    # Reorder time variables
    time = factor(time, levels = c('May', 'July')),
    year = factor(year, levels = c(2022, 2023, 2024)),
    #!!!  orchard = factor(orchard, levels = ...),
    
    # Reorder cultivar levels
    cultivar = factor(cultivar, levels = c("Cortland"  , "Honeycrisp" ,"Liberty" , "Spartan", "Paulared")),
    site = ifelse(site == "PMP", "PMB", site),
    # Create Location variable
    Location = case_when(
      site %in% c('PMB', 'COM') ~ 'Compton',
      site %in% c('MIC', 'MIB') ~ 'Milton',
      site == 'ASB' ~ 'Saint-Benoît',
      site == 'VBS' ~ 'Windsor'
    ),
    MANAGEMENT = factor(
      recode(MANAGEMENT, CONV = "Conventional", ORG = "Organic"),
      levels = c('Conventional', 'Organic'))) %>% 
  # flush useless variables
  dplyr::select(-sample, -seq, -replicate, -type, -code) %>% 
  # Consistent variable naming scheme
  dplyr::rename_with(~ stringr::str_to_sentence(.x))

# add meteo data -----------------------------------------------------------

meteo_raw <- read_rds("data/meteo/processed/Meteo_indexes.rds") 

meteo <- meteo_raw %>% 
  separate_wider_delim(Time, delim = "_", names = c('Month', 'Year')) %>% 
  mutate(Time = case_when(Month == '5' ~ 'May', Month == '7' ~ 'July'),
         .keep = 'unused') %>% 
  select(-Location)

meta_out <- meta_formatted %>% 
  left_join(meteo, by = c('Year', 'Time', 'Site')) %>% 
  rename(Sample = Unique)

#  subsets by barcode
ps.ls.in <- read_rds('data/ps_objects_preproc.rds')
ps_bact <- ps.ls.in$Bacteria$filt
ps_fung <- ps.ls.in$Fungi$filt

setdiff(sample_names(ps_bact), meta_out$Sample)
setdiff(sample_names(ps_fung), meta_out$Sample)

meta_fung <- meta_out %>% 
  filter(Sample %in% sample_names(ps_fung)) %>% 
  as.data.frame() %>% column_to_rownames("Sample")

meta_bact <- meta_out %>% 
  filter(Sample %in% sample_names(ps_bact)) %>% 
  as.data.frame() %>% column_to_rownames("Sample")

sample_data(ps_bact) <- meta_bact
sample_data(ps_fung) <- meta_fung

ps.ls.out <- list(
  Bacteria = ps_bact,
  Fungi = ps_fung
); ps.ls.out

# write out ---------------------------------------------------------------

write_rds(
  ps.ls.out, 
  'data/ps_objects_full.rds', 
  compress = 'xz')


# Visualise sample count per metadata combinations:

dat <- rbind(
  samdat_as_tibble(ps.ls$Fungi) %>% mutate(Barcode = 'Fungi'),
  samdat_as_tibble(ps.ls$Bacteria) %>% mutate(Barcode = 'Bacteria')
) 


dat %>% 
  count(Barcode, Year, Time, Site, Cultivar, Orchard, Management) %>%
  rename(N_samples = n) %>% 
  ggplot(aes(x = Time, y = N_samples, fill = Orchard)) +
  geom_col(position = "dodge") +
  ggh4x::facet_nested(Barcode+Year ~ Management + Cultivar) +  # Facet by 2 variables
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_light()  +
  theme(
    legend.position = 'bottom',
    panel.grid = element_blank())


ggsave('out/summaries/sample_count_by_metadata.pdf',
       bg = 'white', width = 2200, height = 2000, 
       units = 'px', dpi = 220)

