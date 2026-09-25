# Parse metadata and add it to the phyloseq objects that have trees

# Author: Jonathan Rondeau-Leclaire

pacman::p_load(readxl, tidyverse, phyloseq, mgx.tools)
meta_raw <- read_xlsx("data/Meta_interannual.xlsx", sheet = "Combined")

# fix metadata naming and typos ------------------------------------------
meta_formatted <- meta_raw %>% 
  dplyr::mutate(
    # Reorder time variables
    time = factor(time, levels = c('May', 'July')),
    year = factor(year, levels = c(2022, 2023, 2024)),
    #!!!  orchard = factor(orchard, levels = ...),
    TREE_ID = paste(year, site, cultivar, replicate, sep = "-"),
    # Reorder cultivar levels
    cultivar = factor(cultivar, levels = c("Cortland"  , "Honeycrisp" ,"Liberty" , "Spartan", "Paulared")),
    # Create Location variable
    Location = as.factor(str_extract(code, "^.")),
    code = as.factor(code),
    MANAGEMENT = factor(
      recode(MANAGEMENT, CONV = "Conventional", ORG = "Organic"),
      levels = c('Conventional', 'Organic')),
    Dataset = case_when(
      cultivar %in% c('Honeycrisp', 'Spartan') ~ '2-year',
      cultivar %in% c('Cortland', 'Liberty', 'Paulared') ~ '3-year'
    )) %>% 
  # flush useless variables
  dplyr::select(-sample, -seq, -replicate, -type) %>% 
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
  rename(Sample = Unique) %>% 
  select(-Site)

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

# add sequencing depth ----------------------

meta_bact$Seq_depth <- rowSums(otu_table(ps_bact))[rownames(meta_bact)]
meta_fung$Seq_depth <- rowSums(otu_table(ps_fung))[rownames(meta_fung)]

# add centered/scaled covariates (hiermod models, src/hiermod/0_SETUP.R) -----
# deg_h_z: centered WITHIN Time (July is reliably warmer than May every
# year, a real seasonal identity worth keeping separate from the
# covariate). precip_72h_z/seq_depth_z: centered globally -- precip has no
# reliable May-vs-July direction (e.g. 2022 reverses it), and seq_depth's
# Management correlation is a confound we want removed, not preserved.
# Keep this logic in sync with 0_SETUP.R's own if it ever changes.

add_centered_covariates <- function(dat){
  deg_h_season_mean <- tapply(dat$deg_h, dat$Time, mean)
  deg_h_z <- (dat$deg_h - deg_h_season_mean[dat$Time])
  dat$deg_h_z <- deg_h_z / sd(deg_h_z)

  dat$precip_72h_z <- (dat$precip_72h - mean(dat$precip_72h)) / sd(dat$precip_72h)

  log_seq_depth   <- log(dat$Seq_depth)
  dat$seq_depth_z <- (log_seq_depth - mean(log_seq_depth)) / sd(log_seq_depth)

  dat
}

meta_bact <- add_centered_covariates(meta_bact)
meta_fung <- add_centered_covariates(meta_fung)

# build final objects ------------------------

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

