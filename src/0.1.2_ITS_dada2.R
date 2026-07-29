# Script to process 16S sequencing reads from 2026_AppleMicrobiome on ip34/Mammouth
# Author: Anja Werz

######################
###Set up R in ip34###
######################

# nice R
path_cutadapt <- '/cvmfs/soft.mugqic/CentOS6/software/cutadapt/cutadapt-2.10/bin/cutadapt'
ncores <- 24

# install.packages(c("devtools","BiocManager", "tidyverse"))
# library(devtools)
# library(BiocManager)
# BiocManager::install(c("Biostrings", "ShortRead"))
# devtools::install_github(c("trinker/pacman", "jorondo1/mgx.tools"))

pacman::p_load(dada2, tidyverse, mgx.tools, Biostrings, ShortRead, parallel, update = FALSE)


# Define paths -----------------------------------------------------------------

path_data <- normalizePath('/jbod2/def-ilafores/analysis/2026_AppleMicrobiome/ITS')

path_raw <- file.path(path_data, '0_raw')
if (!dir.exists(path_raw))     dir.create(path_raw, recursive = TRUE)

path_filt <- file.path(path_data, '1_filtN')
path_cut <- file.path(path_data, '2_cutadapt')
path_out <- file.path(path_data, '4_out')
path_summary <- file.path(path_data, 'X_DADA2_Summary')

if (!dir.exists(path_filt))    dir.create(path_filt)
if (!dir.exists(path_cut))     dir.create(path_cut)
if (!dir.exists(path_out))     dir.create(path_out)
if (!dir.exists(path_summary)) dir.create(path_summary)

###########################
### DADA2 WORKFLOW ###
###########################

# add primers, taken from Sophie ITS file

FWD <- "CTTGGTCATTTAGAGGAAGTAA" # ITS1F
REV <- "GCTGCGTTCTTCATCGATGC" # ITS2 

# Get files and sample names ---------------------------------------------------

fnFs <- sort(list.files(path_raw, pattern = "_R1\\.fastq\\.gz$", full.names = TRUE))
fnRs <- sort(list.files(path_raw, pattern = "_R2\\.fastq\\.gz$", full.names = TRUE))

(sample.names <- gsub(
  pattern = ".*\\.([^.]+_*)_R1\\.fastq\\.gz", 
  replacement = "\\1", 
  x = basename(fnFs)))

write_delim(data.frame(sample.names), file.path(path_summary, 'sample_names.tsv'))

# 1. N-filtering ------------------------------------------------------------------

fnFs.filtN <- file.path(path_data, "1_filtN", basename(fnFs))
fnRs.filtN <- file.path(path_data, "1_filtN", basename(fnRs))

out.N <- filterAndTrim(
  fnFs, fnFs.filtN,
  fnRs, fnRs.filtN,
  maxN = 0,         # reject any read containing an N
  rm.lowcomplex = TRUE,    # remove low-complexity reads (see note below)
  multithread = ncores
)

names(fnFs.filtN) <- sample.names
names(fnRs.filtN) <- sample.names
rownames(out.N) <- sample.names

head(out.N)

# 2. Primer removal ---------------------------------------------------------------

mgx.tools::primer_occurence(fnFs.filtN, fnRs.filtN, FWD, REV, ncores = ncores)

fnFs.cut <- file.path(path_cut, basename(fnFs))
fnRs.cut <- file.path(path_cut, basename(fnRs))

mclapply(
  seq_along(fnFs),
  mgx.tools::run_cutadapt,
  cutadapt_path = path_cutadapt,
  # Remove FWD primer and reverse-complement of REV from R1
  R1.flags = paste("-g", FWD, "-a", dada2:::rc(REV)),
  # Remove REV primer and reverse-complement of FWD from R2
  R2.flags = paste("-G", REV, "-A", dada2:::rc(FWD)),
  mc.cores = ncores
)

mgx.tools::primer_occurence(fnFs.cut, fnRs.cut, FWD, REV, ncores = ncores)

cutFs <- sort(list.files(path_cut, pattern = "_R1.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path_cut, pattern = "_R2.fastq.gz", full.names = TRUE))

plot_list <- mgx.tools::gen_qplots(Fs = cutFs, Rs = cutRs, nsam = 4)
mgx.tools::save_qplots(plot_list = plot_list,
                       out_dir = path_summary,
                       step_id = 'raw')

# 3. Quality filtering ------------------------------------------------------------

filtFs <- file.path(path_data, "3_filtered_maxEE45", basename(cutFs))
filtRs <- file.path(path_data, "3_filtered_maxEE45", basename(cutRs))

names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(
  cutFs, filtFs, cutRs, filtRs,
  maxEE = c(4, 5), 
  minLen = 50, 
  rm.phix = TRUE, 
  compress = TRUE,
  multithread = ncores
)

rownames(out) <- sample.names

#### >>> tests ---
reads_dropped45 <- out %>% 
  as.data.frame() %>% 
  rownames_to_column('Sample') %>% 
  mutate(
    change = (reads.in-reads.out)/reads.in
  ) %>% 
  summarise(
    min = min(change),
    mean = mean(change),
    median = median(change),
    sd = sd(change),
    max = max(change)
  ) %>% 
  mutate(
    test = 'maxEE45', .before = everything()
  )
rbind(reads_dropped, reads_dropped44, reads_dropped45)
#### /// tests ---

plot_list <- gen_qplots(Fs = filtFs, Rs = filtRs, nsam = 4)
save_qplots(plot_list = plot_list,
            out_dir = path_summary,
            step_id = 'clean_maxEE45')

filtFs_survived <- mgx.tools::dropped_samples(filtFs)
filtRs_survived <- mgx.tools::dropped_samples(filtRs)

# 4. Error model and denoising ----------------------------------------------------

errF <- learnErrors(filtFs_survived, multithread = ncores)
errR <- learnErrors(filtRs_survived, multithread = ncores)

# Visualise the fitted error model — this is an important QC step
err_plot <- plotErrors(errF, nominalQ = TRUE)

ggsave(filename = file.path(path_summary, "error_plot_maxE45.pdf"),
       plot = err_plot, bg = 'white',
       width = 2000, height = 2000, units = 'px', dpi = 180)

dadaFs <- dada(filtFs_survived, err = errF, pool = 'pseudo', multithread = ncores)
dadaRs <- dada(filtRs_survived, err = errR, pool = 'pseudo', multithread = ncores)

# 5. Create sequence table --------------------------------------------------------

mergers_pooled <- mergePairs(
  dadaFs, filtFs_survived,
  dadaRs, filtRs_survived,
  minOverlap = 12, # minimum overlap in bp; increase if you have long amplicons
  maxMismatch = 0 # the default; ASV with overlap mismatchs will be discarded
)

seqtab <- makeSequenceTable(mergers_pooled)

seqtab.nochim <- removeBimeraDenovo(
  seqtab, method = "consensus", multithread = ncores, verbose = TRUE
)

sink(file.path(path_summary, "chimera_report.txt"), append = FALSE, split = TRUE)
mgx.tools::chimera_report(seqtab, seqtab.nochim)
sink() # sink writes a report while displaying it on the console

# Show distribution:
table(nchar(getSequences(seqtab.nochim)))

# inspect abundance of lengths
tibble(
  abund=colSums(seqtab.nochim),
  insert_size = nchar(getSequences(seqtab.nochim))
) %>% 
  group_by(insert_size) %>% 
  summarise(total_abund = sum(abund)) %>%
  mutate(total_abund = round(100*total_abund/sum(total_abund),2)) %>% 
  arrange(insert_size) %>% 
  mutate(insert_size = paste0("L",insert_size)) %>% 
  deframe()
# shows % of ASVs of any length L, *weighted by abundance*

# Subset within target range:
seqtab.nochim_filtered <- seqtab.nochim[, nchar(colnames(seqtab.nochim)) %in% 199:282] 

table(nchar(getSequences(seqtab.nochim)))  # length distribution after chimera removal

# Drop ASVs with fewer than 2 total reads across all samples:
seqtab.filt <- mgx.tools::drop_rare_asvs(seqtab.nochim_filtered, at_least_n = 2)
dim(seqtab.filt)

write_rds(seqtab.filt, file.path(path_out, 'seqtab_ITS.RDS'), compress = 'gz')

# 6. Read tracking ----------------------------------------------------------------

track_change <- mgx.tools::track_dada(
  out.N = out.N,
  out = out,
  sample.names = sample.names,
  dadaFs = dadaFs,
  dadaRs = dadaRs,
  mergers = mergers_pooled,
  seqtab.nochim = seqtab.nochim_filtered
)

p_load(patchwork)

mgx.tools::plot_track_change(track_change) %>%
  ggsave(filename = file.path(path_summary, "track_changes_maxEE45.pdf"),
         plot = ., bg = 'white',
         width = 1600, height = 2200, units = 'px', dpi = 180)

# 7. Taxonomic assignment (DECIPHER) ----------------------------------------------

ref_db_dir <- file.path(path_data, 'ref_taxonomy')
if (!dir.exists(ref_db_dir)) {dir.create(ref_db_dir, recursive = TRUE)}

# name the file appropriately, ending with .RData , for example:
ref_db_path <- file.path(ref_db_dir, "UNITE_v2025.RData")

# This calls the bash command "wget" which downloads files from a url:
system(paste(
  'wget -O',
  ref_db_path, 
  # then REPLACE this path with the address you copied!
  '"https://drive.usercontent.google.com/download?id=1wFob94wNna6RAdYeQudoMJILXPELxWPM&export=download&authuser=0&confirm=t&uuid=fbd67f4d-69b1-4448-af76-d1daafdfdbab&at=AAINaII7Lu2hYh3OElhhu5Mc0z4V%3A1781275116620"')) 

# check the filename:
list.files(ref_db_dir)

# Extract ASVs
dna <- DNAStringSet(getSequences(seqtab.filt)) 
load(ref_db_path)

# Determine taxonomy
ids <- DECIPHER::IdTaxa(
  test = dna,
  trainingSet = trainingSet,
  strand = "both", 
  processors = ncores,
  verbose=TRUE)

# Build dada2-formatted taxonomy table (see function above)
taxonomy <- mgx.tools::format_DECIPHER_for_dada2(
  ids, seqtab.filt,
  ranks =  c('rootrank', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus')
)

write_rds(taxonomy, file.path(path_out, 'taxonomy_ITS_DECIPHER.RDS'), compress = 'gz')









