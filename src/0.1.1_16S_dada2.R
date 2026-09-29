# Script to process 16S sequencing reads from 2026_AppleMicrobiome on ip34/Mammouth
# Author: ANJA WERTZ

# following tutorial: https://jorondo1.github.io/mgx.tutorials/dada2_16S_tutorial.html 

# Set up R in ip34 -------------------------------------------------------

# terminal (bash) commands, run manually on ip34/Mammouth, not R code
# ssh <username>@ip34.ccs.usherbrooke.ca
# newgrp def-ilafores

# Start tmux session -----------------------------------------------------------

# terminal (bash) command, run manually, not R code
# tmux new -s DADA2

# terminal (bash) commands, run manually before launching R, not R code
# module load StdEnv/2023 r/4.4.0 mugqic/cutadapt/2.10
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

path_data <- normalizePath('/net/nfs-bio/jbod2/def-ilafores/analysis/2026_AppleMicrobiome/data/16S')

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


# DADA2 workflow ----------------------------------------------------------

FWD <- "AACMGGATTAGATACCCKG"  # 799F
REV <- "AGGGTTGCGCTCGTTG"     # 1115R

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

# Reprint sample names for ENA
(sample.names <- gsub(
  pattern = ".*\\.([^.]+_*)_R1\\.fastq\\.gz", 
  replacement = "\\1", 
  x = basename(cutFs)))

write_delim(data.frame(sample.names), file.path(path_summary, 'sample_names.tsv'))


# 3. Quality filtering ------------------------------------------------------------

plot_list <- mgx.tools::gen_qplots(Fs = cutFs, Rs = cutRs, nsam = 4)
mgx.tools::save_qplots(plot_list = plot_list,
                       out_dir = path_summary,
                       step_id = 'raw')

filtFs <- file.path(path_data, "3_filtered", basename(cutFs))
filtRs <- file.path(path_data, "3_filtered", basename(cutRs))

names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(
  cutFs, filtFs, cutRs, filtRs,
  maxEE = c(2, 2), 
  truncLen = c(190, 140),
  minLen = 50, 
  rm.phix = TRUE, 
  compress = TRUE,
  multithread = ncores
)

rownames(out) <- sample.names

plot_list <- mgx.tools::gen_qplots(Fs = filtFs, Rs = filtRs, nsam = 4)
mgx.tools::save_qplots(plot_list = plot_list,
                       out_dir = path_summary,
                       step_id = 'clean')

filtFs_survived <- mgx.tools::dropped_samples(filtFs)
filtRs_survived <- mgx.tools::dropped_samples(filtRs)

# 4. Error model and denoising ----------------------------------------------------

errF <- learnErrors(filtFs_survived, multithread = ncores)
errR <- learnErrors(filtRs_survived, multithread = ncores)

err_plot <- plotErrors(errF, nominalQ = TRUE)

ggsave(filename = file.path(path_summary, "error_plot.pdf"),
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

table(nchar(getSequences(seqtab)))

seqtab2 <- seqtab[, nchar(colnames(seqtab)) %in% 290:310]

seqtab.nochim <- removeBimeraDenovo(
  seqtab2, method = "consensus", multithread = ncores, verbose = TRUE
)

sink(file.path(path_summary, "chimera_report.txt"), append = FALSE, split = TRUE)
mgx.tools::chimera_report(seqtab, seqtab.nochim)
sink() # sink writes a report while displaying it on the console

table(nchar(getSequences(seqtab.nochim)))

seqtab.filt <- mgx.tools::drop_rare_asvs(seqtab.nochim, at_least_n = 2)
dim(seqtab.filt)

write_rds(seqtab.filt, file.path(path_out, 'seqtab.RDS'), compress = 'gz')

# 6. Read tracking ----------------------------------------------------------------

track_change <- mgx.tools::track_dada(
  out.N = out.N,
  out = out,
  sample.names = sample.names,
  dadaFs = dadaFs,
  dadaRs = dadaRs, # this was missing in tutorial
  mergers = mergers_pooled,
  seqtab.nochim = seqtab.nochim
)

p_load(patchwork)

mgx.tools::plot_track_change(track_change) %>%
  ggsave(filename = file.path(path_summary, "track_changes.pdf"),
         plot = ., bg = 'white',
         width = 1600, height = 2200, units = 'px', dpi = 180)


# Copy from ip34 to local 

# terminal (bash) command, run manually, not R code
# scp -r <username>@ip34.ccs.usherbrooke.ca:/jbod2/def-ilafores/analysis/2026_AppleMicrobiome/16S/X_DADA2_Summary C:/Users/anjaw/Documents/Canada_UdeS/PhD/A_Apple_microbiome/DADA2

# Copy reads table and taxonomy from ip34 to local -----------------------------

# terminal (bash) command, run manually, not R code
# scp -r <username>@ip34.ccs.usherbrooke.ca:/jbod2/def-ilafores/analysis/2026_AppleMicrobiome/16S/4_out C:/Users/anjaw/Documents/Canada_UdeS/PhD/A_Apple_microbiome/DADA2

# 7. Taxonomic assignment (DECIPHER) ----------------------------------------------

# BiocManager::install("DECIPHER")

# define paths
path_data <- normalizePath('/jbod2/def-ilafores/analysis/2026_AppleMicrobiome/16S')

path_out <- file.path(path_data, '4_out')

seqtab.filt <- read_rds(paste0(path_out, "/seqtab.RDS"))


ref_db_dir <- file.path(path_data, 'ref_taxonomy_DECIPHER')
dir.create(ref_db_dir)

ref_db_path <- file.path(ref_db_dir, "SILVA_SSU_r138_2.RData")


# link to download SILVA SSU r138.2 database

system(paste(
  'wget -O',
  ref_db_path,'"https://drive.usercontent.google.com/download?id=1w3wdSCpSihntWkbP_zvXz7r3s-tNB8DV&export=download&confirm=t&uuid=3d0c7314-9dea-42db-b064-c4d480f47b78"')) 

list.files(ref_db_dir)

# use training set to do assignment
dna <- DNAStringSet(getSequences(seqtab.filt)) 
load(ref_db_path)

ids <- IdTaxa(
  test = dna,
  trainingSet = trainingSet,
  strand = "both", 
  processors = ncores,
  verbose=TRUE)

# Build dada2-formatted taxonomy table (see function above)
taxonomy <- mgx.tools::format_DECIPHER_for_dada2(
  ids, seqtab.filt,
  ranks =  c('rootrank', 'domain', 'phylum', 'class', 'order', 'family', 'genus')
)

write_rds(taxonomy, file.path(path_out, 'taxonomy_DECIPHER.RDS'), compress = 'gz')

# terminal (bash) command, run manually, not R code
# scp -r <username>@ip34.ccs.usherbrooke.ca:/jbod2/def-ilafores/analysis/2026_AppleMicrobiome/16S/4_out/taxonomy_DECIPHER.RDS S:/LaforestLapointeI/ANJA_WERZ/PhD/A_Apple_microbiome/Data/4_out
