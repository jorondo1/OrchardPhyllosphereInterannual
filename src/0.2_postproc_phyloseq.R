# Build phyloseq objects (16S & ITS) from DADA2 outputs
# 
# Authors: Amy Heim, Anja Werz, Jonathan Rondeau-Leclaire, 
#
# For each dataset we only keep two objects:
#   ps_raw_*  - straight from DADA2, nothing removed
#   ps_filt_* - taxonomy-filtered, rare ASVs dropped, low-depth samples dropped

pacman::p_load(tidyverse, phyloseq, magrittr, vegan)

# Paths -----------------------------------------------------------------

path_in  <- "data/dada2_out"
path_out <- "data/phyloseq"
path_summary <- "out/summaries"
dir.create(path_summary, recursive = TRUE, showWarnings = FALSE)

# Helper functions ------------------------------------------------------

# print+log a mean/median/sd summary line for a numeric vector
report_stats <- function(x, label) {
  cat(sprintf("%-45s mean=%.1f  median=%.1f  sd=%.1f\n", label, mean(x), median(x), sd(x)))
}

# rarecurve() on the ASV table, saved as rds + plot; returns the tidy curve
plot_rarecurve <- function(seqtab, thresholds, dataset) {
  rc <- rarecurve(seqtab, tidy = TRUE, step= 100)
  
  rc_exp <- rc %>%
    separate(Site, into = c("Year", "Timepoint", "Site", "Variety", "Tissue", "Replicate")) %>%
    mutate(Unique = paste0(Year, Timepoint, Site, Variety, Tissue, Replicate))
  
  p <- ggplot(rc_exp, aes(x = Sample, y = Species, color = Unique)) +
    geom_line() +
    geom_vline(xintercept = thresholds) +
    facet_wrap(~Site) +
    theme(legend.position = "none")
  print(p)
  ggsave(file.path(path_summary, paste0("rarecurve_", dataset, ".png")), p, dpi = 300)
  invisible(rc)
}

# open/close the per-dataset console+file summary sink (mirrors the sink()
# same pattern as chimera report in DADA2 scripts)
sink_16S <- function(append = TRUE) sink(file.path(path_summary, "16S_summary.txt"), append = append, split = TRUE)
sink_ITS <- function(append = TRUE) sink(file.path(path_summary, "ITS_summary.txt"), append = append, split = TRUE)

# reset both summary files
sink_16S(append = FALSE); cat("== 16S: DADA2 -> phyloseq summary ==\n"); sink()
sink_ITS(append = FALSE); cat("== ITS: DADA2 -> phyloseq summary ==\n"); sink()

##############################################################################
# 1. Load DADA2 outputs -------------------------------------------------
##############################################################################

# --- 16S ---
seqtab_16S <- read_rds(file.path(path_in, "seqtab_16S.RDS")) %>% as.data.frame()
rownames(seqtab_16S) <- rownames(seqtab_16S) %>%
  #gsub("PMP", "PMB", .) %>% # we keep the errors in sample names, we'll fix metadata
  gsub("_16S$", "", .)  # remove the 16S id

tax_16S <- read_rds(file.path(path_in, "taxonomy_16S_DECIPHER.RDS")) %>% 
  as.data.frame() %>% 
  select(-rootrank)

# dirty fix
colnames(tax_16S) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")

# 
tax_16S %<>% 
  mutate(across(everything(), ~ case_when(
    . == "Unassigned" | str_detect(., 'Incertae') | is.na(.) ~ "Unclassified", 
    TRUE ~ .))) 



# --- ITS ---
seqtab_ITS <- read_rds(file.path(path_in, "seqtab_ITS.RDS")) %>% as.data.frame()
rownames(seqtab_ITS) <- rownames(seqtab_ITS) %>% gsub("_ITS$", "", .)

tax_ITS <- read_rds(file.path(path_in, "taxonomy_ITS_DECIPHER.RDS")) %>% as.data.frame()
tax_ITS %<>% mutate(across(everything(), ~ case_when(
  . == "Unassigned" | str_detect(., 'Incertae') | is.na(.) ~ "Unclassified", 
  TRUE ~ .))) 

##############################################################################
# 2. Raw phyloseq objects (nothing filtered/removed) --------------------
##############################################################################

# --- 16S ---
ps_raw_16S <- phyloseq(
  otu_table(as.matrix(seqtab_16S), taxa_are_rows = FALSE),
  tax_table(as.matrix(tax_16S))
)

# --- ITS ---
ps_raw_ITS <- phyloseq(
  otu_table(as.matrix(seqtab_ITS), taxa_are_rows = FALSE),
  tax_table(as.matrix(tax_ITS))
)

##############################################################################
# 3. Taxonomy-based & rare-ASV filtering ---------------------------------
##############################################################################

n_rare <- 10 # minimum total reads for an ASV to be kept

# --- 16S: keep Bacteria, drop Mitochondria, drop ASVs with < n_rare reads ---

sink_16S()
cat("\n--- Step 3: ASV filtering ---\n")

tax_16S_filt <- tax_16S %>%
  filter(Kingdom == "Bacteria", Family != "Mitochondria")
seqtab_16S_filt <- seqtab_16S[, colnames(seqtab_16S) %in% rownames(tax_16S_filt)]

cat(sprintf("Total reads (bacteria only): %d\n", sum(seqtab_16S_filt)))
report_stats(colSums(seqtab_16S_filt), "ASV read count (bacteria only)")
report_stats(rowSums(seqtab_16S_filt != 0), "ASVs per sample")
report_stats(colSums(seqtab_16S_filt != 0), "Samples per ASV")

# negative control check: which ASVs show up there, and how abundant are they elsewhere
neg_ctrl_asv <- colnames(seqtab_16S_filt)[seqtab_16S_filt["ctrl_PCR_neg_LaforestP10", ] != 0]
cat(sprintf("ASVs present in negative control: %d\n", length(neg_ctrl_asv)))
print(colSums(seqtab_16S_filt[, neg_ctrl_asv, drop = FALSE]))

keep_asv <- names(which(colSums(seqtab_16S_filt) >= n_rare))
cat(sprintf("Kept %d / %d ASVs with >= %d reads\n", length(keep_asv), ncol(seqtab_16S_filt), n_rare))
seqtab_16S_filt <- seqtab_16S_filt[, keep_asv]
tax_16S_filt <- tax_16S_filt[keep_asv, ]

sink()

# --- ITS: keep Fungi with >= n_rare reads in one step ---

sink_ITS()
cat("\n--- Step 3: ASV length distribution (pre-filter) ---\n")

length_dist <- tibble(
  length = nchar(sort(rownames(tax_ITS))),
  sequences = sort(colSums(seqtab_ITS))
)
p_len <- ggplot(length_dist, aes(x = length)) + geom_histogram() + theme_minimal()
print(p_len)
ggsave(file.path(path_summary, "asv_length_dist_ITS.png"), p_len, dpi = 300)

p_len_abd <- length_dist %>%
  group_by(length) %>%
  summarise(abundance = sum(sequences)) %>%
  ggplot(aes(x = length, y = abundance)) + geom_col() + theme_minimal()
print(p_len_abd)
ggsave(file.path(path_summary, "asv_abd_length_dist_ITS.png"), p_len_abd, dpi = 300)

cat("\n--- Step 3: ASV filtering ---\n")

asv_reads <- colSums(seqtab_ITS)
cat(sprintf("Total reads: %d\n", sum(seqtab_ITS)))
report_stats(asv_reads, "ASV read count (pre-filter)")
report_stats(rowSums(seqtab_ITS != 0), "ASVs per sample (pre-filter)")
report_stats(colSums(seqtab_ITS != 0), "Samples per ASV (pre-filter)")

keep_asv <- rownames(tax_ITS)[tax_ITS$Kingdom == "Fungi" & asv_reads >= n_rare]
cat(sprintf("Kept %d / %d ASVs (Fungi, >= %d reads)\n", length(keep_asv), ncol(seqtab_ITS), n_rare))

tax_ITS_filt <- tax_ITS[keep_asv, ]
seqtab_ITS_filt <- seqtab_ITS[, keep_asv]

sink()

##############################################################################
# 4. Rarefaction curves ----------------------------------------------------
##############################################################################

plot_rarecurve(seqtab_16S_filt, thresholds = c(7000, 8000), dataset = "16S")
plot_rarecurve(seqtab_ITS_filt, thresholds = c(2500, 4000), dataset = "ITS")

##############################################################################
# 5. Drop low-depth samples ----------------------------------------------
##############################################################################

# --- 16S ---
sink_16S()
cat("\n--- Step 5: sequencing depth ---\n")

seq_depth <- rowSums(seqtab_16S_filt)
report_stats(seq_depth, "Sequencing depth (pre sample-filter)")

seq_min_16S <- 7000
keep_samples <- names(which(seq_depth >= seq_min_16S))
cat(sprintf("Kept %d / %d samples with >= %d reads\n", length(keep_samples), nrow(seqtab_16S_filt), seq_min_16S))
seqtab_16S_filt <- seqtab_16S_filt[keep_samples, ]

sink()

# --- ITS ---
sink_ITS()
cat("\n--- Step 5: sequencing depth ---\n")

seq_depth <- rowSums(seqtab_ITS_filt)
report_stats(seq_depth, "Sequencing depth (pre sample-filter)")

seq_min_ITS <- 3000
keep_samples <- names(which(seq_depth >= seq_min_ITS))
cat(sprintf("Kept %d / %d samples with >= %d reads\n", length(keep_samples), nrow(seqtab_ITS_filt), seq_min_ITS))
seqtab_ITS_filt <- seqtab_ITS_filt[keep_samples, ]

sink()

##############################################################################
# 6. Filtered phyloseq objects -------------------------------------------
##############################################################################

# --- 16S ---
sink_16S()
cat("\n--- Step 6: filtered phyloseq object ---\n")

seqtab_16S_filt <- seqtab_16S_filt[, colSums(seqtab_16S_filt) != 0] # drop now-empty ASVs
tax_16S_filt <- tax_16S_filt[colnames(seqtab_16S_filt), ]

report_stats(rowSums(seqtab_16S_filt), "Sequencing depth (final)")
report_stats(colSums(seqtab_16S_filt != 0), "Samples per ASV (final)")
cat(sprintf("Final table: %d samples x %d ASVs\n", nrow(seqtab_16S_filt), ncol(seqtab_16S_filt)))

ps_filt_16S <- phyloseq(
  otu_table(as.matrix(seqtab_16S_filt), taxa_are_rows = FALSE),
  tax_table(as.matrix(tax_16S_filt))
)

sink()

# --- ITS ---
sink_ITS()
cat("\n--- Step 6: filtered phyloseq object ---\n")

seqtab_ITS_filt <- seqtab_ITS_filt[, colSums(seqtab_ITS_filt) != 0] # drop now-empty ASVs
tax_ITS_filt <- tax_ITS_filt[colnames(seqtab_ITS_filt), ]

report_stats(rowSums(seqtab_ITS_filt), "Sequencing depth (final)")
report_stats(colSums(seqtab_ITS_filt != 0), "Samples per ASV (final)")
cat(sprintf("Final table: %d samples x %d ASVs\n", nrow(seqtab_ITS_filt), ncol(seqtab_ITS_filt)))

ps_filt_ITS <- phyloseq(
  otu_table(as.matrix(seqtab_ITS_filt), taxa_are_rows = FALSE),
  tax_table(as.matrix(tax_ITS_filt))
)

sink()

##############################################################################
# 7. Trees and species clusters ----------------------------------------------
##############################################################################

# Build NJ trees
ps.16S.tree <- mgx.tools::ASV_tree_for_physeq(ps_filt_16S, ncores = 7)
ps.ITS.tree <- mgx.tools::ASV_tree_for_physeq(ps_filt_ITS, ncores = 7)

# Species-level clusters
ps.16S.clust <- mgx.tools::cluster_ASVs_physeq(ps.16S.tree, threshold = 0.03)
ps.ITS.clust <- mgx.tools::cluster_ASVs_physeq(ps.ITS.tree, threshold = 0.0295)

# can't fully resolve 16S inconsistencies in clusters. 
# ARe they important in terms of reelative abundance across samples? 
# 
# ps.16S.clust %>% 
#   #tax_glom2(taxrank = "Species_cluster") %>% 
#   psflashmelt() %>% 
#   select(Sample, Abundance, Species_cluster, Genus) %>% 
#   group_by(Sample) %>% 
#   mutate(relAb = Abundance / sum(Abundance)) %>% 
#   ungroup() %>% 
#   filter(Species_cluster %in% `_problem_clusters`) %>% 
#   group_by(Species_cluster, Genus) %>% 
#   filter(relAb>0) %>% 
#   summarise(mean = mean(relAb), n = n(), max = max(relAb)) %>% 
#   arrange(desc(Species_cluster)) %>% 
#   print(n=100)

# The only "important" ones (according to mean, n or max) are Unclassified
# Let's keep them, with inconsistencies

##############################################################################
# 8. Export objects ---------------------------------------------------------
##############################################################################

list.out <- list(
  Bacteria = list(
    raw = ps_raw_16S,
    filt = ps.16S.clust
  ),
  Fungi = list(
    raw = ps_raw_ITS,
    filt = ps.ITS.clust
  )
)

write_rds(list.out, file.path(path_out, "ps_objects_preproc.rds"), compress = "xz")
