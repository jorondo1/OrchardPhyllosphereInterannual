### ANCOM-BC2 PIPELINE               

pacman::p_load(
  mgx.tools,
  tidyverse,
  ANCOMBC,
  parallel,
  phyloseq,
  update = FALSE
)

source('src/0.0_Config.R')

# 1. PREPARE DATA -----------

# 2yr/3yr subsets already carried in ps_objects_full.rds (src/0.3.2.Metadata_phyloseq.R)
ps.ls  <- readRDS('data/ps_objects_full.rds')
Years3 <- list(Bacteria = ps.ls$Bacteria_3y, Fungi = ps.ls$Fungi_3y)
Years2 <- list(Bacteria = ps.ls$Bacteria_2y, Fungi = ps.ls$Fungi_2y)

ps_all <- list(
  Years3 = Years3,
  Years2 = Years2
)


# Aggregate to species clusters

ps_final <- lapply(
  ps_all,
  function(ps_list) {
    lapply(ps_list, tax_glom2, taxrank = "Species_cluster")
  }
)


# Relab filter
rel_abund_cutoff <- 0.0001 

# Keep taxa whose max relative abundance in any sample >= cutoff
filter_relative_abundance <- function(ps, cutoff = 0.0001) {
  ps <- prune_taxa(taxa_sums(ps) > 0, ps)
  ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))
  otu <- otu_table(ps_rel)
  if (!taxa_are_rows(ps_rel)) otu <- t(otu)
  max_rel_abund <- apply(otu, 1, max, na.rm = TRUE)
  keep <- max_rel_abund >= cutoff
  prune_taxa(keep, ps)
}
ps_final <- lapply(
  ps_final,
  function(ps_list) {
    lapply(ps_list, filter_relative_abundance, cutoff = rel_abund_cutoff)
  }
)


# 2. MODEL SPECIFICATION -------------------


ModRun    <- "Time + Management + seq_depth_z"
RandomRun <- "(1|Tree_id)"

# Dose not effect model, its for table and Figuer 
group_variables <- c("Time", "Management")


# 5. ANCOM-BC2 RUNNER

run_ancom <- function(ps) {
  sample_df <- data.frame(sample_data(ps))
  model_variables <- all.vars(as.formula(paste("~", ModRun)))
  for (var in model_variables) {
    if (is.character(sample_df[[var]]) || is.factor(sample_df[[var]])) {
      sample_df[[var]] <- factor(sample_df[[var]])
    }
  }
  sample_df$Seq_depth <- as.numeric(scale(sample_df$Seq_depth))
  sample_df$Location  <- factor(sample_df$Location)
  sample_df$Year      <- factor(sample_df$Year)
  sample_data(ps) <- sample_df
  ancombc2(
    data        = ps,
    prv_cut     = 0.20,
    fix_formula = ModRun,
    rand_formula = RandomRun,
    struc_zero  = FALSE,
    neg_lb      = FALSE,
    dunnet      = FALSE,
    alpha       = 0.01,
    verbose     = TRUE,
    n_cl        = parallel::detectCores()
  )
}


# 6. RUN ANCOM-BC2 (loop for Dataset and micro)

ancom_results <- lapply(
  names(ps_final),
  function(year_dataset) {
    current_ps <- ps_final[[year_dataset]]
    dataset_results <- lapply(
      names(current_ps),
      function(dataset_name) {
        cat(
          "\n####", year_dataset, "-", dataset_name,
          "| fix_formula:", ModRun,
          "| rand_formula:", RandomRun, "####\n"
        )
        run_ancom(current_ps[[dataset_name]])
      }
    )
    names(dataset_results) <- names(current_ps)
    dataset_results
  }
)
names(ancom_results) <- names(ps_final)

# 7. SAVE RAW RESULTS


ANCOMResults <- list(
  ancom_results   = ancom_results,
  ps_final        = ps_final,
  ps_all          = ps_all,
  ModRun          = ModRun,
  RandomRun       = RandomRun,
  group_variables = group_variables,
  prv_cut         = 0.20,
  struc_zero      = FALSE,
  neg_lb          = FALSE,
  dunnet          = FALSE,
  alpha           = 0.01,
  rel_abund_cutoff = rel_abund_cutoff
)

saveRDS(ANCOMResults, file = "data/ANCOMResults.rds")
