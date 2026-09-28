# ==========================================================
# BLAST TAXONOMY BY SPECIES CLUSTER
# ==========================================================
 
# ==========================================================
# 1. LOAD PACKAGES
# ==========================================================

library(pacman)

p_load(
  tidyverse,
  phyloseq,
  Biostrings
)

library(conflicted)

conflicts_prefer(
  dplyr::slice,
  dplyr::rename,
  dplyr::filter,
  dplyr::select,
  dplyr::count,
  dplyr::first,
  dplyr::mutate,
  dplyr::arrange
)


# ==========================================================
# SELECT FUNGI SPECIES CLUSTERS FOR BLAST
# ==========================================================

ps <- ps_objects_full$Fungi

# Species cluster numbers to export
clusters_to_blast <- c(
  4, 5, 7, 14, 15, 17
)

# ==========================================================
# EXTRACT FUNGAL TAXONOMY
# ==========================================================

tax_df <- as.data.frame(
  tax_table(ps)
) %>%
  rownames_to_column("ASV")

# ==========================================================
# IDENTIFY SPECIES CLUSTER NUMBER
# ==========================================================

tax_df <- tax_df %>%
  mutate(
    Species_cluster_number =
      as.numeric(Species_cluster)
  )

# ==========================================================
# KEEP ONLY THE SELECTED SPECIES CLUSTERS
# ==========================================================

selected_taxa <- tax_df %>%
  filter(
    Species_cluster_number %in% clusters_to_blast
  )


cat(
  "Number of ASVs selected:",
  nrow(selected_taxa),
  "\n"
)


cat(
  "Number of Species clusters selected:",
  n_distinct(selected_taxa$Species_cluster),
  "\n"
)

# ==========================================================
# CREATE FASTA SEQUENCES
# ==========================================================

selected_sequences <- DNAStringSet(
  selected_taxa$ASV
)


# Use Species cluster + ASV number as FASTA names
names(selected_sequences) <- paste0(
  selected_taxa$Species_cluster,
  "_ASV",
  seq_len(nrow(selected_taxa))
)


# ==========================================================
# WRITE FASTA FILE
# ==========================================================

writeXStringSet(
  selected_sequences,
  filepath = "Fungi_Blast.fasta"
)

# ==========================================================
# 2. LOAD FILES AFTER BLAST
# ==========================================================

blast_file <- "Fungi_Fin_Blast.txt"

fasta_file <- "Fungi_Blast.fasta"

# ==========================================================
# 3. READ FASTA
# ==========================================================

seqs <- readDNAStringSet(
  fasta_file
)

fasta_df <- tibble(
  FASTA_ID = names(seqs),
  Sequence = as.character(seqs)
)

cat(
  "Number of FASTA sequences:",
  nrow(fasta_df),
  "\n"
)

# ==========================================================
# 4. EXTRACT PHYLOSEQ TAXONOMY
# ==========================================================

tax_df <- as.data.frame(
  tax_table(ps)
) %>%
  rownames_to_column(
    "ASV"
  )

tax_df$Sequence <- tax_df$ASV

cat(
  "Number of phyloseq taxa:",
  nrow(tax_df),
  "\n"
)


# ==========================================================
# 5. READ BLAST OUTPUT
# ==========================================================

blast_lines <- readLines(
  blast_file,
  warn = FALSE
)

cat(
  "Number of BLAST lines:",
  length(blast_lines),
  "\n"
)

# ==========================================================
# 6. FIND ALL BLAST QUERIES
# ==========================================================

query_lines <- grep(
  "^Query #[0-9]+:",
  blast_lines
)

cat(
  "Number of BLAST queries:",
  length(query_lines),
  "\n"
)


# ==========================================================
# 8. PARSE ALL BLAST HITS
# ==========================================================

all_hits <- list()

hit_counter <- 1


for (i in seq_along(query_lines)) {
  
  
  # --------------------------------------------------------
  # Determine boundaries of this query
  # --------------------------------------------------------
  
  start <- query_lines[i]
  
  
  if (i < length(query_lines)) {
    
    end <- query_lines[i + 1] - 1
    
  } else {
    
    end <- length(blast_lines)
    
  }
  
  
  block <- blast_lines[
    start:end
  ]
  
  
  # --------------------------------------------------------
  # Extract query sequence
  # --------------------------------------------------------
  
  query_line <- block[1]
  
  
  query_sequence <- sub(
    "^Query #[0-9]+: ",
    "",
    query_line
  )
  
  
  query_sequence <- sub(
    " Query ID:.*$",
    "",
    query_sequence
  )
  
  
  query_sequence <- trimws(
    query_sequence
  )
  
  
  # --------------------------------------------------------
  # Extract Query ID
  # --------------------------------------------------------
  
  query_id <- NA_character_
  
  
  if (
    grepl(
      "Query ID:",
      query_line
    )
  ) {
    
    query_id <- sub(
      ".*Query ID:\\s*([^ ]+).*",
      "\\1",
      query_line
    )
    
  }
  
  
  # --------------------------------------------------------
  # Find BLAST summary table
  # --------------------------------------------------------
  
  sig_line <- grep(
    "^Sequences producing significant alignments:",
    block
  )
  
  
  if (
    length(sig_line) == 0
  ) {
    
    next
    
  }
  
  
  sig_line <- sig_line[1]
  
  
  # --------------------------------------------------------
  # Find end of summary table
  # --------------------------------------------------------
  
  alignment_line <- grep(
    "^Alignments:",
    block
  )
  
  
  if (
    length(alignment_line) > 0
  ) {
    
    alignment_line <- alignment_line[1]
    
    
    hit_lines <- block[
      (sig_line + 1):(alignment_line - 1)
    ]
    
  } else {
    
    hit_lines <- block[
      (sig_line + 1):length(block)
    ]
    
  }
  
  
  # --------------------------------------------------------
  # Remove empty lines
  # --------------------------------------------------------
  
  hit_lines <- hit_lines[
    trimws(hit_lines) != ""
  ]
  
  
  # --------------------------------------------------------
  # Remove table headers
  # --------------------------------------------------------
  
  hit_lines <- hit_lines[
    !grepl(
      "^Description|^Scientific|^Name|^Common|^Max Score|^Total Score",
      trimws(hit_lines),
      ignore.case = TRUE
    )
  ]
  
  
  # --------------------------------------------------------
  # Parse every BLAST hit
  # --------------------------------------------------------
  
  for (
    hit_line in hit_lines
  ) {
    
    
    hit_line <- trimws(
      hit_line
    )
    
    
    # ------------------------------------------------------
    # Extract the numeric BLAST fields from the END
    # of each summary-table row.
    #
    # Example:
    #
    # Description
    # Max Score
    # Total Score
    # Query Cover
    # E-value
    # Percent Identity
    # Length
    # Accession
    # ------------------------------------------------------
    
    parsed <- str_match(
      
      hit_line,
      
      paste0(
        "^(.+?)\\s+",
        "([0-9.]+)\\s+",
        "([0-9.]+)\\s+",
        "([0-9]+)%\\s+",
        "([0-9.eE+-]+)\\s+",
        "([0-9.]+)\\s+",
        "([0-9]+)\\s+",
        "(\\S+)$"
      )
      
    )
    
    
    # ------------------------------------------------------
    # Skip rows that do not match the BLAST format
    # ------------------------------------------------------
    
    if (
      is.na(parsed[1, 1])
    ) {
      
      next
      
    }
    
    
    # ------------------------------------------------------
    # Extract fields
    # ------------------------------------------------------
    
    description <- parsed[1, 2]
    
    max_score <- as.numeric(
      parsed[1, 3]
    )
    
    total_score <- as.numeric(
      parsed[1, 4]
    )
    
    query_coverage <- as.numeric(
      parsed[1, 5]
    )
    
    e_value <- parsed[1, 6]
    
    percent_identity <- as.numeric(
      parsed[1, 7]
    )
    
    hit_length <- as.numeric(
      parsed[1, 8]
    )
    
    accession <- parsed[1, 9]
    
    
    # ------------------------------------------------------
    # Extract first word as potential genus
    # ------------------------------------------------------
    
    potential_genus <- word(
      description,
      1
    )
    
    
    # ------------------------------------------------------
    # Store ALL BLAST information
    # ------------------------------------------------------
    
    all_hits[[hit_counter]] <- tibble(
      
      Query_number =
        i,
      
      Query_ID =
        query_id,
      
      Sequence =
        query_sequence,
      
      BLAST_description =
        description,
      
      Potential_genus =
        potential_genus,
      
      Query_coverage =
        query_coverage,
      
      E_value =
        e_value,
      
      Percent_identity =
        percent_identity,
      
      BLAST_length =
        hit_length,
      
      BLAST_accession =
        accession,
      
      Max_score =
        max_score,
      
      Total_score =
        total_score
      
    )
    
    
    hit_counter <- hit_counter + 1
    
  }
  
}


# ==========================================================
# 9. COMBINE ALL BLAST HITS
# ==========================================================

blast_df <- bind_rows(
  all_hits
)


cat(
  "Total BLAST hits parsed:",
  nrow(blast_df),
  "\n"
)


# ==========================================================
# 10. CHECK THAT BLAST SEQUENCE COLUMN EXISTS
# ==========================================================

if (
  !"Sequence" %in% names(blast_df)
) {
  
  stop(
    "The BLAST parser did not create a Sequence column."
  )
  
}

# ==========================================================
# 11. MATCH BLAST QUERIES TO FASTA SEQUENCES
# ==========================================================

blast_df <- blast_df %>%
  left_join(
    fasta_df %>%
      select(
        FASTA_ID,
        Sequence
      ) %>%
      rename(
        FASTA_Sequence = Sequence
      ),
    by = c("Sequence" = "FASTA_ID")
  )


# ==========================================================
# REPLACE FASTA ID WITH ACTUAL DNA SEQUENCE
# ==========================================================

blast_df <- blast_df %>%
  mutate(
    Sequence = FASTA_Sequence
  ) %>%
  select(
    -FASTA_Sequence
  )


# ==========================================================
# MATCH TO PHYLOSEQ TAXONOMY
# ==========================================================

blast_df <- blast_df %>%
  left_join(
    tax_df %>%
      select(
        Sequence,
        ASV,
        Species_cluster
      ),
    by = "Sequence"
  )


# ==========================================================
# CHECK MATCHING
# ==========================================================

matched_hits <- sum(
  !is.na(blast_df$Species_cluster)
)

cat(
  "BLAST hits matched to Species_cluster:",
  matched_hits,
  "of",
  nrow(blast_df),
  "\n"
)


# ==========================================================
# 12. CHECK SPECIES CLUSTER MATCHING
# ==========================================================

matched_hits <- sum(
  !is.na(
    blast_df$Species_cluster
  )
)

cat(
  "BLAST hits matched to Species_cluster:",
  matched_hits,
  "of",
  nrow(blast_df),
  "\n"
)


# ==========================================================
# 13. SHOW UNMATCHED HITS
# ==========================================================

unmatched_hits <- blast_df %>%
  
  filter(
    is.na(Species_cluster)
  )


cat(
  "Unmatched BLAST hits:",
  nrow(unmatched_hits),
  "\n"
)


# ==========================================================
# 14. APPLY BLAST QUALITY THRESHOLDS
# ==========================================================

blast_df <- blast_df %>%
  
  mutate(
    
    E_value_numeric =
      suppressWarnings(
        as.numeric(E_value)
      ),
    
    Pass_identity =
      Percent_identity >= 97,
    
    Pass_coverage =
      Query_coverage >= 90,
    
    Pass_Evalue =
      E_value_numeric <= 1e-50,
    
    Pass_all =
      Pass_identity &
      Pass_coverage &
      Pass_Evalue
    
  )


# ==========================================================
# 15. DEFINE RECOGNIZABLE GENUS
# ==========================================================
#
# A BLAST description does not always begin with a genus.
#
# Examples that should NOT be treated as a genus:
#
#   uncultured
#   unclassified
#   environmental
#   fungus
#   fungi
#   Ascomycota
#   Basidiomycota
#   Fungi
#
# Everything else is retained as the potential genus.
# ==========================================================

blast_df <- blast_df %>%
  
  mutate(
    
    BLAST_genus = case_when(
      
      is.na(Potential_genus) |
        Potential_genus == "" ~
        "No recognizable genus",
      
      str_to_lower(Potential_genus) %in% c(
        "uncultured",
        "unclassified",
        "environmental",
        "fungus",
        "fungi",
        "ascomycota",
        "basidiomycota",
        "na"
      ) ~
        "No recognizable genus",
      
      TRUE ~
        Potential_genus
      
    )
    
  )


# ==========================================================
# 16. KEEP ALL QUALIFYING BLAST HITS
# ==========================================================

blast_filtered <- blast_df %>%
  
  filter(
    
    Pass_all,
    
    !is.na(Species_cluster),
    
    Species_cluster != ""
    
  )


cat(
  "Qualifying BLAST hits:",
  nrow(blast_filtered),
  "\n"
)


# ==========================================================
# 17. TOTAL BLAST HITS AND QUALITY SUMMARY
# ==========================================================

blast_counts <- blast_df %>%
  filter(
    !is.na(Species_cluster),
    Species_cluster != ""
  ) %>%
  group_by(Species_cluster) %>%
  summarise(
    
    Total_BLAST_hits =
      n(),
    
    Qualifying_BLAST_hits =
      sum(Pass_all, na.rm = TRUE),
    
    .groups = "drop"
  )


# ==========================================================
# 18. AVERAGE BLAST QUALITY FOR EACH SPECIES CLUSTER
# ==========================================================

blast_quality_summary <- blast_filtered %>%
  group_by(Species_cluster) %>%
  summarise(
    
    Mean_percent_identity =
      mean(
        Percent_identity,
        na.rm = TRUE
      ),
    
    Mean_query_coverage =
      mean(
        Query_coverage,
        na.rm = TRUE
      ),
    
    Mean_E_value =
      mean(
        E_value_numeric,
        na.rm = TRUE
      ),
    
    Mean_log10_E_value =
      mean(
        log10(E_value_numeric),
        na.rm = TRUE
      ),
    
    .groups = "drop"
  )


# ==========================================================
# 19. COUNT GENUS CATEGORIES
# ==========================================================

genus_category_counts <- blast_filtered %>%
  count(
    Species_cluster,
    BLAST_genus,
    name = "Genus_hits"
  ) %>%
  group_by(Species_cluster) %>%
  mutate(
    
    Genus_support =
      Genus_hits / sum(Genus_hits),
    
    Genus_support_percent =
      Genus_support * 100
  ) %>%
  ungroup()


# ==========================================================
# 20. CREATE GENUS BREAKDOWN
# ==========================================================

genus_breakdown <- genus_category_counts %>%
  arrange(
    Species_cluster,
    desc(Genus_hits)
  ) %>%
  group_by(Species_cluster) %>%
  summarise(
    
    Genus_breakdown =
      paste0(
        BLAST_genus,
        ": ",
        Genus_hits,
        " (",
        round(Genus_support_percent, 1),
        "%)",
        collapse = "; "
      ),
    
    .groups = "drop"
  )


# ==========================================================
# 21. RECOGNIZABLE VS NON-RECOGNIZABLE GENUS
# ==========================================================

genus_recognition_summary <- genus_category_counts %>%
  group_by(Species_cluster) %>%
  summarise(
    
    Recognizable_genus_hits =
      sum(
        Genus_hits[
          BLAST_genus != "No recognizable genus"
        ]
      ),
    
    No_recognizable_genus_hits =
      sum(
        Genus_hits[
          BLAST_genus == "No recognizable genus"
        ]
      ),
    
    .groups = "drop"
  ) %>%
  mutate(
    
    No_recognizable_genus_percent =
      if_else(
        Recognizable_genus_hits +
          No_recognizable_genus_hits > 0,
        
        No_recognizable_genus_hits /
          (
            Recognizable_genus_hits +
              No_recognizable_genus_hits
          ) * 100,
        
        0
      )
  )


# ==========================================================
# 22. FIND MOST COMMON RECOGNIZABLE GENUS
# ==========================================================

best_genus <- genus_category_counts %>%
  filter(
    BLAST_genus != "No recognizable genus"
  ) %>%
  arrange(
    Species_cluster,
    desc(Genus_hits)
  ) %>%
  group_by(Species_cluster) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  transmute(
    
    Species_cluster,
    
    Most_common_genus =
      BLAST_genus,
    
    Most_common_genus_hits =
      Genus_hits,
    
    Most_common_genus_support_percent =
      Genus_support_percent
  )


# ==========================================================
# 23. CREATE FINAL IDENTIFICATION
# ==========================================================

final_identification <- best_genus %>%
  mutate(
    
    Taxonomic_level =
      if_else(
        Most_common_genus_support_percent >= 50,
        "Genus",
        "Neither"
      ),
    
    Final_identification =
      if_else(
        Most_common_genus_support_percent >= 50,
        Most_common_genus,
        "No genus meeting threshold"
      )
  )


# ==========================================================
# 24. COMBINE EVERYTHING INTO ONE TABLE
# ==========================================================

Species_cluster_BLAST_summary <- blast_counts %>%
  
  left_join(
    blast_quality_summary,
    by = "Species_cluster"
  ) %>%
  
  left_join(
    genus_recognition_summary,
    by = "Species_cluster"
  ) %>%
  
  left_join(
    best_genus,
    by = "Species_cluster"
  ) %>%
  
  left_join(
    genus_breakdown,
    by = "Species_cluster"
  ) %>%
  
  left_join(
    final_identification %>%
      select(
        Species_cluster,
        Taxonomic_level,
        Final_identification
      ),
    by = "Species_cluster"
  ) %>%
  
  mutate(
    
    Recognizable_genus_hits =
      replace_na(
        Recognizable_genus_hits,
        0
      ),
    
    No_recognizable_genus_hits =
      replace_na(
        No_recognizable_genus_hits,
        0
      ),
    
    No_recognizable_genus_percent =
      replace_na(
        No_recognizable_genus_percent,
        0
      ),
    
    Genus_breakdown =
      replace_na(
        Genus_breakdown,
        "No qualifying BLAST hits"
      ),
    
    Taxonomic_level =
      replace_na(
        Taxonomic_level,
        "Neither"
      ),
    
    Final_identification =
      replace_na(
        Final_identification,
        "No genus meeting threshold"
      )
  ) %>%
  
  select(
    Species_cluster,
    
    Total_BLAST_hits,
    Qualifying_BLAST_hits,
    
    Mean_percent_identity,
    Mean_query_coverage,
    Mean_E_value,
    Mean_log10_E_value,
    
    Recognizable_genus_hits,
    No_recognizable_genus_hits,
    No_recognizable_genus_percent,
    
    Most_common_genus,
    Most_common_genus_hits,
    Most_common_genus_support_percent,
    
    Genus_breakdown,
    
    Taxonomic_level,
    Final_identification
  ) %>%
  
  arrange(
    Species_cluster
  )


# ==========================================================
# 25. ROUND NUMERIC VALUES
# ==========================================================

Species_cluster_BLAST_summary <- Species_cluster_BLAST_summary %>%
  mutate(
    
    Mean_percent_identity =
      round(
        Mean_percent_identity,
        2
      ),
    
    Mean_query_coverage =
      round(
        Mean_query_coverage,
        2
      ),
    
    Mean_E_value =
      signif(
        Mean_E_value,
        3
      ),
    
    Mean_log10_E_value =
      round(
        Mean_log10_E_value,
        2
      ),
    
    No_recognizable_genus_percent =
      round(
        No_recognizable_genus_percent,
        1
      ),
    
    Most_common_genus_support_percent =
      round(
        Most_common_genus_support_percent,
        1
      )
  )


# ==========================================================
# 26. VIEW FINAL TABLE
# ==========================================================

print(
  Species_cluster_BLAST_summary,
  n = Inf
)

View(Species_cluster_BLAST_summary)
# ==========================================================
# 27. SAVE FINAL TABLE
# ==========================================================

write.csv(
  Species_cluster_BLAST_summary,
  "Species_cluster_BLAST_summary.csv",
  row.names = FALSE
)
