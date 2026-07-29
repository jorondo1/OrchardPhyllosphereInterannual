Currently, the only non-code files in the repository are the final RDS objects, essentially:

`data/ps_objects_full.rds` which is a list containing only two ps objects now: $Bacteria and $Fungi, each of which has phylo trees and species clusters. These are the filtered datasets, meaning without ASVs having <10 reads across all samples. There is another dataset with phyloseq objects including all ASVs (pre-filtering), it's not on the repo yet because I don't think you'll need it; let me know if you do.

`data/diversity_data.rds` are alpha and beta diversity data for each datasests, computed at two levels: ASV and (phylogenetic) species clusters. It is a nested list, containing:
- at the first level, four lists: `$Bacteria`, `$Fungi`, `$Bacteria_sp_clust` and `$Fungi_sp_clust`
- at the second level, `$alpha` or `$beta` sublists:
  - `$alpha` is the phyloseq object's sample_data as a tibble, with added column for the following indices: Richness, Shannon, Simpson, Hill order 1, Hill order 2, Faith's PD and Tail statistic
  - `$beta` is a sublist of distance/dissimilarity matrices, namely Bray-Curtis, Robust Aitchison, Weighted Unifrac and Unweighted Unifrac.
