pacman::p_load(tidyverse, phyloseq, mgx.tools, kableExtra)

source('src/0.0_Config.R')

ps.ls <- readRDS('data/ps_objects_full.rds')[c('Bacteria', 'Fungi')]

# Core microbiome

relab.ls <- map(ps.ls, function(ps) { 
  tax_glom2(ps, 'Species_cluster') %>% 
    psflashmelt() %>% 
    group_by(Sample) %>% 
    mutate(relAb = Abundance/sum(Abundance))
})

prev <- imap(relab.ls, function(relab.tibble, Kingdom) {
  
  relab.tibble%>% 
    group_by(Kingdom, Phylum, Class, Order, Family, Genus, Species_cluster) %>% 
    summarise(
      prev = sum(relAb>0.001)/n(),
      mean_relab = mean(relAb),
      sd_relab = sd(relAb),
      .groups = 'drop'
    ) %>% 
    filter(prev>0.7) %>% 
    mutate(Species_cluster = if_else(
      Species_cluster %in% names(fungi_label_overrides),
      fungi_label_overrides[Species_cluster],
      Species_cluster
    ))
}) %>% 
  list_rbind() %>% 
  arrange(desc(prev))

print(prev, n = 100)

kableExtra::kable(prev, "html", align = "l") %>%
  kableExtra::kable_styling(full_width = FALSE, bootstrap_options = c("striped", "hover")) %>%
  kableExtra::save_kable(file = 'out/manuscript/core_species.html')
