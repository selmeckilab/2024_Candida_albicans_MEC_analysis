## Purpose: Input, manipulate and plot phylogenetic tree data
## Author: Nancy Scott
## Email: scot0854@umn.edu

options(scipen = 999)

## Load packages----
library(tidyverse)
library(ggtree)
library(paletteer)
library(phangorn)
library(readxl)
library(tidytree)
library(treeio)
library(castor)

## Variables----
species <- "Calbicans"
sample_number <- 302
raxml_file <- "data/phylogeny/Calbicans/RAxML_bipartitions.Calbicans_302"

metadata_file <- "data/metadata/2024_Calbicans_sorted_patient_info.xlsx"

# Clade colors
color_file <- read.table("scripts/figures/plotting_data/clade_colors.txt")
clade_colors <- color_file[,2]
names(clade_colors) <- color_file[,1]

# Phylogeny file
Candida_raxml <- read.tree(raxml_file) 
Candida_raxml_midpoint <- midpoint(Candida_raxml, node.labels = "support") # keeps branch support

#ggtree(Candida_raxml_midpoint) + geom_label2(aes(subset = !isTip, label = label), size =2)

## Metadata----
Candida_metadata <- read_excel(metadata_file)
Candida_metadata[Candida_metadata=="NA"] <- NA
Candida_metadata <- rename(Candida_metadata, label = sample)
Candida_metadata$mec_pt_code <- as.character(Candida_metadata$mec_pt_code)
Candida_metadata$cluster <- as.character(Candida_metadata$pt_cluster)
Candida_metadata$ST <- as.character(Candida_metadata$ST)

## Clade determination----
# Selecting a branch length (or relative branch length)
longest_dist <- castor::find_farthest_tip_pair(Candida_raxml_midpoint)
Candida_clades <- as_tibble(Candida_raxml_midpoint) %>% 
  mutate(rel_branch = branch.length/longest_dist$distance) %>% # normalize branch length to longest pairwise dist.
  filter(rel_branch >= 0.057) %>% # set a minimum distance
  filter(!node %in% parent, !grepl("^[[:alpha:]]", label)) %>% # keep only inner nodes
  rename(clade_node = node) %>% 
  mutate(Clade = c("1", "8B","8A", "4B", "4A", "D", "16", "13", "18", "3", "10",
                   "12", "B", "2", "E", "6", "8C", "11", "9", "7"))
  #Calbicans 310 = mutate(Clade = c("1", "7", "9", "11", "12", "3", "18", "10", "D", "16", "13", "B",
  #                 "2", "8B", "8A", "8*", "4", "4*", "6", "E"))

## Tree object combining phylo and metadata----
Candida_tree <- Candida_raxml_midpoint %>% full_join(Candida_metadata, by = "label")

#Candida_distances <- castor::get_all_pairwise_distances(Candida_raxml_midpoint)

# Plot the midpoint-rooted tree with patient codes and clade highlighting , x=0, y=0 , align = TRUE, linetype = "dotted", linesize = 0.2
# option for labeling a subset with geom_tiplab: aes(subset=(label %in% studies)) aes(subset=(label %in% long_reads))
midpoint_plot <- ggtree(Candida_tree,
                        aes(size = (as.numeric(label) < 95 | is.na(as.numeric(label)))))  +
  scale_size_manual(values = c(1, 0.3), 
                    guide = "none") +
  geom_tippoint(aes(subset = !is.na(singletons)),
                color = "firebrick3", 
                shape = 17, 
                size=1.5) +
  geom_tiplab(aes(subset = !is.na(singletons),
                  label = mec_isolate_code),
              size=2, 
              align = TRUE, 
              offset = -0.008,
              linetype = "dotted", 
              linesize = 0.2) + 
  geom_tiplab(aes(subset = !is.na(ref_genome),
                  label = mec_isolate_code),
              size = 2,
              align = TRUE,
              offset = -0.008,
              linetype = "dotted",
              linesize = 0.2) + 
  geom_rootedge() +
  scale_fill_manual(values = clade_colors, 
                    guide = "none") + 
  geom_treescale(x=0, 
                 y=-0.3) +
  geom_cladelab(data = Candida_clades, 
                mapping = aes(node = clade_node, label = Clade),
                align = TRUE,
                offset = 0.001,
               fontsize = 3
                ) +
  geom_hilight(data = Candida_clades, aes(node = clade_node, fill = Clade)) 

## Save tree----
ggsave(paste0(Sys.Date(),"_",species,"_midpoint_gatk_",sample_number,".pdf"),
       midpoint_plot, bg="white", width = 7, height = 9, units = "in")

## Write out sample list in order of tips----
write.csv(get_taxa_name(midpoint_plot), 
          paste0(Sys.Date(),"_",species,"_midpoint_gatk_",sample_number,".csv"),
          quote = FALSE,
          row.names = FALSE)
