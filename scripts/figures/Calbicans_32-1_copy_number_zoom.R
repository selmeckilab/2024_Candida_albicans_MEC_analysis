## Purpose: Replot genome-view of LOH and CNV from saved summary data
## Author: Nancy Scott
## Email: scot0854@umn.edu

# Load packages----
library(readxl)
library(RcppRoll)
library(tidyverse)

# Variables----
read_depth_file <- "data/genome_plots/Calbicans/500bp/2024-03-27_MEC218.xlsx"
sample_id <- "Isolate_32-1"
feature_file <- "scripts/figures/plotting_data/Calbicans_SC5314_A21_plotting_features.txt"
label_file <- "scripts/figures/plotting_data/Calbicans_SC5314_A21_chr_labels.txt"

window <- 500 # size of window used for rolling mean and snp density
ploidy <- 2
chrom_of_interest <- 8
y_max_limit <- 5
region_start <- 1000000 
region_stop <- 2290000
x_breaks <- c(1000000, 1300000, 1900000, 2280000)
x_labels <- c("1.0","1.3", "1.9", "2.28")
bp1_start <- 1158902
bp1_stop <- 1158997


feature_colors <- c("#0072B2", "#56B4E9")
repeat_color <- "#D55E00"
feature_shapes <- c(23,21)
feature_sizes <- c(3, 2)

save_dir <- "" # path with trailing slash, or  "" to save locally

# X-axis labels overwrite input scaffold names in final plot
chr_ids <- scan(label_file, what = character())

# Read saved dataframe and feature file----
genome_depth <- read_xlsx(read_depth_file, sheet=1)

features <- read_tsv(feature_file, show_col_types = FALSE)

features <- features %>%
  group_by(chr, index=consecutive_id(chr)) %>%
  mutate(ymin=0, ymax=y_max_limit)

chr_feature <- features %>%
  filter(index == chrom_of_interest, Feature %in% c("Centromere", "MTL"))

# Plot region of interest----
region_depth <- genome_depth %>%
  filter(index == chrom_of_interest,
         pos>= region_start,
         pos<= region_stop) %>%
  ggplot(aes(x = pos, y=copy_number)) +
  geom_point(alpha = 0.6, size=1) +
  theme_minimal() +
  xlab(paste0("Chr",
              chrom_of_interest,
              " position (Mb)")) +
  ylab("Relative copy number") +
  theme(axis.title = element_text(family = "Helvetica")) +
  scale_y_continuous(limits = c(0, y_max_limit)) +
  scale_x_continuous(breaks = x_breaks,labels = x_labels) +
  annotate("rect", 
           xmin = bp1_start,
           xmax = bp1_stop, 
           ymin = 0, 
           ymax = 5,
           alpha = 0.6,
           color = repeat_color,
           fill = repeat_color) +
  geom_point(data = chr_feature, 
             aes(x = start, 
                 y = 0, 
                 shape = Feature,
                 fill = Feature,
                 size = Feature)) +
  scale_shape_manual(values = feature_shapes) +
  scale_fill_manual(values = feature_colors) +
  scale_size_manual(values = feature_sizes) +
  theme(legend.position = "none")

ggsave(paste0(sample_id, "_zoom_copy_number.pdf"),
       region_depth,
       device = pdf,
       width = 3.5,
       height = 3.6,
       units = "in"
       )
