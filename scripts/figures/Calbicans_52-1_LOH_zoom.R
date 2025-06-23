## Purpose: Replot LOHfrom saved summary data
## Author: Nancy Scott
## Email: scot0854@umn.edu

# Load packages----
library(readxl)
library(RcppRoll)
library(tidyverse)

# Variables----
read_depth_file <- "data/genome_plots/Calbicans/2023-11-29_MEC318.xlsx"
sample_id <- "Isolate 52-1"
feature_file <- "scripts/figures/plotting_data/Calbicans_SC5314_A21_plotting_features.txt"
label_file <- "scripts/figures/plotting_data/Calbicans_SC5314_A21_chr_labels.txt"

window <- 5000 # size of window used for rolling mean and snp density
ploidy <- 2
chrom_of_interest <- 5
y_max_limit <- 5
region_start <- 0
region_stop <- 760000
x_breaks <- c(0, 250000, 500000, 750000)
x_labels <- c("0", "0.25", "0.5", "0.75")

snp_low <- "white"
snp_high <- "black"
feature_colors <- c("#0072B2", "#56B4E9")
feature_shapes <- c(23,21)
feature_sizes <- c(3, 2)

save_dir <- "" # path with trailing slash, or  "" to save locally
ref <- "SC5314" # short label for file name or "" to leave out

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
region_loh <- genome_depth %>%
  filter(index == chrom_of_interest,
         pos>= region_start,
         pos<= region_stop) %>%
  ggplot() +
  geom_segment(aes(x = pos, 
                   xend = pos, 
                   y = 0, 
                   yend = y_max_limit,
                   color = snp_count)) +
  theme_minimal() +
  xlab(paste0("Chr",
              chrom_of_interest,
              " position (Mb)")) +
  ylab("Heterozygous SNP Density") +
  theme(axis.title = element_text(family = "Helvetica")) +
  scale_color_gradient(low=snp_low,
                       high=snp_high, 
                       na.value = "white", 
                       guide = "none") +
  scale_y_continuous(limits = c(0, y_max_limit), breaks = NULL) +
  scale_x_continuous(breaks = x_breaks,labels = x_labels) +
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

ggsave(paste0(sample_id, "_zoom_LOH.pdf"),
       region_loh,
       device = pdf,
       width = 3.5,
       height = 2.2,
       units = "in"
       )
