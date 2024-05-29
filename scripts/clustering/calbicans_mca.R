## ---------------------------
## Script name: calbicans_mca.R
##
## Purpose of script: Perform MCA and plot resulting clusters from normalized
## SNP genotype data in tab-delimited table format
##
## Author: Nancy Scott
##
## Date Created: 2023-04-25
##
## Email: scot0854@umn.edu
## ---------------------------
#!/usr/bin/env Rscript

options(scipen = 999)

# Load packages
library(data.table)
library(ade4)
library(FactoMineR)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(paletteer)
#library(ggsci)
#library(ggtext)
#library(glue)

##############################################################################
# Data preparation
species <- "Calbicans"
genotype_table <-"Calbicans_MEC_bwa_genotypes.txt"
patient_data <- "data/metadata/2022_Calbicans_sorted_patient_pop.csv"
mcaviz <- paste0("images/Calbicans/",Sys.Date(),"_",species,"_SNP_mca_scatterplot.png")
patient_colors <- c(paletteer_d("ggsci::default_igv"),
                    "#838b8b",
                    "black",
                    "#cdcdb4",
                    "#155F83")

# Load vcf table, transpose to samples as rows,
genotypes <- read.table(genotype_table,
                   header = TRUE,
                   sep="\t")
# Use info rows (chr, pos, snp, etc) for column names
# Convert strings to factors and proceed
gt <- data.table::transpose(genotypes)
genotypes_id <- paste(gt[1,],gt[2,])
colnames(gt) <- genotypes_id
gt <- gt[-c(1:2),]
gt <- as.data.frame(unclass(gt), stringsAsFactors = TRUE)
gt <- gt[, sapply(gt, nlevels) !=1]
num_variants <- ncol(gt)

# Load patient info, use sample IDs as row names in snp dataframe
pop.data <- read.table(patient_data,
                       sep = ",",
                       header = TRUE)
row.names(gt) <- pop.data$sample
# Add patient factor to genotype table
gt$pop <- as.factor(pop.data$population)
# Optional, add additional factors to snp dataframe for coloring individuals
gt$group <- as.factor(pop.data$group)
###############################################################################
# Multiple correspondence analysis of genotypes
# following example at http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/114-mca-multiple-correspondence-analysis-in-r-essentials/
# assuming only SNPs in dataset, with fixed SNPs and missing genotypes are removed
# specify active individuals and active variables (variant sites)
gt.active <- gt[,1:num_variants]
var.mca <- MCA(gt.active, graph = FALSE)
save(var.mca, file=paste0(Sys.Date(),"_Calbicans_MEC_SNP_mca.rda"))

# var.mca object is a list of objects including a matrix of eigenvalues with percentage of variance
# and var/ind - matrices of all results for active variables and individuals (including coordinates)

eig_vals <- data.frame(var.mca$eig)
ind_coords <- data.frame(var.mca$ind$coord)
ind_coords <- ind_coords %>%
  mutate(sample = rownames(ind_coords)) %>%
  left_join(pop.data, by= "sample")

################################################################################
# Override too many overlapping labels (for session)
options(ggrepel.max.overlaps = Inf)
mca_plot <- ggplot(ind_coords, aes(x=Dim.1, y=Dim.2, color=as.factor(population), label=sample)) +
  geom_point() +
  scale_color_manual(values=patient_colors, guide = "none") +
  geom_text_repel(show.legend = FALSE, size = 4) +
  theme_bw() +
  theme(axis.title.x = element_text(family = "Helvetica", color = "black", size=12),
        axis.title.y = element_text(family = "Helvetica", color = "black", size = 12)) +
  xlab(paste0("Dimension 1, ",sprintf("%.2f",eig_vals$percentage.of.variance[1] ), "%")) +
  ylab(paste0("Dimension 2, ",sprintf("%.2f",eig_vals$percentage.of.variance[2]),"%"))

ggsave(mcaviz, mca_plot,
       width = 8, height = 6, units = "in",
       device = png, bg = "white", dpi = 300)
