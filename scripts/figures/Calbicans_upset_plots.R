## Purpose: Plot intersection of SNPs in series 
## Author: Nancy Scott
## Email: scot0854@umn.edu

## Load packages----
library(UpSetR)
library(tidyverse)

in_file <- "data/variants/Calbicans/manual_review/manual_N_unfixed.txt"
all_snps <- read.delim(in_file,
                 header = T)
#samples <- c("52-1", "52-2", "52-5", "52-3", "52-4")
samples <- c("14-1","14-2","14-3","14-4","14-5")
colnames(all_snps)[-c(1,2)] <- samples

series <- sub('\\.txt$', '', basename(in_file))

genotypes <- unique(as.vector(unlist(all_snps[,-c(1:2)])))
genotypes

all_snps[all_snps=="0/0"] <- 0
all_snps[all_snps=="0/1"] <- 1
all_snps[all_snps=="1/1"] <- 1

all_snps[,-c(1:2)] <- sapply(all_snps[, -c(1:2)], as.numeric)

p <- upset(all_snps[,-c(1:2)],
       order.by = "degree", 
       mainbar.y.label = "Inter-isolate SNPs",
       main.bar.color = "#499894",
       sets.x.label = "SNPs per isolate")

pdf(paste0(series,"_upset.pdf"), width = 4, height = 3.33, bg = "white")
print(p)
dev.off()


gt <- as.matrix(all_snps)
headers <- gt[,"POS"]
gt_sample_rows <- t(gt)
gt_sample_rows <- gt_sample_rows[3:7,]
colnames(gt_sample_rows) <- headers
dist(gt_sample_rows)
hclust(dist(gt_sample_rows))
case_cluster <- hclust(dist(gt_sample_rows))

pdf(paste0(series,"_dendrogram.pdf"), width = 4, height = 4, bg = "white")
print(plot(case_cluster))
dev.off()
