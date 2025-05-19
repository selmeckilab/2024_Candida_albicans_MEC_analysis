## Purpose: Overlap GC-corrected read depth with gene locations to find candidate
## copy number changes
## Author: Nancy Scott
## Email: scot0854@umn.edu

# Load packages----
library(readxl)
library(tidyverse)
library(GenomicRanges)
library(regioneR)
library(rtracklayer)

# Sample inputs----
spreadsheet_list <- "~/umn/data/metadata/2024_Calbicans_MEC_500bp_depth.txt"
copy_number_files <- scan(spreadsheet_list, what=character())

# Genome inputs----
candida_gff <- "~/umn/data/reference_genome_coords/Calbicans/C_albicans_SC5314_A21_current_features.gff"
region_file <- "~/umn/data/reference_genome_coords/SC5314_A21_region_file.bed"

# Window size, CNV size, gaps and depth----
window_size <- 499 # rolling mean window size minus 1
cnv_size <- 15000 # how big of a run to search for
gap_width <- 2002 # merge adjacent runs with gaps less than this
low_cov <- 0.6 # relative depth cutoff for deletions
high_cov <- 1.3 # relative depth cutoff for amplifications
overlap_size <- 350

## Get features----
features <- import.gff(candida_gff)
exclude_features <- c("repeat_region", "rRNA", "long_terminal_repeat")
features <- features[features$type %in% exclude_features]

## Get ranges of interest
regions <- import.bed(region_file)
regions <- setdiff(regions, features, ignore.strand=TRUE)

## Function to get candidate genes----
# Returns GRange of gene candidates, either above or below relative_depth threshold, with minimum specified overlap
find_cnv_candidates <- function(depth_range, op, threshold){
  cov1 <-depth_range[get(op)(depth_range$relative_depth, threshold)]
  cov2 <- reduce(cov1, min.gapwidth = 2)
  cov2 <- cov2[width(cov2) > cnv_size]
  cov1 <- subsetByOverlaps(cov1, cov2)
  return(cov1)
}


## Process all files from list----
for(i in 1:length(copy_number_files)){
  sample_id <- str_extract(copy_number_files[i], "AMS[:digit:]+|MEC[:digit:]+")

  # Make Granges of depth data
  genome_depth <- read_xlsx(copy_number_files[i]) 
  colnames(genome_depth)[2] <- "begin"

  genome_depth <- genome_depth %>% 
    mutate(end = begin + window_size) %>% 
    select(chr,begin,end,relative_depth) %>% 
    mutate(relative_depth = case_when(is.na(relative_depth) ~ 0,
                                      .default = relative_depth))
  
  genome_mean_depth <- genome_depth %>% 
    group_by(chr) %>% 
    summarise(avg_depth = mean(relative_depth))
  
  genome_mean_depth_keep <- genome_mean_depth %>% 
    filter(avg_depth < 1.28)
  
  genome_cov <- toGRanges(as.data.frame(genome_depth))

  # Get low coverage regions
  deletions <- genome_cov[genome_cov$relative_depth < low_cov,]
  deletions <- subsetByOverlaps(deletions, regions, ignore.strand = TRUE, minoverlap = overlap_size)
  del_runs <- reduce(deletions, min.gapwidth = gap_width)
  del_runs <- del_runs[width(del_runs)>cnv_size]
  deletions <- subsetByOverlaps(deletions, del_runs, ignore.strand = TRUE, minoverlap = overlap_size)
  
  # Get high coverage regions
  amplifications <- genome_cov[genome_cov$relative_depth > high_cov,]
  amplifications <- subsetByOverlaps(amplifications, regions, ignore.strand = TRUE, minoverlap = overlap_size)
  amp_runs <- reduce(amplifications, min.gapwidth = gap_width)
  amp_runs <- amp_runs[width(amp_runs)>cnv_size]
  amplifications <- subsetByOverlaps(amplifications, amp_runs, ignore.strand = TRUE, minoverlap = overlap_size)
  
  # data frames and write to excel
  del_details <- as.data.frame(deletions)
  del_summary <- as.data.frame(del_runs)
  amp_details <- as.data.frame(amplifications)
  amp_summary <- as.data.frame(amp_runs)
  
  write_xlsx(list(del_summary = del_summary,
                  del_details = del_details, 
                  amp_summary = amp_summary,
                  amp_details = amp_details), 
                  
             paste0("cnv_files/",sample_id, "_cnv_candidates.xlsx"))
}
