## Purpose: Combine het SNP counts and generate heatmap across chrs for a group of isolates
## Author: Nancy Scott
## Email: scot0854@umn.edu

## Load packages----
library(readxl)
library(tidyverse)
library(paletteer)
library(writexl)

## Variables----
spreadsheet_list <- "data/metadata/2024_Calbicans_snp_depth_paths.txt"
ordered_patient_data <- "data/metadata/2025_Calbicans_midpoint_gatk_302.csv"
in_patient_data <- "data/metadata/2024_Calbicans_sorted_patient_info.xlsx"
save_dir <- "images/Calbicans/"

binned_color_scale <- c("white", paletteer_c("grDevices::Turku", 30))

# Manually determined cen positions
cens <- data.frame(x = c(314, 1024, 1251, 1644, 1860, 2202, 2298, 2751),
                  y = rep(0.6,8), z = rep(102.5, 8))

## Read in SNP counts for all samples----
snp_files <- scan(spreadsheet_list, what=character())

genome_snp <- read_xlsx(snp_files[1]) %>%
  select(index, pos, snp_count)
names(genome_snp)[names(genome_snp)=="snp_count"] <-str_extract(snp_files[1], "AMS[:digit:]+|MEC[:digit:]+")

for(i in 2:length(snp_files)){
  new_snp <- read_xlsx(snp_files[i]) %>%
    select(index, pos, snp_count)
  names(new_snp)[names(new_snp)=="snp_count"] <-str_extract(snp_files[i], "AMS[:digit:]+|MEC[:digit:]+")

  genome_snp <- genome_snp %>%
    left_join(new_snp, by = join_by(index,pos))
}

## Plotting details----

# Regular intervals for x-axis
genome_snp$x_num <- seq.int(nrow(genome_snp))

# Finish tidying data
snp_again <- genome_snp %>%
  pivot_longer(names_to = "sample", values_to = "snp_count", cols=-c(index,pos, x_num))

snp_again <- snp_again %>% 
  mutate(sample = replace(sample, sample=="MEC103", "MEC103-2"))

# Chr outlines and tick marks
chrs <- genome_snp %>%
  group_by(index) %>%
  summarise(border_start=min(x_num),
            border_stop=max(x_num),
            tick=min(x_num) + (max(x_num)-min(x_num))/2
            )

# Metadata, sort samples by tree order
pop.data <- read_excel(in_patient_data)
pop.data[pop.data=="NA"] <- NA

pt_order <- read.csv(ordered_patient_data, header = TRUE)
pt_order <- pt_order %>%
  right_join(pop.data %>%
               filter(study=="MEC", sample !="AMS6466") %>%
               select(sample, mec_pt_code, mec_isolate_code)) %>% 
  mutate(mec_isolate_code = replace(mec_isolate_code, mec_isolate_code=="SC5314-1", "SC5314"))

## Plot with manually annotated clade breaks----
p <- snp_again %>%
  filter(sample != "AMS6466") %>% 
  mutate(clustered_samples = fct_relevel(sample, rev(pt_order$sample))) %>%
  ggplot(aes(x=x_num, y=clustered_samples)) +
  geom_tile(aes(fill=snp_count, colour = snp_count), linewidth=0.27) +
  scale_color_stepsn(na.value = "white", 
                     limits = c(0,300),
                     n.breaks = 31,
                     colors = binned_color_scale,
                     guide = "none") +
  scale_fill_stepsn(na.value = "white", 
                    limits = c(0,300),
                    n.breaks = 31,
                    labels = c("0", rep("",9), "100", rep("", 9), "200", rep("", 9), "300"),
                    colors = binned_color_scale) +
  scale_y_discrete(labels = rev(pt_order$mec_isolate_code)) +
  theme_minimal() +
  geom_rect(data=chrs, aes(group=index, xmin=border_start, xmax=border_stop, ymin=0.5, ymax=Inf),
                            fill=NA, inherit.aes=FALSE,colour = "grey26", linejoin = "round") +
  geom_point(data = cens, aes(x =x, y = y), size = 2, shape = 23, fill = "black") +
  geom_point(data = cens, aes(x = x, y = z), size = 2, shape = 23, fill = "black") +
  scale_x_continuous(expand = c(0,0),
                      breaks = chrs$tick,
                     position = "top",
                      labels=c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5", "Chr6", "Chr7", "ChrR")) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=5.4, color = "black"),
        axis.text.x = element_text(size = 8, color = "black"),
        plot.margin = margin(0,10.5,5.5,5.5,"pt"),
        legend.text = element_text(size = 5, color = "black"),
        legend.position = "right", 
        legend.box.margin = margin(0,0,0,0),
        legend.margin = margin(0,0,0,0),
        legend.justification = "bottom",
        legend.key.size = unit(0.5, 'lines'),
        legend.title = element_blank()) +
  annotate("segment", x=0, xend = 2860, y= 20.5, yend = 20.5, linewidth = 0.2, colour = "grey39") +
  annotate("segment", x=0, xend = 2860, y= 21.5, yend = 21.5, linewidth = 0.2, colour = "grey39") +
  annotate("segment", x=0, xend = 2860, y= 22.5, yend = 22.5, linewidth = 0.2, colour = "grey39") +
  annotate("segment", x=0, xend = 2860, y= 23.5, yend = 23.5, linewidth = 0.2, colour = "grey39") +
  annotate("segment", x=0, xend = 2860, y= 38.5, yend = 38.5, linewidth = 0.2, colour = "grey39") +
  annotate("segment", x=0, xend = 2860, y= 42.5, yend = 42.5, linewidth = 0.2, colour = "grey39") +
  annotate("segment", x=0, xend = 2860, y= 47.5, yend = 47.5, linewidth = 0.2, colour = "grey39") +
  annotate("segment", x=0, xend = 2860, y= 55.5, yend = 55.5, linewidth = 0.2, colour = "grey39") +
  annotate("segment", x=0, xend = 2860, y= 60.5, yend = 60.5, linewidth = 0.2, colour = "grey39") +
  annotate("segment", x=0, xend = 2860, y= 64.5, yend = 64.5, linewidth = 0.2, colour = "grey39") 

# Troubleshooting: what are the axis positions per sample?
#details <- ggplot_build(p)

## Save plot----
ggsave(paste0(save_dir,Sys.Date(),"_MEC_Calbicans_LOH_heatmap.pdf"), 
       p, 
       #device=png, 
       #dpi=300, 
       #bg="white",
       width = 6.8, 
       height = 6, 
       units="in")

