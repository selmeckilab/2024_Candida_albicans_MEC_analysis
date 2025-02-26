## ---------------------------
## Purpose: Summary clinical and phenotypic data 
## Author: Nancy Scott
## Email: scot0854@umn.edu
## ---------------------------
options(scipen = 999)

## Load packages
library(tidyverse)
library(readxl)
library(writexl)

# Local antifungal data
af_spreadsheet <- "~/umn/data/metadata/MEC_all_isolates_antifungal_history.xlsx"

# manual review of CHEF and CNVs
interesting_isolates <- c("MEC352", "MEC079", "MEC218", "MEC219", "MEC279", "MEC293", 
                          "MEC318", "MEC319", "MEC322", "MEC135", "MEC157", "MEC185",
                          "MEC131", "MEC324", "MEC246", "MEC172", "MEC257", "MEC195",
                          "MEC085", "MEC080", "MEC268", "MEC247", "MEC198", "MEC297",
                          "MEC174")

# Redcap report IDs
samples <- '58043'
mic_results <- '58044'
growth_curves <- '58045'

# Function to import report from redcap
import_report <- function(report_number) {
  url <- "https://redcap.ahc.umn.edu/redcap/api/"
  formData <- list("token"=Sys.getenv("redcap_api_key"),
                   content='report',
                   format='csv',
                   report_id=report_number,
                   csvDelimiter='',
                   rawOrLabel='label',
                   rawOrLabelHeaders='raw',
                   exportCheckboxLabel='true',
                   returnFormat='csv'
  )
  response <- httr::POST(url, body = formData, encode = "form")
  result <- httr::content(response, show_col_types = FALSE)
}
################################################################################
# Read in reports.

# Sample ID, species, series and cluster IDs
sample_info <- import_report(samples) %>%
    select(-c(starts_with('redcap_repeat'))) %>%
    filter(isolate_type == "clinical", genus_species == "Candida albicans")

sample_info$genus_species <- str_replace(sample_info$genus_species, "Candida", "C.")

sample_naming <- sample_info %>% 
  mutate(join_id = case_when(!is.na(secondary_id) ~ secondary_id, .default = primary_id)) %>%
  select(primary_id, join_id)

# MIC and SMG results
mic_info <- import_report(mic_results) %>%
    filter(redcap_repeat_instrument != "NA") %>%
    select(primary_id, 
           redcap_repeat_instance, 
           drug, 
           mic_media,
           mic_date, 
           mic50,
           mean_mic50_relative_growth,
           mic90,
           mean_mic90_relative_growth,
           eucast_breakpoint, 
           mean_smg, 
           sem_smg,
           qc_ok) %>% 
    left_join((sample_info %>% 
                   select(primary_id, genus_species, series_id, patient_code)), 
              by=join_by(primary_id))

# Subset to RPMI data with valid control results
# Slice head removes repeated assays, assuming most recent is best (for now, filtering isn't working)
mic_info <- mic_info %>% 
    filter(mic_media %in% c("RPMI liquid", "RPMI agar, Etest"), 
           !(primary_id %in% c("AMS5122","AMS5123","AMS2401")),
           !(mic_date %in% c(as.Date("2024-03-07"))),
           qc_ok=="Yes") %>% 
  group_by(primary_id, drug) %>% 
  arrange(desc(mic_date)) %>% 
  slice_head()

# For ordering species and drug levels
species_count <- sample_info %>%
    group_by(genus_species) %>%
    summarize(species_count=n()) %>%
    arrange(desc(species_count))

mic_info$mic50 <- factor(mic_info$mic50, levels=c("0.016", "0.023", "0.032", "0.047", "0.064", "0.125",
                                                  "0.256", "0.38", "0.5", "1", ">1", "2", 
                                                  "4", "8", "16", "32", ">32", "64", 
                                                  "96" ,"128","160", "256", ">256"))

mic_info$mic90 <- factor(mic_info$mic90, levels=c("0.016", "0.023", "0.032", "0.047", "0.064", "0.125",
                                                  "0.256", "0.38", "0.5", "1", ">1", "2", 
                                                  "4", "8", "16", "32", ">32", "64", 
                                                  "96" ,"128","160", "256", ">256"))

# Adjust Etest MICs to bin with broth MIC levels
mic_info <- mic_info %>% 
  mutate(mic90=replace(mic90, mic90=="0.023", "0.016"))

mic_info <- mic_info %>% 
  mutate(mic50 = replace(mic50, mic50=="0.047", "0.032"))

# Growth curve
gc <- import_report(growth_curves) %>%
  filter(redcap_repeat_instrument != "NA", 
         !primary_id %in% c("AMS5122", "AMS5123")) %>%
  select(primary_id, redcap_repeat_instance, 
         gc_date, drug_used, 
         gc_temp, gc_time, 
         k, r, t_gen, auc_l, auc_e) 

gc <- gc %>% 
  group_by(primary_id) %>% 
  filter(drug_used=="None", 
         gc_date==max(gc_date), 
         gc_time < 24.1,
         primary_id %in% mic_info$primary_id)  


# No-drug well OD vals
control_od <- mic_info %>% 
  group_by(primary_id) %>% 
  summarise(mean_stationary_k = mean(mean_no_drug_stationary_k, na.rm = TRUE)) %>%
  left_join(sample_info)

################################################################################
# For subsetting SMG to isolates with adequate 24-hour growth
carrying_cap <- mic_info %>% 
  group_by(primary_id) %>% 
  filter(mean_no_drug_stationary_k < 0.2)

flc_carrying_cap <- carrying_cap %>% 
  filter(drug=="fluconazole")

# Summarise all FLC SMG by species
smg <- mic_info %>% 
  filter(!primary_id %in% flc_carrying_cap$primary_id, drug=="fluconazole") %>% 
  group_by(drug) %>% 
  filter(!is.na(mean_smg)) %>% 
  summarise(number_smg = n(),
            min_smg = round(min(mean_smg), digits = 2),
            max_smg = round(max(mean_smg), digits = 2),
            overall_mean_smg = round(mean(mean_smg), digits = 2), 
            median_smg = round(median(mean_smg), digits = 2),
            IQR_smg = round(IQR(mean_smg), digits = 2)
            )

# Summarise changes in FLC SMG within each series
serial_flc_smgs <- mic_info %>% 
  filter(!primary_id %in% flc_carrying_cap$primary_id, drug=="fluconazole") %>% 
  group_by(drug, series_id) %>% 
  filter(!is.na(series_id), !is.na(mean_smg)) %>% 
  summarise(number_smg = n(),
            min_smg = round(min(mean_smg), digits =2),
            max_smg = round(max(mean_smg), digits = 2),
            diff_smg = max(mean_smg) - min(mean_smg))

################################################################################
# Summary tables for resistant isolates

# All isolates done for each drug? No duplicates?
isolate_counts <- mic_info %>%
    group_by(drug) %>%
    count(genus_species)

# Differences within series
serial_mic50 <- mic_info %>% 
    filter(!is.na(series_id), drug %in% c("fluconazole", "micafungin")) %>% 
    group_by(series_id, drug, mic_media, mic50) %>% 
    count()

serial_mic90 <- mic_info %>% 
  filter(!is.na(series_id), drug %in% c("amphotericin B")) %>% 
  group_by(series_id, drug, mic_media, mic90) %>% 
  count()

changed_amb_series <- serial_mic90 %>% 
  filter(drug %in% c("amphotericin B")) %>% 
  group_by(series_id, mic_media) %>% 
  filter(n() >1)

changed_mic50_series <- serial_mic50 %>% 
  filter(drug %in% c("fluconazole", "micafungin")) %>% 
  group_by(series_id, mic_media) %>% 
  filter(n() >1)  

changes_by_drug <- changed_mic50_series %>% 
  rbind(changed_amb_series) %>% 
  group_by(series_id,drug) %>% 
  count() %>% 
  filter (n>1) %>% 
  arrange(series_id)

################################################################################
# Drug exposure data
antifungals <- read_xlsx(af_spreadsheet)

no_results <- sample_info %>% filter(!(primary_id %in% antifungals$primary_id)) %>% 
  select(primary_id, genus_species, relative_days, patient_code, series_id, secondary_id)

calbicans_remote_exposure <- antifungals %>% 
  filter(primary_id %in% sample_info$primary_id) %>% 
  filter(relative_start_days <= 0, relative_collection_day == 0)

# remote exposure in isolates with CNA, ploidy change, large scale karyotype changes
interesting_remote <- calbicans_remote_exposure %>% 
  filter(primary_id %in% interesting_isolates)

topical_remotes <- interesting_remote[str_detect(interesting_remote$DRUG_NAME_ORIG, "NYS|CREA|POWD|OINT"),]
systemic_remotes <- interesting_remote[!str_detect(interesting_remote$DRUG_NAME_ORIG, "NYS|CREA|POWD|OINT"),]

calbicans_all_exposure <- antifungals %>% 
  filter(primary_id %in% sample_info$primary_id)

calbicans_topical_exposure <- calbicans_all_exposure[str_detect(calbicans_all_exposure$DRUG_NAME_ORIG, "CLOTRIMAZOLE|NYS|CREA|POWD|OINT"),]
calbicans_systemic_exposure <- calbicans_all_exposure[!str_detect(calbicans_all_exposure$DRUG_NAME_ORIG, "CLOTRIMAZOLE|NYS|CREA|POWD|OINT"),]
