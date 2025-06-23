## Purpose: Summary clinical and phenotypic data 
## Author: Nancy Scott
## Email: scot0854@umn.edu

## Get cleaned and processed data from REDCap----
# Use script from candidemia phenotyping manuscript for function, sample IDs and MICs
source("~/umn/Candida_clinical_isolate_data/redcap_reports/MIC_data_summary.R")

## Variables----

# Redcap report ID
growth_curves <- "58045"

# Local antifungal data
af_spreadsheet <- "data/metadata/2024_MEC_antifungal_history.xlsx"

# Isolates of interest from manual review of CHEF and CNV data
interesting_isolates <- scan("data/metadata/Calbicans_MEC_aneuploids_CNVs.txt",
                             what = character())

# Subset to species of interest
mic_info <- mic_info %>% 
  filter(genus_species=="C. albicans")

sample_info <- sample_info %>% 
  filter(genus_species=="C. albicans")

# Get growth curve data
gc <- import_report(growth_curves) %>%
  filter(
    redcap_repeat_instrument != "NA",
    !primary_id %in% c("AMS5122", "AMS5123")
  ) %>%
  select(
    primary_id, redcap_repeat_instance,
    gc_date, drug_used,
    gc_temp, gc_time,
    k, sem_k,
    r, sem_r,
    t_gen, sem_tgen,
    auc_l, sem_auc_l,
    auc_e, sem_auc_e
  )

gc <- gc %>%
  group_by(primary_id) %>%
  filter(drug_used == "None", gc_date == max(gc_date), gc_time < 24.1)

gc <- gc %>%
  left_join(mic_info, by = join_by(primary_id)) %>% 
  filter(!is.na(genus_species))

## Drug exposure history----
antifungals <- read_xlsx(af_spreadsheet)

no_results <- sample_info %>% filter(!(primary_id %in% antifungals$primary_id)) %>% 
  select(primary_id, genus_species, relative_days, patient_code, series_id, secondary_id)

remote_exposure <- antifungals %>% 
  filter(primary_id %in% sample_info$primary_id) %>% 
  filter(relative_start_days <= 0, relative_collection_day == 0)

# remote exposure in isolates with CNA, ploidy change, large scale karyotype changes
interesting_remote <- remote_exposure %>% 
  filter(primary_id %in% interesting_isolates)

topical_remotes <- interesting_remote[str_detect(interesting_remote$DRUG_NAME_ORIG, "NYS|CREA|POWD|OINT"),]
systemic_remotes <- interesting_remote[!str_detect(interesting_remote$DRUG_NAME_ORIG, "NYS|CREA|POWD|OINT"),]

calbicans_all_exposure <- antifungals %>% 
  filter(primary_id %in% sample_info$primary_id)

calbicans_topical_exposure <- calbicans_all_exposure[str_detect(calbicans_all_exposure$DRUG_NAME_ORIG, "CLOTRIMAZOLE|NYS|CREA|POWD|OINT"),]
calbicans_systemic_exposure <- calbicans_all_exposure[!str_detect(calbicans_all_exposure$DRUG_NAME_ORIG, "CLOTRIMAZOLE|NYS|CREA|POWD|OINT"),]
