#######################################
###### eDNA Amphibian Detection #######
#######################################
# Script 03: Reformat replicate level data 


###### ABOUT: #########################
# Script to reformat qPCR and metabarcoding data
# Want detection at the replicate-level
# This script will get detection/non-detection data from PCR-level
# Reformat for hierarchical occupancy analyses
rm(list=ls())
##### LOAD LIBRARIES #######
library(tidyverse)
library(eDNAoccupancy)
library(lubridate)
library(lme4)



###### 1: READ IN EDNA SAMPLE COVARIATES ######
sample_covariates <- read.csv(file = "05_data/Raw/eDNA_sample_covariates.csv")
sample_covariates %>% select(Purified, Filter_volume_mL, Area, Day, Salinity, Turbidity) %>% cor()
###### 2: READ IN qPCR DATA #####
# samples collected and run in 2018
# qPCR assays on California red legged frog (RADR) and California tiger salamander (AMCA)
# data were in triplicate and then re run if they were not consistent 
# therefore there are varying numbers of assays per sample
# inhibited samples are purified (do not keep these runs)

qpcr <- read.csv("05_data/Raw/qPCR_raw_data.csv")

# clean up field code names
qpcr %>% mutate(Field_Sample = toupper(Field_Sample)) %>% 
  mutate(Field_Sample = str_replace_all(Field_Sample, "GRAMPS", "GRAMP")) -> qpcr 

# reformat so it can link to co-variate data
qpcr %>% separate(Field_Sample, into = c("SITE_CODE", "Collection_Date", "Replicate"), remove = FALSE, sep = "_") %>% 
  mutate(Collection_Date = as.Date(Collection_Date, "%Y%m%d"),
         SITE_CODE = toupper(SITE_CODE),
         Field_Sample = toupper(Field_Sample)) %>% 
  group_by(SITE_CODE) %>% mutate(Visit = dense_rank(Collection_Date)) %>% 
  mutate(SITE_VISIT = paste(SITE_CODE, Visit, sep = "_V0")) %>% 
  select(SITE_CODE, Visit, SITE_VISIT, Collection_Date, Replicate, Sample_ID = Field_Sample, Inhibited, Quant = Starting.Quantity..SQ., Plate,
         Target) %>% ungroup() -> qpcr

# check for false positives
# false positive amplification in one plate but not the other
qpcr %>% group_by(Sample_ID, Plate) %>% filter(Target == "AMCA3" | Target=="CTS") %>%  filter(is.na(Inhibited)) %>% summarise(POS = sum(Quant>0), RUNS = sum(is.na(Inhibited)), FRAC = POS/RUNS)  %>% group_by(Sample_ID) %>% arrange(Plate) %>% mutate(Plate_no = rank(Plate)) %>% 
  pivot_wider(id_cols = Sample_ID, names_from = Plate_no, names_prefix = "Plate_", values_from = FRAC) %>% 
  filter(Plate_1 < 1) %>% arrange(Sample_ID) %>% filter(!is.na(Plate_2)) %>% filter(Plate_1 ==0|Plate_2 == 0) %>% 
  filter(!str_detect(Sample_ID, "NEG")) %>% pull(Sample_ID) -> false.pos.amca


qpcr %>% filter(!(Sample_ID %in% false.pos.amca & (Target == "AMCA3"|Target=="CTS"))) -> qpcr # remove this sample

qpcr %>% group_by(Sample_ID, Plate) %>% filter(Target == "RADR") %>%  filter(is.na(Inhibited)) %>% summarise(POS = sum(Quant>0), RUNS = sum(is.na(Inhibited)), FRAC = POS/RUNS)  %>% group_by(Sample_ID) %>% arrange(Plate) %>% mutate(Plate_no = rank(Plate)) %>% 
  pivot_wider(id_cols = Sample_ID, names_from = Plate_no, names_prefix = "Plate_", values_from = FRAC) %>% 
  filter(Plate_1 < 1) %>% arrange(Sample_ID) %>% filter(!is.na(Plate_2)) %>% filter(Plate_1 ==0|Plate_2 == 0) %>% 
  filter(!str_detect(Sample_ID, "NEG")) %>% pull(Sample_ID) -> false.pos.radr
# there are no false positives 

# add a species column
qpcr %>% mutate(
  Species = case_when(Target == "AMCA3" | Target=="CTS" ~ "AMCA",
                      Target == "RADR" ~ "RADR")) -> qpcr

# filter out inhibited samples 
qpcr %>% filter(is.na(Inhibited)) -> qpcr

# filter out negative controls after verifying they do not contain DNA
qpcr %>% filter(str_detect(SITE_CODE, "NEG")) %>% pull(Quant) %>% quantile()
qpcr %>% filter(!str_detect(SITE_CODE, "NEG")) -> qpcr


###### 3: READ IN METABARCODING DATA ######
mb <- read.csv("05_data/metabar_12s_postClean.csv")

mb %>% filter(str_detect(Original_Sample, "PCR0")) %>% # just get unpooled samples
  mutate(Species = case_when(
    str_detect(Taxa, "Taricha")~"TATO",
    str_detect(Taxa, "Anaxyrus")~"BUBO",
    str_detect(Taxa, "Pseudacris")~"PSRE",
    str_detect(Taxa, "draytonii")~"RADR",
    str_detect(Taxa, "Ambystoma")~"AMCA",
    str_detect(Taxa, "catesbeiana")~"RACA",
    TRUE~"Other"
  )) %>% 
  group_by(Species, Original_Sample, SiteCode, SampleName) %>% 
  summarise(ReadCount = sum(ReadCount)) %>% 
  mutate(Field_Sample = toupper(Original_Sample)) %>% 
  mutate(Field_Sample = str_replace_all(Field_Sample, "CA\\.", "CA-")) %>% 
  mutate(Field_Sample = str_replace_all(Field_Sample, pattern = "Ca.P", "CA-P")) %>% 
  mutate(Field_Sample = str_replace_all(Field_Sample, "GRAMPS", "Gramp")) %>% 
  mutate(Field_Sample = str_replace_all(toupper(Field_Sample), "\\.", "_"),
         ReadCount = case_when(ReadCount>0~1, ReadCount==0~0)) %>% 
  separate(Field_Sample, into = c("SITE_CODE", "Collection_Date", "Replicate"), remove = FALSE, sep = "_")  %>% 
  mutate(Collection_Date = as.Date(Collection_Date, "%Y%m%d")) %>% 
  group_by(SITE_CODE) %>% mutate(Visit = dense_rank(Collection_Date)) %>% 
  mutate(SITE_VISIT = paste(SITE_CODE, Visit, sep = "_V0")) %>% 
  select(SITE_CODE, Visit, SITE_VISIT, Collection_Date, Replicate, Sample_ID = Field_Sample, ReadCount, Species) %>% ungroup() %>% 
  rename(Original_Sample_MB = Sample_ID) %>% 
  mutate(Sample_ID = str_replace_all(Original_Sample_MB, "_PCR.*", "")) -> mb

mb$SITE_CODE %>% unique() # five sites
mb$SITE_VISIT %>% unique() # just five site-visits
mb$Sample_ID %>% unique() %>% length() # 10 samples

mb %>% mutate(
  Visit=as.integer(Visit),
  ReadCount = as.numeric(ReadCount)
) -> mb


###### 4: REFORMAT ALL qPCR DATA FOR OCCUPANCY ######
# for all 20 sites
qpcr %>% filter(Species=="AMCA") %>% mutate(Quant = case_when(
  Quant > 0 ~ 1,
  Quant == 0 ~0))  %>% 
  group_by(Sample_ID) %>%
  mutate(qPCR_replicate = rank(Quant, ties.method = "random"),
         site = SITE_CODE) %>% 
  pivot_wider(id_cols = c(site, Sample_ID),values_from = Quant, names_from = qPCR_replicate, names_sort = TRUE,
              names_prefix = "pcr") %>% group_by(site) %>%  mutate(sample = rank(Sample_ID)) %>% arrange(site, sample) %>% mutate(sample = as.integer(sample)) %>% 
  ungroup()-> amca_det_q

qpcr %>% filter(Species=="RADR") %>% mutate(Quant = case_when(
  Quant > 0 ~ 1,
  Quant == 0 ~0))  %>% 
  group_by(Sample_ID) %>%
  mutate(qPCR_replicate = rank(Quant, ties.method = "random"),
         site = SITE_CODE) %>% 
  pivot_wider(id_cols = c(site, Sample_ID),values_from = Quant, names_from = qPCR_replicate, names_sort = TRUE,
              names_prefix = "pcr") %>% group_by(site) %>%  mutate(sample = rank(Sample_ID)) %>% arrange(site, sample) %>% mutate(sample = as.integer(sample))%>% 
  ungroup()-> radr_det_q

#### get covariates
sample_covariates %>% ungroup() %>% mutate(log_Turbidity = log(Turbidity)) %>% 
  mutate(sArea = scale(log(Area)), sDay = scale(Day),
         sTurbidity = scale(log_Turbidity),
         sVolume = scale(Filter_volume_mL), 
         sSalinity = scale(Salinity),
         filter_type_coarse = as.numeric(case_when(Filter_type == "Normal"~0, TRUE ~1))) %>% 
  select(site = SITE_CODE, Sample_ID,  Purified, sArea, sTurbidity,  sVolume, sSalinity,filter_type_coarse, sDay) %>% data.frame() -> survey_covs

radr_det_q %>% left_join(survey_covs, by = c("site", "Sample_ID"))  %>% data.frame() -> radr_qcovars
amca_det_q %>% left_join(survey_covs, by = c("site", "Sample_ID"))  %>% data.frame() -> amca_qcovars

# format with occupancy package
amcaDetections_q = occData(as.data.frame(amca_det_q%>% select(-Sample_ID) %>% select(site, sample, everything(.))), siteColName = "site", sampleColName = "sample")
radrDetections_q = occData(as.data.frame(radr_det_q%>% select(-Sample_ID) %>% select(site, sample, everything(.))), siteColName = "site", sampleColName = "sample")


###### 5: REFORMAT ALL COMBINED qPCR + METABAR DATA FOR OCCUPANCY ######
# for 5 sites where there is replicate level mb data
qpcr %>% mutate(Sample_ID_Method = paste(Sample_ID, "qPCR", sep = "_")) -> qpcr
mb %>% mutate(Sample_ID_Method = paste(Sample_ID, "MB", sep = "_")) -> mb


#### AMCA ####
mb %>% filter(Species == "AMCA") %>% mutate(Quant = case_when(ReadCount > 0 ~ 1,ReadCount == 0 ~0)) %>% select(site = SITE_CODE, Sample_ID, Sample_ID_Method, Quant) -> mb_amca
qpcr %>% filter(Species == "AMCA") %>% mutate(Quant = case_when(Quant > 0 ~ 1,Quant == 0 ~0)) %>% select(site = SITE_CODE, Sample_ID, Sample_ID_Method, Quant) -> qpcr_amca

# combine metabarcoding and qPCR data for AMCA
rbind(mb_amca, qpcr_amca) -> amca_edna_comb

# now reformat to just detection data
amca_edna_comb %>% group_by(Sample_ID_Method) %>%
  mutate(qPCR_replicate = rank(Quant, ties.method = "random")) %>% 
  pivot_wider(id_cols = c(site, Sample_ID_Method, Sample_ID),values_from = Quant, names_from = qPCR_replicate, names_sort = TRUE,
              names_prefix = "pcr") %>% group_by(site) %>%  mutate(sample = rank(Sample_ID, ties.method = "random")) %>% arrange(site, sample) %>% mutate(sample = as.integer(sample)) %>% 
  ungroup()-> amca_det_eDNA

# get covariates
amca_det_eDNA %>% filter(site %in% mb_amca$site) %>% # just filter to 5 sites
  select(site, Sample_ID, sample, Sample_ID_Method) %>% 
  left_join(survey_covs, by = c("site", "Sample_ID")) %>%
  mutate(Method_MB = case_when(str_detect(Sample_ID_Method, "MB")~1, TRUE~0)) %>% data.frame() -> amca_edna_covars_combined

# get detection data
amca_det_eDNA %>% filter(site %in% mb_amca$site) %>% 
  select(-Sample_ID, -Sample_ID_Method) %>% 
  select(site, sample, everything(.)) -> amca_edna_obs_combined

# format for eDNA occupancy package:
amcaDetections_eDNA_combined= occData(as.data.frame(amca_edna_obs_combined), siteColName = "site", sampleColName = "sample")

#### RADR ####
mb %>% filter(Species == "RADR") %>% mutate(Quant = case_when(ReadCount > 0 ~ 1,ReadCount == 0 ~0)) %>% select(site = SITE_CODE, Sample_ID, Sample_ID_Method, Quant) -> mb_radr
qpcr %>% filter(Species == "RADR") %>% mutate(Quant = case_when(Quant > 0 ~ 1,Quant == 0 ~0)) %>% select(site = SITE_CODE, Sample_ID, Sample_ID_Method, Quant) -> qpcr_radr

rbind(mb_radr, qpcr_radr) -> radr_edna_comb

# now reformat to just detection data
radr_edna_comb %>% group_by(Sample_ID_Method) %>%
  mutate(qPCR_replicate = rank(Quant, ties.method = "random")) %>% 
  pivot_wider(id_cols = c(site, Sample_ID_Method, Sample_ID),values_from = Quant, names_from = qPCR_replicate, names_sort = TRUE,
              names_prefix = "pcr") %>% group_by(site) %>%  mutate(sample = rank(Sample_ID, ties.method = "random")) %>% arrange(site, sample) %>% mutate(sample = as.integer(sample)) %>% 
  ungroup()-> radr_det_eDNA

# get covariates
radr_det_eDNA %>% filter(site %in% mb_radr$site) %>% 
  select(site, Sample_ID, sample, Sample_ID_Method) %>%
  left_join(survey_covs, by = c("site", "Sample_ID")) %>%
  mutate(Method_MB = case_when(str_detect(Sample_ID_Method, "MB")~1, TRUE~0)) %>% data.frame() -> radr_edna_covars_combined

# get detection data
radr_det_eDNA %>% filter(site %in% mb_radr$site) %>% 
  select(-Sample_ID, -Sample_ID_Method) %>% 
  select(site, sample, everything(.)) -> radr_edna_obs_combined

# format for eDNA occupancy package:
radrDetections_eDNA_combined = occData(as.data.frame(radr_edna_obs_combined), siteColName = "site", sampleColName = "sample")

###### EXPORT DATA ######
saveRDS(amcaDetections_q, "05_data/Occupancy_objects_PCR_replicate/amcaDetections_qPCR_replicate.Rds") # full qPCR replicate data - detections (AMCA)
saveRDS(radrDetections_q, "05_data/Occupancy_objects_PCR_replicate/radrDetections_qPCR_replicate.Rds") # full qPCR replicate data - detections (RADR)

covars <- mget(ls()[grep("covars", ls())])
saveRDS(covars, "05_data/Occupancy_objects_PCR_replicate/covariates_qPCR_replicate.Rds") # full qPCR replicate data - covariates (both species)


saveRDS(amcaDetections_eDNA_combined, "05_data/Occupancy_objects_PCR_replicate/amcaDetections_combined_replicate.Rds") # combined qPCR/metabar data - detections
saveRDS(amca_edna_covars_combined, "05_data/Occupancy_objects_PCR_replicate/amcaCovariates_combined_replicate.Rds")# combined qPCR/metabar data - dcovariates

saveRDS(radrDetections_eDNA_combined, "05_data/Occupancy_objects_PCR_replicate/radrDetections_combined_replicate.Rds")
saveRDS(radr_edna_covars_combined, "05_data/Occupancy_objects_PCR_replicate/radrCovariates_combined_replicate.Rds")


