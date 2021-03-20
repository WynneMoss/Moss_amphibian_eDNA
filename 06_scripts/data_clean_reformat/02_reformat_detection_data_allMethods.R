#######################################
###### eDNA Amphibian Detection #######
#######################################
# Script 01: Reformat data 


###### ABOUT: #########################
# Project: compare methods for surveying amphibians
# Data: Amphibian species detection data collected in 2018 from 20 ponds within California Bay Area
# 3 different types of data: 
## (1) field surveys: dipnets, visual surveys, seines
## (2) qPCR: assays for R. draytonii and A. californiense (Goldberg lab)
## (3) metabarcoding using 12S primers (Davis lab)
# This script will combine detection/non-detection data from methods
# Will get data into a format suitable for occupancy analyses

###### LOAD LIBRARIES ######
library(lubridate) # fix dates
library(unmarked) # occupancy formatting
library(fastDummies) # dummy code variables
library(tidyverse) # manipulate data

###### 1: METABARCODING DATA (MICRODECON PACKAGE) #######
mb<- read.csv("05_data/metabar_12s_postClean.csv") 
# fix taxonomic formatting
mb$Taxa <- str_replace_all(mb$Taxa, " ", "_")

# remove the unpooled replicate data 
mb %>% filter(!str_detect(Original_Sample, "PCR0")) -> mb
# only extract amphibian data
# 6 possible amphibian species:
# Ambystoma californinse (AMCA), Anaxyrus boreas (BUBO), Rana draytonii (RADR), 
# Rana catesbeiana (RACA), Pseudacris regilla (PSRE), Taricha torosa (TATO)
mb %>% 
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
  mutate(Field_Sample = str_replace_all(Field_Sample, "CA\\.", "CA-")) %>%  # rename to match field data
  mutate(Field_Sample = str_replace_all(Field_Sample, pattern = "Ca.P", "CA-P")) %>% # rename to match field data
  mutate(Field_Sample = str_replace_all(Field_Sample, "GRAMPS", "Gramp")) %>% # rename to match field data
  mutate(Field_Sample = str_replace_all(toupper(Field_Sample), "\\.", "_"), # rename to match field data
         ReadCount = case_when(ReadCount>0~1, ReadCount==0~0)) %>% # standardize to detection (1) nondetection (0)
  separate(Field_Sample, into = c("SITE_CODE", "Collection_Date", "Replicate"), remove = FALSE, sep = "_")  %>% 
  mutate(Collection_Date = as.Date(Collection_Date, "%Y%m%d")) %>% 
  group_by(SITE_CODE) %>% mutate(Visit = dense_rank(Collection_Date)) %>% # number the visits by date
  mutate(SITE_VISIT = paste(SITE_CODE, Visit, sep = "_V0")) %>% 
  dplyr::select(SITE_CODE, Visit, SITE_VISIT, Collection_Date, Replicate, Sample_ID = Field_Sample, ReadCount, Species) %>% ungroup() %>% 
  rename(Original_Sample_MB = Sample_ID) %>% 
  mutate(Sample_ID = str_replace_all(Original_Sample_MB, "_PCR.*", "")) -> mb

mb %>% group_by(SITE_CODE, Collection_Date, Species) %>% summarise(ReadCount = max(ReadCount)) %>% 
  mutate(Method = "MB_md", SITE_CODE = toupper(SITE_CODE)) %>% # add method in (metabarcoding_microdecon)
  data.frame() %>% 
  filter(Species!="Other")-> mb.amphib 

# Clean up column names
mb.amphib %>%
  dplyr::select(SITE_CODE, DATE = Collection_Date, SPECIES = Species, PRESENCE=ReadCount, METHOD=Method) %>%  
  ungroup() %>% data.frame()-> mb.amphib

# wide format (each row a site/date, each column a species detection)
mb.amphib %>% pivot_wider(id_cols = c(SITE_CODE, DATE, METHOD), names_from = SPECIES, values_from = PRESENCE) %>% data.frame()-> mb.amphib

# 1 row per water sample (39)

###### 2: FIELD DATA (SEINES/DIPNETS/VES) ########
# detection/non-detection method for each field method
# Dip Net: Detection (1) or nondetection (0) across all 10-15 dipnets per site-visit
# Seine: Detection (1) or nondetection (0) across all 3-5 seines 
# VES: Detection (1) or nondetection (0) from single visual encounter survey

# Dipnet surveys conducted on visit 1 and seines and VES conducted on visit 1 and visit 2

field.amphib <- read.csv(file ="05_data/Raw/Field_survey_detections.csv")

###### 3: QPCR DATA ######
# assays run on only for A. californiense and R. draytonii
# negative controls had no amplification for AMCA, RADR, or BD so no corrections needed
qpcr <- read.csv("05_data/Raw/qPCR_raw_data.csv")


# clean up site codes/columns
qpcr %>% mutate(Field_Sample = toupper(Field_Sample)) %>% 
  mutate(Field_Sample = str_replace_all(Field_Sample, "GRAMPS", "GRAMP")) %>%  # fix site name typo
  separate(Field_Sample, into = c("SITE_CODE", "DATE", "Replicate"), remove = FALSE, sep = "_") %>% 
  mutate(DATE = as.Date(DATE, "%Y%m%d"),
         Species = case_when(Target == "AMCA3" | Target=="CTS" ~ "AMCA", # group together same assay that was named two diff things
                             Target == "RADR" ~ "RADR")) %>% 
  rename(Quant = Starting.Quantity..SQ.) %>% 
  mutate(PRESENCE = case_when(Quant > 0~1, Quant == 0~0)) %>%  # reformat to detection/nondetection data if Quant > 0
  filter(is.na(Inhibited),!str_detect(SITE_CODE, "NEG")) -> qpcr

# check for false positives: amplification in one plate but not the other
qpcr %>% filter(Species == "AMCA") %>% group_by(Field_Sample, Plate) %>%
  summarise(POS = sum(Quant>0), RUNS = sum(is.na(Inhibited)), FRAC = POS/RUNS)  %>% 
  group_by(Field_Sample) %>% arrange(Plate) %>% mutate(Plate_no = rank(Plate)) %>% 
  pivot_wider(id_cols = Field_Sample, names_from = Plate_no, names_prefix = "Plate_", values_from = FRAC) %>% 
  filter(Plate_1 < 1) %>% arrange(Field_Sample) %>%
  filter(!is.na(Plate_2)) %>% filter(Plate_1 ==0|Plate_2 == 0) %>% # check for samples that didn't re-amplify
  filter(!str_detect(Field_Sample, "NEG")) %>% pull(Field_Sample) -> false.pos.amca

# change 1 to a 0 because this is technically a "false positive" (even though RONJRGCP site has AMCA in other samples)

# re-code the false pos
false.pos.amca
qpcr$PRESENCE[which(qpcr$Field_Sample == false.pos.amca)] = 0
# in this case it doesn't matter because the other DNA extract from that site visit contained it

# check for RADR false positives
qpcr %>% filter(Species == "RADR") %>% group_by(Field_Sample, Plate) %>%
  summarise(POS = sum(Quant>0), RUNS = sum(is.na(Inhibited)), FRAC = POS/RUNS)  %>% 
  group_by(Field_Sample) %>% arrange(Plate) %>% mutate(Plate_no = rank(Plate)) %>% 
  pivot_wider(id_cols = Field_Sample, names_from = Plate_no, names_prefix = "Plate_", values_from = FRAC) %>% 
  filter(Plate_1 < 1) %>% arrange(Field_Sample) %>%
  filter(!is.na(Plate_2)) %>% filter(Plate_1 ==0|Plate_2 == 0) %>% # check for samples that didn't re-amplify
  filter(!str_detect(Field_Sample, "NEG")) %>% pull(Field_Sample) -> false.pos.radr # none



# now group across PCR replicates to get survey level detection
qpcr %>%  
  group_by(SITE_CODE, DATE, Species) %>% summarise(PRESENCE = max(PRESENCE)) %>%   # overall presence across runs
  pivot_wider(id_cols = c(SITE_CODE, DATE), names_from = Species, values_from = PRESENCE) %>% 
  ungroup() %>% 
  mutate(METHOD = "qPCR") %>% 
  dplyr::select(one_of(colnames(mb.amphib))) %>% 
  data.frame() -> qpcr.amphib


###### 4: JOIN ALL DETECTION DATA  #####
detection.data <- gtools::smartbind(field.amphib, mb.amphib, qpcr.amphib)

detection.data %>% dplyr::select(SITE_CODE, DATE, METHOD, AMCA, BUBO, PSRE, RACA, RADR, TATO) %>% 
  mutate(DATE = as.Date(DATE, "%Y-%m-%d")) -> detection.data


# how many surveys per site?
# should have 9: 2 eDNA, 2 qPCR, 2 ves, 2 seine, 1 dipnet
detection.data %>% group_by(SITE_CODE) %>% mutate(n=n())

# only the site CA-EDWD doesn't have 9, because it went dry
detection.data %>% group_by(SITE_CODE) %>% mutate(n=n()) %>% filter(n<9) %>% arrange(SITE_CODE, n, DATE) %>% data.frame()

###### 5: COVARIATES FOR DETECTION ######
field.covariates <- read.csv("05_data/Raw/Field_survey_covariates.csv")
head(field.covariates)
# Turbidity, Salinity, Area were only measured during "field" visits
# For eDNA, these values come from the nearest field visit by date

# scale variables 
field.covariates %>% mutate( slog_Turb = scale(log(Turbidity)), 
                         slog_Area = scale(log(Area)),
                         sDay = scale(jDay), # center = 166, scale = 22.8
                         s_week = (jDay - mean(jDay))/7,
                         DATE = as.Date(DATE)) -> covariate.data

rm(mb, qpcr, mb.amphib, qpcr.amphib, field.amphib, field.covariates)



all.data.surv <- left_join(detection.data, covariate.data, by = c("SITE_CODE", "DATE", "METHOD"))

# remove the second CA-EDWD visit where it was dry (violates closure)
all.data.surv <- all.data.surv %>% filter(!(SITE_CODE=="CA-EDWD" &  DATE == "2018-07-17"))

###### 6: REFORMAT FOR OCCUPANCY ########
species <- c("AMCA", "BUBO", "PSRE", "RACA", "RADR", "TATO")
# define function to get detections in a format where each row is a site each column a visit
# separate table for each species
det.func <- function(SPECIES){
  if(SPECIES %in% c("AMCA", "RADR")){
    all.data.surv %>% group_by(SITE_CODE) %>% arrange(DATE) %>%
      mutate(Visit = rank(DATE, ties.method = c("first"))) %>%
      dplyr::select(SITE_CODE, Visit, SPECIES) %>% pivot_wider(id_cols = SITE_CODE, names_from = Visit, names_prefix = "y.", values_from = SPECIES) %>% 
      arrange(SITE_CODE) %>% 
      column_to_rownames("SITE_CODE")->y
    return(y)
  } else 
    all.data.surv %>% filter(METHOD != "qPCR") %>% group_by(SITE_CODE) %>% arrange(DATE) %>%
    mutate(Visit = rank(DATE, ties.method = c("first"))) %>%
    dplyr::select(SITE_CODE, Visit, SPECIES) %>% pivot_wider(id_cols = SITE_CODE, names_from = Visit, names_prefix = "y.", values_from = SPECIES) %>% 
    arrange(SITE_CODE) %>% 
    column_to_rownames("SITE_CODE")->y
  return(y)
}
# detections data frames
AMCA.y <- det.func("AMCA")
BUBO.y <- det.func("BUBO")
PSRE.y <- det.func("PSRE")
RACA.y <- det.func("RACA")
RADR.y <- det.func("RADR")
TATO.y <- det.func("TATO")

# get obs level data 
obs.covar <- c("sDay", "METHOD", "slog_Area")

# function to reformat observation level data
obs_func_wQpcr <- function(covar){
  all.data.surv %>%  group_by(SITE_CODE) %>% arrange(DATE) %>%
    mutate(Visit = rank(DATE, ties.method = c("first"))) %>% dplyr::select(SITE_CODE, Visit, covar) %>% 
    pivot_wider(id_cols = SITE_CODE, names_from = Visit, names_prefix = "y.", values_from = covar) %>% 
    arrange(SITE_CODE) %>% 
    column_to_rownames("SITE_CODE") -> obs_table
assign(covar, obs_table,envir = globalenv())
}

# apply func to each object in the list
for(i in 1:length(obs.covar)){
  obs_func_wQpcr(obs.covar[i]) # leave out pc
}
# now add all those objects to one list
obs.level<- mget(obs.covar)
# delete separate table objects
rm(list=obs.covar)

### create unmarked frames
amca.occu <- unmarkedFrameOccu(y = AMCA.y, obsCovs = obs.level)
radr.occu <- unmarkedFrameOccu(y = RADR.y, obsCovs = obs.level)

# now do other species (no qPCR data)
obs_func <- function(covar){
  all.data.surv %>% filter(METHOD!="qPCR") %>%  group_by(SITE_CODE) %>% arrange(DATE) %>%
    mutate(Visit = rank(DATE, ties.method = c("first"))) %>% dplyr::select(SITE_CODE, Visit, covar) %>% 
    pivot_wider(id_cols = SITE_CODE, names_from = Visit, names_prefix = "y.", values_from = covar) %>% 
    arrange(SITE_CODE) %>% 
    column_to_rownames("SITE_CODE") -> obs_table
  assign(covar, obs_table,envir = globalenv())
}

# apply func to each object in the list
for(i in 1:length(obs.covar)){
  obs_func(obs.covar[i]) # leave out pc
}
# now add all those objects to one list
obs.level<- mget(obs.covar)
# delete separate table objects
rm(list=obs.covar)

bubo.occu <- unmarkedFrameOccu(y = BUBO.y, obsCovs = obs.level)
raca.occu <- unmarkedFrameOccu(y = RACA.y, obsCovs = obs.level)
tato.occu <- unmarkedFrameOccu(y = TATO.y, obsCovs = obs.level)
psre.occu <- unmarkedFrameOccu(y = PSRE.y, obsCovs = obs.level)

###### 7: AMCA DATA ######
# AMCA: don't use dipnets and ves due to perfect separation issue (no detections)
all.data.surv %>% filter(!METHOD %in% c("VES", "Dip Net")) %>% group_by(SITE_CODE) %>% arrange(DATE) %>% 
  mutate(Visit = rank(DATE, ties.method = c("first"))) %>% dplyr::select(SITE_CODE, AMCA, Visit) %>% pivot_wider(
    id_cols = SITE_CODE, names_from = Visit, names_prefix = "y.", values_from = AMCA) %>% arrange(SITE_CODE) %>% 
  column_to_rownames("SITE_CODE") -> AMCA.y.subset

obs.covar <- c("slog_Area",  "sDay",  "METHOD")

obs_func_amca <- function(covar){
  all.data.surv %>% filter(!METHOD %in% c("VES", "Dip Net")) %>% group_by(SITE_CODE) %>% arrange(DATE) %>%
    mutate(Visit = rank(DATE, ties.method = c("first"))) %>% dplyr::select(SITE_CODE, Visit, covar) %>% 
    pivot_wider(id_cols = SITE_CODE, names_from = Visit, names_prefix = "y.", values_from = covar) %>% arrange(SITE_CODE) %>% 
    column_to_rownames("SITE_CODE") -> obs_table
  assign(covar, obs_table,envir = globalenv())
}

# apply func to each object in the list
for(i in 1:length(obs.covar)){
  obs_func_amca(obs.covar[i]) 
}
# now add all those objects to one list
obs.level_amca_subset<- mget(obs.covar)
# delete separate table objects
rm(list=obs.covar)

amca.occu.subset <- unmarkedFrameOccu(y = AMCA.y.subset, obsCovs = obs.level_amca_subset)

###### 8: EXPORT DATA ######
saveRDS(amca.occu, "05_data/Occupancy_objects_combined_methods/amca.occu.2021.Rdata")
saveRDS(amca.occu.subset, "05_data/Occupancy_objects_combined_methods/amca.occu.subset.2021.Rdata")
saveRDS(bubo.occu, "05_data/Occupancy_objects_combined_methods/bubo.occu.2021.Rdata")
saveRDS(psre.occu, "05_data/Occupancy_objects_combined_methods/psre.occu.2021.Rdata")
saveRDS(radr.occu, "05_data/Occupancy_objects_combined_methods/radr.occu.2021.Rdata")
saveRDS(raca.occu, "05_data/Occupancy_objects_combined_methods/raca.occu.2021.Rdata")
saveRDS(tato.occu, "05_data/Occupancy_objects_combined_methods/tato.occu.2021.Rdata")


detection.data %>% write.csv("05_data/Detection_all_methods.csv", row.names = FALSE)



  