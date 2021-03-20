#######################################
###### eDNA Amphibian Detection #######
#######################################
# Script 4: GLMMs for metabarcoding sample level data

##### Analyze probability of detection at the sample level
# using GLMMs
# metabarcoding detection data
# evaluate the probability of detection as a function of pond attributes and sample attributes
# each row of the data = detection of a species from one DNA extract (usually a total of 4 extracts per pond)


#### load libraries
library(tidyverse)
library(MuMIn)
library(glmmTMB)
rm(list=ls())
###########################################
###### PREPARE DATA #######################
###########################################
###### READ IN DATA #######################
# sample level metabarcoding data
mb <- read.csv("05_data/metabar_12s_postClean.csv")
# collection/survey covariates
sample_covars <- read.csv("05_data/Raw/eDNA_sample_covariates.csv")
# true occupancy
detection.all <- read.csv("05_data/Detection_all_methods.csv")

###### REFORMAT METABARCODING  DATA #######
# Want a detection for each DNA extract (across all PCR replicates/sequences)
mb$Taxa <- str_replace_all(mb$Taxa, " ", "_")

# remove the unpooled replicate data 
mb %>% filter(!str_detect(Original_Sample, "PCR0")) -> mb

# clean up columns/sites names
mb %>% 
  mutate(Species = case_when(
    str_detect(Taxa, "Taricha")~"TATO", str_detect(Taxa, "Anaxyrus")~"BUBO", str_detect(Taxa, "Pseudacris")~"PSRE",
    str_detect(Taxa, "draytonii")~"RADR",str_detect(Taxa, "Ambystoma")~"AMCA",str_detect(Taxa, "catesbeiana")~"RACA",
    TRUE~"Other"
  )) %>% 
  group_by(Species, Original_Sample, SiteCode, SampleName) %>% 
  summarise(ReadCount = sum(ReadCount)) %>% 
  mutate(Field_Sample = toupper(Original_Sample)) %>% 
  mutate(Field_Sample = str_replace_all(Field_Sample, "CA\\.", "CA-")) %>%  # rename to match field data
  mutate(Field_Sample = str_replace_all(Field_Sample, pattern = "Ca.P", "CA-P")) %>% # rename to match field data
  mutate(Field_Sample = str_replace_all(Field_Sample, "GRAMPS", "Gramp")) %>% # rename to match field data
  mutate(Field_Sample = str_replace_all(toupper(Field_Sample), "\\.", "_"), # rename to match field data
         MB_Pos = case_when(ReadCount>0~1, ReadCount==0~0)) %>% # standardize to detection (1) nondetection (0)
  separate(Field_Sample, into = c("SITE_CODE", "Collection_Date", "Replicate"), remove = FALSE, sep = "_")  %>% 
  mutate(Collection_Date = as.Date(Collection_Date, "%Y%m%d")) %>% 
  group_by(SITE_CODE) %>% mutate(Visit = dense_rank(Collection_Date)) %>% # number the visits by date
  mutate(SITE_VISIT = paste(SITE_CODE, Visit, sep = "_V0")) %>% ungroup() %>% 
  filter(Species!="Other") -> mb_clean

mb_clean %>% select(SITE_CODE, Sample_ID = Field_Sample, Species, MB_Pos) -> mb_clean


###### JOIN WITH COVARIATE DATA ###########
mb_clean %>% left_join(sample_covars, by = c("Sample_ID", "SITE_CODE")) -> joined_mb
joined_mb %>% select(Sample_ID, SITE_CODE, Species, Filter_volume_mL, Filter_type, Area, Turbidity, MB_Pos, Day) -> joined_mb

# true detection
detection.all %>% pivot_longer(AMCA:TATO, names_to = "Species", values_to = "DETECTION") %>% 
  group_by(SITE_CODE, Species) %>% summarise(True.Occupancy = max(DETECTION, na.rm=TRUE)) %>% 
  right_join(joined_mb, by = c("SITE_CODE", "Species")) -> joined_mb

joined_mb %>% mutate(log_Area = log(Area), log_Turbidity = log(Turbidity)) -> joined_mb
###########################################
###### FIT MODELS #########################
###########################################
###### CALIFORNIA TIGER SALAMANDER  ######
# filter to sites where true occupancy (based on all methods pooled) is 1
joined_mb %>% filter(Species == "AMCA", True.Occupancy==1) %>% data.frame() -> true.amca
# correlation btw predictor variables

true.amca %>% mutate(scale_vol = scale(Filter_volume_mL), scale_log_area = scale(log_Area),
                     scale_log_turb = scale(log_Turbidity), scale_day = scale(Day)) %>% 
  select(scale_vol, scale_log_area, scale_log_turb, scale_day) %>% cor()

# turbidity and area are correlated; (-0.66)
# remove area and just use turbidity

# fit full model
amca.m.full <- glmmTMB(MB_Pos~scale(Filter_volume_mL)  + Filter_type + 
                       scale(log_Area) + scale(Day) + 
                       (1|SITE_CODE), true.amca, family = "binomial", na.action = "na.fail")
# fit all combinations of predictors
amca_dredge <- dredge(amca.m.full, rank = "AICc")
# get models within 2 AICc of best ranking model 
amca_mods <- get.models(amca_dredge, subset = delta <= 2)

# view model objects
amca_mods

###### CALIFORNIA RED-LEGGED FROG  ########
# repeat above 
joined_mb %>% filter(Species == "RADR",True.Occupancy==1) %>% data.frame() -> true.radr
true.radr %>% mutate(scale_vol = scale(Filter_volume_mL), scale_log_area = scale(log_Area), 
                     scale_log_turb = scale(log_Turbidity), scale_day = scale(Day)) %>% 
  select(scale_vol, scale_log_area, scale_log_turb, scale_day) %>% cor()
# volume and turbidity are correlated (-0.67)
# take volume out and keep turbidity
radr.m.full <- glmmTMB(MB_Pos~ Filter_type  +scale(log_Area) +
                       scale(log_Turbidity) + scale(Day) +
                       (1|SITE_CODE), true.radr, family = "binomial", na.action = "na.fail")
radr_dredge <- dredge(radr.m.full, rank = "AICc")
radr_mods <- get.models(radr_dredge, subset = delta <= 2)


###### WESTERN TOAD #######################
joined_mb %>% filter(Species == "BUBO") %>% filter(True.Occupancy==1) %>% data.frame() -> true.bubo

true.bubo %>% mutate(scale_vol = scale(Filter_volume_mL), scale_log_area = scale(log_Area), 
                     scale_log_turb = scale(log_Turbidity), scale_day = scale(Day)) %>% 
  select(scale_vol, scale_log_area, scale_log_turb, scale_day) %>% cor()
# no correlations

bubo.m.full<- glmmTMB(MB_Pos~scale(Filter_volume_mL)  + Filter_type  +
                      scale(log_Turbidity) + scale(log_Area)+ scale(Day)+
                      (1|SITE_CODE),true.bubo, family = "binomial", na.action = "na.fail")
bubo_dredge<- dredge(bubo.m.full, rank = "AICc") 
bubo_mods <- get.models(bubo_dredge, subset = delta <= 2) 

###### PACIFIC CHORUS FROG ################
joined_mb %>% filter(Species == "PSRE") %>% filter(True.Occupancy==1) %>% data.frame() -> true.psre
true.psre %>% mutate(scale_vol = scale(Filter_volume_mL), scale_log_area = scale(log_Area), 
                     scale_log_turb = scale(log_Turbidity), scale_day = scale(Day)) %>% 
  select(scale_vol, scale_log_area, scale_log_turb, scale_day) %>% cor()
# turbidity and area r = -0.55; remove area 
psre.m.full <- glmmTMB(MB_Pos~scale(Filter_volume_mL)  + Filter_type  + scale(log_Turbidity) + 
                        + scale(Day) +
                       (1|SITE_CODE), true.psre, family = "binomial", na.action = "na.fail")
psre_dredge <- dredge(psre.m.full, rank = "AICc") # 
psre_mods <- get.models(psre_dredge, subset = delta <= 2)
# example of getting coeff estimates for best model
summary(psre_mods$`11`)
confint(psre_mods$`11`)

###### AMERICAN BULLFROG ##################
joined_mb %>% filter(Species == "RACA", True.Occupancy==1) %>% data.frame() -> true.raca

# samples where metabarcoding missed the species:
true.raca%>% filter(MB_Pos == 0) %>% nrow() # only 4 samples
# not going to run a model because sample size is too low


###### CALIFORNIA NEWT #####################
joined_mb %>% filter(Species == "TATO") %>% filter(True.Occupancy==1) %>% data.frame() -> true.tato

true.tato %>% mutate(scale_vol = scale(Filter_volume_mL), scale_log_area = scale(log_Area), 
                     scale_log_turb = scale(log_Turbidity), scale_day = scale(Day)) %>% 
  select(scale_vol, scale_log_area, scale_log_turb, scale_day) %>% cor()
# turbidity and area are correlated

tato.m.full <- glmmTMB(MB_Pos~scale(Filter_volume_mL)  + Filter_type  + scale(log_Turbidity) + scale(Day) +
                     (1|SITE_CODE), true.tato, family = "binomial", na.action = "na.fail")

tato_dredge <- dredge(tato.m.full, rank = "AICc")  
tato_mods <- get.models(tato_dredge, subset = delta <=2)




