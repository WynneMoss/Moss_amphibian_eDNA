#### Script for cleaning  eDNA metabarcoding data

#### About the data
## eDNA metabarcoding run on samples from 20 ponds in Bay Area
## Collections done in summer 2018
## Riaz metabarcoding primers used
## BLAST setting of 95 and CLUSTER of 100% 

#### In this script:
## Remove contamination

#### Load Libraries ####
# library(devtools)
# install_github("https://github.com/donaldtmcknight/microDecon")
library(microDecon)
library(plyr)
library(tidyverse)
library(gtools)
##### Read in metabarcoding data and re-format ######

# read in data from metaBEAT pipeline

# assigned reads
mba <- read.csv("05_data/Raw/Riaz_12S-trim30_crop114_min90_merge-forwonly_nonchimera_c1cov3_blast95-by-taxonomy-readcounts.blast.csv")

# originally unassigned reads
mbu <- read.csv("05_data/Raw//Riaz_12S-trim30_crop114_min90_merge-forwonly_nonchimera_c1cov3_blast95_unassigned-by-taxonomy-readcounts.blast.csv")

# Remove '-nc.blast' from column names
colnames(mba) <- gsub('.nc.blast', '', colnames(mba))
colnames(mbu) <- gsub('.nc.blast.blast', '', colnames(mbu))

colnames(mba)[1] <- "Assignment"
colnames(mbu)[1] <- "Assignment"

# remove taxonomic level info
mba$Assignment <- gsub("s__|g__|f__|o__|c__|__|p__|sk__", "", 
                       mba$Assignment)
mbu$Assignment <- gsub("s__|g__|f__|o__|c__|p__|__|p__|sk__", "", 
                       mbu$Assignment)


###### MICRODECON #######
# first remove the PCR negatives and positive controls \
# microDecon suggests using blanks, e.g. the field collected blanks that were carried throughout the process

# names:
## PCR.Negatives = PCR negatives
## EN.MES/EN.MJS = Extraction negatives
## PCR.Positives = PCR positives
## NegCtrl = Field negative controls 


mba <- mba %>%select(-contains("EN."), -contains("Plate"))
mbu <- mbu %>%select(-contains("EN."), -contains("Plate"))


# Remove non vertebrates
mba %>% filter(str_detect(taxomomy, "Chordata")) -> mba
mbu %>% filter(str_detect(taxomomy, "Chordata")) -> mbu

# Remove unassigned reads
mba %>% filter(Assignment!="unassigned") %>% select(-taxomomy) %>% smartbind(mbu %>% select(-taxomomy)) -> mb

# Replace NAs with 0
mb[is.na(mb)] <- 0

# Merge read counts by taxonomic assignment - anything with the same 
# Name will be merged and anything unique will be retained
mb[,2:ncol(mb)] <- lapply(mb[,2:ncol(mb)], function(x) as.numeric(as.character(x)))
mb <- plyr::ddply(mb, plyr::.(Assignment), plyr::numcolwise(sum))

# remove empty taxonomic assignments
mb <- mb[rowSums(mb[, -1]) > 0,]


# change genus/family level assignments when only one exists within the study area


# *Taricha granulosa* -> *Taricha torosa* 
# (only Taricha torosa are observed in ponds in the study area over the past decade; Taricha granulosa not present)

## Rename the genus/family assignments 
mb$Assignment <- as.character(mb$Assignment)
mb %>% mutate(Assignment = recode(Assignment,  "Taricha granulosa" = "Taricha torosa")) -> mb
# collapse reads if this created duplicates
mb <- plyr::ddply(mb, plyr::.(Assignment), plyr::numcolwise(sum))

# Rename columns
mb %>% rowid_to_column("OTU_ID") %>% rename(Taxa = Assignment) %>% transform(OTU_ID = as.character(OTU_ID))  -> mb

colSums(mb %>% select_if(is.numeric))

# Need to remove any samples with no reads otherwise package won't work
mb %>% select(-names(which(colSums(mb %>% select_if(is.numeric)) == 0))) -> mb

sum(colSums(mb %>% select_if(is.numeric))==0)

# Make better sample names

# MicroDecon wants to know how many samples are in a "population" 
# Population = pond
sample.names <- data.frame(Original_Sample = colnames(mb))

sample.names %>% mutate(Sample = str_replace_all(Original_Sample, pattern = "CA.", replacement = "CA-")) %>%
  mutate(Sample = str_replace_all(Sample, pattern = "Ca.P", replacement = "CA-P")) %>%
  mutate(Name = sub("\\..*", "", Sample)) %>% group_by(Name) %>% dplyr::mutate(rank = row_number()) %>%
  mutate(SampleName = case_when(
    str_detect(Name, "NegCtrl")~paste("Blank", rank, sep = ""),
    str_detect(Name, "Assignment")~Name,
    str_detect(Name, "OTU_ID")~Name,
    str_detect(Name, "Taxa")~Name,
    TRUE~paste(Name, "_Sample", rank, sep = ""))) -> sample.names

sample.names$SampleName
colnames(mb) <- sample.names$SampleName

# how many samples per "population" 
sample.names %>% group_by(Name) %>% dplyr::summarise(n = n()) %>% filter(!str_detect(Name, "NegCtrl")) %>% filter(!str_detect(Name, "Taxa")) %>% 
  filter(!str_detect(Name, "OTU_ID")) %>% pull(n) -> ns

# order the data frame
mb <- mb[, order(names(mb))]
mb %>% dplyr::select(OTU_ID, contains("Blank"), everything(.)) %>% dplyr::select(-Taxa, everything(.)) -> mb
colnames(mb)

##### Run microDecon #####
dcontam.eDNA <- decon(data = mb, numb.blanks = 5, numb.ind = ns, taxa = T)
dcontam.eDNA$reads.removed  # reads removed from each sample
# which OTUs are getting removed?

mb %>% dplyr::filter(OTU_ID %in% dcontam.eDNA$reads.removed$OTU_ID) %>% pull(Taxa) # only removing bos_taurus, homo_sapiens, and sus_scrofa


dcontam.eDNA$decon.table # final dataset
names(dcontam.eDNA)

mb.clean <- dcontam.eDNA$decon.table
mb.clean %>% pivot_longer(cols = BassGCP_Sample1:YBBA_Sample4, names_to = "SampleName", values_to = "ReadCount") %>% 
  left_join(sample.names %>% select(Original_Sample, Name, SampleName), by = "SampleName") %>% 
  rename(SiteCode = Name) -> mb.clean.long

# save output
write.csv(mb.clean.long, "05_data/metabar_12s_postClean.csv",row.names = FALSE)
