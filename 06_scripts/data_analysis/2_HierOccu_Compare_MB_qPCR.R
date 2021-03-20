#######################################
###### eDNA Amphibian Detection #######
#######################################
# Script 2: Run hierarchical occupancy models for eDNA methods


##### Analyze probability of detection with metabarcoding and qCPR
# hierarchical occupancy analyses 
# using PCR replicate-level data
# evaluate probability of detection in a given PCR replicate and a given DNA extract
# compare metabarcoding and qPCR 
# separate models for California tiger salamander and California red legged frog

###### LOAD LIBRARIES AND DATA ######
library(tidyverse)
library(eDNAoccupancy)

# load detection data, tiger salamanders
amcaDetections <- readRDS("05_data/Occupancy_objects_PCR_replicate/amcaDetections_combined_replicate.Rds")
amcaDetections 
# y is observation history ; K is number of PCRs for each sample
# each column is a sample-by-method combination
# only 5 sites because only 5 sites have unpooled metabarcoding data 

# load detection data, red legged frogs
radrDetections <- readRDS("05_data/Occupancy_objects_PCR_replicate/radrDetections_combined_replicate.Rds")
radrDetections

# load covariate data
# in this case we are only interested in method (qPCR vs. metabar)
surveyCovars_amcaboth <- readRDS("05_data/Occupancy_objects_PCR_replicate/amcaCovariates_combined_replicate.Rds")
surveyCovars_radrboth <- readRDS("05_data/Occupancy_objects_PCR_replicate/radrCovariates_combined_replicate.Rds")

##### FIT HIERARCHICAL OCCUPANCY MODELS #######
# compare qPCR and metabarcoding using hierarchical occupancy models

##### California tiger salamanders ####
# model fitting will take several minutes 
# iterations
niter = 10000
burnin = 5000
amca_both_1 = occModel(formulaSite= ~1,
                       formulaSiteAndSample = ~ Method_MB, # predictor variable: whether MB or qPCR was used
                       formulaReplicate = ~Method_MB,
                       detectionMats = amcaDetections,
                       siteAndSampleData = surveyCovars_amcaboth,
                       niter=niter,
                       niterInterval = 1000,
                       siteColName = "site", sampleColName = "sample")

plotTrace(amca_both_1, paramName = c("alpha.(Intercept)", "beta.(Intercept)", "delta.(Intercept)")) 
plotTrace(amca_both_1, paramName = c("alpha.Method_MB","delta.Method_MB" )) 

# posterior coeff estimates
amca_both_1sum  = posteriorSummary(amca_both_1, burnin = burnin, mcError=TRUE, outputSummary = TRUE)

theta.amca_both_1 = posteriorSummaryOfSampleOccupancy(amca_both_1, burnin = burnin) # sample specific occurrence probabilities

p.amca_both_1 = posteriorSummaryOfDetection(amca_both_1, burnin = burnin)
  

##### California red legged frog #####
# both methods are identical! 
radrDetections$y/radrDetections$K # all are either 1 (all PCRs amplified) or 0 (none amplified) 

surveyCovars_radrboth %>% select(site, sample, Method_MB) %>% pivot_wider(id_cols = site, names_from = sample, values_from = Method_MB) # 1s are using metabar, 0 are qPCR

# can try to fit a model but it won't converge because both methods generated same data

radr_both_1 = occModel(formulaSite= ~1,
                       formulaSiteAndSample = ~ Method_MB,
                       formulaReplicate = ~Method_MB,
                       detectionMats =radrDetections,
                       siteAndSampleData = surveyCovars_radrboth,
                       niter=niter,
                       niterInterval = 1000,
                       siteColName = "site", sampleColName = "sample")
plotTrace(radr_both_1, paramName = c("alpha.(Intercept)", "beta.(Intercept)", "delta.(Intercept)")) # yikes
plotTrace(radr_both_1, paramName = c("alpha.Method_MB","delta.Method_MB" )) 
posteriorSummaryOfDetection(radr_both_1, burnin = burnin)






