#######################################
###### eDNA Amphibian Detection #######
#######################################
# Script 3: Run hierarchical occupancy models for qPCR replicate data

##### Analyze probability of detection with qCPR
# hierarchical occupancy analyses 
# using PCR replicate-level data (n = 20 sites)
# evaluate probability of detection in a given PCR replicate and a given DNA extract
# separate models for California tiger salamander and California red legged frog

#############################################################
###### LOAD LIBRARIES AND DATA ######
#############################################################
library(tidyverse)
library(eDNAoccupancy)

# covariate data for each DNA extract
survey_covars <- readRDS("05_data/Occupancy_objects_PCR_replicate/covariates_qPCR_replicate.Rds") # covariate data
list2env(survey_covars, envir = .GlobalEnv)
rm(survey_covars)

# detections from individual qPCR reactions
amcaDetections_q <- readRDS("05_data/Occupancy_objects_PCR_replicate/amcaDetections_qPCR_replicate.Rds")
radrDetections_q <- readRDS("05_data/Occupancy_objects_PCR_replicate/radrDetections_qPCR_replicate.Rds")

#############################################################
########## FIT HIERARCHICAL OCCUPANCY MODELS ###############
#############################################################

##### California tiger salamanders ######

# possible covariates for detection probability
varlist <- c("filter_type_coarse", "sArea", "sVolume", "sDay",
             "sSalinity", "sTurbidity", "sVolume", "sVolume.rep", "Purified")

# volume of water filtered will be in both levels of model; need distinct names
amca_qcovars$sVolume.rep <- amca_qcovars$sVolume

# examine correlation ####
amca_qcovars %>% select(varlist) %>% cor()
# sArea and turbidity are correlated (-0.56) but are not in the same component of the model

# assign predictors to one of the hierarchical levels

# prob of detection at the sample (DNA extract level): theta
sample.level <- c("filter_type_coarse", "sArea", "sVolume", "sDay")

# prob of detection at the PCR replicate level: p 
rep.level <- c("sSalinity", "sTurbidity", "sVolume.rep", "Purified")

# all possible combinations of covariates
combos.sample.level <- NULL
for(n in 0:length(sample.level)){
  test <- combn(sample.level, n, simplify = F)
  combos.sample.level <- append(combos.sample.level, test)
}
combos.replicate.level <- NULL
for(n in 0:length(rep.level)){
  test <- combn(rep.level, n, simplify = F)
  combos.replicate.level <- append(combos.replicate.level, test)
}

# Create the output data frame
model.summary.amca <- data.frame(model.number = numeric(),  
                            model.sample.vars = character(),
                            model.replicate.vars = character(),
                            PPLC = numeric(),
                            WAIC = numeric())
model.number = 0
iterations = 10000  # 
burnin.no = 5000 # 

# Now run a for loop to perform all model combinations #

##### WARNING!!!!! 
# This will take a long time # 
# Change following to TRUE if you want to fit all the models 
run_models = FALSE
if(run_models == TRUE){
  for (r in 1:length(combos.replicate.level)){
    for(s in 1:length(combos.sample.level)){
      model.number <- model.number + 1  # for output
      print(paste0("Starting model number ",model.number))
      samp.vars.vec <- combos.sample.level[[s]]    # pick out the vars to be used
      rep.vars.vec <- combos.replicate.level[[r]]   # and which are replicate-level
      
      # prepare lists of variables - or replacements if no variables of that type
      if (length(samp.vars.vec) > 0){  # if using any sample variable, write it out...
        samp.vars <- paste("~ ", paste(samp.vars.vec, collapse = ' + '))
      } else {                         # otherwise write "~ 1"
        samp.vars <- "~ 1"
      }
      
      if (length(rep.vars.vec) > 0){   # as above for replicate variables
        rep.vars <- paste("~ ", paste(rep.vars.vec, collapse = ' + '))
      } else {
        rep.vars <- "~ 1"
      }
      # construct the model
      testm = occModel(formulaSite = ~ 1,
                       formulaSiteAndSample = samp.vars,
                       formulaReplicate = rep.vars,
                       detectionMats = amcaDetections_q,
                       siteAndSampleData=amca_qcovars,
                       niter=iterations,
                       niterInterval=iterations,
                       siteColName = "site", sampleColName = "sample")
      
      ## Evaluate fit of sample model
      PPLC <- posteriorPredictiveLoss(testm, burnin=burnin.no)[[1]] 
      WAIC <- WAIC(testm, burnin=burnin.no)[[1]]
      
      # prepare output    
      out <- cbind(model.number, samp.vars,rep.vars, PPLC, WAIC)
      model.summary.amca <- rbind(model.summary.amca,out)
    }
  }
  model.summary.amca %>% mutate(WAIC = as.numeric(WAIC), PPLC = as.numeric(PPLC)) -> model.summary.amca
  model.summary.amca %>% arrange(WAIC)
  model.summary.amca %>% arrange(PPLC)
}

# re fit best model
bestm = occModel(formulaSite = ~ 1,
                 formulaSiteAndSample = ~ sDay,
                 formulaReplicate = ~ sVolume.rep,
                 siteAndSampleData  = amca_covar,
                 detectionMats=amcaDetections_q,
                 niter=iterations,
                 niterInterval=iterations,
                 siteColName = 'site', sampleColName = "sample")
posteriorSummary(bestm, burnin=burnin.no, mcError=TRUE)

# assess convergence
plotTrace(bestm, c('beta.(Intercept)', 'alpha.(Intercept)', 'delta.(Intercept)'), 
          burnin=burnin.no) 
# Autocorrelation plots of the parameters
plotACF(bestm, c('beta.(Intercept)', 'alpha.(Intercept)', 'delta.(Intercept)'), 
        burnin=burnin.no)
psi = posteriorSummaryOfSiteOccupancy(bestm, burnin=burnin.no)
theta = posteriorSummaryOfSampleOccupancy(bestm, burnin=burnin.no)
p = posteriorSummaryOfDetection(bestm, burnin=burnin.no)
preds.best <- cbind(psi=psi$median, theta=theta$median[,1], p=p$median[,1])

## Evaluate fit of best model
PPLC.best <- posteriorPredictiveLoss(bestm, burnin=burnin.no)
WAIC.best <- WAIC(bestm, burnin=burnin.no)
AUC.best <- posteriorSummaryOfAUC(bestm, burnin=burnin.no)

# compare to a null model
nullm = occModel(formulaSite = ~ 1,
                 formulaSiteAndSample = ~ 1,
                 formulaReplicate = ~ 1,
                 siteAndSampleData  = amca_covar,
                 detectionMats=amcaDetections_q,
                 niter=iterations,
                 niterInterval=iterations,
                 siteColName = 'site', sampleColName = "sample")
null.table <- posteriorSummary(nullm, burnin=burnin.no, mcError=TRUE)
plotTrace(nullm, c('beta.(Intercept)', 'alpha.(Intercept)', 'delta.(Intercept)'), 
          burnin=burnin.no) 
plotACF(nullm, c('beta.(Intercept)', 'alpha.(Intercept)', 'delta.(Intercept)'), 
        burnin=burnin.no)
psi.n = posteriorSummaryOfSiteOccupancy(nullm, burnin=burnin.no)
theta.n = posteriorSummaryOfSampleOccupancy(nullm, burnin=burnin.no)
p.n = posteriorSummaryOfDetection(nullm, burnin=burnin.no)
preds.null <- cbind(psi=psi$median, theta=theta$median[,1], p=p$median[,1])

PPLC.null <- posteriorPredictiveLoss(nullm, burnin=burnin.no)
WAIC.null <- WAIC(nullm, burnin=burnin.no)
AUC.null <- posteriorSummaryOfAUC(nullm, burnin=burnin.no)

##### California red-legged frogs ######
# possible covariates for detection probability
sample.level <- c("filter_type_coarse", "sArea", "sVolume", "sDay")
# no replicate level data because of high detection

combos.sample.level <- NULL
for(n in 0:length(sample.level)){
  test <- combn(sample.level, n, simplify = F)
  combos.sample.level <- append(combos.sample.level, test)
}

# Create the output data frame 
model.summary.radr <- data.frame(model.number = numeric(),  
                            model.sample.vars = character(),
                            PPLC = numeric(),
                            WAIC = numeric())

model.number = 0
iterations = 10000  
burnin.no = 5000 

# Now run for loop to fit all model combinations 
##### WARNING!!!!! 
# This will take a long time # 
# Change following to TRUE if you want to fit all the models 
run.models = FALSE
if(run.models==TRUE){
  for(s in 1:length(combos.sample.level)){
    model.number <- model.number + 1  # for output
    print(paste0("Starting model number ",model.number))
    samp.vars.vec <- combos.sample.level[[s]]    # pick out the vars to be used
    # prepare lists of variables - or replacements if no variables of that type
    if (length(samp.vars.vec) > 0){  # if using any sample variable, write it out...
      samp.vars <- paste("~ ", paste(samp.vars.vec, collapse = ' + '))
    } else {                         # otherwise write "~ 1"
      samp.vars <- "~ 1"
    }
    
    # construct the model
    testm = occModel(formulaSite = ~ 1,
                     formulaSiteAndSample = samp.vars,
                     formulaReplicate = ~1,
                     detectionMats = radrDetections_q,
                     siteAndSampleData=radr_qcovars,
                     niter=iterations,
                     niterInterval=iterations,
                     siteColName = "site", sampleColName = "sample")
    
    ## Evaluate fit of sample model
    PPLC <- posteriorPredictiveLoss(testm, burnin=burnin.no)[[1]]
    WAIC <- WAIC(testm, burnin=burnin.no)[[1]]
    
    # prepare output    
    out <- cbind(model.number, samp.vars, PPLC, WAIC)
    model.summary.radr <- rbind(model.summary.radr,out)
  }
  
  model.summary.radr %>% mutate(WAIC = as.numeric(WAIC), PPLC = as.numeric(PPLC)) -> model.summary.radr
  model.summary.radr %>% arrange(WAIC)
  model.summary.radr %>% arrange(PPLC)
  
}

# re fit best model
bestm = occModel(formulaSite = ~ 1,
                 formulaSiteAndSample = ~ 1,
                 formulaReplicate = ~ sVolume,
                 siteAndSampleData  = radr_covar,
                 detectionMats=radrDetections_q,
                 niter=iterations,
                 niterInterval=iterations,
                 siteColName = 'site', sampleColName = "sample")
bestm.table <- posteriorSummary(bestm, burnin=burnin.no, mcError=TRUE)


# diagnostics 
plotTrace(bestm, c('beta.(Intercept)', 'alpha.(Intercept)', 'delta.(Intercept)', "alpha.sVolume"), 
          burnin=burnin.no) 
plotACF(bestm, c('beta.(Intercept)', 'alpha.(Intercept)', 'delta.(Intercept)', "alpha.sVolume"), 
        burnin=burnin.no)

psi = posteriorSummaryOfSiteOccupancy(bestm, burnin=burnin.no)
theta = posteriorSummaryOfSampleOccupancy(bestm, burnin=burnin.no)
p = posteriorSummaryOfDetection(bestm, burnin=burnin.no)
preds.best <- cbind(psi=psi$median, theta=theta$median[,1], p=p$median[,1])

## Evaluate fit of best model
PPLC.best <- posteriorPredictiveLoss(bestm, burnin=burnin.no)
WAIC.best <- WAIC(bestm, burnin=burnin.no)
AUC.best <- posteriorSummaryOfAUC(bestm, burnin=burnin.no)

nullm = occModel(formulaSite = ~ 1,
                 formulaSiteAndSample = ~ 1,
                 formulaReplicate = ~ 1,
                 siteAndSampleData  = radr_covar,
                 detectionMats=radrDetections_q,
                 niter=iterations,
                 niterInterval=burnin.no,
                 siteColName = 'site', sampleColName = "sample")

null.table <- posteriorSummary(nullm, burnin=burnin.no, mcError=TRUE)
plotTrace(nullm, c('beta.(Intercept)', 'alpha.(Intercept)', 'delta.(Intercept)'), 
          burnin=burnin.no)
plotACF(nullm, c('beta.(Intercept)', 'alpha.(Intercept)', 'delta.(Intercept)'), 
        burnin=burnin.no)

psi.n = posteriorSummaryOfSiteOccupancy(nullm, burnin=burnin.no)
theta.n = posteriorSummaryOfSampleOccupancy(nullm, burnin=burnin.no)
p.n = posteriorSummaryOfDetection(nullm, burnin=burnin.no)
preds.null <- cbind(psi=psi$median, theta=theta$median[,1], p=p$median[,1])

PPLC.null <- posteriorPredictiveLoss(nullm, burnin=burnin.no)
WAIC.null <- WAIC(nullm, burnin=burnin.no)
AUC.null <- posteriorSummaryOfAUC(nullm, burnin=burnin.no)
