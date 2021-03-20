#######################################
###### eDNA Amphibian Detection #######
#######################################
# Script 1: Run occupancy models for all methods


##### Analyze probability of detection with all survey methods
# occupancy analyses 
# using survey-level data
# evaluate probability of detection for a given survey at a given visit
# compare metabarcoding, qPCR, dipnet, seine, VES methods
# separate models for all six species


#######################################
###### LOAD LIBRARIES AND DATA ######
#######################################
library(unmarked)
library(tidyverse)
library(MuMIn)


amca.occu <- readRDS("05_data/Occupancy_objects_combined_methods/amca.occu.subset.2021.Rdata") 
bubo.occu <- readRDS("05_data/Occupancy_objects_combined_methods/bubo.occu.2021.Rdata") 
psre.occu <- readRDS("05_data/Occupancy_objects_combined_methods/psre.occu.2021.Rdata") 
raca.occu <- readRDS("05_data/Occupancy_objects_combined_methods/raca.occu.2021.Rdata") 
radr.occu <- readRDS("05_data/Occupancy_objects_combined_methods/radr.occu.2021.Rdata") 
tato.occu <- readRDS("05_data/Occupancy_objects_combined_methods/tato.occu.2021.Rdata") 

#######################################
##### FIT OCCUPANCY MODELS ###########
#######################################
##### California tiger salamander ###########
# we removed observations from dipnets and VES because all were 0s
# cant fit a model with method if we include these due to perfect separation

# fit full model with all covariates
# default reference level is metabarcoding (1st alphabetically)
amca_occ_full <- occu(~METHOD + sDay + METHOD*sDay + slog_Area~1, amca.occu)
# fit all other combinations
amca_occ_dredge <- dredge(amca_occ_full)
# extract  models within 2 AICc values
amca_occ_dredge_delta <- get.models(amca_occ_dredge, subset = delta <=2)
amca_occ_dredge_delta %>% length() # 

summary(amca_occ_dredge_delta[[1]])# coefficient estimates from best ranking model
confint(amca_occ_dredge_delta[[1]], type = "det") # CI of coefficients

# set a new reference level to obtain other pairwise contrasts 
# new reference level is qPCR
amca.occu@obsCovs %>% mutate(METHOD2 = relevel(METHOD, ref = "qPCR")) %>% data.frame()->amca.occu@obsCovs 
amca_occ_full_qREF <- occu(~METHOD2 + sDay + METHOD2*sDay + slog_Area ~1, amca.occu)
amca_occ_dredge_qREF <- dredge(amca_occ_full_qREF)
amca_occ_dredge_delta_qREF <- get.models(amca_occ_dredge_qREF, subset = delta <= 2)

summary(amca_occ_dredge_delta_qREF[[1]])
confint(amca_occ_dredge_delta_qREF[[1]], type = "det")

# make predictions
amca_pred_df <- expand.grid(METHOD = levels(amca.occu@obsCovs$METHOD), 
                            sDay = seq(from = -1.6, to = 1.6, by = .2), 
                            slog_Area = seq(from = -2, to = 2, by = .2))

amca_occ_pred <- predict(fitList(fits = amca_occ_dredge_delta), # model averaged predictors
                         newdata = amca_pred_df, 
                         type = "det", append = TRUE)

# now, manually add in the 0s for ves and dipnets 
# we don't know the true probability/uncertainty but we want to represent them on the plot
# this is really just for plotting purposes so they have the same categories as other species
add_empty_amca <- data.frame(Predicted = c(0,0), SE = c(0,0), 
                             lower = c(0,0), upper = c(0,0), 
                             METHOD = c("VES", "Dip Net"), sDay = c(0,0), slog_Area = c(0,0))

amca_occ_pred <- rbind(amca_occ_pred, add_empty_amca) 

##### California red legged frog  #####
# includes all methods: qPCR, metabarcoding, seine, dipnet, VES
# dipnet is the reference level
radr_occ_full <- occu(~METHOD + sDay + METHOD*sDay + slog_Area ~1, radr.occu)
radr_occ_dredge <- dredge(radr_occ_full)
radr_occ_dredge_delta <- get.models(radr_occ_dredge, subset = delta <= 2)
length(radr_occ_dredge_delta) # 3 competitive models

summary(radr_occ_dredge_delta[[1]])
confint(radr_occ_dredge_delta[[1]], type = "det")

summary(radr_occ_dredge_delta[[2]])
confint(radr_occ_dredge_delta[[2]], type = "det")

# now set qPCR as reference level
radr.occu@obsCovs %>% mutate(METHOD2 = relevel(METHOD, ref = "qPCR")) %>% data.frame()->radr.occu@obsCovs 
radr_occ_full_qREF <- occu(~METHOD2 + sDay + METHOD2*sDay + slog_Area ~1, radr.occu)
radr_occ_dredge_qREF <- dredge(radr_occ_full_qREF)
radr_occ_dredge_delta_qREF <- get.models(radr_occ_dredge_qREF, subset = delta <= 2)
summary(radr_occ_dredge_delta_qREF[[1]])
confint(radr_occ_dredge_delta_qREF[[1]], type = "det")


# make predictions (model averaged)
radr_pred_df <- expand.grid(METHOD = levels(radr.occu@obsCovs$METHOD), 
                            sDay = seq(from = -1.6, to = 1.6, by = .2), slog_Area = seq(from = -2, to = 2, by = .2))
radr_occ_pred <- predict(fitList(fits = radr_occ_dredge_delta), # model averaged predictors
                         newdata = radr_pred_df, 
                         type = "det", append = TRUE)

##### Pacific chorus frog  #####
psre_occ_full <- occu(~METHOD + sDay + METHOD*sDay + slog_Area~1, psre.occu)
psre_occ_dredge <- dredge(psre_occ_full)
psre_occ_dredge_delta <- get.models(psre_occ_dredge, subset = delta <= 2)

length(psre_occ_dredge_delta) # 2 models

summary(psre_occ_dredge_delta[[1]]) 
confint(psre_occ_dredge_delta[[1]], type = "det")

summary(psre_occ_dredge_delta[[2]]) #
confint(psre_occ_dredge_delta[[2]], type = "det")

# method was not significant so we do not need to relevel
# make predictions
psre_pred_df <- expand.grid(METHOD = levels(psre.occu@obsCovs$METHOD), 
                            sDay = seq(from = -1.6, to = 1.6, by = .2), slog_Area = seq(from = -2, to = 2, by = .2))
psre_occ_pred <- predict(fitList(fits = psre_occ_dredge_delta), # model averaged predictions
                         newdata = psre_pred_df, 
                         type = "det", append = TRUE)

##### American bullfrog  #####
raca_occ_full <- occu(~METHOD + sDay + METHOD*sDay + slog_Area~1,raca.occu)
raca_occ_dredge <- dredge(raca_occ_full)
raca_occ_dredge_delta <- get.models(raca_occ_dredge, subset = delta <= 2)
length(raca_occ_dredge_delta) 
summary(raca_occ_dredge_delta[[1]])
confint(raca_occ_dredge_delta[[1]], type = "det")
summary(raca_occ_dredge_delta[[2]])
confint(raca_occ_dredge_delta[[2]], type = "det")

# set seine as ref level
raca.occu@obsCovs %>% mutate(METHOD2 = relevel(METHOD, ref = "Seine")) %>% data.frame()->raca.occu@obsCovs 
raca_occ_full_sREF <- occu(~METHOD2 + sDay + METHOD2*sDay + slog_Area ~1, raca.occu)
raca_occ_dredge_sREF <- dredge(raca_occ_full_sREF)
raca_occ_dredge_delta_sREF <- get.models(raca_occ_dredge_sREF, subset = delta <= 2)

summary(raca_occ_dredge_delta_sREF[[1]])
confint(raca_occ_dredge_delta_sREF[[1]],type = "det")

# make predictions
raca_pred_df <- expand.grid(METHOD = levels(raca.occu@obsCovs$METHOD), 
                            sDay = seq(from = -1.6, to = 1.6, by = .2), slog_Area = seq(from = -2, to = 2, by = .2))
raca_occ_pred <- predict(fitList(fits = raca_occ_dredge_delta), 
                         newdata = raca_pred_df, 
                         type = "det", append = TRUE)


##### Western toad ######
bubo_occ_full <- occu(~METHOD + sDay + METHOD*sDay + slog_Area ~1,bubo.occu)
bubo_occ_dredge <- dredge(bubo_occ_full)
bubo_occ_dredge_delta <- get.models(bubo_occ_dredge, subset = delta <= 2)
length(bubo_occ_dredge_delta) # only one model is competitive
summary(bubo_occ_dredge_delta[[1]]) # reference level is dipnets
confint(bubo_occ_dredge_delta[[1]], type = "det") # reference level is dipnets

# set VES as reference level
bubo.occu@obsCovs %>% mutate(METHOD2 = relevel(METHOD, ref = "VES")) %>% data.frame()->bubo.occu@obsCovs 
bubo_occ_full_vREF <- occu(~METHOD2 + sDay + METHOD2*sDay + slog_Area ~1, bubo.occu)
bubo_occ_dredge_vREF <- dredge(bubo_occ_full_vREF)
bubo_occ_dredge_delta_vREF <- get.models(bubo_occ_dredge_vREF, subset = delta <= 2)
summary(bubo_occ_dredge_delta_vREF[[1]]) 

# set metabarcoding as ref level
bubo.occu@obsCovs %>% mutate(METHOD2 = relevel(METHOD, ref = "MB_md")) %>% data.frame()->bubo.occu@obsCovs 
bubo_occ_full_mbREF <- occu(~METHOD2 + sDay + METHOD2*sDay + slog_Area ~1, bubo.occu)
bubo_occ_dredge_mbREF <- dredge(bubo_occ_full_mbREF)
bubo_occ_dredge_delta_mbREF <- get.models(bubo_occ_dredge_mbREF, subset = delta <= 2)
summary(bubo_occ_dredge_delta_mbREF[[1]]) 

# make predictions
bubo_pred_df <- expand.grid(METHOD = levels(bubo.occu@obsCovs$METHOD), 
                            sDay = seq(from = -1.6, to = 1.6, by = .2), slog_Area = seq(from = -2, to = 2, by = .2))
bubo_occ_pred <- predict(fitList(fits = bubo_occ_dredge_delta), # model averaged predictors
                         newdata = bubo_pred_df, 
                         type = "det", append = TRUE)

##### California newt ######
tato_occ_full <- occu(~METHOD + sDay + METHOD*sDay + slog_Area ~1,tato.occu)
tato_occ_dredge <- dredge(tato_occ_full)
tato_occ_dredge_delta <- get.models(tato_occ_dredge, subset = delta <= 2)
length(tato_occ_dredge_delta)
summary(tato_occ_dredge_delta[[1]])
confint(tato_occ_dredge_delta[[1]], type = "det")
summary(tato_occ_dredge_delta[[2]])
confint(tato_occ_dredge_delta[[2]], type = "det")

# now set  VES as ref level
tato.occu@obsCovs %>% mutate(METHOD2 = relevel(METHOD, ref = "VES")) %>% data.frame()->tato.occu@obsCovs 
tato_occ_full_vREF <- occu(~METHOD2 + sDay + METHOD2*sDay + slog_Area ~1, tato.occu)
tato_occ_dredge_vREF <- dredge(tato_occ_full_vREF)
tato_occ_dredge_delta_vREF <- get.models(tato_occ_dredge_vREF, subset = delta <= 2)
summary(tato_occ_dredge_delta_vREF[[1]])
confint(tato_occ_dredge_delta_vREF[[1]], type = "det")

## predict:
tato_pred_df <- expand.grid(METHOD = levels(tato.occu@obsCovs$METHOD), 
                            sDay = seq(from = -1.6, to = 1.6, by = .2), slog_Area = seq(from = -2, to = 2, by = .2))
tato_occ_pred <- predict(fitList(fits = tato_occ_dredge_delta), # model averaged predictors
                         newdata = tato_pred_df, 
                         type = "det", append = TRUE)

#######################################
####### PLOT DETECTION PROBS #######
#######################################
#### Figure 3 (detection prob; all species) #####
# bind predictions into a dataframe
rbind(
  radr_occ_pred %>% mutate(SPECIES = "California red-legged frog"),
  amca_occ_pred %>% mutate(SPECIES = "California tiger salamander"),
  bubo_occ_pred %>% mutate(SPECIES = "Western toad"),
  psre_occ_pred %>% mutate(SPECIES = "Pacific chorus frog"),
  raca_occ_pred %>% mutate(SPECIES = "American bullfrog"),
  tato_occ_pred %>% mutate(SPECIES = "California newt")) -> all.preds
all.preds %>% mutate(METHOD = fct_reorder(METHOD, Predicted)) -> all.preds

# we want to predict detection when date and area are at the average (scaled value = 0)
all.preds %>% filter(sDay == 0, slog_Area == 0) ->all.preds.0area0day 

# manually add in the group annotations to denote significance
all.preds.0area0day %>% mutate(Significance = case_when(
  SPECIES == "American bullfrog" & METHOD %in% c("Dip Net", "Seine") ~ "a",
  SPECIES == "American bullfrog" & METHOD %in% c("MB_md", "VES") ~ "b",
  
  SPECIES == "California newt" & METHOD %in% c("VES") ~ "a",
  SPECIES == "California newt" & METHOD %in% c("Dip Net", "Seine", "MB_md") ~ "b",
  
  SPECIES == "California red-legged frog" & METHOD %in% c("Dip Net") ~ "a",
  SPECIES == "California red-legged frog" & METHOD %in% c("VES", "Seine") ~ "a, b",
  SPECIES == "California red-legged frog" & METHOD %in% c("MB_md", "qPCR") ~ "b",
  
  SPECIES == "Pacific chorus frog" ~ "a",
  
  SPECIES == "California tiger salamander" & METHOD %in% c("Seine") ~ "a",
  SPECIES == "California tiger salamander" & METHOD %in% c("MB_md") ~ "a, b",
  SPECIES == "California tiger salamander" & METHOD %in% c("qPCR") ~ "b",
  
  SPECIES == "Western toad" & METHOD %in% c("Seine") ~ "a*",
  SPECIES == "Western toad" & METHOD %in% c("Dip Net", "MB_md") ~ "a, b*",
  SPECIES == "Western toad" & METHOD %in% c("VES") ~ "b*"
)) -> all.preds.0area0day

# mean predicted probability per method
all.preds.0area0day %>% group_by(METHOD) %>% summarise(mean_p = mean(Predicted))

# relevel for plotting
all.preds.0area0day<- all.preds.0area0day %>%
  mutate(METHOD = factor(METHOD, levels = c("Dip Net", "VES", "Seine", "MB_md", "qPCR")))

# plot all species detection 
all.preds.0area0day %>% ggplot(aes(x=METHOD, y = Predicted, group = SPECIES, col = SPECIES)) +
  geom_point(size = 3)+
  geom_text(aes(label = Significance, y = upper + .12), col = "black", size =4)+
  geom_errorbar(aes(ymin = lower, ymax = upper), width = .1, size = 1)+
  facet_wrap(~SPECIES, scales = "free_x", ncol = 2)+
  scale_x_discrete(labels = c("VES", "Dipnet", "Seine", "eDNA\n(MB)", "eDNA\n(qPCR)")) +
  xlab("") +
  ylab("")+
  ylim(0,1.15)+
  scale_y_continuous(breaks = c(0, .25, 0.5, 0.75, 1))+
  scale_color_manual(values = c("#698B22", "orange",  "red", "#104E8B", "#66CDAA","#8B3E2F")) +
  theme_bw()+
  theme(legend.position = "none", panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        strip.background = element_rect(color = "white", size = 3),
        strip.text = element_text(size = 12, vjust =1))

#### Figure 5 (detection prob; data) #####
# for unscaling date and undoing julian date:
# center = 166, scale = 22.8
# 135 : May 15
# 152: June 1
# 166: June 15
# 182: Jul 1
# 196: Jul 15
jDays <- c(135, 152, 166, 182, 196)
breaks <- (jDays - 166)/22.8 # 166 is mean sampling day and how week was scaled in models

# date matters for Pacific chorus frog, Western toad 
all.preds%>% filter(slog_Area ==0,SPECIES %in% c("Pacific chorus frog", "Western toad")) %>% mutate(
  METHOD = as.character(METHOD)) -> preds.date

# combine methods that don't significantly differ
# we don't have repeated dipnet data so we won't plot this method
preds.date %>% mutate(METHOD_POOL = case_when(
  SPECIES == "Western toad" & METHOD %in% c("Seine", "MB_md") ~ "Seine/MB",
  SPECIES == "Western toad" & METHOD == "VES" ~ "VES",
  SPECIES == "Pacific chorus frog" ~ "VES/Seine/eDNA(MB)")) -> preds.date

# group by identical methods and take mean prediction
preds.date %>% filter(!is.na(METHOD_POOL)) %>% group_by(METHOD_POOL, SPECIES, sDay) %>% summarise(Predicted = mean(Predicted), lower = mean(lower), upper = mean(upper)) -> preds.date 

preds.date %>%
  ggplot(aes(x=sDay, y = Predicted, group=METHOD_POOL, col = METHOD_POOL, fill = METHOD_POOL)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .2, col = NA)+
  facet_wrap(~SPECIES, nrow = 1) +
  scale_x_continuous(breaks = breaks, labels = c("15 May", "1 June", "15 June", "1 July", "15 July"))+
  xlab("Date") +
  ylab("Probability detection (p)") +
  geom_vline(xintercept = (148-166)/22.8) + # eDNA visit 1 line (julian date = 148)
  geom_vline(xintercept = (196-166)/22.8, linetype = "dashed")+ # eDNA visit 2
  scale_fill_manual(values = c("#2F4F4F", "#CD5C5C", "#2F4F4F", "#CD5C5C", "#2F4F4F", "#CD5C5C"))+
  scale_color_manual(values = c("#2F4F4F", "#CD5C5C", "#2F4F4F", "#CD5C5C", "#2F4F4F", "#CD5C5C"))+
  theme_bw()+
  theme(legend.position = "none", panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),strip.background = element_blank(),
        strip.text = element_text(size = 14, vjust =1))



