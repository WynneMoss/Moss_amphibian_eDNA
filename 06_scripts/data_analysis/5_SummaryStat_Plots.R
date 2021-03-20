#######################################
###### eDNA Amphibian Detection #######
#######################################

##### Create plots and examine detection histories
# using raw detection data
# also incorporate amphibian density data from net sweeps
#######################################
###### LOAD LIBRARIES AND DATA ######
#######################################
library(colorspace)
library(tidyverse)
library(glmmTMB)
# detection data: each row a site-visit-method combination
# covariates included

detections <- read.csv("05_data/Detection_all_methods.csv")

# data frame for species and colors
species <- data.frame(SPP = c("AMCA", "BUBO", "PSRE", "RACA", "RADR", "TATO"),
                      SPP_NAME = c("California tiger salamander", "Western toad", 
                                   "Pacific chorus frog", "American bullfrog", "California red-legged frog",
                                   "California newt"),
                      color = c("#104E8B", "#8B3E2F", "#66CDAA", "#698B22", "red", "orange"))

# amphibian density data (see Appendix)
net_densities <- read.csv("05_data/Raw/Net_sweep_larvae.csv")


###############################################
##### SUMMARY STATISTICS  #####################
###############################################
#### Naive occupancy #####
# based on pooled detections
detections %>% select(SITE_CODE, METHOD, AMCA:TATO) %>% pivot_longer(AMCA:TATO, names_to = "SPP_CODE", values_to  = "DETECTION") %>% group_by(SITE_CODE, SPP_CODE) %>% mutate(TRUE_OCCUPANCY = max(DETECTION, na.rm=TRUE)) -> true.occupancy

# summarise detections across species/methods
true.occupancy %>% group_by(SITE_CODE, SPP_CODE, METHOD) %>%  summarise(DETECTED = max(DETECTION)) %>% 
  group_by(SPP_CODE, METHOD) %>% summarise(N_SITES_DETECTED = sum(DETECTED, na.rm=TRUE)) %>% 
  left_join(true.occupancy %>% distinct(SITE_CODE, SPP_CODE, TRUE_OCCUPANCY) %>% 
              group_by(SPP_CODE) %>% summarise(N_SITES_OCCUPIED = sum(TRUE_OCCUPANCY))) -> true.occupancy

###############################################
##### LARVAL DENSITY #####################
###############################################
density_model <- glmmTMB(TOTAL_COUNT ~ scale(jDay)*SPECIES +I(scale(jDay)^2) + offset(log(TOTAL_DIST)) + (1|SITE_CODE), data = net_densities, family = "poisson")
summary(density_model)

# predict density for each species and across time
predict_densities <- expand.grid(jDay = seq(from = min(net_densities$jDay), to = max(net_densities$jDay), by = .1),
                                 TOTAL_DIST = 1, SPECIES = unique(net_densities$SPECIES), SITE_CODE="fake")
predict_densities$predict <- predict(density_model,newdata = predict_densities,allow.new.levels = TRUE) 
predict_densities$se <- predict(density_model,newdata = predict_densities,allow.new.levels = TRUE, se.fit = TRUE,
                                re.form = NA) %>% pluck("se.fit") 
###############################################
##### PLOT DATA  ##############################
###############################################
##### Figure 4: heat map, detection by site #####
# first group field and eDNA methods together
# get a dataframe for heat map
detections %>% mutate(method_broad = case_when(METHOD %in% c("VES", "Seine", "Dip Net") ~ "Field", TRUE ~"eDNA")) %>% 
  group_by(SITE_CODE, DATE, method_broad) %>% 
  summarise(AMCA = max(AMCA, na.rm=TRUE), BUBO = max(BUBO, na.rm=TRUE), PSRE = max(PSRE, na.rm=TRUE), RACA = max(RACA, na.rm=TRUE),
            RADR = max(RADR, na.rm=TRUE), TATO = max(TATO, na.rm=TRUE)) %>% 
  group_by(SITE_CODE, method_broad) %>% mutate(Visit = rank(DATE)) %>% pivot_longer(AMCA:TATO, names_to = "SPECIES", values_to = "DETECTION") %>% 
  mutate(METHOD_SPP = paste(SPECIES, method_broad, sep = "_")) %>% 
  pivot_wider(id_cols = c(SITE_CODE, Visit), names_from = METHOD_SPP, values_from = DETECTION, names_sort = TRUE ) %>% 
  pivot_longer(AMCA_eDNA:TATO_Field, names_to = "SPP_SURVEY", values_to = "DETECTION") %>% 
  separate(SPP_SURVEY, into = c("SPP", "SURVEY"), remove = FALSE, sep = "_") %>% 
  left_join(species) %>%
  mutate(color_plot = case_when(DETECTION == 0 ~ lighten(color, amount = .5), TRUE ~ color)) -> survey.level

species <- arrange(species, SPP_NAME)

survey.level %>% arrange(desc(SPP_NAME))%>% mutate(SPP_SURVEY =factor(SPP_SURVEY, levels = unique(SPP_SURVEY))) %>% mutate(VISIT = paste("Visit", Visit)) -> survey.level

(survey.level %>% 
  ggplot(aes(x=SITE_CODE, y =SPP_SURVEY, fill = SPP_NAME, alpha = DETECTION)) +
  geom_tile() +
  scale_y_discrete(labels =rep(c("eDNA", "Field"), 6), expand = c(0,0), breaks = levels(survey.level$SPP_SURVEY)) +
  scale_x_discrete(expand = c(0,0)) +
  facet_wrap(~VISIT) +
  theme_classic()+
  scale_alpha_continuous(range = c(0.4,1))+
  scale_fill_manual(values = species$color, name = "Species")+
  xlab("Site") +
  ylab("Survey Method") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 7, hjust = 1), 
        legend.position = "none",
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title =element_text(size = 12),
        panel.spacing = unit(2, "lines"),
        strip.text = element_text(size = 12),
        panel.border = element_rect(fill = NA),
        strip.background = element_blank()) +
  guides(alpha = FALSE) -> detection.histories)



##### Figure 2: overall detection summary #####
# manually add "true" occupancy as another method (pooled)
# group detections among field and eDNA based methods
detections %>% 
  select(SITE_CODE, METHOD, AMCA:TATO) %>% 
  pivot_longer(AMCA:TATO, names_to = "SPP_CODE", values_to  = "DETECTION") %>%
  mutate(METHOD_POOLED = case_when(METHOD %in% c("Dip Net", "Seine", "VES") ~ "Field", 
                                   METHOD == "MB_md" ~ "Metabarcoding",
                                   TRUE~METHOD)) %>% 
  group_by(SITE_CODE, METHOD_POOLED, SPP_CODE) %>% summarise(DETECTION = max(DETECTION, na.rm=TRUE)) %>% ungroup() %>% 
  group_by(SITE_CODE, SPP_CODE) %>% mutate(TRUE_OCCUPANCY = max(DETECTION, na.rm=TRUE)) %>% ungroup() %>% 
  pivot_wider(names_from = METHOD_POOLED, values_from = DETECTION) %>% 
  pivot_longer(cols = TRUE_OCCUPANCY:qPCR, values_to = "DETECTION", names_to = "METHOD_POOLED") %>% 
  mutate(METHOD_POOLED = str_replace(METHOD_POOLED, "TRUE_OCCUPANCY", "Pooled")) %>% 
  group_by(SPP_CODE, METHOD_POOLED) %>% summarise(N_SITES = sum(DETECTION, na.rm=TRUE)) %>% ungroup() %>% 
  group_by(SPP_CODE) %>% mutate(N_SITES_OCCUPIED = max(N_SITES, na.rm = TRUE)) %>% ungroup() %>% 
  mutate(r_prop_detected = N_SITES/N_SITES_OCCUPIED) %>% 
  mutate(Label = case_when(is.finite(N_SITES) ~paste(N_SITES, "/", 20, sep = ""))) %>% 
  mutate(Method = fct_relevel(METHOD_POOLED, c("Pooled", "Field", "Metabarcoding", "qPCR"))) %>% 
  left_join(species, by = c("SPP_CODE" = "SPP")) %>% 
  mutate(SPP_NAME = fct_relevel(SPP_NAME, 
                                c("California red-legged frog", "American bullfrog", "California tiger salamander",
                                  "California newt", "Western toad", "Pacific chorus frog"))) %>% 
  ggplot() +
  geom_tile(aes(x=Method, y = SPP_NAME, fill = r_prop_detected), color = "white") +
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0))+
  scale_fill_gradient(low = "#7AC5CD", high = "#8B668B", na.value = "white")+
  geom_text(aes(x=Method, y = SPP_NAME, label = Label), size =6) +
  annotate(geom = "text", x=-0.75, y = 7, label = "More common", hjust = 0, size = 6)+
  annotate(geom = "text", x=-0.75, y = 0, label = "Less common", hjust = 0, size = 6)+
  theme_bw()+
  ylab("")+
  labs(fill = expression(paste("Relative detection:   ", over("# detections with method i", "# pooled detections"))))+
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 16),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.position = "top") +
  coord_cartesian(clip = "off", xlim = c(0.5,4.5), ylim = c(0.5,6.5))








##### Figure 5: densities across time ######
colors <- c("#66CDAA", "#8B3E2F", "#104E8B", "#698B22", "red", "orange")

predict_densities %>% mutate(SPECIES = factor(SPECIES, levels = c("PSRE", "BUBO", "AMCA", "RACA", "RADR", "TATO"))) %>% 
  ggplot(aes(x=jDay, y = exp(predict), group = SPECIES)) +
  geom_ribbon(aes(ymin =exp(predict-se), ymax = exp(predict+se), fill = SPECIES), alpha = .5)+
  geom_line(aes(col = SPECIES), size = 1.5) +
  geom_vline(xintercept = 148) +
  geom_vline(xintercept = 196, linetype = "dashed")+
  scale_color_manual(values = colors, name = "",
                     labels = c("Pacific chorus frog", "Western toad","California tiger salamander", 
                                "American bullfrog", "California red-legged frog", "California newt"))+  
  scale_fill_manual(values = colors, name = "",
                    labels = c("Pacific chorus frog", "Western toad","California tiger salamander", 
                               "American bullfrog", "California red-legged frog", "California newt")) +
  scale_y_sqrt(breaks = c(0.01,.1,.5)) +
  xlab("Date") +
  ylab("Larval Density") +
  scale_x_continuous(breaks = c(135, 152, 166, 182, 196), labels = c("15 May", "1 June", "15 June", "1 July", "15 July"))+
  theme_bw() +
  guides(fill=guide_legend(ncol=1), col = guide_legend(ncol = 1))+
  theme(legend.position = "right",  legend.text = element_text(size = 12), legend.title = element_text(size = 14),
        axis.text = element_text(size = 12), axis.title = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
