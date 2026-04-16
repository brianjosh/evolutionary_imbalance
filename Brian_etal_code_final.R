#Code for 'Plants that evolved under high phylogenetic diversity have higher invasion success, particularly in undisturbed communities'
#Joshua Brian, Mark van Kleunen, Wayne Dawson, Anne Kempel, Weihan Zhao, Jane Catford

#Code prepared by Joshua Brian; joshua.brian@kcl.ac.uk OR jshbrian@gmail.com

#This code was run using R v.4.1.1 and the following packages:
library(tidyverse) #version 1.3.2
library(glmmTMB) #version 1.1.4
library(DHARMa) #version 0.4.6
library(car) #version 3.1-1
library(MuMIn) #version 1.46.0
library(bbmle) #version 1.0.25
library(visreg) #version 2.7.0
library(lme4) #version 1.1-31
library(performance) #version 0.12.3
library(interactions) #version 1.2.0
library(viridis) #version 0.6.2
library(ggpubr) #version 0.4.0
library(cowplot) #version 1.1.1
library(rtry) #version 1.1.0
library(piecewiseSEM) #version 2.1.2
library(arm) #version 1.13-1
library(sf) #version 1.0-12
library(patchwork) #version 1.3.2.9000

#This code provides all details required to perform the statistical analyses and obtain the 
#results described in the text. It also provides code to produce all main text figures 
#and supplementary results figures (Figs 1-6 and S9-S17).

#Models for all four phylogenetic diversity metrics are included (mean alpha PD, max alpha PD, gamma PD and 
#median alpha PD). Results for mean alpha PD correspond to those reported and 
#discussed in the main text, results for the other three metrics are reported in the
#Supporting Information (Tables S2-S9 and Figs S9-S17).

#In all cases, models described below are the final selected models 
#(e.g. only interactions found to be significant were maintained in models).

#Main analyses for all three research questions are presented first. Supplementary analyses
#as described in the Supporting Information then follow. See the Supporting Information file (.pdf) for 
#full details. Note the 'Species origin' and 'Hnd' supplementary analyses are included in Main text
#Question 1 and 2 analyses.

#Please note this code contains lots of robustness tests i.e. the same tests carried out/same models run on 
#different subsets of the data. As the statistical analysis is extensive, this means that sometimes the same
#model names and dataset names are recycled and reassigned through the code. Each section in the contents should
#be able to be run alone (i.e. one section does not rely on another section, all datasets are loaded fresh and all 
#models constructed from scratch), but if possible the code should be run in order, and running only latter parts
#of a given section should not be attempted. 

##################################################################################################

#CONTENTS

##################################################################################################

#Construction of Methods figure: Line 67 
#Main text Question 1: Line 120
#Main text Question 2: Line 1632
#Main text Question 3: Line 4144
#Supplementary analysis K+M models: Line 4409
#Supplementary analysis K+H models: Line 5084
#Supplementary analysis SC models: Line 5867
#Supplementary analysis Structural Equation Models: Line 7438

##################################################################################################

#METHODS FIGURE: Mean Alpha PD and number of experimental species by region

##################################################################################################

#Read in the shapefile of regions
sf_tdwg3_ea <- st_read("level3.shp") |> #note the .shx, .dbf and .prj files also need to be in the working directory
  st_make_valid() |>
  st_wrap_dateline() |>
  st_transform("+proj=eqearth")
plot(sf_tdwg3_ea["geometry"])
st_crs(sf_tdwg3_ea)

#Read in phylogenetic diversity by region
PDbyregion <- read.csv("TDWG3_PD.csv", header=T)

sf_tdwg3_ea <- left_join(sf_tdwg3_ea, PDbyregion, by="LEVEL3_COD")

#Plot the Mean Alpha PD for each region
geomPD <- ggplot(sf_tdwg3_ea, aes(fill=AlphaPhy)) + 
  geom_sf() + 
  coord_sf(
    xlim = c(-14.5e6, 16e6),   
    ylim = c(-6.8e6,  8.5e6),    
    expand = FALSE
  ) +
  theme_bw() + 
  scale_fill_gradient2(name = "Mean Alpha \nPhylogenetic \nDiversity", midpoint=0, low="blue", mid="white", high="red") 
geomPD

#Plot the number of species in the experiments native to each region
geomcount <- ggplot(sf_tdwg3_ea, aes(fill=Species)) + 
  geom_sf() + 
  coord_sf(
    xlim = c(-14.5e6, 16e6),   
    ylim = c(-6.8e6,  8.5e6),    
    expand = FALSE
  ) +
  theme_bw() + 
  scale_fill_gradient(name = "Experimental \nspecies", low="white", high="#1bb260ff") 
geomcount

geomfull <- ggarrange(geomPD, geomcount, nrow=2, labels = c("(c)", "(d)"))
geomfull

ggsave(geomfull, 
       filename = "AlphaPhymap.svg",
       height = 280, width = 180, units = "mm")

#This was then joined with the vector artwork showing how PD was calculated
#(Fig. 1 panels (a) and (b)) using the program Inkscape

##################################################################################################

#QUESTION ONE: Do sown species with higher home-range-PD have higher rates of colonisation and 
#survival than species with low home-range-PD?

##################################################################################################

#Read in all four data files 

germination <- readRDS("germinationrevised.RDS") %>% #Referred to as "Colonisation" in main text
  mutate(suitabilitypc1V2abs = abs(suitabilitypc1V2)) %>%
  mutate(suitabilitypc2V2abs = abs(suitabilitypc2V2))
str(germination)

surviveyr1 <- readRDS("surviveyr1revised.RDS") %>% #Referred to as "First year" in main text
  mutate(suitabilitypc1V2abs = abs(suitabilitypc1V2)) %>%
  mutate(suitabilitypc2V2abs = abs(suitabilitypc2V2))
str(surviveyr1)

surviveyr2early <- readRDS("surviveyr2earlyrevised.RDS") %>% #Referred to as "Overwinter" in main text
  mutate(suitabilitypc1V2abs = abs(suitabilitypc1V2)) %>%
  mutate(suitabilitypc2V2abs = abs(suitabilitypc2V2))
str(surviveyr2early)

surviveyr2late <- readRDS("surviveyr2laterevised.RDS") %>% #Referred to as "Second year" in main text
  mutate(suitabilitypc1V2abs = abs(suitabilitypc1V2)) %>%
  mutate(suitabilitypc2V2abs = abs(suitabilitypc2V2))
str(surviveyr2late)

################################################################################

#DO SURVIVAL Y/N FOR ALL DATASETS

#For each model, two versions were run: 
# -version with standardised continuous variables (rescale() function following Gelman 2008)
# -version with unstandardised continuous variables (allows for direct interpretation of phylogenetic diversity measure effect size)
#For this reason, additional tests where the strength of effect size of phylogenetic diversity is compared when including biogeographic
#origin or native/non-native status, the unstandardised model is used, but statistically it is the same outcome.

#COLONISATION YES NO

germination$Disturbance <- relevel(germination$Disturbance, ref="0") #to ensure that 'Undisturbed' is the reference level

#Mean alpha diversity 

germ_YN_Mean_AlphaPhy <- glmmTMB(germ_YN ~ rescale(Mean_AlphaPhy) + Disturbance + rescale(propagule_pressure) +
                                   rescale(suitabilitypc1V2abs) + rescale(suitabilitypc2V2abs) +
                                   (1|Family) + (1|POWO.name) +
                                   (1|Site) + offset(log(Density)), 
                                 family=binomial(link="logit"), 
                                 data=germination)

germ_YN_Mean_AlphaPhy <- glmmTMB(germ_YN ~ Mean_AlphaPhy + Disturbance + propagule_pressure +
                                   suitabilitypc1V2abs + suitabilitypc2V2abs +
                                   (1|Family) + (1|POWO.name) +
                                   (1|Site) + offset(log(Density)), 
                                 family=binomial(link="logit"), 
                                 data=germination)
summary(germ_YN_Mean_AlphaPhy)
Anova(germ_YN_Mean_AlphaPhy, type="III")
confint(germ_YN_Mean_AlphaPhy)
r.squaredGLMM(germ_YN_Mean_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = germ_YN_Mean_AlphaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(germ_YN_Mean_AlphaPhy)
check_singularity(germ_YN_Mean_AlphaPhy)

#Testing whether worth adding native status - all results in Table S10

germ_YN_Mean_AlphaPhyN <- glmmTMB(germ_YN ~ Mean_AlphaPhy + Disturbance + propagule_pressure +
                                    suitabilitypc1V2abs + suitabilitypc2V2abs + Status +
                                    (1|Family) + (1|POWO.name) +
                                    (1|Site) + offset(log(Density)), 
                                  family=binomial(link="logit"), 
                                  data=germination)

AICtab(germ_YN_Mean_AlphaPhyN, germ_YN_Mean_AlphaPhy)

#Testing whether worth adding biogeographic origin - all results in Table S10

BG <- read.csv("biogeographic_ID.csv", header=T) %>% dplyr::select(plant_name_id, BG_ID)
BG$BG_ID <- as.factor(BG$BG_ID)

germinationBG <- left_join(germination, BG, by="plant_name_id")

germ_YN_Mean_AlphaPhyBG <- glmmTMB(germ_YN ~ Mean_AlphaPhy + Disturbance + propagule_pressure +
                                     suitabilitypc1V2abs + suitabilitypc2V2abs + 
                                     (1|Family) + (1|POWO.name) + (1|BG_ID) +
                                     (1|Site) + offset(log(Density)), 
                                   family=binomial(link="logit"), 
                                   data=germinationBG)
summary(germ_YN_Mean_AlphaPhyBG)

AICtab(germ_YN_Mean_AlphaPhyBG, germ_YN_Mean_AlphaPhy)

#Maximum alpha diversity

germ_YN_Max_AlphaPhy <- glmmTMB(germ_YN ~ rescale(Max_AlphaPhy) + Disturbance + rescale(propagule_pressure) +
                                  rescale(suitabilitypc1V2abs) + rescale(suitabilitypc2V2abs) +
                                  (1|Family) + (1|POWO.name) +
                                  (1|Site) + offset(log(Density)), 
                                family=binomial(link="logit"), 
                                data=germination)

germ_YN_Max_AlphaPhy <- glmmTMB(germ_YN ~ Max_AlphaPhy + Disturbance + propagule_pressure +
                                  suitabilitypc1V2abs + suitabilitypc2V2abs +
                                  (1|Family) + (1|POWO.name) +
                                  (1|Site) + offset(log(Density)), 
                                family=binomial(link="logit"), 
                                data=germination)
summary(germ_YN_Max_AlphaPhy)
Anova(germ_YN_Max_AlphaPhy, type="III")
confint(germ_YN_Max_AlphaPhy)
r.squaredGLMM(germ_YN_Max_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = germ_YN_Max_AlphaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(germ_YN_Max_AlphaPhy)
check_singularity(germ_YN_Max_AlphaPhy)

#Gamma diversity

germ_YN_GammaPhy <- glmmTMB(germ_YN ~ rescale(log(GammaPhy)) + Disturbance + rescale(propagule_pressure) +
                              rescale(suitabilitypc1V2abs) + rescale(suitabilitypc2V2abs) +
                              (1|Family) + (1|POWO.name) +
                              (1|Site) + offset(log(Density)), 
                            family=binomial(link="logit"), 
                            data=germination)

germ_YN_GammaPhy <- glmmTMB(germ_YN ~ log(GammaPhy) + Disturbance + propagule_pressure +
                              suitabilitypc1V2abs + suitabilitypc2V2abs +
                              (1|Family) + (1|POWO.name) +
                              (1|Site) + offset(log(Density)), 
                            family=binomial(link="logit"), 
                            data=germination)
summary(germ_YN_GammaPhy)
Anova(germ_YN_GammaPhy, type="III")
confint(germ_YN_GammaPhy)
r.squaredGLMM(germ_YN_GammaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = germ_YN_GammaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(germ_YN_GammaPhy)
check_singularity(germ_YN_GammaPhy)

#Median alpha diversity 

germ_YN_MedPhy <- glmmTMB(germ_YN ~ rescale(Med_AlphaPhy)*Disturbance + rescale(propagule_pressure) +
                            rescale(suitabilitypc1V2abs) + rescale(suitabilitypc2V2abs) +
                            (1|Family) + (1|POWO.name) +
                            (1|Site) + offset(log(Density)), 
                          family=binomial(link="logit"), 
                          data=germination)

germ_YN_MedPhy <- glmmTMB(germ_YN ~ Med_AlphaPhy*Disturbance + propagule_pressure +
                            suitabilitypc1V2abs + suitabilitypc2V2abs +
                            (1|Family) + (1|POWO.name) +
                            (1|Site) + offset(log(Density)), 
                          family=binomial(link="logit"), 
                          data=germination)
summary(germ_YN_MedPhy)
Anova(germ_YN_MedPhy, type="III")
confint(germ_YN_MedPhy)
r.squaredGLMM(germ_YN_MedPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = germ_YN_MedPhy, plot=F)
plot(simulationOutput)
check_overdispersion(germ_YN_MedPhy)
check_singularity(germ_YN_MedPhy)

#FIRST YEAR SURVIVAL YES NO

surviveyr1$Disturbance <- relevel(surviveyr1$Disturbance, ref="0") 

#Mean alpha diversity 

surviveyr1_YN_Mean_AlphaPhy <- glmmTMB(surviveyr1_YN ~ rescale(Mean_AlphaPhy) + Disturbance + rescale(propagule_pressure) +
                                         rescale(suitabilitypc1V2abs) + rescale(suitabilitypc2V2abs) + Herbivory + 
                                         (1|Family) + (1|POWO.name) + 
                                         (1|Site) + offset(log(Density)), 
                                       family=binomial(link="logit"), 
                                       data=surviveyr1)

surviveyr1_YN_Mean_AlphaPhy <- glmmTMB(surviveyr1_YN ~ Mean_AlphaPhy + Disturbance + propagule_pressure +
                                         suitabilitypc1V2abs + suitabilitypc2V2abs + Herbivory + 
                                         (1|Family) + (1|POWO.name) + 
                                         (1|Site) + offset(log(Density)), 
                                       family=binomial(link="logit"), 
                                       data=surviveyr1)
summary(surviveyr1_YN_Mean_AlphaPhy)
Anova(surviveyr1_YN_Mean_AlphaPhy, type="III")
confint(surviveyr1_YN_Mean_AlphaPhy)
r.squaredGLMM(surviveyr1_YN_Mean_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr1_YN_Mean_AlphaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr1_YN_Mean_AlphaPhy)
check_singularity(surviveyr1_YN_Mean_AlphaPhy)

#### Testing whether native explains anything else

surviveyr1_YN_Mean_AlphaPhyN <- glmmTMB(surviveyr1_YN ~ Mean_AlphaPhy + Disturbance + propagule_pressure +
                                          suitabilitypc1V2abs + suitabilitypc2V2abs + Herbivory + Status +
                                          (1|Family) + (1|POWO.name) + 
                                          (1|Site) + offset(log(Density)), 
                                        family=binomial(link="logit"), 
                                        data=surviveyr1)

AICtab(surviveyr1_YN_Mean_AlphaPhy, surviveyr1_YN_Mean_AlphaPhyN)

#Testing whether worth adding biogeographic origin

surviveyr1BG <- left_join(surviveyr1, BG, by="plant_name_id")

surviveyr1_YN_Mean_AlphaPhyBG <- glmmTMB(surviveyr1_YN ~ Mean_AlphaPhy + Disturbance + propagule_pressure +
                                           suitabilitypc1V2abs + suitabilitypc2V2abs + Herbivory + 
                                           (1|Family) + (1|POWO.name) + (1|BG_ID) +
                                           (1|Site) + offset(log(Density)), 
                                         family=binomial(link="logit"), 
                                         data=surviveyr1BG)
summary(surviveyr1_YN_Mean_AlphaPhyBG)

AICtab(surviveyr1_YN_Mean_AlphaPhy, surviveyr1_YN_Mean_AlphaPhyBG)

#Maximum alpha diversity

surviveyr1_YN_Max_AlphaPhy <- glmmTMB(surviveyr1_YN ~ rescale(Max_AlphaPhy)*Disturbance + rescale(propagule_pressure) +
                                        rescale(suitabilitypc1V2abs) + rescale(suitabilitypc2V2abs) + Herbivory +
                                        (1|Family) + (1|POWO.name) + 
                                        (1|Site) + offset(log(Density)), 
                                      family=binomial(link="logit"), 
                                      data=surviveyr1)

surviveyr1_YN_Max_AlphaPhy <- glmmTMB(surviveyr1_YN ~ Max_AlphaPhy*Disturbance + propagule_pressure +
                                        suitabilitypc1V2abs + suitabilitypc2V2abs + Herbivory +
                                        (1|Family) + (1|POWO.name) + 
                                        (1|Site) + offset(log(Density)), 
                                      family=binomial(link="logit"), 
                                      data=surviveyr1)
summary(surviveyr1_YN_Max_AlphaPhy)
Anova(surviveyr1_YN_Max_AlphaPhy, type="III")
confint(surviveyr1_YN_Max_AlphaPhy)
r.squaredGLMM(surviveyr1_YN_Max_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr1_YN_Max_AlphaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr1_YN_Max_AlphaPhy)
check_singularity(surviveyr1_YN_Max_AlphaPhy)

#Gamma diversity 

surviveyr1_YN_GammaPhy <- glmmTMB(surviveyr1_YN ~ rescale(log(GammaPhy))*rescale(suitabilitypc1V2abs) + Disturbance + rescale(propagule_pressure) +
                                    rescale(suitabilitypc2V2abs) + Herbivory +
                                    (1|Family) + (1|POWO.name) + 
                                    (1|Site) + offset(log(Density)), 
                                  family=binomial(link="logit"), 
                                  data=surviveyr1)

surviveyr1_YN_GammaPhy <- glmmTMB(surviveyr1_YN ~ log(GammaPhy)*suitabilitypc1V2abs + Disturbance + propagule_pressure +
                                    suitabilitypc2V2abs + Herbivory +
                                    (1|Family) + (1|POWO.name) + 
                                    (1|Site) + offset(log(Density)), 
                                  family=binomial(link="logit"), 
                                  data=surviveyr1)
summary(surviveyr1_YN_GammaPhy)
Anova(surviveyr1_YN_GammaPhy, type="III")
confint(surviveyr1_YN_GammaPhy)
r.squaredGLMM(surviveyr1_YN_GammaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr1_YN_GammaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr1_YN_GammaPhy)
check_singularity(surviveyr1_YN_GammaPhy)

#Median alpha diversity 

surviveyr1_YN_MedPhy <- glmmTMB(surviveyr1_YN ~ recsale(Med_AlphaPhy) + Disturbance + recsale(propagule_pressure) +
                                  recsale(suitabilitypc1V2abs) + recsale(suitabilitypc2V2abs) + Herbivory +
                                  (1|Family) + (1|POWO.name) + 
                                  (1|Site) + offset(log(Density)), 
                                family=binomial(link="logit"), 
                                data=surviveyr1)

surviveyr1_YN_MedPhy <- glmmTMB(surviveyr1_YN ~ Med_AlphaPhy + Disturbance + propagule_pressure +
                                  suitabilitypc1V2abs + suitabilitypc2V2abs + Herbivory +
                                  (1|Family) + (1|POWO.name) + 
                                  (1|Site) + offset(log(Density)), 
                                family=binomial(link="logit"), 
                                data=surviveyr1)
summary(surviveyr1_YN_MedPhy)
Anova(surviveyr1_YN_MedPhy, type="III")
confint(surviveyr1_YN_MedPhy)
r.squaredGLMM(surviveyr1_YN_MedPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr1_YN_MedPhy, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr1_YN_MedPhy)
check_singularity(surviveyr1_YN_MedPhy)

#SECOND YEAR EARLY SURVIVAL YES NO

surviveyr2early$Disturbance <- relevel(surviveyr2early$Disturbance, ref="0") 

#Mean alpha diversity 

surviveyr2early_YN_Mean_AlphaPhy <- glmmTMB(surviveyr2early_YN ~ rescale(Mean_AlphaPhy) + Disturbance + rescale(propagule_pressure) +
                                              rescale(suitabilitypc1V2abs) + rescale(suitabilitypc2V2abs) + Herbivory +
                                              (1|Family) + (1|POWO.name) +  
                                              (1|Site) + offset(log(Density)), 
                                            family=binomial(link="logit"), 
                                            data=surviveyr2early)

surviveyr2early_YN_Mean_AlphaPhy <- glmmTMB(surviveyr2early_YN ~ Mean_AlphaPhy + Disturbance + propagule_pressure +
                                              suitabilitypc1V2abs + suitabilitypc2V2abs + Herbivory +
                                              (1|Family) + (1|POWO.name) +  
                                              (1|Site) + offset(log(Density)), 
                                            family=binomial(link="logit"), 
                                            data=surviveyr2early)
summary(surviveyr2early_YN_Mean_AlphaPhy)
Anova(surviveyr2early_YN_Mean_AlphaPhy, type="III")
confint(surviveyr2early_YN_Mean_AlphaPhy)
r.squaredGLMM(surviveyr2early_YN_Mean_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2early_YN_Mean_AlphaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr2early_YN_Mean_AlphaPhy)
check_singularity(surviveyr2early_YN_Mean_AlphaPhy)

#Testing to see if native explains anything else

surviveyr2early_YN_Mean_AlphaPhyN <- glmmTMB(surviveyr2early_YN ~ Mean_AlphaPhy + Disturbance + propagule_pressure +
                                               suitabilitypc1V2abs + suitabilitypc2V2abs + Herbivory + Status +
                                               (1|Family) + (1|POWO.name) +  
                                               (1|Site) + offset(log(Density)), 
                                             family=binomial(link="logit"), 
                                             data=surviveyr2early)
summary(surviveyr2early_YN_Mean_AlphaPhyN)

AICtab(surviveyr2early_YN_Mean_AlphaPhyN, surviveyr2early_YN_Mean_AlphaPhy)

#Testing whether worth adding biogeographic origin

surviveyr2earlyBG <- left_join(surviveyr2early, BG, by="plant_name_id")

surviveyr2early_YN_Mean_AlphaPhyBG <- glmmTMB(surviveyr2early_YN ~ Mean_AlphaPhy + Disturbance + propagule_pressure +
                                                suitabilitypc1V2abs + suitabilitypc2V2abs + Herbivory +
                                                (1|Family) + (1|POWO.name) + (1|BG_ID) +  
                                                (1|Site) + offset(log(Density)), 
                                              family=binomial(link="logit"), 
                                              data=surviveyr2earlyBG)
summary(surviveyr2early_YN_Mean_AlphaPhyBG)

AICtab(surviveyr2early_YN_Mean_AlphaPhyBG, surviveyr2early_YN_Mean_AlphaPhy)

# Maximum alpha diversity

surviveyr2early_YN_Max_AlphaPhy <- glmmTMB(surviveyr2early_YN ~ rescale(Max_AlphaPhy)*rescale(suitabilitypc2V2abs) +  Disturbance + rescale(propagule_pressure) +
                                             rescale(suitabilitypc1V2abs) + Herbivory +
                                             (1|Family) + (1|POWO.name) + 
                                             (1|Site) + offset(log(Density)), 
                                           family=binomial(link="logit"), 
                                           data=surviveyr2early)

surviveyr2early_YN_Max_AlphaPhy <- glmmTMB(surviveyr2early_YN ~ Max_AlphaPhy*suitabilitypc2V2abs +  Disturbance + propagule_pressure +
                                             suitabilitypc1V2abs + Herbivory +
                                             (1|Family) + (1|POWO.name) + 
                                             (1|Site) + offset(log(Density)), 
                                           family=binomial(link="logit"), 
                                           data=surviveyr2early)
summary(surviveyr2early_YN_Max_AlphaPhy)
Anova(surviveyr2early_YN_Max_AlphaPhy, type="III")
confint(surviveyr2early_YN_Max_AlphaPhy)
r.squaredGLMM(surviveyr2early_YN_Max_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2early_YN_Max_AlphaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr2early_YN_Max_AlphaPhy)
check_singularity(surviveyr2early_YN_Max_AlphaPhy)

#Gamma diversity 

surviveyr2early_YN_GammaPhy <- glmmTMB(surviveyr2early_YN ~ rescale(log(GammaPhy)) + Disturbance + rescale(propagule_pressure) +
                                         rescale(suitabilitypc1V2abs) + rescale(suitabilitypc2V2abs) + Herbivory +
                                         (1|Family) + (1|POWO.name) + 
                                         (1|Site) + offset(log(Density)), 
                                       family=binomial(link="logit"), 
                                       data=surviveyr2early)

surviveyr2early_YN_GammaPhy <- glmmTMB(surviveyr2early_YN ~ log(GammaPhy) + Disturbance + propagule_pressure +
                                         suitabilitypc1V2abs + suitabilitypc2V2abs + Herbivory +
                                         (1|Family) + (1|POWO.name) + 
                                         (1|Site) + offset(log(Density)), 
                                       family=binomial(link="logit"), 
                                       data=surviveyr2early)
summary(surviveyr2early_YN_GammaPhy)
Anova(surviveyr2early_YN_GammaPhy, type="III")
confint(surviveyr2early_YN_GammaPhy)
r.squaredGLMM(surviveyr2early_YN_GammaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2early_YN_GammaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr2early_YN_GammaPhy)
check_singularity(surviveyr2early_YN_GammaPhy)

#Median alpha diversity 

surviveyr2early_YN_MedPhy <- glmmTMB(surviveyr2early_YN ~ rescale(Med_AlphaPhy)*rescale(suitabilitypc1V2abs) + Disturbance + rescale(propagule_pressure) +
                                       rescale(suitabilitypc2V2abs) + Herbivory +
                                       (1|Family) + (1|POWO.name) + 
                                       (1|Site) + offset(log(Density)), 
                                     family=binomial(link="logit"), 
                                     data=surviveyr2early)

surviveyr2early_YN_MedPhy <- glmmTMB(surviveyr2early_YN ~ Med_AlphaPhy*suitabilitypc1V2abs + Disturbance + propagule_pressure +
                                       suitabilitypc2V2abs + Herbivory +
                                       (1|Family) + (1|POWO.name) + 
                                       (1|Site) + offset(log(Density)), 
                                     family=binomial(link="logit"), 
                                     data=surviveyr2early)
summary(surviveyr2early_YN_MedPhy)
Anova(surviveyr2early_YN_MedPhy, type="III")
confint(surviveyr2early_YN_MedPhy)
r.squaredGLMM(surviveyr2early_YN_MedPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2early_YN_MedPhy, plot=F)
plot(simulationOutput)
testOutliers(simulationOutput)
check_overdispersion(surviveyr2early_YN_MedPhy)
check_singularity(surviveyr2early_YN_MedPhy)

#SECOND YEAR LATE SURVIVAL YES NO

surviveyr2late$Disturbance <- relevel(surviveyr2late$Disturbance, ref="0")

#Mean alpha diversity 

surviveyr2late_YN_Mean_AlphaPhy <- glmmTMB(surviveyr2late_YN ~ rescale(Mean_AlphaPhy)*rescale(suitabilitypc1V2abs) + Disturbance + rescale(propagule_pressure) +
                                             rescale(suitabilitypc2V2abs) + 
                                             (1|Family) + (1|POWO.name) + 
                                             (1|Site) + offset(log(Density)), 
                                           family=binomial(link="logit"), 
                                           data=surviveyr2late)

surviveyr2late_YN_Mean_AlphaPhy <- glmmTMB(surviveyr2late_YN ~ Mean_AlphaPhy*suitabilitypc1V2abs + Disturbance +  propagule_pressure +
                                             suitabilitypc2V2abs + 
                                             (1|Family) + (1|POWO.name) + 
                                             (1|Site) + offset(log(Density)), 
                                           family=binomial(link="logit"), 
                                           data=surviveyr2late)
summary(surviveyr2late_YN_Mean_AlphaPhy)
Anova(surviveyr2late_YN_Mean_AlphaPhy, type="III")
confint(surviveyr2late_YN_Mean_AlphaPhy)
r.squaredGLMM(surviveyr2late_YN_Mean_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2late_YN_Mean_AlphaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr2late_YN_Mean_AlphaPhy)
check_singularity(surviveyr2late_YN_Mean_AlphaPhy)

#Test if native status explains anything

surviveyr2late_YN_Mean_AlphaPhyN <- glmmTMB(surviveyr2late_YN ~ Mean_AlphaPhy*suitabilitypc1V2abs + Disturbance + propagule_pressure +
                                              suitabilitypc2V2abs + Status +
                                              (1|Family) + (1|POWO.name) + 
                                              (1|Site) + offset(log(Density)), 
                                            family=binomial(link="logit"), 
                                            data=surviveyr2late)
summary(surviveyr2late_YN_Mean_AlphaPhyN)

AICtab(surviveyr2late_YN_Mean_AlphaPhyN, surviveyr2late_YN_Mean_AlphaPhy)

#Testing whether worth adding biogeographic origin

surviveyr2lateBG <- left_join(surviveyr2late, BG, by="plant_name_id")

surviveyr2late_YN_Mean_AlphaPhyBG <- glmmTMB(surviveyr2late_YN ~ Mean_AlphaPhy*suitabilitypc1V2abs + Disturbance +  propagule_pressure +
                                               suitabilitypc2V2abs + 
                                               (1|Family) + (1|POWO.name) + (1|BG_ID) +
                                               (1|Site) + offset(log(Density)), 
                                             family=binomial(link="logit"), 
                                             data=surviveyr2lateBG)
summary(surviveyr2late_YN_Mean_AlphaPhyBG)

AICtab(surviveyr2late_YN_Mean_AlphaPhyBG, surviveyr2late_YN_Mean_AlphaPhy)

#Maximum alpha diversity

surviveyr2late_YN_Max_AlphaPhy <- glmmTMB(surviveyr2late_YN ~ rescale(Max_AlphaPhy) + Disturbance + rescale(propagule_pressure) +
                                            rescale(suitabilitypc1V2abs) + rescale(suitabilitypc2V2abs) + 
                                            (1|Family) + (1|POWO.name) + 
                                            (1|Site) + offset(log(Density)), 
                                          family=binomial(link="logit"), 
                                          data=surviveyr2late)

surviveyr2late_YN_Max_AlphaPhy <- glmmTMB(surviveyr2late_YN ~ Max_AlphaPhy + Disturbance + propagule_pressure +
                                            suitabilitypc1V2abs + suitabilitypc2V2abs + 
                                            (1|Family) + (1|POWO.name) + 
                                            (1|Site) + offset(log(Density)), 
                                          family=binomial(link="logit"), 
                                          data=surviveyr2late)
summary(surviveyr2late_YN_Max_AlphaPhy)
Anova(surviveyr2late_YN_Max_AlphaPhy, type="III")
confint(surviveyr2late_YN_Max_AlphaPhy)
r.squaredGLMM(surviveyr2late_YN_Max_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2late_YN_Max_AlphaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr2late_YN_Max_AlphaPhy)
check_singularity(surviveyr2late_YN_Max_AlphaPhy)

#Gamma diversity 

surviveyr2late_YN_GammaPhy <- glmmTMB(surviveyr2late_YN ~ rescale(log(GammaPhy))*rescale(suitabilitypc1V2abs) + Disturbance + rescale(propagule_pressure) +
                                        rescale(suitabilitypc2V2abs) + 
                                        (1|Family) + (1|POWO.name) + 
                                        (1|Site) + offset(log(Density)), 
                                      family=binomial(link="logit"), 
                                      data=surviveyr2late)

surviveyr2late_YN_GammaPhy <- glmmTMB(surviveyr2late_YN ~ log(GammaPhy)*suitabilitypc1V2abs + Disturbance + propagule_pressure +
                                        suitabilitypc2V2abs + 
                                        (1|Family) + (1|POWO.name) + 
                                        (1|Site) + offset(log(Density)), 
                                      family=binomial(link="logit"), 
                                      data=surviveyr2late)
summary(surviveyr2late_YN_GammaPhy)
Anova(surviveyr2late_YN_GammaPhy, type="III")
confint(surviveyr2late_YN_GammaPhy)
r.squaredGLMM(surviveyr2late_YN_GammaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2late_YN_GammaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr2late_YN_GammaPhy)
check_singularity(surviveyr2late_YN_GammaPhy)

#Median alpha diversity 

surviveyr2late_YN_MedPhy <- glmmTMB(surviveyr2late_YN ~ rescale(Med_AlphaPhy)*rescale(suitabilitypc1V2abs) + Disturbance + rescale(propagule_pressure) +
                                      rescale(suitabilitypc2V2abs) + 
                                      (1|Family) + (1|POWO.name) + 
                                      (1|Site) + offset(log(Density)), 
                                    family=binomial(link="logit"), 
                                    data=surviveyr2late)

surviveyr2late_YN_MedPhy <- glmmTMB(surviveyr2late_YN ~ Med_AlphaPhy*suitabilitypc1V2abs + Disturbance + propagule_pressure +
                                      suitabilitypc2V2abs + 
                                      (1|Family) + (1|POWO.name) + 
                                      (1|Site) + offset(log(Density)), 
                                    family=binomial(link="logit"), 
                                    data=surviveyr2late)
summary(surviveyr2late_YN_MedPhy)
Anova(surviveyr2late_YN_MedPhy, type="III")
confint(surviveyr2late_YN_MedPhy)
r.squaredGLMM(surviveyr2late_YN_MedPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2late_YN_MedPhy, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr2late_YN_MedPhy)
check_singularity(surviveyr2late_YN_MedPhy)

#######################################################################################

#DO NUMBER SURVIVING FOR ALL DATASETS

#GERMINATION COUNT

germinationcount <- filter(germination, germ_num>0)

#Mean alpha diversity 

germ_count_Mean_AlphaPhy <- glmmTMB(germ_num ~ rescale(Mean_AlphaPhy)*Disturbance + rescale(propagule_pressure) +
                                      rescale(suitabilitypc1V2abs) + rescale(suitabilitypc2V2abs) + rescale(Density) +
                                      (1|Family) + (1|POWO.name) + (1|Study) + (1|Site),
                                    family=nbinom2(link = "log"), 
                                    data=germinationcount)

germ_count_Mean_AlphaPhy <- glmmTMB(germ_num ~ Mean_AlphaPhy*Disturbance + propagule_pressure +
                                      suitabilitypc1V2abs + suitabilitypc2V2abs + Density +
                                      (1|Family) + (1|POWO.name) + (1|Study) + (1|Site),
                                    family=nbinom2(link = "log"), 
                                    data=germinationcount)
summary(germ_count_Mean_AlphaPhy)
Anova(germ_count_Mean_AlphaPhy, type="III")
confint(germ_count_Mean_AlphaPhy)
r.squaredGLMM(germ_count_Mean_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = germ_count_Mean_AlphaPhy, plot=F)
plot(simulationOutput)

check_overdispersion(germ_count_Mean_AlphaPhy)
check_singularity(germ_count_Mean_AlphaPhy)
check_model(germ_count_Mean_AlphaPhy)

#Testing if native status explains anything

germ_count_Mean_AlphaPhyN <- glmmTMB(germ_num ~ Mean_AlphaPhy*Disturbance + propagule_pressure +
                                       suitabilitypc1V2abs + suitabilitypc2V2abs + Density + Status +
                                       (1|Family) + (1|POWO.name) + (1|Study) + (1|Site),
                                     family=nbinom2(link = "log"), 
                                     data=germinationcount)
summary(germ_count_Mean_AlphaPhyN)

AICtab(germ_count_Mean_AlphaPhyN, germ_count_Mean_AlphaPhy)

#Test if biogeographic origin explains anything

germinationcountBG <- filter(germinationBG, germ_num>0)

germ_count_Mean_AlphaPhyBG <- glmmTMB(germ_num ~ Mean_AlphaPhy*Disturbance + propagule_pressure +
                                        suitabilitypc1V2abs + suitabilitypc2V2abs + Density +
                                        (1|Family) + (1|POWO.name) + (1|Study) + (1|Site) + (1|BG_ID),
                                      family=nbinom2(link = "log"), 
                                      data=germinationcountBG)
summary(germ_count_Mean_AlphaPhyBG)

AICtab(germ_count_Mean_AlphaPhyBG, germ_count_Mean_AlphaPhy)

#Maximum alpha diversity

germ_count_Max_AlphaPhy <- glmmTMB(germ_num ~ rescale(Max_AlphaPhy)*Disturbance + rescale(propagule_pressure) +
                                     rescale(suitabilitypc1V2abs) + rescale(suitabilitypc2V2abs) + rescale(Density) +
                                     (1|Family) + (1|POWO.name) + (1|Study) + (1|Site),
                                   family=nbinom2(link = "log"), 
                                   data=germinationcount)

germ_count_Max_AlphaPhy <- glmmTMB(germ_num ~ Max_AlphaPhy*Disturbance + propagule_pressure +
                                     suitabilitypc1V2abs + suitabilitypc2V2abs + Density +
                                     (1|Family) + (1|POWO.name) + (1|Study) + (1|Site),
                                   family=nbinom2(link = "log"), 
                                   data=germinationcount)
summary(germ_count_Max_AlphaPhy)
Anova(germ_count_Max_AlphaPhy, type="III")
confint(germ_count_Max_AlphaPhy)
r.squaredGLMM(germ_count_Max_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = germ_count_Max_AlphaPhy, plot=F)
plot(simulationOutput)

check_overdispersion(germ_count_Max_AlphaPhy)
check_singularity(germ_count_Max_AlphaPhy)

#Gamma diversity 

germ_count_GammaPhy <- glmmTMB(germ_num ~ recsale(log(GammaPhy)) + Disturbance + recsale(propagule_pressure) +
                                 recsale(suitabilitypc1V2abs) + recsale(suitabilitypc2V2abs) + recsale(Density) +
                                 (1|Family) + (1|POWO.name) + (1|Study) + (1|Site),
                               family=nbinom2(link = "log"), 
                               data=germinationcount)

germ_count_GammaPhy <- glmmTMB(germ_num ~ log(GammaPhy) + Disturbance + propagule_pressure +
                                 suitabilitypc1V2abs + suitabilitypc2V2abs + Density +
                                 (1|Family) + (1|POWO.name) + (1|Study) + (1|Site),
                               family=nbinom2(link = "log"), 
                               data=germinationcount)
summary(germ_count_GammaPhy)
Anova(germ_count_GammaPhy, type="III")
confint(germ_count_GammaPhy)
r.squaredGLMM(germ_count_GammaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = germ_count_GammaPhy, plot=F)
plot(simulationOutput)

check_overdispersion(germ_count_GammaPhy)
check_singularity(germ_count_GammaPhy)

#Median Alpha diversity

germ_count_MedPhy <- glmmTMB(germ_num ~ rescale(Med_AlphaPhy)*rescale(suitabilitypc1V2abs) + rescale(Med_AlphaPhy)*Disturbance + rescale(propagule_pressure) +
                               rescale(suitabilitypc2V2abs) + rescale(Density) +
                               (1|Family) + (1|POWO.name) + (1|Study) + (1|Site),
                             family=nbinom2(link = "log"), 
                             data=germinationcount)

germ_count_MedPhy <- glmmTMB(germ_num ~ Med_AlphaPhy*suitabilitypc1V2abs + Med_AlphaPhy*Disturbance + propagule_pressure +
                               suitabilitypc2V2abs + Density +
                               (1|Family) + (1|POWO.name) + (1|Study) + (1|Site),
                             family=nbinom2(link = "log"), 
                             data=germinationcount)
summary(germ_count_MedPhy)
Anova(germ_count_MedPhy, type="III")
confint(germ_count_MedPhy)
r.squaredGLMM(germ_count_MedPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = germ_count_MedPhy, plot=F)
plot(simulationOutput)

check_overdispersion(germ_count_MedPhy)
check_singularity(germ_count_MedPhy)

#FIRST YEAR SURVIVAL COUNT

surviveyr1count <- filter(surviveyr1, surviveyr1_num>0)

#Mean alpha diversity 

surviveyr1_count_Mean_AlphaPhy <- glmmTMB(surviveyr1_num ~ rescale(Mean_AlphaPhy)*Disturbance + rescale(propagule_pressure) +
                                            rescale(suitabilitypc1V2abs) + rescale(suitabilitypc2V2abs) + rescale(Density) + Herbivory +
                                            (1|Family) + (1|POWO.name) + (1|Site) + (1|Study),
                                          family=nbinom2(link = "log"), 
                                          data=surviveyr1count)

surviveyr1_count_Mean_AlphaPhy <- glmmTMB(surviveyr1_num ~ Mean_AlphaPhy*Disturbance + propagule_pressure +
                                            suitabilitypc1V2abs + suitabilitypc2V2abs + Density + Herbivory +
                                            (1|Family) + (1|POWO.name) + (1|Site) + (1|Study),
                                          family=nbinom2(link = "log"), 
                                          data=surviveyr1count)
summary(surviveyr1_count_Mean_AlphaPhy)
Anova(surviveyr1_count_Mean_AlphaPhy, type="III")
confint(surviveyr1_count_Mean_AlphaPhy)
r.squaredGLMM(surviveyr1_count_Mean_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr1_count_Mean_AlphaPhy, plot=F)
plot(simulationOutput)

check_overdispersion(surviveyr1_count_Mean_AlphaPhy)
check_singularity(surviveyr1_count_Mean_AlphaPhy)

#Testing if native status explains anything

surviveyr1_count_Mean_AlphaPhyN <- glmmTMB(surviveyr1_num ~ Mean_AlphaPhy*Disturbance + propagule_pressure +
                                             suitabilitypc1V2abs + suitabilitypc2V2abs + Density + Herbivory + Status +
                                             (1|Family) + (1|POWO.name) + (1|Site) + (1|Study),
                                           family=nbinom2(link = "log"), 
                                           data=surviveyr1count)
summary(surviveyr1_count_Mean_AlphaPhyN)

AICtab(surviveyr1_count_Mean_AlphaPhyN, surviveyr1_count_Mean_AlphaPhy)

#Test if biogeographic origin explains anything

surviveyr1countBG <- filter(surviveyr1BG, surviveyr1_num>0)

surviveyr1_count_Mean_AlphaPhyBG <- glmmTMB(surviveyr1_num ~ Mean_AlphaPhy + Disturbance + propagule_pressure +
                                              suitabilitypc1V2abs + suitabilitypc2V2abs + Density + Herbivory +
                                              (1|Family) + (1|POWO.name) + (1|Site) + (1|Study) + (1|BG_ID),
                                            family=nbinom2(link = "log"), 
                                            data=surviveyr1countBG)
summary(surviveyr1_count_Mean_AlphaPhyBG)

AICtab(surviveyr1_count_Mean_AlphaPhyBG, surviveyr1_count_Mean_AlphaPhy)

#Maximum alpha diversity

surviveyr1_count_Max_AlphaPhy <- glmmTMB(surviveyr1_num ~ rescale(Max_AlphaPhy)*Disturbance + rescale(propagule_pressure) +
                                           rescale(suitabilitypc1V2abs) + rescale(suitabilitypc2V2abs) + rescale(Density) + Herbivory +
                                           (1|Family) + (1|POWO.name) + (1|Site) + (1|Study),
                                         family=nbinom2(link = "log"), 
                                         data=surviveyr1count)

surviveyr1_count_Max_AlphaPhy <- glmmTMB(surviveyr1_num ~ Max_AlphaPhy*Disturbance + propagule_pressure +
                                           suitabilitypc1V2abs + suitabilitypc2V2abs + Density + Herbivory +
                                           (1|Family) + (1|POWO.name) + (1|Site) + (1|Study),
                                         family=nbinom2(link = "log"), 
                                         data=surviveyr1count)
summary(surviveyr1_count_Max_AlphaPhy)
Anova(surviveyr1_count_Max_AlphaPhy, type="III")
confint(surviveyr1_count_Max_AlphaPhy)
r.squaredGLMM(surviveyr1_count_Max_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr1_count_Max_AlphaPhy, plot=F)
plot(simulationOutput)

check_overdispersion(surviveyr1_count_Max_AlphaPhy)
check_singularity(surviveyr1_count_Max_AlphaPhy)

#Gamma diversity

surviveyr1_count_GammaPhy <- glmmTMB(surviveyr1_num ~ rescale(log(GammaPhy)) + Disturbance + rescale(propagule_pressure) +
                                       rescale(suitabilitypc1V2abs) + rescale(suitabilitypc2V2abs) + rescale(Density) + Herbivory +
                                       (1|Family) + (1|POWO.name) + (1|Site) + (1|Study),
                                     family=nbinom2(link = "log"), 
                                     data=surviveyr1count)

surviveyr1_count_GammaPhy <- glmmTMB(surviveyr1_num ~ log(GammaPhy) + Disturbance + propagule_pressure +
                                       suitabilitypc1V2abs + suitabilitypc2V2abs + Density + Herbivory +
                                       (1|Family) + (1|POWO.name) + (1|Site) + (1|Study),
                                     family=nbinom2(link = "log"), 
                                     data=surviveyr1count)
summary(surviveyr1_count_GammaPhy)
Anova(surviveyr1_count_GammaPhy, type="III")
confint(surviveyr1_count_GammaPhy)
r.squaredGLMM(surviveyr1_count_GammaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr1_count_GammaPhy, plot=F)
plot(simulationOutput)

check_overdispersion(surviveyr1_count_GammaPhy)
check_singularity(surviveyr1_count_GammaPhy)

#Median alpha diversity

surviveyr1_count_MedPhy <- glmmTMB(surviveyr1_num ~ rescale(Med_AlphaPhy)*Disturbance + rescale(propagule_pressure) +
                                     rescale(suitabilitypc1V2abs) + rescale(suitabilitypc2V2abs) + rescale(Density) + Herbivory +
                                     (1|Family) + (1|POWO.name) + (1|Site) + (1|Study),
                                   family=nbinom2(link = "log"), 
                                   data=surviveyr1count)

surviveyr1_count_MedPhy <- glmmTMB(surviveyr1_num ~ Med_AlphaPhy*Disturbance + propagule_pressure +
                                     suitabilitypc1V2abs + suitabilitypc2V2abs + Density + Herbivory +
                                     (1|Family) + (1|POWO.name) + (1|Site) + (1|Study),
                                   family=nbinom2(link = "log"), 
                                   data=surviveyr1count)
summary(surviveyr1_count_MedPhy)
Anova(surviveyr1_count_MedPhy, type="III")
confint(surviveyr1_count_MedPhy)
r.squaredGLMM(surviveyr1_count_MedPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr1_count_MedPhy, plot=F)
plot(simulationOutput)

check_overdispersion(surviveyr1_count_MedPhy)
check_singularity(surviveyr1_count_MedPhy)

#EARLY SECOND YEAR SURVIVAL COUNT

surviveyr2earlycount <- filter(surviveyr2early, surviveyr2early_num>0)

#Mean alpha diversity 

surviveyr2early_count_Mean_AlphaPhy <- glmmTMB(surviveyr2early_num ~ rescale(Mean_AlphaPhy) + Disturbance + rescale(propagule_pressure) +
                                                 rescale(suitabilitypc1V2abs) + rescale(suitabilitypc2V2abs) + rescale(Density) + Herbivory +
                                                 (1|Family) + (1|POWO.name) + (1|Site),
                                               family=nbinom2(link = "log"), 
                                               data=surviveyr2earlycount)

surviveyr2early_count_Mean_AlphaPhy <- glmmTMB(surviveyr2early_num ~ Mean_AlphaPhy + Disturbance + propagule_pressure +
                                                 suitabilitypc1V2abs + suitabilitypc2V2abs + Density + Herbivory +
                                                 (1|Family) + (1|POWO.name) + (1|Site),
                                               family=nbinom2(link = "log"), 
                                               data=surviveyr2earlycount)
summary(surviveyr2early_count_Mean_AlphaPhy)
Anova(surviveyr2early_count_Mean_AlphaPhy, type="III")
confint(surviveyr2early_count_Mean_AlphaPhy)
r.squaredGLMM(surviveyr2early_count_Mean_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2early_count_Mean_AlphaPhy, plot=F)
plot(simulationOutput)

check_overdispersion(surviveyr2early_count_Mean_AlphaPhy)
check_singularity(surviveyr2early_count_Mean_AlphaPhy)

#Testing if native status explains anything

surviveyr2early_count_Mean_AlphaPhyN <- glmmTMB(surviveyr2early_num ~ Mean_AlphaPhy + Disturbance + propagule_pressure +
                                                  suitabilitypc1V2abs + suitabilitypc2V2abs + Density + Herbivory + Status +
                                                  (1|Family) + (1|POWO.name) + (1|Site),
                                                family=nbinom2(link = "log"), 
                                                data=surviveyr2earlycount)
summary(surviveyr2early_count_Mean_AlphaPhyN)

AICtab(surviveyr2early_count_Mean_AlphaPhyN, surviveyr2early_count_Mean_AlphaPhy)

#Test if biogeographic origin explains anything

surviveyr2earlycountBG <- filter(surviveyr2earlyBG, surviveyr2early_num>0)

surviveyr2early_count_Mean_AlphaPhyBG <- glmmTMB(surviveyr2early_num ~ Mean_AlphaPhy + Disturbance + propagule_pressure +
                                                   suitabilitypc1V2abs + suitabilitypc2V2abs + Density + Herbivory +
                                                   (1|Family) + (1|POWO.name) + (1|Site) + (1|BG_ID),
                                                 family=nbinom2(link = "log"), 
                                                 data=surviveyr2earlycountBG)
summary(surviveyr2early_count_Mean_AlphaPhyBG)

AICtab(surviveyr2early_count_Mean_AlphaPhyBG, surviveyr2early_count_Mean_AlphaPhy)

#Max alpha diversity 

surviveyr2early_count_Max_AlphaPhy <- glmmTMB(surviveyr2early_num ~ rescale(Max_AlphaPhy) + Disturbance + rescale(propagule_pressure) +
                                                rescale(suitabilitypc1V2abs) + rescale(suitabilitypc2V2abs) + rescale(Density) + Herbivory +
                                                (1|Family) + (1|POWO.name) + (1|Site),
                                              family=nbinom2(link = "log"), 
                                              data=surviveyr2earlycount)

surviveyr2early_count_Max_AlphaPhy <- glmmTMB(surviveyr2early_num ~ Max_AlphaPhy + Disturbance + propagule_pressure +
                                                suitabilitypc1V2abs + suitabilitypc2V2abs + Density + Herbivory +
                                                (1|Family) + (1|POWO.name) + (1|Site),
                                              family=nbinom2(link = "log"), 
                                              data=surviveyr2earlycount)
summary(surviveyr2early_count_Max_AlphaPhy)
Anova(surviveyr2early_count_Max_AlphaPhy, type="III")
confint(surviveyr2early_count_Max_AlphaPhy)
r.squaredGLMM(surviveyr2early_count_Max_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2early_count_Max_AlphaPhy, plot=F)
plot(simulationOutput)

check_overdispersion(surviveyr2early_count_Max_AlphaPhy)
check_singularity(surviveyr2early_count_Max_AlphaPhy)

#Gamma diversity 

surviveyr2early_count_GammaPhy <- glmmTMB(surviveyr2early_num ~ rescale(log(GammaPhy))*Disturbance + rescale(propagule_pressure) +
                                            rescale(suitabilitypc1V2abs) + rescale(suitabilitypc2V2abs) + rescale(Density) + Herbivory +
                                            (1|Family) + (1|POWO.name) + (1|Site),
                                          family=nbinom2(link = "log"), 
                                          data=surviveyr2earlycount)

surviveyr2early_count_GammaPhy <- glmmTMB(surviveyr2early_num ~ log(GammaPhy)*Disturbance + propagule_pressure +
                                            suitabilitypc1V2abs + suitabilitypc2V2abs + Density + Herbivory +
                                            (1|Family) + (1|POWO.name) + (1|Site),
                                          family=nbinom2(link = "log"), 
                                          data=surviveyr2earlycount)
summary(surviveyr2early_count_GammaPhy)
Anova(surviveyr2early_count_GammaPhy, type="III")
confint(surviveyr2early_count_GammaPhy)
r.squaredGLMM(surviveyr2early_count_GammaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2early_count_GammaPhy, plot=F)
plot(simulationOutput)

check_overdispersion(surviveyr2early_count_GammaPhy)
check_singularity(surviveyr2early_count_GammaPhy)

#Median alpha diversity 

surviveyr2early_count_MedPhy <- glmmTMB(surviveyr2early_num ~ rescale(Med_AlphaPhy) + Disturbance + rescale(propagule_pressure) +
                                          rescale(suitabilitypc1V2abs) + rescale(suitabilitypc2V2abs) + rescale(Density) + Herbivory +
                                          (1|Family) +(1|POWO.name) + (1|Site),
                                        family=nbinom2(link = "log"), 
                                        data=surviveyr2earlycount)

surviveyr2early_count_MedPhy <- glmmTMB(surviveyr2early_num ~ Med_AlphaPhy + Disturbance + propagule_pressure +
                                          suitabilitypc1V2abs + suitabilitypc2V2abs + Density + Herbivory +
                                          (1|Family) +(1|POWO.name) + (1|Site),
                                        family=nbinom2(link = "log"), 
                                        data=surviveyr2earlycount)
summary(surviveyr2early_count_MedPhy)
Anova(surviveyr2early_count_MedPhy, type="III")
confint(surviveyr2early_count_MedPhy)
r.squaredGLMM(surviveyr2early_count_MedPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2early_count_MedPhy, plot=F)
plot(simulationOutput)

check_overdispersion(surviveyr2early_count_MedPhy)
check_singularity(surviveyr2early_count_MedPhy)

#LATE SECOND YEAR SURVIVAL COUNT

surviveyr2latecount <- filter(surviveyr2late, surviveyr2late_num>0)

#Mean alpha diversity 

surviveyr2late_count_Mean_AlphaPhy <- glmmTMB(surviveyr2late_num ~ rescale(Mean_AlphaPhy) + Disturbance + rescale(propagule_pressure) +
                                                rescale(suitabilitypc1V2abs) + rescale(suitabilitypc2V2abs) + rescale(Density) +
                                                (1|POWO.name),
                                              family=nbinom2(link = "log"), 
                                              data=surviveyr2latecount)

surviveyr2late_count_Mean_AlphaPhy <- glmmTMB(surviveyr2late_num ~ Mean_AlphaPhy + Disturbance + propagule_pressure +
                                                suitabilitypc1V2abs + suitabilitypc2V2abs + Density +
                                                (1|POWO.name),
                                              family=nbinom2(link = "log"), 
                                              data=surviveyr2latecount)
summary(surviveyr2late_count_Mean_AlphaPhy)
Anova(surviveyr2late_count_Mean_AlphaPhy, type="III")
confint(surviveyr2late_count_Mean_AlphaPhy)
r.squaredGLMM(surviveyr2late_count_Mean_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2late_count_Mean_AlphaPhy, plot=F)
plot(simulationOutput)

check_overdispersion(surviveyr2late_count_Mean_AlphaPhy)
check_singularity(surviveyr2late_count_Mean_AlphaPhy)

#Testing whether native status explains anything

surviveyr2late_count_Mean_AlphaPhyN <- glmmTMB(surviveyr2late_num ~ Mean_AlphaPhy + Disturbance + propagule_pressure +
                                                 suitabilitypc1V2abs + suitabilitypc2V2abs + Density + Status +
                                                 (1|POWO.name),
                                               family=nbinom2(link = "log"), 
                                               data=surviveyr2latecount)
summary(surviveyr2late_count_Mean_AlphaPhyN)

AICtab(surviveyr2late_count_Mean_AlphaPhyN, surviveyr2late_count_Mean_AlphaPhy)

#Test if biogeographic origin explains anything

surviveyr2latecountBG <- filter(surviveyr2lateBG, surviveyr2late_num>0)

surviveyr2late_count_Mean_AlphaPhyBG <- glmmTMB(surviveyr2late_num ~ Mean_AlphaPhy + Disturbance + propagule_pressure +
                                                  suitabilitypc1V2abs + suitabilitypc2V2abs + Density +
                                                  (1|POWO.name) + (1|BG_ID),
                                                family=nbinom2(link = "log"), 
                                                data=surviveyr2latecountBG)
summary(surviveyr2late_count_Mean_AlphaPhyBG)

AICtab(surviveyr2late_count_Mean_AlphaPhyBG, surviveyr2late_count_Mean_AlphaPhy)

#Max alpha diversity 

surviveyr2late_count_Max_AlphaPhy <- glmmTMB(surviveyr2late_num ~ rescale(Max_AlphaPhy) + Disturbance + rescale(propagule_pressure) +
                                               rescale(suitabilitypc1V2abs) + rescale(suitabilitypc2V2abs) + rescale(Density) +
                                               (1|POWO.name),
                                             family=nbinom2(link = "log"), 
                                             data=surviveyr2latecount)

surviveyr2late_count_Max_AlphaPhy <- glmmTMB(surviveyr2late_num ~ Max_AlphaPhy + Disturbance + propagule_pressure +
                                               suitabilitypc1V2abs + suitabilitypc2V2abs + Density +
                                               (1|POWO.name),
                                             family=nbinom2(link = "log"), 
                                             data=surviveyr2latecount)
summary(surviveyr2late_count_Max_AlphaPhy)
Anova(surviveyr2late_count_Max_AlphaPhy, type="III")
confint(surviveyr2late_count_Max_AlphaPhy)
r.squaredGLMM(surviveyr2late_count_Max_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2late_count_Max_AlphaPhy, plot=F)
plot(simulationOutput)

check_overdispersion(surviveyr2late_count_Max_AlphaPhy)
check_singularity(surviveyr2late_count_Max_AlphaPhy)

#Gamma diversity 

surviveyr2late_count_GammaPhy <- glmmTMB(surviveyr2late_num ~ rescale(log(GammaPhy)) + Disturbance + rescale(propagule_pressure) +
                                           rescale(suitabilitypc1V2abs) + rescale(suitabilitypc2V2abs) + rescale(Density) +
                                           (1|POWO.name),
                                         family=nbinom2(link = "log"), 
                                         data=surviveyr2latecount)

surviveyr2late_count_GammaPhy <- glmmTMB(surviveyr2late_num ~ log(GammaPhy) + Disturbance + propagule_pressure +
                                           suitabilitypc1V2abs + suitabilitypc2V2abs + Density +
                                           (1|POWO.name),
                                         family=nbinom2(link = "log"), 
                                         data=surviveyr2latecount)
summary(surviveyr2late_count_GammaPhy)
Anova(surviveyr2late_count_GammaPhy, type="III")
confint(surviveyr2late_count_GammaPhy)
r.squaredGLMM(surviveyr2late_count_GammaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2late_count_GammaPhy, plot=F)
plot(simulationOutput)

check_overdispersion(surviveyr2late_count_GammaPhy)
check_singularity(surviveyr2late_count_GammaPhy)

#Median alpha diversity 

surviveyr2late_count_MedPhy <- glmmTMB(surviveyr2late_num ~ rescale(Med_AlphaPhy)*rescale(suitabilitypc2V2abs) + Disturbance + rescale(propagule_pressure) +
                                         rescale(suitabilitypc1V2abs) + rescale(Density) +
                                         (1|POWO.name),
                                       family=nbinom2(link = "log"), 
                                       data=surviveyr2latecount)

surviveyr2late_count_MedPhy <- glmmTMB(surviveyr2late_num ~ Med_AlphaPhy*suitabilitypc2V2abs + Disturbance + propagule_pressure +
                                         suitabilitypc1V2abs + Density +
                                         (1|POWO.name),
                                       family=nbinom2(link = "log"), 
                                       data=surviveyr2latecount)
summary(surviveyr2late_count_MedPhy)
Anova(surviveyr2late_count_MedPhy, type="III")
confint(surviveyr2late_count_MedPhy)
r.squaredGLMM(surviveyr2late_count_MedPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2late_count_MedPhy, plot=F)
plot(simulationOutput)

check_overdispersion(surviveyr2late_count_MedPhy)
check_singularity(surviveyr2late_count_MedPhy)

#########################################################################
#########################################################################

#Plotting all the results

#Plot all effect sizes for Mean_AlphaPhy - Figure 2

effectsizes <- read.csv("all3results_effectsize_revised.csv", header=T, stringsAsFactors = T)
effectsizes$Stage = factor(effectsizes$Stage, levels = c("Second year","Overwinter","First year","Colonisation"))
effectsizes$Variable = factor(effectsizes$Variable, levels = c("Alpha PD x Climate (temp.)", "Alpha PD x Disturbed","Herbivory", "Climate dissimilarity (temp.)",
                                                               "Climate dissimilarity (precip.)", "Disturbed","Mean Alpha PD"))

effectplot <- ggplot(effectsizes, aes(x=effect, y=Variable, color=Stage, shape=significant)) +
  geom_point(size=2.3, position=position_dodge(0.5)) +
  scale_shape_manual(values=c(2, 1, 16)) +
  geom_errorbar(data=effectsizes, aes(y=Variable, xmin=effect_low, xmax=effect_high), 
                width=0, size=1, position=position_dodge(0.5)) +
  facet_grid(. ~ Metric, scales="free") +
  scale_color_viridis(discrete = TRUE, guide = guide_legend(reverse = TRUE)) +
  theme_bw() + geom_vline(xintercept=0, linetype="solid", color="grey", size=0.9) +
  xlab("Standardised effect size ± 95% C.I.") + guides(shape = "none") + theme(legend.position = "top") +
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=10), 
        legend.text=element_text(size=8), 
        legend.title=element_text(size=10))
effectplot

ggsave(effectplot, 
       filename = "Fig_2.svg",
       height = 130, width = 180, units = "mm")

#PLOT ALL RESULTS FOR MEAN ALPHA PHY - Figure 3

#Fig 3a

it_vr <- visreg(germ_YN_Mean_AlphaPhy, "Mean_AlphaPhy", ylab="Germination Success",
                xlab="Mean_AlphaPhy", scale="response", partial=T)
it_vr_df <- as.data.frame(cbind(it_vr$fit$Mean_AlphaPhy, it_vr$fit$visregFit,
                                it_vr$fit$visregLwr, it_vr$fit$visregUpr))
colnames(it_vr_df) <- c("it", "fit", "lwr", "upr")
it_vr_pr <- as.data.frame(cbind(it_vr$res$Mean_AlphaPhy,
                                it_vr$res$visregRes))
colnames(it_vr_pr) <- c("it", "resids")

fig3a <- ggplot() +
  geom_jitter(aes(y=resids, x= it), alpha = 0.8, data = it_vr_pr, width=0.1, height=0.0001, colour="darkgray") +
  #geom_line(aes(x=it, y=fit), data=it_vr_df, size = 1.2, linetype="solid") +
  theme_bw() +
  ylab("Colonisation") +
  xlab("Mean Alpha Phylogenetic Diversity") + 
  geom_vline(xintercept = 0, linetype="dashed", colour="grey", size=1) +
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=10), 
        legend.text=element_text(size=8), 
        legend.title=element_text(size=10)) +
  xlim(-1.5, 1.09) + xlab(NULL)
fig3a

#Fig 3b

germmeanalphacountinteraction <- interact_plot(germ_count_Mean_AlphaPhy, pred="Mean_AlphaPhy",
                                               modx = "Disturbance",
                                               plot.points=TRUE, partial.residuals = TRUE) +
  scale_y_continuous(trans='log10') + 
  xlab("Mean Alpha Phylogenetic Diversity") + ylab("Number of plants colonising") +
  scale_color_manual(values=c("#999999", "#E69F00")) +
  scale_linetype_manual(values=c("solid", "solid")) +
  theme_bw() + guides(linetype = "none") + theme(legend.position = "top")  +
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=10), 
        legend.text=element_text(size=8), 
        legend.title=element_text(size=10))
germmeanalphacountinteraction

#Fig 3c

it_vr <- visreg(surviveyr1_YN_Mean_AlphaPhy, "Mean_AlphaPhy", ylab="Germination Success",
                xlab="Mean_AlphaPhy", scale="response", partial=T)
it_vr_df <- as.data.frame(cbind(it_vr$fit$Mean_AlphaPhy, it_vr$fit$visregFit,
                                it_vr$fit$visregLwr, it_vr$fit$visregUpr))
colnames(it_vr_df) <- c("it", "fit", "lwr", "upr")
it_vr_pr <- as.data.frame(cbind(it_vr$res$Mean_AlphaPhy,
                                it_vr$res$visregRes))
colnames(it_vr_pr) <- c("it", "resids")

fig3c <- ggplot() +
  geom_jitter(aes(y=resids, x= it), alpha = 0.8, data = it_vr_pr, width=0.1, height=0.0001, colour="darkgray") +
  geom_line(aes(x=it, y=fit), data=it_vr_df, size = 1.2, linetype="solid") +
  theme_bw() +
  ylab("First year") +
  xlab("Mean Alpha Phylogenetic Diversity") +
  geom_vline(xintercept = 0, linetype="dashed", colour="grey", size=1) +
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=10), 
        legend.text=element_text(size=8), 
        legend.title=element_text(size=10))  +
  xlim(-1.5, 1.09) + xlab(NULL) 
fig3c

#Fig 3d

yr1meanalphacountinteraction <- interact_plot(surviveyr1_count_Mean_AlphaPhy, pred="Mean_AlphaPhy", modx = "Disturbance",
                                              plot.points=TRUE, partial.residuals = TRUE) +
  scale_y_continuous(trans='log10') + 
  xlab("Mean Alpha Phylogenetic Diversity") + ylab("Number of plants surviving to end of first year") +
  scale_color_manual(values=c("#999999", "#E69F00")) +
  scale_linetype_manual(values=c("solid", "solid")) +
  theme_bw() + guides(linetype = "none") + theme(legend.position = "top")  +
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=10), 
        legend.text=element_text(size=9), 
        legend.title=element_text(size=10))
yr1meanalphacountinteraction

#Fig 3e

it_vr <- visreg(surviveyr2early_YN_Mean_AlphaPhy, "Mean_AlphaPhy", ylab="Germination Success",
                xlab="Mean_AlphaPhy", scale="response", partial=T)
it_vr_df <- as.data.frame(cbind(it_vr$fit$Mean_AlphaPhy, it_vr$fit$visregFit,
                                it_vr$fit$visregLwr, it_vr$fit$visregUpr))
colnames(it_vr_df) <- c("it", "fit", "lwr", "upr")
it_vr_pr <- as.data.frame(cbind(it_vr$res$Mean_AlphaPhy,
                                it_vr$res$visregRes))
colnames(it_vr_pr) <- c("it", "resids")

fig3e <- ggplot() +
  geom_jitter(aes(y=resids, x= it), alpha = 0.8, data = it_vr_pr, width=0.1, height=0.0001, colour="darkgray") +
  #geom_line(aes(x=it, y=fit), data=it_vr_df, size = 1.2, linetype="solid") +
  theme_bw() +
  ylab("Overwinter") +
  xlab("Mean Alpha Phylogenetic Diversity") +
  geom_vline(xintercept = 0, linetype="dashed", colour="grey", size=1) +
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=10), 
        legend.text=element_text(size=8), 
        legend.title=element_text(size=10))  +
  xlim(-1.5, 1.09) + xlab(NULL)
fig3e

#Fig 3f

it_vr <- visreg(surviveyr2early_count_Mean_AlphaPhy, "Mean_AlphaPhy", ylab="Germination Success",
                xlab="Mean_AlphaPhy",  partial=T) #scale="response",
it_vr_df <- as.data.frame(cbind(it_vr$fit$Mean_AlphaPhy, it_vr$fit$visregFit,
                                it_vr$fit$visregLwr, it_vr$fit$visregUpr))
colnames(it_vr_df) <- c("it", "fit", "lwr", "upr")
it_vr_pr <- as.data.frame(cbind(it_vr$res$Mean_AlphaPhy,
                                it_vr$res$visregRes))
colnames(it_vr_pr) <- c("it", "resids")

fig3f <- ggplot() +
  geom_jitter(aes(y=resids, x= it), alpha = 0.8, data = it_vr_pr, width=0.1, colour="darkgray") +
  #geom_line(aes(x=it, y=fit), data=it_vr_df, size = 1.2, linetype="solid") +
  theme_bw() +
  #scale_y_continuous(trans='log10') +
  ylab("Overwinter") +
  xlab("Mean Alpha Phylogenetic Diversity") +
  geom_vline(xintercept = 0, linetype="dashed", colour="grey", size=1) +
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=10), 
        legend.text=element_text(size=8), 
        legend.title=element_text(size=10))  +
  xlim(-1.5, 1.09) + xlab(NULL) + ylab(NULL)
fig3f

#Fig 3g

yr2YNtempinteraction <- interact_plot(surviveyr2late_YN_Mean_AlphaPhy, pred="Mean_AlphaPhy", modx = "suitabilitypc1V2abs",
                                      plot.points=TRUE, partial.residuals = TRUE, legend.main = "Climate dissimilarity (temp.)") +
  scale_y_continuous(trans='log10') +
  xlab("Mean Alpha Phylogenetic Diversity") + ylab("Likelihood of any plant in plot surviving to end of second year") +
  scale_linetype_manual(values=c("solid", "solid", "solid")) +
  theme_bw() + guides(linetype = "none") + theme(legend.position = "top")  +
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=10), 
        legend.text=element_text(size=8), 
        legend.title=element_text(size=10))
yr2YNtempinteraction

#Fig 3h

it_vr <- visreg(surviveyr2late_count_Mean_AlphaPhy, "Mean_AlphaPhy", ylab="Germination Success",
                xlab="Mean_AlphaPhy",  partial=T) #scale="response",
it_vr_df <- as.data.frame(cbind(it_vr$fit$Mean_AlphaPhy, it_vr$fit$visregFit,
                                it_vr$fit$visregLwr, it_vr$fit$visregUpr))
colnames(it_vr_df) <- c("it", "fit", "lwr", "upr")
it_vr_pr <- as.data.frame(cbind(it_vr$res$Mean_AlphaPhy,
                                it_vr$res$visregRes))
colnames(it_vr_pr) <- c("it", "resids")

fig3h <- ggplot() +
  geom_jitter(aes(y=resids, x= it), alpha = 0.8, data = it_vr_pr, width=0.1, colour="darkgray") +
  #geom_line(aes(x=it, y=fit), data=it_vr_df, size = 1.2, linetype="solid") +
  theme_bw() +
  #scale_y_continuous(trans='log10') +
  ylab("Second year") +
  xlab("Mean Alpha Phylogenetic Diversity") +
  geom_vline(xintercept = 0, linetype="dashed", colour="grey", size=1) +
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=10), 
        legend.text=element_text(size=8), 
        legend.title=element_text(size=10))  +
  xlim(-1.5, 1.09) + ylab(NULL)
fig3h

Fig3 <- (fig3a | germmeanalphacountinteraction)/
  (fig3c | yr1meanalphacountinteraction)/
  (fig3e | fig3f)/
  (yr2YNtempinteraction | fig3h)
Fig3
ggsave(Fig3, 
       filename = "Figure_3.svg",
       height = 240, width = 180, units = "mm")

#MAX ALPHA PHY - INTERACTIONS - Figure S9

germmaxalphacountinteraction <- interact_plot(germ_count_Max_AlphaPhy, pred="Max_AlphaPhy", modx = "Disturbance",
                                              plot.points=TRUE, partial.residuals = TRUE) +
  scale_y_continuous(trans='log10') + 
  xlab("Maximum Alpha Phylogenetic Diversity") + ylab("Number of plants colonising") +
  scale_color_manual(values=c("#999999", "#E69F00")) +
  scale_linetype_manual(values=c("solid", "solid")) +
  theme_bw() + guides(linetype = "none") + theme(legend.position = "top")  +
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=10), 
        legend.text=element_text(size=8), 
        legend.title=element_text(size=10))
germmaxalphacountinteraction

yr1YNmaxinteraction <- interact_plot(surviveyr1_YN_Max_AlphaPhy, pred="Max_AlphaPhy", modx = "Disturbance",
                                     plot.points=TRUE, partial.residuals = TRUE) +
  scale_y_continuous(trans='log10') + 
  xlab("Maximum Alpha Phylogenetic Diversity") + ylab("Likelihood of any plant in plot surviving to end of first year") +
  scale_color_manual(values=c("#999999", "#E69F00")) +
  scale_linetype_manual(values=c("solid", "solid")) +
  theme_bw() + guides(linetype = "none") + theme(legend.position = "top")  +
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=10), 
        legend.text=element_text(size=8), 
        legend.title=element_text(size=10))
yr1YNmaxinteraction

yr1maxalphacountinteraction <- interact_plot(surviveyr1_count_Max_AlphaPhy, pred="Max_AlphaPhy", modx = "Disturbance",
                                             plot.points=TRUE, partial.residuals = TRUE) +
  scale_y_continuous(trans='log10') + 
  xlab("Maximum Alpha Phylogenetic Diversity") + ylab("Number of plants surviving to end of first year") +
  scale_color_manual(values=c("#999999", "#E69F00")) +
  scale_linetype_manual(values=c("solid", "solid")) +
  theme_bw() + guides(linetype = "none") + theme(legend.position = "top")  +
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=10), 
        legend.text=element_text(size=8), 
        legend.title=element_text(size=10))
yr1maxalphacountinteraction

yr2earlyYNmaxtempinteraction <- interact_plot(surviveyr2early_YN_Max_AlphaPhy, pred="Max_AlphaPhy", modx = "suitabilitypc2V2abs",
                                              plot.points=TRUE, partial.residuals = TRUE, legend.main = "Climate dissimilarity (precip.)") +
  scale_y_continuous(trans='log10') +
  xlab("Maximum Alpha Phylogenetic Diversity") + ylab("Likelihood of any plant in plot surviving overwinter") +
  scale_linetype_manual(values=c("solid", "solid", "solid")) +
  theme_bw() + guides(linetype = "none") + theme(legend.position = "top")  +
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=10), 
        legend.text=element_text(size=8), 
        legend.title=element_text(size=10))
yr2earlyYNmaxtempinteraction

maxalphafig <- ggarrange(germmaxalphacountinteraction, yr1YNmaxinteraction, 
                         yr1maxalphacountinteraction, yr2earlyYNmaxtempinteraction, nrow=2, ncol=2)
maxalphafig
ggsave(maxalphafig, 
       filename = "Fig_S9.svg",
       height = 210, width = 210, units = "mm")

#MEDIAN ALPHA PHY - INTERACTIONS - Figure S10

germYNmedinteraction <- interact_plot(germ_YN_MedPhy, pred="Med_AlphaPhy", modx = "Disturbance",
                                      plot.points=TRUE, partial.residuals = TRUE) +
  scale_y_continuous(trans='log10') +
  xlab("Median Alpha Phylogenetic Diversity") + ylab("Likelihood of any plant in plot colonising") +
  scale_color_manual(values=c("#999999", "#E69F00")) +
  scale_linetype_manual(values=c("solid", "solid")) +
  theme_bw() + guides(linetype = "none") + theme(legend.position = "top")  +
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=10), 
        legend.text=element_text(size=8), 
        legend.title=element_text(size=10))
germYNmedinteraction

germmedalphacountinteraction <- interact_plot(germ_count_MedPhy, pred="Med_AlphaPhy", modx = "Disturbance",
                                              plot.points=TRUE, partial.residuals = TRUE) +
  scale_y_continuous(trans='log10') + 
  xlab("Median Alpha Phylogenetic Diversity") + ylab("Number of plants colonising") +
  scale_color_manual(values=c("#999999", "#E69F00")) +
  scale_linetype_manual(values=c("solid", "solid")) +
  theme_bw() + guides(linetype = "none") + theme(legend.position = "top")  +
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=10), 
        legend.text=element_text(size=8), 
        legend.title=element_text(size=10))
germmedalphacountinteraction

germcountmedtempinteraction <- interact_plot(germ_count_MedPhy, pred="Med_AlphaPhy", modx = "suitabilitypc1V2abs",
                                             plot.points=TRUE, partial.residuals = TRUE, legend.main = "Climate dissimilarity (temp.)") +
  scale_y_continuous(trans='log10') +
  xlab("Median Alpha Phylogenetic Diversity") + ylab("Number of plants colonising") +
  scale_linetype_manual(values=c("solid", "solid", "solid")) +
  theme_bw() + guides(linetype = "none") + theme(legend.position = "top")  +
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=10), 
        legend.text=element_text(size=8), 
        legend.title=element_text(size=10))
germcountmedtempinteraction

yr1medalphacountinteraction <- interact_plot(surviveyr1_count_MedPhy, pred="Med_AlphaPhy", modx = "Disturbance",
                                             plot.points=TRUE, partial.residuals = TRUE) +
  scale_y_continuous(trans='log10') + 
  xlab("Median Alpha Phylogenetic Diversity") + ylab("Number of plants surviving to end of first year") +
  scale_color_manual(values=c("#999999", "#E69F00")) +
  scale_linetype_manual(values=c("solid", "solid")) +
  theme_bw() + guides(linetype = "none") + theme(legend.position = "top")  +
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=10), 
        legend.text=element_text(size=8), 
        legend.title=element_text(size=10))
yr1medalphacountinteraction

yr2earlyYNmedtempinteraction <- interact_plot(surviveyr2early_YN_MedPhy, pred="Med_AlphaPhy", modx = "suitabilitypc1V2abs",
                                              plot.points=TRUE, partial.residuals = TRUE, legend.main = "Climate dissimilarity (temp.)") +
  scale_y_continuous(trans='log10') +
  xlab("Median Alpha Phylogenetic Diversity") + ylab("Number of plants surviving ovewinter") +
  scale_linetype_manual(values=c("solid", "solid", "solid")) +
  theme_bw() + guides(linetype = "none") + theme(legend.position = "top")  +
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=10), 
        legend.text=element_text(size=8), 
        legend.title=element_text(size=10))
yr2earlyYNmedtempinteraction

yr2lateYNmedtempinteraction <- interact_plot(surviveyr2late_YN_MedPhy, pred="Med_AlphaPhy", modx = "suitabilitypc1V2abs",
                                             plot.points=TRUE, partial.residuals = TRUE, legend.main = "Climate dissimilarity (temp.)") +
  scale_y_continuous(trans='log10') +
  xlab("Median Alpha Phylogenetic Diversity") + ylab("Number of plants surviving to end of second year") +
  scale_linetype_manual(values=c("solid", "solid", "solid")) +
  theme_bw() + guides(linetype = "none") + theme(legend.position = "top")  +
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=10), 
        legend.text=element_text(size=8), 
        legend.title=element_text(size=10))
yr2lateYNmedtempinteraction

yr2latecountmedprecipinteraction <- interact_plot(surviveyr2late_count_MedPhy, pred="Med_AlphaPhy", modx = "suitabilitypc2V2abs",
                                                  plot.points=TRUE, partial.residuals = TRUE, legend.main = "Climate dissimilarity (precip.)") +
  scale_y_continuous(trans='log10') +
  xlab("Median Alpha Phylogenetic Diversity") + ylab("Number of plants surviving to end of second year") +
  scale_linetype_manual(values=c("solid", "solid", "solid")) +
  theme_bw() + guides(linetype = "none") + theme(legend.position = "top")  +
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=10), 
        legend.text=element_text(size=8), 
        legend.title=element_text(size=10))
yr2latecountmedprecipinteraction

medianalphafig <- ggarrange(germYNmedinteraction, germmedalphacountinteraction,
                            yr1medalphacountinteraction, yr2latecountmedprecipinteraction,
                            germcountmedtempinteraction, yr2earlyYNmedtempinteraction,
                            yr2lateYNmedtempinteraction, nrow=4, ncol=2)
medianalphafig
ggsave(medianalphafig, 
       filename = "Fig_S10.svg",
       height = 420, width = 210, units = "mm")

#GAMMAPHY - INTERACTIONS - Figure S11

yr1YNgammatempinteraction <- interact_plot(surviveyr1_YN_GammaPhy, pred="GammaPhy", modx = "suitabilitypc1V2abs",
                                           plot.points=TRUE, partial.residuals = TRUE, legend.main = "Climate dissimilarity (temp.)") +
  scale_y_continuous(trans='log10') +
  xlab("Gamma Phylogenetic Diversity") + ylab("Likelihood of any plant in plot surviving to end of first year") +
  scale_linetype_manual(values=c("solid", "solid", "solid")) +
  theme_bw() + guides(linetype = "none") + theme(legend.position = "top")  +
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=10), 
        legend.text=element_text(size=8), 
        legend.title=element_text(size=10))
yr1YNgammatempinteraction

yr2earlygammacountinteraction <- interact_plot(surviveyr2early_count_GammaPhy, pred="GammaPhy", modx = "Disturbance",
                                               plot.points=TRUE, partial.residuals = TRUE) +
  scale_y_continuous(trans='log10') + 
  xlab("Gamma Phylogenetic Diversity") + ylab("Number of plants surviving overwinter") +
  scale_color_manual(values=c("#999999", "#E69F00")) +
  scale_linetype_manual(values=c("solid", "solid")) +
  theme_bw() + guides(linetype = "none") + theme(legend.position = "top")  +
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=10), 
        legend.text=element_text(size=8), 
        legend.title=element_text(size=10))
yr2earlygammacountinteraction

yr2YNgammatempinteraction <- interact_plot(surviveyr2late_YN_GammaPhy, pred="GammaPhy", modx = "suitabilitypc1V2abs",
                                           plot.points=TRUE, partial.residuals = TRUE, legend.main = "Climate dissimilarity (temp.)") +
  scale_y_continuous(trans='log10') +
  xlab("Gamma Phylogenetic Diversity") + ylab("Likelihood of any plant in plot surviving to end of second year") +
  scale_linetype_manual(values=c("solid", "solid", "solid")) +
  theme_bw() + guides(linetype = "none") + theme(legend.position = "top")  +
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=10), 
        legend.text=element_text(size=8), 
        legend.title=element_text(size=10))
yr2YNgammatempinteraction

gammafig <- ggarrange(yr1YNgammatempinteraction, yr2earlygammacountinteraction, 
                      yr2YNgammatempinteraction, nrow=2, ncol=2)
gammafig
ggsave(gammafig, 
       filename = "Fig_S11.svg",
       height = 210, width = 210, units = "mm")

##################################################################################################

#QUESTION TWO: Do sown species with higher relative home-range-PD have higher rates of 
#colonisation and survival than species with lower relative home-range-PD?

##################################################################################################

haueserdata <- readRDS("haueserdata.RDS") %>%
  mutate(suitabilitypc1V2abs = abs(suitabilitypc1V2)) %>%
  mutate(suitabilitypc2V2abs = abs(suitabilitypc2V2))
#this contains data for germination, 1st year survival, 2nd year survival
#when analysing, will need to filter when analysing 1st and 2nd year data to only allow plants that survived previous stage
haueserdata$awmpd <- as.numeric(haueserdata$awmpd)
haueserdata$Species <- as.factor(haueserdata$Species)
haueserdata$Plot <- as.factor(haueserdata$Plot)
haueserdata$Disturbance <- relevel(haueserdata$Disturbance, ref="0")
str(haueserdata)

#models checked with 'DHARMa' and 'performance' packages 

####################################################################

#BINARY DATA

####################################################################

#For all of Q2, looked at both raw PD and PD-difference. The raw PD metrics are just for the supplementary, so the two models
#(with continuous variables scaled or not scaled) are just done for the '_diff' models

#GERMINATION YES/NO 

#Raw Mean AlphaPhy

germ_YN_Mean_AlphaPhy <- glmmTMB(germ_YN ~ Mean_AlphaPhy + Disturbance + suitabilitypc1V2abs + 
                                   suitabilitypc2V2abs + Heating + scale(awmpd) + OptGermRt +
                                   (1|Species) + (1|Plot), 
                                 family=binomial(link="logit"), 
                                 data=haueserdata)
summary(germ_YN_Mean_AlphaPhy)
Anova(germ_YN_Mean_AlphaPhy, type="III")
confint(germ_YN_Mean_AlphaPhy)
r.squaredGLMM(germ_YN_Mean_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = germ_YN_Mean_AlphaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(germ_YN_Mean_AlphaPhy)
check_singularity(germ_YN_Mean_AlphaPhy)

#Difference in Mean AlphaPhy 

germ_YN_Mean_AlphaPhy_diff <- glmmTMB(germ_YN ~ rescale(Mean_AlphaPhy_diff) + Disturbance + rescale(suitabilitypc1V2abs) +
                                        rescale(suitabilitypc2V2abs) + Heating + rescale(awmpd) + rescale(OptGermRt) +
                                        (1|Species) + (1|Plot), 
                                      family=binomial(link="logit"), 
                                      data=haueserdata)

germ_YN_Mean_AlphaPhy_diff <- glmmTMB(germ_YN ~ Mean_AlphaPhy_diff + Disturbance + suitabilitypc1V2abs +
                                        suitabilitypc2V2abs + Heating + scale(awmpd) + OptGermRt +
                                        (1|Species) + (1|Plot), 
                                      family=binomial(link="logit"), 
                                      data=haueserdata)
summary(germ_YN_Mean_AlphaPhy_diff)
Anova(germ_YN_Mean_AlphaPhy_diff, type="III")
confint(germ_YN_Mean_AlphaPhy_diff)
r.squaredGLMM(germ_YN_Mean_AlphaPhy_diff)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = germ_YN_Mean_AlphaPhy_diff, plot=F)
plot(simulationOutput)
check_overdispersion(germ_YN_Mean_AlphaPhy_diff)
check_singularity(germ_YN_Mean_AlphaPhy_diff)

#Test if an effect of native status

germ_YN_Mean_AlphaPhy_diffN <- glmmTMB(germ_YN ~ Mean_AlphaPhy_diff + Disturbance + suitabilitypc1V2abs +
                                         suitabilitypc2V2abs + Heating + scale(awmpd) + OptGermRt + Status +
                                         (1|Species) + (1|Plot), 
                                       family=binomial(link="logit"), 
                                       data=haueserdata)
summary(germ_YN_Mean_AlphaPhy_diffN)

AICtab(germ_YN_Mean_AlphaPhy_diffN, germ_YN_Mean_AlphaPhy_diff)

#Test effect of biogeographic origin

BG <- read.csv("biogeographic_ID.csv", header=T) %>% dplyr::select(plant_name_id, BG_ID)
BG$BG_ID <- as.factor(BG$BG_ID)

haueserdataBG <- left_join(haueserdata, BG, by="plant_name_id")

germ_YN_Mean_AlphaPhy_diffBG <- glmmTMB(germ_YN ~ Mean_AlphaPhy_diff + Disturbance + suitabilitypc1V2abs +
                                          suitabilitypc2V2abs + Heating + scale(awmpd) + OptGermRt +
                                          (1|Species) + (1|Plot) + (1|BG_ID), 
                                        family=binomial(link="logit"), 
                                        data=haueserdataBG)
summary(germ_YN_Mean_AlphaPhy_diffBG)

AICtab(germ_YN_Mean_AlphaPhy_diffBG, germ_YN_Mean_AlphaPhy_diff)

##

#Raw Max AlphaPhy  

germ_YN_Max_AlphaPhy <- glmmTMB(germ_YN ~ Max_AlphaPhy*Disturbance + suitabilitypc1V2abs + 
                                  suitabilitypc2V2abs + Heating + awmpd + OptGermRt +
                                  (1|Species) + (1|Plot), 
                                family=binomial(link="logit"), 
                                data=haueserdata)
summary(germ_YN_Max_AlphaPhy)
Anova(germ_YN_Max_AlphaPhy, type="III")
confint(germ_YN_Max_AlphaPhy)
r.squaredGLMM(germ_YN_Max_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = germ_YN_Max_AlphaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(germ_YN_Max_AlphaPhy)
check_singularity(germ_YN_Max_AlphaPhy)

#Difference in Max AlphaPhy 

germ_YN_Max_AlphaPhy_diff <- glmmTMB(germ_YN ~ rescale(Max_AlphaPhy_diff)*Disturbance + rescale(suitabilitypc1V2abs) +
                                       rescale(suitabilitypc2V2abs) + Heating + rescale(awmpd) + rescale(OptGermRt) +
                                       (1|Species) + (1|Plot), 
                                     family=binomial(link="logit"), 
                                     data=haueserdata)

germ_YN_Max_AlphaPhy_diff <- glmmTMB(germ_YN ~ Max_AlphaPhy_diff*Disturbance + suitabilitypc1V2abs +
                                       suitabilitypc2V2abs + Heating + awmpd + OptGermRt +
                                       (1|Species) + (1|Plot), 
                                     family=binomial(link="logit"), 
                                     data=haueserdata)
summary(germ_YN_Max_AlphaPhy_diff)
Anova(germ_YN_Max_AlphaPhy_diff, type="III")
confint(germ_YN_Max_AlphaPhy_diff)
r.squaredGLMM(germ_YN_Max_AlphaPhy_diff)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = germ_YN_Max_AlphaPhy_diff, plot=F)
plot(simulationOutput)
check_overdispersion(germ_YN_Max_AlphaPhy_diff)
check_singularity(germ_YN_Max_AlphaPhy_diff)

#Raw GammaPhy

germ_YN_GammaPhy <- glmmTMB(germ_YN ~ scale(GammaPhy)*Disturbance + suitabilitypc1V2abs +   
                              suitabilitypc2V2abs + Heating + awmpd + OptGermRt +
                              (1|Species) + (1|Plot), 
                            family=binomial(link="logit"), 
                            data=haueserdata)
summary(germ_YN_GammaPhy)
Anova(germ_YN_GammaPhy, type="III")
confint(germ_YN_GammaPhy)
r.squaredGLMM(germ_YN_GammaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = germ_YN_GammaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(germ_YN_GammaPhy)
check_singularity(germ_YN_GammaPhy)

#Difference in GammaPhy 

germ_YN_GammaPhy_diff <- glmmTMB(germ_YN ~ rescale(suitabilitypc1V2abs) + scale(GammaPhy_diff)*Disturbance + 
                                   rescale(suitabilitypc2V2abs) + Heating + rescale(awmpd) + rescale(OptGermRt) +
                                   (1|Species) + (1|Plot), 
                                 family=binomial(link="logit"), 
                                 data=haueserdata)

germ_YN_GammaPhy_diff <- glmmTMB(germ_YN ~ suitabilitypc1V2abs + scale(GammaPhy_diff)*Disturbance + 
                                   suitabilitypc2V2abs + Heating + awmpd + OptGermRt +
                                   (1|Species) + (1|Plot), 
                                 family=binomial(link="logit"), 
                                 data=haueserdata)
summary(germ_YN_GammaPhy_diff)
Anova(germ_YN_GammaPhy_diff, type="III")
confint(germ_YN_GammaPhy_diff)
r.squaredGLMM(germ_YN_GammaPhy_diff)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = germ_YN_GammaPhy_diff, plot=F)
plot(simulationOutput)
check_overdispersion(germ_YN_GammaPhy_diff)
check_singularity(germ_YN_GammaPhy_diff)

#Raw Median AlphaPhy

germ_YN_Med_AlphaPhy <- glmmTMB(germ_YN ~ Med_AlphaPhy*Disturbance + suitabilitypc1V2abs +   
                                  suitabilitypc2V2abs + Heating + awmpd + OptGermRt +
                                  (1|Species) + (1|Plot), 
                                family=binomial(link="logit"), 
                                data=haueserdata)
summary(germ_YN_Med_AlphaPhy)
Anova(germ_YN_Med_AlphaPhy, type="III")
confint(germ_YN_Med_AlphaPhy)
r.squaredGLMM(germ_YN_Med_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = germ_YN_Med_AlphaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(germ_YN_Med_AlphaPhy)
check_singularity(germ_YN_Med_AlphaPhy)

#Difference in Median AlphaPhy  

germ_YN_Med_AlphaPhy_diff <- glmmTMB(germ_YN ~ rescale(suitabilitypc2V2abs) + rescale(Med_AlphaPhy_diff)*Disturbance + 
                                       rescale(suitabilitypc1V2abs) + Heating + rescale(awmpd) + rescale(OptGermRt) +
                                       (1|Species) + (1|Plot), 
                                     family=binomial(link="logit"), 
                                     data=haueserdata)

germ_YN_Med_AlphaPhy_diff <- glmmTMB(germ_YN ~ suitabilitypc2V2abs + Med_AlphaPhy_diff*Disturbance + 
                                       suitabilitypc1V2abs + Heating + awmpd + OptGermRt +
                                       (1|Species) + (1|Plot), 
                                     family=binomial(link="logit"), 
                                     data=haueserdata)
summary(germ_YN_Med_AlphaPhy_diff)
Anova(germ_YN_Med_AlphaPhy_diff, type="III")
confint(germ_YN_Med_AlphaPhy_diff)
r.squaredGLMM(germ_YN_Med_AlphaPhy_diff)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = germ_YN_Med_AlphaPhy_diff, plot=F)
plot(simulationOutput)
check_overdispersion(germ_YN_Med_AlphaPhy_diff)
check_singularity(germ_YN_Med_AlphaPhy_diff)

##########################################################

#FIRST YEAR SURVIVAL

haueserdatafirstyear <- haueserdata %>% filter(germ_YN>0 | surviveyr1_YN>0)

#Raw Mean AlphaPhy

surviveyr1_YN_Mean_AlphaPhy <- glmmTMB(surviveyr1_YN ~ Mean_AlphaPhy + suitabilitypc2V2abs + Disturbance +
                                         suitabilitypc1V2abs + Heating + awmpd + OptGermRt +
                                         (1|Species) + (1|Plot), 
                                       family=binomial(link="logit"), 
                                       data=haueserdatafirstyear)
summary(surviveyr1_YN_Mean_AlphaPhy)
Anova(surviveyr1_YN_Mean_AlphaPhy, type="III")
confint(surviveyr1_YN_Mean_AlphaPhy)
r.squaredGLMM(surviveyr1_YN_Mean_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr1_YN_Mean_AlphaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr1_YN_Mean_AlphaPhy)
check_singularity(surviveyr1_YN_Mean_AlphaPhy)
#slight issue with singularity, but this likely due to the plot random variable, which is still 
#important to have in there

#Difference in Mean AlphaPhy 

surviveyr1_YN_Mean_AlphaPhy_diff <- glmmTMB(surviveyr1_YN ~ rescale(Mean_AlphaPhy_diff) + rescale(suitabilitypc1V2abs) + Disturbance + 
                                              rescale(suitabilitypc2V2abs) + Heating + rescale(awmpd) + rescale(OptGermRt) +
                                              (1|Species) + (1|Plot), 
                                            family=binomial(link="logit"), 
                                            data=haueserdatafirstyear)

surviveyr1_YN_Mean_AlphaPhy_diff <- glmmTMB(surviveyr1_YN ~ Mean_AlphaPhy_diff + suitabilitypc1V2abs + Disturbance + 
                                              suitabilitypc2V2abs + Heating + scale(awmpd) + OptGermRt +
                                              (1|Species) + (1|Plot), 
                                            family=binomial(link="logit"), 
                                            data=haueserdatafirstyear)
summary(surviveyr1_YN_Mean_AlphaPhy_diff)
Anova(surviveyr1_YN_Mean_AlphaPhy_diff, type="III")
confint(surviveyr1_YN_Mean_AlphaPhy_diff)
r.squaredGLMM(surviveyr1_YN_Mean_AlphaPhy_diff)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr1_YN_Mean_AlphaPhy_diff, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr1_YN_Mean_AlphaPhy_diff)
check_singularity(surviveyr1_YN_Mean_AlphaPhy_diff)
#slight issue with singularity, but this likely due to the plot random variable, which is still 
#important to have in there

#Testing whether native status has an effect

surviveyr1_YN_Mean_AlphaPhy_diffN <- glmmTMB(surviveyr1_YN ~ Mean_AlphaPhy_diff + suitabilitypc1V2abs + Disturbance + 
                                               suitabilitypc2V2abs + Heating + scale(awmpd) + OptGermRt + Status +
                                               (1|Species) + (1|Plot), 
                                             family=binomial(link="logit"), 
                                             data=haueserdatafirstyear)
summary(surviveyr1_YN_Mean_AlphaPhy_diffN)

AICtab(surviveyr1_YN_Mean_AlphaPhy_diffN, surviveyr1_YN_Mean_AlphaPhy_diff)

#Test effect of biogeographic origin

haueserdatafirstyearBG <- haueserdataBG %>% filter(germ_YN>0 | surviveyr1_YN>0)

surviveyr1_YN_Mean_AlphaPhy_diffBG <- glmmTMB(surviveyr1_YN ~ Mean_AlphaPhy_diff + suitabilitypc1V2abs + Disturbance + 
                                                suitabilitypc2V2abs + Heating + scale(awmpd) + OptGermRt +
                                                (1|Species) + (1|Plot) + (1|BG_ID), 
                                              family=binomial(link="logit"), 
                                              data=haueserdatafirstyearBG)
summary(surviveyr1_YN_Mean_AlphaPhy_diffBG)

AICtab(surviveyr1_YN_Mean_AlphaPhy_diffBG, surviveyr1_YN_Mean_AlphaPhy_diff)

##

#Raw Max AlphaPhy 

surviveyr1_YN_Max_AlphaPhy <- glmmTMB(surviveyr1_YN ~ Max_AlphaPhy*Disturbance + suitabilitypc1V2abs + 
                                        suitabilitypc2V2abs + Heating + scale(awmpd) + OptGermRt +
                                        (1|Species) + (1|Plot), 
                                      family=binomial(link="logit"), 
                                      data=haueserdatafirstyear)
summary(surviveyr1_YN_Max_AlphaPhy)
Anova(surviveyr1_YN_Max_AlphaPhy, type="III")
confint(surviveyr1_YN_Max_AlphaPhy)
r.squaredGLMM(surviveyr1_YN_Max_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr1_YN_Max_AlphaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr1_YN_Max_AlphaPhy)
check_singularity(surviveyr1_YN_Max_AlphaPhy)
#slight issue with singularity, but this likely due to the plot random variable, which is still 
#important to have in there

#Difference in Max AlphaPhy 

surviveyr1_YN_Max_AlphaPhy_diff <- glmmTMB(surviveyr1_YN ~ rescale(Max_AlphaPhy_diff)*Disturbance + rescale(suitabilitypc2V2abs) +    
                                             rescale(suitabilitypc1V2abs) + Heating + rescale(awmpd) + rescale(OptGermRt) +
                                             (1|Species) + (1|Plot), 
                                           family=binomial(link="logit"), 
                                           data=haueserdatafirstyear)

surviveyr1_YN_Max_AlphaPhy_diff <- glmmTMB(surviveyr1_YN ~ Max_AlphaPhy_diff*Disturbance + suitabilitypc2V2abs +    
                                             suitabilitypc1V2abs + Heating + scale(awmpd) + OptGermRt +
                                             (1|Species) + (1|Plot), 
                                           family=binomial(link="logit"), 
                                           data=haueserdatafirstyear)
summary(surviveyr1_YN_Max_AlphaPhy_diff)
Anova(surviveyr1_YN_Max_AlphaPhy_diff, type="III")
confint(surviveyr1_YN_Max_AlphaPhy_diff)
r.squaredGLMM(surviveyr1_YN_Max_AlphaPhy_diff)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr1_YN_Max_AlphaPhy_diff, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr1_YN_Max_AlphaPhy_diff)
check_singularity(surviveyr1_YN_Max_AlphaPhy_diff)
#slight issue with singularity, but this likely due to the plot random variable, which is still 
#important to have in there

#Raw GammaPhy 

surviveyr1_YN_GammaPhy <- glmmTMB(surviveyr1_YN ~ scale(GammaPhy) + suitabilitypc2V2abs +  Disturbance + 
                                    suitabilitypc1V2abs + Heating + awmpd + OptGermRt +
                                    (1|Species) + (1|Plot), 
                                  family=binomial(link="logit"), 
                                  data=haueserdatafirstyear)
summary(surviveyr1_YN_GammaPhy)
Anova(surviveyr1_YN_GammaPhy, type="III")
confint(surviveyr1_YN_GammaPhy)
r.squaredGLMM(surviveyr1_YN_GammaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr1_YN_GammaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr1_YN_GammaPhy)
check_singularity(surviveyr1_YN_GammaPhy)
#slight issue with singularity, but this likely due to the plot random variable, which is still 
#important to have in there

#Difference in GammaPhy 

surviveyr1_YN_GammaPhy_diff <- glmmTMB(surviveyr1_YN ~ scale(GammaPhy_diff) + rescale(suitabilitypc2V2abs) + Disturbance +   
                                         rescale(suitabilitypc1V2abs) + Heating + rescale(awmpd) + rescale(OptGermRt) +
                                         (1|Species) + (1|Plot), 
                                       family=binomial(link="logit"), 
                                       data=haueserdatafirstyear)

surviveyr1_YN_GammaPhy_diff <- glmmTMB(surviveyr1_YN ~ scale(GammaPhy_diff) + suitabilitypc2V2abs + Disturbance +   
                                         suitabilitypc1V2abs + Heating + scale(awmpd) + OptGermRt +
                                         (1|Species) + (1|Plot), 
                                       family=binomial(link="logit"), 
                                       data=haueserdatafirstyear)
summary(surviveyr1_YN_GammaPhy_diff)
Anova(surviveyr1_YN_GammaPhy_diff, type="III")
confint(surviveyr1_YN_GammaPhy_diff)
r.squaredGLMM(surviveyr1_YN_GammaPhy_diff)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr1_YN_GammaPhy_diff, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr1_YN_GammaPhy_diff)
check_singularity(surviveyr1_YN_GammaPhy_diff)
#slight issue with singularity, but this likely due to the plot random variable, which is still 
#important to have in there

#Raw Median AlphaPhy 

surviveyr1_YN_Med_AlphaPhy <- glmmTMB(surviveyr1_YN ~ Med_AlphaPhy*Disturbance +  suitabilitypc1V2abs +  
                                        suitabilitypc2V2abs + Heating + awmpd + OptGermRt +
                                        (1|Species) + (1|Plot), 
                                      family=binomial(link="logit"), 
                                      data=haueserdatafirstyear)
summary(surviveyr1_YN_Med_AlphaPhy)
Anova(surviveyr1_YN_Med_AlphaPhy, type="III")
confint(surviveyr1_YN_Med_AlphaPhy)
r.squaredGLMM(surviveyr1_YN_Med_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr1_YN_Med_AlphaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr1_YN_Med_AlphaPhy)
check_singularity(surviveyr1_YN_Med_AlphaPhy)
#slight issue with singularity, but this likely due to the plot random variable, which is still 
#important to have in there

#Difference in Median AlphaPhy 

surviveyr1_YN_Med_AlphaPhy_diff <- glmmTMB(surviveyr1_YN ~ rescale(Med_AlphaPhy_diff)*Disturbance + rescale(suitabilitypc1V2abs) +    
                                             rescale(suitabilitypc2V2abs) + Heating + rescale(awmpd) + rescale(OptGermRt) +
                                             (1|Species) + (1|Plot), 
                                           family=binomial(link="logit"), 
                                           data=haueserdatafirstyear)

surviveyr1_YN_Med_AlphaPhy_diff <- glmmTMB(surviveyr1_YN ~ Med_AlphaPhy_diff*Disturbance + suitabilitypc1V2abs +    
                                             suitabilitypc2V2abs + Heating + awmpd + OptGermRt +
                                             (1|Species) + (1|Plot), 
                                           family=binomial(link="logit"), 
                                           data=haueserdatafirstyear)
summary(surviveyr1_YN_Med_AlphaPhy_diff)
Anova(surviveyr1_YN_Med_AlphaPhy_diff, type="III")
confint(surviveyr1_YN_Med_AlphaPhy_diff)
r.squaredGLMM(surviveyr1_YN_Med_AlphaPhy_diff)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr1_YN_Med_AlphaPhy_diff, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr1_YN_Med_AlphaPhy_diff)
check_singularity(surviveyr1_YN_Med_AlphaPhy_diff)
#slight issue with singularity, but this likely due to the plot random variable, which is still 
#important to have in there

##########################################################

#OVERWINTER SURVIVAL

haueserdataoverwinter <- haueserdata %>% filter(surviveyr1_YN>0 | surviveyr2early_YN > 0)

#Raw Mean AlphaPhy 

surviveyr2early_YN_Mean_AlphaPhy <- glmmTMB(surviveyr2early_YN ~ Mean_AlphaPhy*suitabilitypc1V2abs + Disturbance + 
                                              suitabilitypc2V2abs + Heating + awmpd + OptGermRt +
                                              (1|Species) + (1|Plot), 
                                            family=binomial(link="logit"), 
                                            data=haueserdataoverwinter)
summary(surviveyr2early_YN_Mean_AlphaPhy)
Anova(surviveyr2early_YN_Mean_AlphaPhy, type="III")
confint(surviveyr2early_YN_Mean_AlphaPhy)
r.squaredGLMM(surviveyr2early_YN_Mean_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2early_YN_Mean_AlphaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr2early_YN_Mean_AlphaPhy)
check_singularity(surviveyr2early_YN_Mean_AlphaPhy)
#slight issue with singularity, but this likely due to the plot random variable, which is still 
#important to have in there

#Difference in Mean AlphaPhy 

surviveyr2early_YN_Mean_AlphaPhy_diff <- glmmTMB(surviveyr2early_YN ~ rescale(Mean_AlphaPhy_diff)*rescale(suitabilitypc1V2abs) + Disturbance + 
                                                   rescale(suitabilitypc2V2abs) + Heating + rescale(awmpd) + rescale(OptGermRt) +
                                                   (1|Species) + (1|Plot), 
                                                 family=binomial(link="logit"), 
                                                 data=haueserdataoverwinter)

surviveyr2early_YN_Mean_AlphaPhy_diff <- glmmTMB(surviveyr2early_YN ~ Mean_AlphaPhy_diff*suitabilitypc1V2abs + Disturbance + 
                                                   suitabilitypc2V2abs + Heating + scale(awmpd) + OptGermRt +
                                                   (1|Species) + (1|Plot), 
                                                 family=binomial(link="logit"), 
                                                 data=haueserdataoverwinter)
summary(surviveyr2early_YN_Mean_AlphaPhy_diff)
Anova(surviveyr2early_YN_Mean_AlphaPhy_diff, type="III")
confint(surviveyr2early_YN_Mean_AlphaPhy_diff)
r.squaredGLMM(surviveyr2early_YN_Mean_AlphaPhy_diff)

simulationOutput <- simulateResiduals(fittedModel = surviveyr2early_YN_Mean_AlphaPhy_diff, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr2early_YN_Mean_AlphaPhy_diff)
check_singularity(surviveyr2early_YN_Mean_AlphaPhy_diff)

###

#Testing whether native status has an effect

surviveyr2early_YN_Mean_AlphaPhy_diffN <- glmmTMB(surviveyr2early_YN ~ Mean_AlphaPhy_diff*suitabilitypc1V2abs + Disturbance + 
                                                    suitabilitypc2V2abs + Heating + scale(awmpd) + OptGermRt + Status +
                                                    (1|Species) + (1|Plot), 
                                                  family=binomial(link="logit"), 
                                                  data=haueserdataoverwinter)
summary(surviveyr2early_YN_Mean_AlphaPhy_diffN)

AICtab(surviveyr2early_YN_Mean_AlphaPhy_diffN, surviveyr2early_YN_Mean_AlphaPhy_diff)

#Test effect of biogeographic origin

haueserdataoverwinterBG <- haueserdataBG %>% filter(surviveyr1_YN>0 | surviveyr2early_YN > 0)

surviveyr2early_YN_Mean_AlphaPhy_diffBG <- glmmTMB(surviveyr2early_YN ~ Mean_AlphaPhy_diff*suitabilitypc1V2abs + Disturbance + 
                                                     suitabilitypc2V2abs + Heating + scale(awmpd) + OptGermRt +
                                                     (1|Species) + (1|Plot) + (1|BG_ID), 
                                                   family=binomial(link="logit"), 
                                                   data=haueserdataoverwinterBG)
summary(surviveyr2early_YN_Mean_AlphaPhy_diffBG)

AICtab(surviveyr2early_YN_Mean_AlphaPhy_diffBG, surviveyr2early_YN_Mean_AlphaPhy_diff)

#Segue: checking climate suitability - showing why we use absolute values; performamce declines when gets hotter OR colder / wetter OR drier

surviveyr2early_YN_Mean_AlphaPhy_diff <- glmmTMB(surviveyr2early_YN ~ Mean_AlphaPhy_diff + suitabilitypc2V2 + Disturbance + 
                                                   suitabilitypc1V2 + Heating + scale(awmpd) + OptGermRt +
                                                   (1|Species) + (1|Plot), 
                                                 family=binomial(link="logit"), 
                                                 data=haueserdataoverwinter)
summary(surviveyr2early_YN_Mean_AlphaPhy_diff)

surviveyr2early_YN_Mean_AlphaPhy_diffsq <- glmmTMB(surviveyr2early_YN ~ Mean_AlphaPhy_diff + I(suitabilitypc2V2^2) + Disturbance + 
                                                     suitabilitypc1V2 + Heating + scale(awmpd) + OptGermRt +
                                                     (1|Species) + (1|Plot), 
                                                   family=binomial(link="logit"), 
                                                   data=haueserdataoverwinter)
summary(surviveyr2early_YN_Mean_AlphaPhy_diffsq)

surviveyr2early_YN_Mean_AlphaPhy_diffsq2 <- glmmTMB(surviveyr2early_YN ~ Mean_AlphaPhy_diff + I(suitabilitypc1V2^2) + Disturbance + 
                                                      suitabilitypc2V2 + Heating + scale(awmpd) + OptGermRt +
                                                      (1|Species) + (1|Plot), 
                                                    family=binomial(link="logit"), 
                                                    data=haueserdataoverwinter)
summary(surviveyr2early_YN_Mean_AlphaPhy_diffsq2)

precipplot <- visreg(surviveyr2early_YN_Mean_AlphaPhy_diffsq, 
                     "suitabilitypc2V2", gg=TRUE, xlab="Untransformed precipitation dissimilarity",
                     ylab="Scaled likelihood of overwinter survival", line=list(col="black"),
                     fill=list(fill=rgb(red=0.8, green=0.8, blue=0.8, alpha=0.6)),
                     points=list(size=2, pch=16, col="black")) + theme_bw() + 
  theme(axis.text=element_text(size=7), axis.title=element_text(size=9),
        legend.text=element_text(size=7), legend.title=element_text(size=9))
precipplot

tempplot <- visreg(surviveyr2early_YN_Mean_AlphaPhy_diffsq2, 
                   "suitabilitypc1V2", gg=TRUE, xlab="Untransformed temperature dissimilarity",
                   ylab="Scaled likelihood of overwinter survival", line=list(col="black"),
                   fill=list(fill=rgb(red=0.8, green=0.8, blue=0.8, alpha=0.6)),
                   points=list(size=2, pch=16, col="black")) + theme_bw() + 
  theme(axis.text=element_text(size=7), axis.title=element_text(size=9),
        legend.text=element_text(size=7), legend.title=element_text(size=9))
tempplot

AICtab(surviveyr2early_YN_Mean_AlphaPhy_diff, surviveyr2early_YN_Mean_AlphaPhy_diffsq)
AICtab(surviveyr2early_YN_Mean_AlphaPhy_diff, surviveyr2early_YN_Mean_AlphaPhy_diffsq2)

###

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2early_YN_Mean_AlphaPhy_diff, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr2early_YN_Mean_AlphaPhy_diff)
check_singularity(surviveyr2early_YN_Mean_AlphaPhy_diff)
#slight issue with singularity, but this likely due to the plot random variable, which is still 
#important to have in there

##

#Raw Max AlphaPhy 

surviveyr2early_YN_Max_AlphaPhy <- glmmTMB(surviveyr2early_YN ~ Max_AlphaPhy + suitabilitypc1V2abs + Disturbance + 
                                             suitabilitypc2V2abs + Heating + awmpd + OptGermRt +
                                             (1|Species) + (1|Plot), 
                                           family=binomial(link="logit"), 
                                           data=haueserdataoverwinter)
summary(surviveyr2early_YN_Max_AlphaPhy)
Anova(surviveyr2early_YN_Max_AlphaPhy, type="III")
confint(surviveyr2early_YN_Max_AlphaPhy)
r.squaredGLMM(surviveyr2early_YN_Max_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2early_YN_Max_AlphaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr2early_YN_Max_AlphaPhy)
check_singularity(surviveyr2early_YN_Max_AlphaPhy)
#slight issue with singularity, but this likely due to the plot random variable, which is still 
#important to have in there

#Difference in Max AlphaPhy 

surviveyr2early_YN_Max_AlphaPhy_diff <- glmmTMB(surviveyr2early_YN ~ rescale(Max_AlphaPhy_diff) + rescale(suitabilitypc1V2abs) + Disturbance + 
                                                  rescale(suitabilitypc2V2abs) + Heating + rescale(awmpd) + rescale(OptGermRt) +
                                                  (1|Species) + (1|Plot), 
                                                family=binomial(link="logit"), 
                                                data=haueserdataoverwinter)

surviveyr2early_YN_Max_AlphaPhy_diff <- glmmTMB(surviveyr2early_YN ~ Max_AlphaPhy_diff + suitabilitypc1V2abs + Disturbance + 
                                                  suitabilitypc2V2abs + Heating + scale(awmpd) + OptGermRt +
                                                  (1|Species) + (1|Plot), 
                                                family=binomial(link="logit"), 
                                                data=haueserdataoverwinter)
summary(surviveyr2early_YN_Max_AlphaPhy_diff)
Anova(surviveyr2early_YN_Max_AlphaPhy_diff, type="III")
confint(surviveyr2early_YN_Max_AlphaPhy_diff)
r.squaredGLMM(surviveyr2early_YN_Max_AlphaPhy_diff)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2early_YN_Max_AlphaPhy_diff, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr2early_YN_Max_AlphaPhy_diff)
check_singularity(surviveyr2early_YN_Max_AlphaPhy_diff)
#slight issue with singularity, but this likely due to the plot random variable, which is still 
#important to have in there

#Raw GammaPhy  

surviveyr2early_YN_GammaPhy <- glmmTMB(surviveyr2early_YN ~ scale(GammaPhy)*suitabilitypc2V2abs + Disturbance + 
                                         suitabilitypc1V2abs + Heating + awmpd + OptGermRt +
                                         (1|Species) + (1|Plot), 
                                       family=binomial(link="logit"), 
                                       data=haueserdataoverwinter)
summary(surviveyr2early_YN_GammaPhy)
Anova(surviveyr2early_YN_GammaPhy, type="III")
confint(surviveyr2early_YN_GammaPhy)
r.squaredGLMM(surviveyr2early_YN_GammaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2early_YN_GammaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr2early_YN_GammaPhy)
check_singularity(surviveyr2early_YN_GammaPhy)
#slight issue with singularity, but this likely due to the plot random variable, which is still 
#important to have in there

#Difference in GammaPhy

surviveyr2early_YN_GammaPhy_diff <- glmmTMB(surviveyr2early_YN ~ scale(GammaPhy_diff)*rescale(suitabilitypc2V2abs) + Disturbance + 
                                              rescale(suitabilitypc1V2abs) + Heating + rescale(awmpd) + rescale(OptGermRt) +
                                              (1|Species) + (1|Plot), 
                                            family=binomial(link="logit"), 
                                            data=haueserdataoverwinter)

surviveyr2early_YN_GammaPhy_diff <- glmmTMB(surviveyr2early_YN ~ scale(GammaPhy_diff)*suitabilitypc2V2abs +  Disturbance + 
                                              suitabilitypc1V2abs + Heating + scale(awmpd) + OptGermRt +
                                              (1|Species) + (1|Plot), 
                                            family=binomial(link="logit"), 
                                            data=haueserdataoverwinter)
summary(surviveyr2early_YN_GammaPhy_diff)
Anova(surviveyr2early_YN_GammaPhy_diff, type="III")
confint(surviveyr2early_YN_GammaPhy_diff)
r.squaredGLMM(surviveyr2early_YN_GammaPhy_diff)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2early_YN_GammaPhy_diff, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr2early_YN_GammaPhy_diff)
check_singularity(surviveyr2early_YN_GammaPhy_diff)
#slight issue with singularity, but this likely due to the plot random variable, which is still 
#important to have in there

#Raw Median AlphaPhy 

surviveyr2early_YN_Med_AlphaPhy <- glmmTMB(surviveyr2early_YN ~ Med_AlphaPhy + suitabilitypc1V2abs +  Disturbance + 
                                             suitabilitypc2V2abs + Heating + awmpd + OptGermRt +
                                             (1|Species) + (1|Plot), 
                                           family=binomial(link="logit"), 
                                           data=haueserdataoverwinter)
summary(surviveyr2early_YN_Med_AlphaPhy)
Anova(surviveyr2early_YN_Med_AlphaPhy, type="III")
confint(surviveyr2early_YN_Med_AlphaPhy)
r.squaredGLMM(surviveyr2early_YN_Med_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2early_YN_Med_AlphaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr2early_YN_Med_AlphaPhy)
check_singularity(surviveyr2early_YN_Med_AlphaPhy)
#slight issue with singularity, but this likely due to the plot random variable, which is still 
#important to have in there

#Difference in Median AlphaPhy 

surviveyr2early_YN_Med_AlphaPhy_diff <- glmmTMB(surviveyr2early_YN ~ rescale(Med_AlphaPhy_diff)*rescale(suitabilitypc1V2abs) + Disturbance +  
                                                  rescale(suitabilitypc2V2abs) + Heating + rescale(awmpd) + rescale(OptGermRt) +
                                                  (1|Species) + (1|Plot), 
                                                family=binomial(link="logit"), 
                                                data=haueserdataoverwinter)

surviveyr2early_YN_Med_AlphaPhy_diff <- glmmTMB(surviveyr2early_YN ~ Med_AlphaPhy_diff*suitabilitypc1V2abs + Disturbance +  
                                                  suitabilitypc2V2abs + Heating + scale(awmpd) + OptGermRt +
                                                  (1|Species) + (1|Plot), 
                                                family=binomial(link="logit"), 
                                                data=haueserdataoverwinter)
summary(surviveyr2early_YN_Med_AlphaPhy_diff)
Anova(surviveyr2early_YN_Med_AlphaPhy_diff, type="III")
confint(surviveyr2early_YN_Med_AlphaPhy_diff)
r.squaredGLMM(surviveyr2early_YN_Med_AlphaPhy_diff)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2early_YN_Med_AlphaPhy_diff, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr2early_YN_Med_AlphaPhy_diff)
check_singularity(surviveyr2early_YN_Med_AlphaPhy_diff)
#slight issue with singularity, but this likely due to the plot random variable, which is still 
#important to have in there

##########################################################

#LATE SECOND YEAR SURVIVAL

haueserdatayr2late <- haueserdata %>% filter(surviveyr2early_YN > 0 | surviveyr2late_YN>0)

#Raw Mean AlphaPhy 

surviveyr2late_YN_Mean_AlphaPhy <- glmmTMB(surviveyr2late_YN ~ Mean_AlphaPhy + Disturbance + suitabilitypc1V2abs +
                                             suitabilitypc2V2abs + Heating + awmpd + OptGermRt +
                                             (1|Species) + (1|Plot), 
                                           family=binomial(link="logit"), 
                                           data=haueserdatayr2late)
summary(surviveyr2late_YN_Mean_AlphaPhy)
Anova(surviveyr2late_YN_Mean_AlphaPhy, type="III")
confint(surviveyr2late_YN_Mean_AlphaPhy)
r.squaredGLMM(surviveyr2late_YN_Mean_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2late_YN_Mean_AlphaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr2late_YN_Mean_AlphaPhy)
check_singularity(surviveyr2late_YN_Mean_AlphaPhy)

#Difference in Mean AlphaPhy 

surviveyr2late_YN_Mean_AlphaPhy_diff <- glmmTMB(surviveyr2late_YN ~ rescale(Mean_AlphaPhy_diff)*Disturbance + rescale(suitabilitypc1V2abs) +
                                                  rescale(suitabilitypc2V2abs) + Heating + rescale(awmpd) + rescale(OptGermRt) +
                                                  (1|Species) + (1|Plot), 
                                                family=binomial(link="logit"), 
                                                data=haueserdatayr2late)

surviveyr2late_YN_Mean_AlphaPhy_diff <- glmmTMB(surviveyr2late_YN ~ Mean_AlphaPhy_diff*Disturbance + suitabilitypc1V2abs +
                                                  suitabilitypc2V2abs + Heating + scale(awmpd) + OptGermRt +
                                                  (1|Species) + (1|Plot), 
                                                family=binomial(link="logit"), 
                                                data=haueserdatayr2late)
summary(surviveyr2late_YN_Mean_AlphaPhy_diff)
Anova(surviveyr2late_YN_Mean_AlphaPhy_diff, type="III")
confint(surviveyr2late_YN_Mean_AlphaPhy_diff)
r.squaredGLMM(surviveyr2late_YN_Mean_AlphaPhy_diff)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2late_YN_Mean_AlphaPhy_diff, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr2late_YN_Mean_AlphaPhy_diff)
check_singularity(surviveyr2late_YN_Mean_AlphaPhy_diff)

#Testing whether native status has an effect

surviveyr2late_YN_Mean_AlphaPhy_diffN <- glmmTMB(surviveyr2late_YN ~ Mean_AlphaPhy_diff*Disturbance + suitabilitypc1V2abs +
                                                   suitabilitypc2V2abs + Heating + scale(awmpd) + OptGermRt + Status +
                                                   (1|Species) + (1|Plot), 
                                                 family=binomial(link="logit"), 
                                                 data=haueserdatayr2late)
summary(surviveyr2late_YN_Mean_AlphaPhy_diffN)

AICtab(surviveyr2late_YN_Mean_AlphaPhy_diffN, surviveyr2late_YN_Mean_AlphaPhy_diff)

#Test effect of biogeographic origin

haueserdatayr2lateBG <- haueserdataBG %>% filter(surviveyr2early_YN > 0 | surviveyr2late_YN>0)

surviveyr2late_YN_Mean_AlphaPhy_diffBG <- glmmTMB(surviveyr2late_YN ~ Mean_AlphaPhy_diff*Disturbance + suitabilitypc1V2abs +
                                                    suitabilitypc2V2abs + Heating + scale(awmpd) + OptGermRt +
                                                    (1|Species) + (1|Plot) + (1|BG_ID), 
                                                  family=binomial(link="logit"), 
                                                  data=haueserdatayr2lateBG)
summary(surviveyr2late_YN_Mean_AlphaPhy_diffBG)

AICtab(surviveyr2late_YN_Mean_AlphaPhy_diffBG, surviveyr2late_YN_Mean_AlphaPhy_diff)

#Segue: checking climate suitability

surviveyr2late_YN_Mean_AlphaPhy_diff <- glmmTMB(surviveyr2late_YN ~ Mean_AlphaPhy_diff + suitabilitypc2V2 + Disturbance + 
                                                  suitabilitypc1V2 + Heating + scale(awmpd) + OptGermRt +
                                                  (1|Species) + (1|Plot), 
                                                family=binomial(link="logit"), 
                                                data=haueserdatayr2late)
summary(surviveyr2late_YN_Mean_AlphaPhy_diff)

surviveyr2late_YN_Mean_AlphaPhy_diffsq <- glmmTMB(surviveyr2late_YN ~ Mean_AlphaPhy_diff + I(suitabilitypc2V2^2) + Disturbance + 
                                                    suitabilitypc1V2 + Heating + scale(awmpd) + OptGermRt +
                                                    (1|Species) + (1|Plot), 
                                                  family=binomial(link="logit"), 
                                                  data=haueserdatayr2late)
summary(surviveyr2late_YN_Mean_AlphaPhy_diffsq)

surviveyr2late_YN_Mean_AlphaPhy_diffsq2 <- glmmTMB(surviveyr2late_YN ~ Mean_AlphaPhy_diff + I(suitabilitypc1V2^2) + Disturbance + 
                                                     suitabilitypc2V2 + Heating + scale(awmpd) + OptGermRt +
                                                     (1|Species) + (1|Plot), 
                                                   family=binomial(link="logit"), 
                                                   data=haueserdatayr2late)
summary(surviveyr2late_YN_Mean_AlphaPhy_diffsq2)

precipplotlate <- visreg(surviveyr2late_YN_Mean_AlphaPhy_diffsq, 
                         "suitabilitypc2V2", gg=TRUE, xlab="Untransformed precipitation dissimilarity",
                         ylab="Scaled likelihood of survival to end of second growing season", line=list(col="black"),
                         fill=list(fill=rgb(red=0.8, green=0.8, blue=0.8, alpha=0.6)),
                         points=list(size=2, pch=16, col="black")) + theme_bw() +
  theme(axis.text=element_text(size=7),
        axis.title=element_text(size=9), 
        legend.text=element_text(size=7), 
        legend.title=element_text(size=9))
precipplotlate

tempplotlate <- visreg(surviveyr2late_YN_Mean_AlphaPhy_diffsq2, 
                       "suitabilitypc1V2", gg=TRUE, xlab="Untransformed temperature dissimilarity",
                       ylab="Scaled likelihood of survival to end of second growing season", line=list(col="black"),
                       fill=list(fill=rgb(red=0.8, green=0.8, blue=0.8, alpha=0.6)),
                       points=list(size=2, pch=16, col="black")) + theme_bw() +
  theme(axis.text=element_text(size=7),
        axis.title=element_text(size=9), 
        legend.text=element_text(size=7), 
        legend.title=element_text(size=9))
tempplotlate

AICtab(surviveyr2late_YN_Mean_AlphaPhy_diff, surviveyr2late_YN_Mean_AlphaPhy_diffsq)
AICtab(surviveyr2late_YN_Mean_AlphaPhy_diff, surviveyr2late_YN_Mean_AlphaPhy_diffsq2)

climplottotalrevised <- ggarrange(precipplot, tempplot, precipplotlate, tempplotlate,
                                  nrow=2, ncol=2)
climplottotalrevised
ggsave(climplottotalrevised, 
       filename = "Fig_S17.svg",
       height = 210, width = 210, units = "mm")
##

#Raw Max AlphaPhy

surviveyr2late_YN_Max_AlphaPhy <- glmmTMB(surviveyr2late_YN ~ Max_AlphaPhy + Disturbance + suitabilitypc1V2abs + 
                                            suitabilitypc2V2abs + Heating + awmpd + OptGermRt +
                                            (1|Species) + (1|Plot), 
                                          family=binomial(link="logit"), 
                                          data=haueserdatayr2late)
summary(surviveyr2late_YN_Max_AlphaPhy)
Anova(surviveyr2late_YN_Max_AlphaPhy, type="III")
confint(surviveyr2late_YN_Max_AlphaPhy)
r.squaredGLMM(surviveyr2late_YN_Max_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2late_YN_Max_AlphaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr2late_YN_Max_AlphaPhy)
check_singularity(surviveyr2late_YN_Max_AlphaPhy)

#Difference in Max AlphaPhy 

surviveyr2late_YN_Max_AlphaPhy_diff <- glmmTMB(surviveyr2late_YN ~ rescale(Max_AlphaPhy_diff) + Disturbance + rescale(suitabilitypc1V2abs) +  
                                                 rescale(suitabilitypc2V2abs) + Heating + rescale(awmpd) + rescale(OptGermRt) +
                                                 (1|Species) + (1|Plot), 
                                               family=binomial(link="logit"), 
                                               data=haueserdatayr2late)

surviveyr2late_YN_Max_AlphaPhy_diff <- glmmTMB(surviveyr2late_YN ~ Max_AlphaPhy_diff + Disturbance + suitabilitypc1V2abs +  
                                                 suitabilitypc2V2abs + Heating + scale(awmpd) + OptGermRt +
                                                 (1|Species) + (1|Plot), 
                                               family=binomial(link="logit"), 
                                               data=haueserdatayr2late)
summary(surviveyr2late_YN_Max_AlphaPhy_diff)
Anova(surviveyr2late_YN_Max_AlphaPhy_diff, type="III")
confint(surviveyr2late_YN_Max_AlphaPhy_diff)
r.squaredGLMM(surviveyr2late_YN_Max_AlphaPhy_diff)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2late_YN_Max_AlphaPhy_diff, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr2late_YN_Max_AlphaPhy_diff)
check_singularity(surviveyr2late_YN_Max_AlphaPhy_diff)

#Raw GammaPhy #no interactions 

surviveyr2late_YN_GammaPhy <- glmmTMB(surviveyr2late_YN ~ scale(GammaPhy) + Disturbance + suitabilitypc1V2abs +
                                        suitabilitypc2V2abs + Heating + awmpd + OptGermRt +
                                        (1|Species) + (1|Plot), 
                                      family=binomial(link="logit"), 
                                      data=haueserdatayr2late)
summary(surviveyr2late_YN_GammaPhy)
Anova(surviveyr2late_YN_GammaPhy, type="III")
confint(surviveyr2late_YN_GammaPhy)
r.squaredGLMM(surviveyr2late_YN_GammaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2late_YN_GammaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr2late_YN_GammaPhy)
check_singularity(surviveyr2late_YN_GammaPhy)

#Difference in GammaPhy 

surviveyr2late_YN_GammaPhy_diff <- glmmTMB(surviveyr2late_YN ~ scale(GammaPhy_diff) + rescale(Disturbance) + rescale(suitabilitypc1V2abs) + 
                                             rescale(suitabilitypc2V2abs) + Heating + rescale(awmpd) + rescale(OptGermRt) +
                                             (1|Species) + (1|Plot), 
                                           family=binomial(link="logit"), 
                                           data=haueserdatayr2late)

surviveyr2late_YN_GammaPhy_diff <- glmmTMB(surviveyr2late_YN ~ scale(GammaPhy_diff) +  Disturbance + suitabilitypc1V2abs + 
                                             suitabilitypc2V2abs + Heating + scale(awmpd) + OptGermRt +
                                             (1|Species) + (1|Plot), 
                                           family=binomial(link="logit"), 
                                           data=haueserdatayr2late)
summary(surviveyr2late_YN_GammaPhy_diff)
Anova(surviveyr2late_YN_GammaPhy_diff, type="III")
confint(surviveyr2late_YN_GammaPhy_diff)
r.squaredGLMM(surviveyr2late_YN_GammaPhy_diff)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2late_YN_GammaPhy_diff, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr2late_YN_GammaPhy_diff)
check_singularity(surviveyr2late_YN_GammaPhy_diff)

#Raw Median AlphaPhy  

surviveyr2late_YN_Med_AlphaPhy <- glmmTMB(surviveyr2late_YN ~ Med_AlphaPhy*Disturbance + suitabilitypc1V2abs +
                                            suitabilitypc2V2abs + Heating + awmpd + OptGermRt +
                                            (1|Species) + (1|Plot), 
                                          family=binomial(link="logit"), 
                                          data=haueserdatayr2late)
summary(surviveyr2late_YN_Med_AlphaPhy)
Anova(surviveyr2late_YN_Med_AlphaPhy, type="III")
confint(surviveyr2late_YN_Med_AlphaPhy)
r.squaredGLMM(surviveyr2late_YN_Med_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2late_YN_Med_AlphaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr2late_YN_Med_AlphaPhy)
check_singularity(surviveyr2late_YN_Med_AlphaPhy)

#Difference in Median AlphaPhy #no interactions

surviveyr2late_YN_Med_AlphaPhy_diff <- glmmTMB(surviveyr2late_YN ~ rescale(Med_AlphaPhy_diff)*Disturbance + rescale(suitabilitypc1V2abs) + 
                                                 rescale(suitabilitypc2V2abs) + Heating + rescale(awmpd) + rescale(OptGermRt) +
                                                 (1|Species) + (1|Plot), 
                                               family=binomial(link="logit"), 
                                               data=haueserdatayr2late)

surviveyr2late_YN_Med_AlphaPhy_diff <- glmmTMB(surviveyr2late_YN ~ Med_AlphaPhy_diff*Disturbance + suitabilitypc1V2abs + 
                                                 suitabilitypc2V2abs + Heating + scale(awmpd) + OptGermRt +
                                                 (1|Species) + (1|Plot), 
                                               family=binomial(link="logit"), 
                                               data=haueserdatayr2late)
summary(surviveyr2late_YN_Med_AlphaPhy_diff)
Anova(surviveyr2late_YN_Med_AlphaPhy_diff, type="III")
confint(surviveyr2late_YN_Med_AlphaPhy_diff)
r.squaredGLMM(surviveyr2late_YN_Med_AlphaPhy_diff)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2late_YN_Med_AlphaPhy_diff, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr2late_YN_Med_AlphaPhy_diff)
check_singularity(surviveyr2late_YN_Med_AlphaPhy_diff)

##########################################################

#FLOWERING

#there was NO flowering in undisturbed plots
haueserdisturb <- filter(haueserdata, Disturbance=="1")

#Raw Mean AlphaPhy

flowers_YN_Mean_AlphaPhy <- glmmTMB(flowers_YN ~ Mean_AlphaPhy + suitabilitypc2V2abs +
                                      suitabilitypc1V2abs + Heating + awmpd + OptGermRt +
                                      (1|Species) + (1|Plot), 
                                    family=binomial(link="logit"), 
                                    data=haueserdisturb)
summary(flowers_YN_Mean_AlphaPhy)
Anova(flowers_YN_Mean_AlphaPhy, type="III")
confint(flowers_YN_Mean_AlphaPhy)
r.squaredGLMM(flowers_YN_Mean_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = flowers_YN_Mean_AlphaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(flowers_YN_Mean_AlphaPhy)
check_singularity(flowers_YN_Mean_AlphaPhy)

#Difference in Mean AlphaPhy 

flowers_YN_Mean_AlphaPhy_diff <- glmmTMB(flowers_YN ~ rescale(Mean_AlphaPhy_diff) + rescale(suitabilitypc2V2abs) +
                                           rescale(suitabilitypc1V2abs) + Heating + rescale(awmpd) + rescale(OptGermRt) +
                                           (1|Species) + (1|Plot), 
                                         family=binomial(link="logit"), 
                                         data=haueserdisturb)

flowers_YN_Mean_AlphaPhy_diff <- glmmTMB(flowers_YN ~ Mean_AlphaPhy_diff + suitabilitypc2V2abs +
                                           suitabilitypc1V2abs + Heating + scale(awmpd) + OptGermRt +
                                           (1|Species) + (1|Plot), 
                                         family=binomial(link="logit"), 
                                         data=haueserdisturb)
summary(flowers_YN_Mean_AlphaPhy_diff)
Anova(flowers_YN_Mean_AlphaPhy_diff, type="III")
confint(flowers_YN_Mean_AlphaPhy_diff)
r.squaredGLMM(flowers_YN_Mean_AlphaPhy_diff)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = flowers_YN_Mean_AlphaPhy_diff, plot=F)
plot(simulationOutput)
check_overdispersion(flowers_YN_Mean_AlphaPhy_diff)
check_singularity(flowers_YN_Mean_AlphaPhy_diff)

#Raw Max AlphaPhy #no interactions

flowers_YN_Max_AlphaPhy <- glmmTMB(flowers_YN ~ Max_AlphaPhy + suitabilitypc2V2abs + 
                                     suitabilitypc1V2abs + Heating + awmpd + OptGermRt +
                                     (1|Species) + (1|Plot), 
                                   family=binomial(link="logit"), 
                                   data=haueserdisturb)
summary(flowers_YN_Max_AlphaPhy)
Anova(flowers_YN_Max_AlphaPhy, type="III")
confint(flowers_YN_Max_AlphaPhy)
r.squaredGLMM(flowers_YN_Max_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = flowers_YN_Max_AlphaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(flowers_YN_Max_AlphaPhy)
check_singularity(flowers_YN_Max_AlphaPhy)

#Difference in Max AlphaPhy 

flowers_YN_Max_AlphaPhy_diff <- glmmTMB(flowers_YN ~ rescale(Max_AlphaPhy_diff) + rescale(suitabilitypc2V2abs) +  
                                          rescale(suitabilitypc1V2abs) + Heating + rescale(awmpd) + rescale(OptGermRt) +
                                          (1|Species) + (1|Plot), 
                                        family=binomial(link="logit"), 
                                        data=haueserdisturb)

flowers_YN_Max_AlphaPhy_diff <- glmmTMB(flowers_YN ~ Max_AlphaPhy_diff + suitabilitypc2V2abs +  
                                          suitabilitypc1V2abs + Heating + scale(awmpd) + OptGermRt +
                                          (1|Species) + (1|Plot), 
                                        family=binomial(link="logit"), 
                                        data=haueserdisturb)
summary(flowers_YN_Max_AlphaPhy_diff)
Anova(flowers_YN_Max_AlphaPhy_diff, type="III")
confint(flowers_YN_Max_AlphaPhy_diff)
r.squaredGLMM(flowers_YN_Max_AlphaPhy_diff)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = flowers_YN_Max_AlphaPhy_diff, plot=F)
plot(simulationOutput)
check_overdispersion(flowers_YN_Max_AlphaPhy_diff)
check_singularity(flowers_YN_Max_AlphaPhy_diff)

#Raw GammaPhy #no interactions 

flowers_YN_GammaPhy <- glmmTMB(flowers_YN ~ scale(GammaPhy) + suitabilitypc2V2abs +
                                 suitabilitypc1V2abs + Heating + awmpd + OptGermRt +
                                 (1|Species) + (1|Plot), 
                               family=binomial(link="logit"), 
                               data=haueserdisturb)
summary(flowers_YN_GammaPhy)
Anova(flowers_YN_GammaPhy, type="III")
confint(flowers_YN_GammaPhy)
r.squaredGLMM(flowers_YN_GammaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = flowers_YN_GammaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(flowers_YN_GammaPhy)
check_singularity(flowers_YN_GammaPhy)

#Difference in GammaPhy

flowers_YN_GammaPhy_diff <- glmmTMB(flowers_YN ~ scale(GammaPhy_diff) + rescale(suitabilitypc2V2abs) + 
                                      rescale(suitabilitypc1V2abs) + Heating + rescale(awmpd) + rescale(OptGermRt) +
                                      (1|Species) + (1|Plot), 
                                    family=binomial(link="logit"), 
                                    data=haueserdisturb)

flowers_YN_GammaPhy_diff <- glmmTMB(flowers_YN ~ scale(GammaPhy_diff) + suitabilitypc2V2abs + 
                                      suitabilitypc1V2abs + Heating + scale(awmpd) + OptGermRt +
                                      (1|Species) + (1|Plot), 
                                    family=binomial(link="logit"), 
                                    data=haueserdisturb)
summary(flowers_YN_GammaPhy_diff)
Anova(flowers_YN_GammaPhy_diff, type="III")
confint(flowers_YN_GammaPhy_diff)
r.squaredGLMM(flowers_YN_GammaPhy_diff)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = flowers_YN_GammaPhy_diff, plot=F)
plot(simulationOutput)
check_overdispersion(flowers_YN_GammaPhy_diff)
check_singularity(flowers_YN_GammaPhy_diff)

#Raw Median AlphaPhy #no interactions

flowers_YN_Med_AlphaPhy <- glmmTMB(flowers_YN ~ Med_AlphaPhy + suitabilitypc1V2abs + 
                                     suitabilitypc2V2abs + Heating + awmpd + OptGermRt +
                                     (1|Species) + (1|Plot), 
                                   family=binomial(link="logit"), 
                                   data=haueserdisturb)
summary(flowers_YN_Med_AlphaPhy)
Anova(flowers_YN_Med_AlphaPhy, type="III")
confint(flowers_YN_Med_AlphaPhy)
r.squaredGLMM(flowers_YN_Med_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = flowers_YN_Med_AlphaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(flowers_YN_Med_AlphaPhy)
check_singularity(flowers_YN_Med_AlphaPhy)

#Difference in Med AlphaPhy #no interactions

flowers_YN_Med_AlphaPhy_diff <- glmmTMB(flowers_YN ~ rescale(Med_AlphaPhy_diff) + rescale(suitabilitypc1V2abs) +  
                                          rescale(suitabilitypc2V2abs) + Heating + rescale(awmpd) + rescale(OptGermRt) +
                                          (1|Species) + (1|Plot), 
                                        family=binomial(link="logit"), 
                                        data=haueserdisturb)

flowers_YN_Med_AlphaPhy_diff <- glmmTMB(flowers_YN ~ Med_AlphaPhy_diff + suitabilitypc1V2abs +  
                                          suitabilitypc2V2abs + Heating + scale(awmpd) + OptGermRt +
                                          (1|Species) + (1|Plot), 
                                        family=binomial(link="logit"), 
                                        data=haueserdisturb)
summary(flowers_YN_Med_AlphaPhy_diff)
Anova(flowers_YN_Med_AlphaPhy_diff, type="III")
confint(flowers_YN_Med_AlphaPhy_diff)
r.squaredGLMM(flowers_YN_Med_AlphaPhy_diff)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = flowers_YN_Med_AlphaPhy_diff, plot=F)
plot(simulationOutput)
check_overdispersion(flowers_YN_Med_AlphaPhy_diff)
check_singularity(flowers_YN_Med_AlphaPhy_diff)

####################################################################

#COUNT DATA

####################################################################

#GERMINATION COUNT

hauesergermcount <- filter(haueserdata, germ_num>0)

#Raw Mean AlphaPhy 

germ_count_Mean_AlphaPhy <- glmmTMB(germ_num ~ Mean_AlphaPhy + Disturbance + suitabilitypc1V2abs + 
                                      suitabilitypc2V2abs + Heating + awmpd + OptGermRt +
                                      (1|Species) + (1|Plot),
                                    family=poisson(link = "log"), 
                                    data=hauesergermcount)
summary(germ_count_Mean_AlphaPhy)
Anova(germ_count_Mean_AlphaPhy, type="III")
confint(germ_count_Mean_AlphaPhy)
r.squaredGLMM(germ_count_Mean_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = germ_count_Mean_AlphaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(germ_count_Mean_AlphaPhy)
check_singularity(germ_count_Mean_AlphaPhy)

#Mean AlphaPhy difference 

germ_count_Mean_AlphaPhy_diff <- glmmTMB(germ_num ~ rescale(Mean_AlphaPhy_diff) + Disturbance + rescale(suitabilitypc1V2abs) + 
                                           rescale(suitabilitypc2V2abs) + Heating + rescale(awmpd) + rescale(OptGermRt) +
                                           (1|Species) + (1|Plot),
                                         family=poisson(link = "log"), 
                                         data=hauesergermcount)

germ_count_Mean_AlphaPhy_diff <- glmmTMB(germ_num ~ Mean_AlphaPhy_diff + Disturbance + suitabilitypc1V2abs + 
                                           suitabilitypc2V2abs + Heating + scale(awmpd) + OptGermRt +
                                           (1|Species) + (1|Plot),
                                         family=poisson(link = "log"), 
                                         data=hauesergermcount)
summary(germ_count_Mean_AlphaPhy_diff)
Anova(germ_count_Mean_AlphaPhy_diff, type="III")
confint(germ_count_Mean_AlphaPhy_diff)
r.squaredGLMM(germ_count_Mean_AlphaPhy_diff)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = germ_count_Mean_AlphaPhy_diff, plot=F)
plot(simulationOutput)
check_overdispersion(germ_count_Mean_AlphaPhy_diff)
check_singularity(germ_count_Mean_AlphaPhy_diff)
check_model(germ_count_Mean_AlphaPhy_diff)

#Testing if native status has an effect

germ_count_Mean_AlphaPhy_diffN <- glmmTMB(germ_num ~ Mean_AlphaPhy_diff + Disturbance + suitabilitypc1V2abs + 
                                            suitabilitypc2V2abs + Heating + scale(awmpd) + OptGermRt + Status +
                                            (1|Species) + (1|Plot),
                                          family=poisson(link = "log"), 
                                          data=hauesergermcount)
summary(germ_count_Mean_AlphaPhy_diffN)

AICtab(germ_count_Mean_AlphaPhy_diffN, germ_count_Mean_AlphaPhy_diff)

#Test effect of biogeographic origin

hauesergermcountBG <- filter(haueserdataBG, germ_num>0)

germ_count_Mean_AlphaPhy_diffBG <- glmmTMB(germ_num ~ Mean_AlphaPhy_diff + Disturbance + suitabilitypc1V2abs + 
                                             suitabilitypc2V2abs + Heating + scale(awmpd) + OptGermRt +
                                             (1|Species) + (1|Plot) + (1|BG_ID),
                                           family=poisson(link = "log"), 
                                           data=hauesergermcountBG)
summary(germ_count_Mean_AlphaPhy_diffBG)

AICtab(germ_count_Mean_AlphaPhy_diffBG, germ_count_Mean_AlphaPhy_diff)

#Raw Max AlphaPhy 

germ_count_Max_AlphaPhy <- glmmTMB(germ_num ~ Max_AlphaPhy*suitabilitypc2V2abs + Disturbance +  
                                     suitabilitypc1V2abs + Heating + awmpd + OptGermRt +
                                     (1|Species) + (1|Plot),
                                   family=poisson(link = "log"), 
                                   data=hauesergermcount)
summary(germ_count_Max_AlphaPhy)
Anova(germ_count_Max_AlphaPhy, type="III")
confint(germ_count_Max_AlphaPhy)
r.squaredGLMM(germ_count_Max_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = germ_count_Max_AlphaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(germ_count_Max_AlphaPhy)
check_singularity(germ_count_Max_AlphaPhy)

#Mean AlphaPhy difference 

germ_count_Max_AlphaPhy_diff <- glmmTMB(germ_num ~ rescale(Max_AlphaPhy_diff)*rescale(suitabilitypc2V2abs) + Disturbance + rescale(suitabilitypc1V2abs) + 
                                          Heating + rescale(awmpd) + rescale(OptGermRt) +
                                          (1|Species) + (1|Plot),
                                        family=poisson(link = "log"), 
                                        data=hauesergermcount)

germ_count_Max_AlphaPhy_diff <- glmmTMB(germ_num ~ Max_AlphaPhy_diff*suitabilitypc2V2abs + Disturbance + suitabilitypc1V2abs + 
                                          Heating + scale(awmpd) + OptGermRt +
                                          (1|Species) + (1|Plot),
                                        family=poisson(link = "log"), 
                                        data=hauesergermcount)
summary(germ_count_Max_AlphaPhy_diff)
Anova(germ_count_Max_AlphaPhy_diff, type="III")
confint(germ_count_Max_AlphaPhy_diff)
r.squaredGLMM(germ_count_Max_AlphaPhy_diff)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = germ_count_Max_AlphaPhy_diff, plot=F)
plot(simulationOutput)
check_overdispersion(germ_count_Max_AlphaPhy_diff)
check_singularity(germ_count_Max_AlphaPhy_diff)

#Raw GammaPhy 

germ_count_GammaPhy <- glmmTMB(germ_num ~ scale(GammaPhy) + Disturbance + suitabilitypc1V2abs + 
                                 suitabilitypc2V2abs + Heating + awmpd + OptGermRt +
                                 (1|Species) + (1|Plot),
                               family=poisson(link = "log"), 
                               data=hauesergermcount)
summary(germ_count_GammaPhy)
Anova(germ_count_GammaPhy, type="III")
confint(germ_count_GammaPhy)
r.squaredGLMM(germ_count_GammaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = germ_count_GammaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(germ_count_GammaPhy)
check_singularity(germ_count_GammaPhy)
check_model(germ_count_GammaPhy)

#GammaPhy difference 

germ_count_GammaPhy_diff <- glmmTMB(germ_num ~ scale(GammaPhy_diff) + rescale(suitabilitypc1V2abs) + Disturbance +  
                                      rescale(suitabilitypc2V2abs) + Heating + rescale(awmpd) + rescale(OptGermRt) +
                                      (1|Species) + (1|Plot),
                                    family=poisson(link = "log"), 
                                    data=hauesergermcount)

germ_count_GammaPhy_diff <- glmmTMB(germ_num ~ scale(GammaPhy_diff) + suitabilitypc1V2abs + Disturbance +  
                                      suitabilitypc2V2abs + Heating + scale(awmpd) + OptGermRt +
                                      (1|Species) + (1|Plot),
                                    family=poisson(link = "log"), 
                                    data=hauesergermcount)
summary(germ_count_GammaPhy_diff)
Anova(germ_count_GammaPhy_diff, type="III")
confint(germ_count_GammaPhy_diff)
r.squaredGLMM(germ_count_GammaPhy_diff)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = germ_count_GammaPhy_diff, plot=F)
plot(simulationOutput)
check_overdispersion(germ_count_GammaPhy_diff)
check_singularity(germ_count_GammaPhy_diff)

#Raw Median AlphaPhy 

germ_count_Med_AlphaPhy <- glmmTMB(germ_num ~ Med_AlphaPhy + suitabilitypc1V2abs +  Disturbance +
                                     suitabilitypc2V2abs + Heating + awmpd + OptGermRt +
                                     (1|Species) + (1|Plot),
                                   family=poisson(link = "log"), 
                                   data=hauesergermcount)
summary(germ_count_Med_AlphaPhy)
Anova(germ_count_Med_AlphaPhy, type="III")
confint(germ_count_Med_AlphaPhy)
r.squaredGLMM(germ_count_Med_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = germ_count_Med_AlphaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(germ_count_Med_AlphaPhy)
check_singularity(germ_count_Med_AlphaPhy)

#Med AlphaPhy difference 

germ_count_Med_AlphaPhy_diff <- glmmTMB(germ_num ~ rescale(Med_AlphaPhy_diff) + rescale(suitabilitypc1V2abs) + Disturbance +  
                                          rescale(suitabilitypc2V2abs) + Heating + rescale(awmpd) + rescale(OptGermRt) +
                                          (1|Species) + (1|Plot),
                                        family=poisson(link = "log"), 
                                        data=hauesergermcount)

germ_count_Med_AlphaPhy_diff <- glmmTMB(germ_num ~ Med_AlphaPhy_diff + suitabilitypc1V2abs + Disturbance +  
                                          suitabilitypc2V2abs + Heating + scale(awmpd) + OptGermRt +
                                          (1|Species) + (1|Plot),
                                        family=poisson(link = "log"), 
                                        data=hauesergermcount)
summary(germ_count_Med_AlphaPhy_diff)
Anova(germ_count_Med_AlphaPhy_diff, type="III")
confint(germ_count_Med_AlphaPhy_diff)
r.squaredGLMM(germ_count_Med_AlphaPhy_diff)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = germ_count_Med_AlphaPhy_diff, plot=F)
plot(simulationOutput)
check_overdispersion(germ_count_Med_AlphaPhy_diff)
check_singularity(germ_count_Med_AlphaPhy_diff)
check_model(germ_count_Med_AlphaPhy_diff)

####################################################################

#FIRST YEAR SURVIVAL COUNT

haueseryr1count <- filter(haueserdata, surviveyr1_num>0)

#Raw Mean AlphaPhy

yr1_count_Mean_AlphaPhy <- glmmTMB(surviveyr1_num ~ Mean_AlphaPhy + Disturbance + suitabilitypc1V2abs + 
                                     suitabilitypc2V2abs + Heating + awmpd + OptGermRt +
                                     (1|Species) + (1|Plot),
                                   family=poisson(link = "log"), 
                                   data=haueseryr1count)
summary(yr1_count_Mean_AlphaPhy)
Anova(yr1_count_Mean_AlphaPhy, type="III")
confint(yr1_count_Mean_AlphaPhy)
r.squaredGLMM(yr1_count_Mean_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = yr1_count_Mean_AlphaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(yr1_count_Mean_AlphaPhy)
check_singularity(yr1_count_Mean_AlphaPhy)

#Mean AlphaPhy difference

yr1_count_Mean_AlphaPhy_diff <- glmmTMB(surviveyr1_num ~ rescale(Mean_AlphaPhy_diff) + Disturbance + rescale(suitabilitypc1V2abs) + 
                                          rescale(suitabilitypc2V2abs) + Heating + rescale(awmpd) + rescale(OptGermRt) +
                                          (1|Species) + (1|Plot),
                                        family=poisson(link = "log"), 
                                        data=haueseryr1count)

yr1_count_Mean_AlphaPhy_diff <- glmmTMB(surviveyr1_num ~ Mean_AlphaPhy_diff + Disturbance + suitabilitypc1V2abs + 
                                          suitabilitypc2V2abs + Heating + scale(awmpd) + OptGermRt +
                                          (1|Species) + (1|Plot),
                                        family=poisson(link = "log"), 
                                        data=haueseryr1count)
summary(yr1_count_Mean_AlphaPhy_diff)
Anova(yr1_count_Mean_AlphaPhy_diff, type="III")
confint(yr1_count_Mean_AlphaPhy_diff)
r.squaredGLMM(yr1_count_Mean_AlphaPhy_diff)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = yr1_count_Mean_AlphaPhy_diff, plot=F)
plot(simulationOutput)
check_overdispersion(yr1_count_Mean_AlphaPhy_diff)
check_singularity(yr1_count_Mean_AlphaPhy_diff)

#Testing if native status has an effect

yr1_count_Mean_AlphaPhy_diffN <- glmmTMB(surviveyr1_num ~ Mean_AlphaPhy_diff + Disturbance + suitabilitypc1V2abs + 
                                           suitabilitypc2V2abs + Heating + scale(awmpd) + OptGermRt + Status +
                                           (1|Species) + (1|Plot),
                                         family=poisson(link = "log"), 
                                         data=haueseryr1count)
summary(yr1_count_Mean_AlphaPhy_diffN)

AICtab(yr1_count_Mean_AlphaPhy_diffN, yr1_count_Mean_AlphaPhy_diff)

#Test effect of biogeographic origin

haueseryr1countBG <- filter(haueserdataBG, surviveyr1_num>0)

yr1_count_Mean_AlphaPhy_diffBG <- glmmTMB(surviveyr1_num ~ Mean_AlphaPhy_diff + Disturbance + suitabilitypc1V2abs + 
                                            suitabilitypc2V2abs + Heating + scale(awmpd) + OptGermRt +
                                            (1|Species) + (1|Plot) + (1|BG_ID),
                                          family=poisson(link = "log"), 
                                          data=haueseryr1countBG)
summary(yr1_count_Mean_AlphaPhy_diffBG)

AICtab(yr1_count_Mean_AlphaPhy_diffBG, yr1_count_Mean_AlphaPhy_diff)

#Raw Max AlphaPhy 

yr1_count_Max_AlphaPhy <- glmmTMB(surviveyr1_num ~ Max_AlphaPhy*suitabilitypc2V2abs + Disturbance + awmpd +   
                                    suitabilitypc1V2abs + Heating + OptGermRt +
                                    (1|Species) + (1|Plot),
                                  family=poisson(link = "log"), 
                                  data=haueseryr1count)
summary(yr1_count_Max_AlphaPhy)
Anova(yr1_count_Max_AlphaPhy, type="III")
confint(yr1_count_Max_AlphaPhy)
r.squaredGLMM(yr1_count_Max_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = yr1_count_Max_AlphaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(yr1_count_Max_AlphaPhy)
check_singularity(yr1_count_Max_AlphaPhy)

#Max AlphaPhy difference 

yr1_count_Max_AlphaPhy_diff <- glmmTMB(surviveyr1_num ~ rescale(Max_AlphaPhy_diff)*rescale(suitabilitypc2V2abs) + Disturbance +  
                                         rescale(suitabilitypc1V2abs) + Heating + rescale(awmpd) + rescale(OptGermRt) +
                                         (1|Species) + (1|Plot),
                                       family=poisson(link = "log"), 
                                       data=haueseryr1count)

yr1_count_Max_AlphaPhy_diff <- glmmTMB(surviveyr1_num ~ Max_AlphaPhy_diff*suitabilitypc2V2abs + Disturbance +  
                                         suitabilitypc1V2abs + Heating + scale(awmpd) + OptGermRt +
                                         (1|Species) + (1|Plot),
                                       family=poisson(link = "log"), 
                                       data=haueseryr1count)
summary(yr1_count_Max_AlphaPhy_diff)
Anova(yr1_count_Max_AlphaPhy_diff, type="III")
confint(yr1_count_Max_AlphaPhy_diff)
r.squaredGLMM(yr1_count_Max_AlphaPhy_diff)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = yr1_count_Max_AlphaPhy_diff, plot=F)
plot(simulationOutput)
check_overdispersion(yr1_count_Max_AlphaPhy_diff)
check_singularity(yr1_count_Max_AlphaPhy_diff)

#Raw GammaPhy 

yr1_count_GammaPhy <- glmmTMB(surviveyr1_num ~ scale(GammaPhy) + Disturbance + suitabilitypc1V2abs + 
                                suitabilitypc2V2abs + Heating + awmpd + OptGermRt +
                                (1|Species) + (1|Plot),
                              family=poisson(link = "log"), 
                              data=haueseryr1count)
summary(yr1_count_GammaPhy)
Anova(yr1_count_GammaPhy, type="III")
confint(yr1_count_GammaPhy)
r.squaredGLMM(yr1_count_GammaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = yr1_count_GammaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(yr1_count_GammaPhy)
check_singularity(yr1_count_GammaPhy)
check_model(yr1_count_GammaPhy)

#GammaPhy difference 

yr1_count_GammaPhy_diff <- glmmTMB(surviveyr1_num ~ scale(GammaPhy_diff) + Disturbance + rescale(suitabilitypc1V2abs) + 
                                     rescale(suitabilitypc2V2abs) + Heating + rescale(awmpd) + rescale(OptGermRt) +
                                     (1|Species) + (1|Plot),
                                   family=poisson(link = "log"), 
                                   data=haueseryr1count)

yr1_count_GammaPhy_diff <- glmmTMB(surviveyr1_num ~ scale(GammaPhy_diff) + Disturbance + suitabilitypc1V2abs + 
                                     suitabilitypc2V2abs + Heating + scale(awmpd) + OptGermRt +
                                     (1|Species) + (1|Plot),
                                   family=poisson(link = "log"), 
                                   data=haueseryr1count)
summary(yr1_count_GammaPhy_diff)
Anova(yr1_count_GammaPhy_diff, type="III")
confint(yr1_count_GammaPhy_diff)
r.squaredGLMM(yr1_count_GammaPhy_diff)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = yr1_count_GammaPhy_diff, plot=F)
plot(simulationOutput)
check_overdispersion(yr1_count_GammaPhy_diff)
check_singularity(yr1_count_GammaPhy_diff)

#Raw Median AlphaPhy

yr1_count_Med_AlphaPhy <- glmmTMB(surviveyr1_num ~ Med_AlphaPhy + suitabilitypc1V2abs + Disturbance +  
                                    suitabilitypc2V2abs + Heating + awmpd + OptGermRt +
                                    (1|Species) + (1|Plot),
                                  family=poisson(link = "log"), 
                                  data=haueseryr1count)
summary(yr1_count_Med_AlphaPhy)
Anova(yr1_count_Med_AlphaPhy, type="III")
confint(yr1_count_Med_AlphaPhy)
r.squaredGLMM(yr1_count_Med_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = yr1_count_Med_AlphaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(yr1_count_Med_AlphaPhy)
check_singularity(yr1_count_Med_AlphaPhy)

#Med AlphaPhy difference 

yr1_count_Med_AlphaPhy_diff <- glmmTMB(surviveyr1_num ~ rescale(Med_AlphaPhy_diff) + rescale(suitabilitypc1V2abs) +  Disturbance + 
                                         rescale(suitabilitypc2V2abs) + Heating + rescale(awmpd) + rescale(OptGermRt) +
                                         (1|Species) + (1|Plot),
                                       family=poisson(link = "log"), 
                                       data=haueseryr1count)

yr1_count_Med_AlphaPhy_diff <- glmmTMB(surviveyr1_num ~ Med_AlphaPhy_diff + suitabilitypc1V2abs +  Disturbance + 
                                         suitabilitypc2V2abs + Heating + scale(awmpd) + OptGermRt +
                                         (1|Species) + (1|Plot),
                                       family=poisson(link = "log"), 
                                       data=haueseryr1count)
summary(yr1_count_Med_AlphaPhy_diff)
Anova(yr1_count_Med_AlphaPhy_diff, type="III")
confint(yr1_count_Med_AlphaPhy_diff)
r.squaredGLMM(yr1_count_Med_AlphaPhy_diff)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = yr1_count_Med_AlphaPhy_diff, plot=F)
plot(simulationOutput)
check_overdispersion(yr1_count_Med_AlphaPhy_diff)
check_singularity(yr1_count_Med_AlphaPhy_diff)

####################################################################

#OVERWINTER SURVIVAL COUNT

haueseryr2earlycount <- filter(haueserdata, surviveyr2early_num>0)

#Raw Mean AlphaPhy 

yr2early_count_Mean_AlphaPhy <- glmmTMB(surviveyr2early_num ~ Mean_AlphaPhy + suitabilitypc2V2abs + Disturbance + suitabilitypc1V2abs + 
                                          + Heating + awmpd + OptGermRt +
                                          (1|Species) + (1|Plot),
                                        family=poisson(link = "log"), 
                                        data=haueseryr2earlycount)
summary(yr2early_count_Mean_AlphaPhy)
Anova(yr2early_count_Mean_AlphaPhy, type="III")
confint(yr2early_count_Mean_AlphaPhy)
r.squaredGLMM(yr2early_count_Mean_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = yr2early_count_Mean_AlphaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(yr2early_count_Mean_AlphaPhy)
check_singularity(yr2early_count_Mean_AlphaPhy)

#Mean AlphaPhy difference 

yr2early_count_Mean_AlphaPhy_diff <- glmmTMB(surviveyr2early_num ~ rescale(Mean_AlphaPhy_diff) + rescale(suitabilitypc2V2abs) + Disturbance + rescale(suitabilitypc1V2abs) + 
                                               Heating + rescale(awmpd) + rescale(OptGermRt) +
                                               (1|Species) + (1|Plot),
                                             family=poisson(link = "log"), 
                                             data=haueseryr2earlycount)

yr2early_count_Mean_AlphaPhy_diff <- glmmTMB(surviveyr2early_num ~ Mean_AlphaPhy_diff + suitabilitypc2V2abs + Disturbance + suitabilitypc1V2abs + 
                                               Heating + scale(awmpd) + OptGermRt +
                                               (1|Species) + (1|Plot),
                                             family=poisson(link = "log"), 
                                             data=haueseryr2earlycount)
summary(yr2early_count_Mean_AlphaPhy_diff)
Anova(yr2early_count_Mean_AlphaPhy_diff, type="III")
confint(yr2early_count_Mean_AlphaPhy_diff)
r.squaredGLMM(yr2early_count_Mean_AlphaPhy_diff)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = yr2early_count_Mean_AlphaPhy_diff, plot=F)
plot(simulationOutput)
check_overdispersion(yr2early_count_Mean_AlphaPhy_diff)
check_singularity(yr2early_count_Mean_AlphaPhy_diff)

#Testing whether native status has an effect

yr2early_count_Mean_AlphaPhy_diffN <- glmmTMB(surviveyr2early_num ~ Mean_AlphaPhy_diff + suitabilitypc2V2abs + Disturbance + suitabilitypc1V2abs + 
                                                Heating + scale(awmpd) + OptGermRt + Status +
                                                (1|Species) + (1|Plot),
                                              family=poisson(link = "log"), 
                                              data=haueseryr2earlycount)
summary(yr2early_count_Mean_AlphaPhy_diffN)

AICtab(yr2early_count_Mean_AlphaPhy_diffN, yr2early_count_Mean_AlphaPhy_diff)

#Test effect of biogeographic origin

haueseryr2earlycountBG <- filter(haueserdataBG, surviveyr2early_num>0)

yr2early_count_Mean_AlphaPhy_diffBG <- glmmTMB(surviveyr2early_num ~ Mean_AlphaPhy_diff + suitabilitypc2V2abs + Disturbance + suitabilitypc1V2abs + 
                                                 Heating + scale(awmpd) + OptGermRt +
                                                 (1|Species) + (1|Plot) + (1|BG_ID),
                                               family=poisson(link = "log"), 
                                               data=haueseryr2earlycountBG)
summary(yr2early_count_Mean_AlphaPhy_diffBG)

AICtab(yr2early_count_Mean_AlphaPhy_diffBG, yr2early_count_Mean_AlphaPhy_diff)

#Raw Max AlphaPhy 

yr2early_count_Max_AlphaPhy <- glmmTMB(surviveyr2early_num ~ Max_AlphaPhy + suitabilitypc2V2abs +  Disturbance +
                                         suitabilitypc1V2abs + Heating + awmpd + OptGermRt +
                                         (1|Species) + (1|Plot),
                                       family=poisson(link = "log"), 
                                       data=haueseryr2earlycount)
summary(yr2early_count_Max_AlphaPhy)
Anova(yr2early_count_Max_AlphaPhy, type="III")
confint(yr2early_count_Max_AlphaPhy)
r.squaredGLMM(yr2early_count_Max_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = yr2early_count_Max_AlphaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(yr2early_count_Max_AlphaPhy)
check_singularity(yr2early_count_Max_AlphaPhy)

#Max AlphaPhy difference 

yr2early_count_Max_AlphaPhy_diff <- glmmTMB(surviveyr2early_num ~ rescale(Max_AlphaPhy_diff) + Disturbance + rescale(suitabilitypc1V2abs) + 
                                              rescale(suitabilitypc2V2abs) + Heating + rescale(awmpd) + rescale(OptGermRt) +
                                              (1|Species) + (1|Plot),
                                            family=poisson(link = "log"), 
                                            data=haueseryr2earlycount)

yr2early_count_Max_AlphaPhy_diff <- glmmTMB(surviveyr2early_num ~ Max_AlphaPhy_diff + Disturbance + suitabilitypc1V2abs + 
                                              suitabilitypc2V2abs + Heating + scale(awmpd) + OptGermRt +
                                              (1|Species) + (1|Plot),
                                            family=poisson(link = "log"), 
                                            data=haueseryr2earlycount)
summary(yr2early_count_Max_AlphaPhy_diff)
Anova(yr2early_count_Max_AlphaPhy_diff, type="III")
confint(yr2early_count_Max_AlphaPhy_diff)
r.squaredGLMM(yr2early_count_Max_AlphaPhy_diff)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = yr2early_count_Max_AlphaPhy_diff, plot=F)
plot(simulationOutput)
check_overdispersion(yr2early_count_Max_AlphaPhy_diff)
check_singularity(yr2early_count_Max_AlphaPhy_diff)

#Raw GammaPhy 

yr2early_count_GammaPhy <- glmmTMB(surviveyr2early_num ~ scale(GammaPhy) + suitabilitypc2V2abs + Disturbance + suitabilitypc1V2abs + 
                                     Heating + awmpd + OptGermRt +
                                     (1|Species) + (1|Plot),
                                   family=poisson(link = "log"), 
                                   data=haueseryr2earlycount)
summary(yr2early_count_GammaPhy)
Anova(yr2early_count_GammaPhy, type="III")
confint(yr2early_count_GammaPhy)
r.squaredGLMM(yr2early_count_GammaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = yr2early_count_GammaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(yr2early_count_GammaPhy)
check_singularity(yr2early_count_GammaPhy)
check_model(yr2early_count_GammaPhy)

#GammaPhy difference 

yr2early_count_GammaPhy_diff <- glmmTMB(surviveyr2early_num ~ scale(GammaPhy_diff) + rescale(suitabilitypc2V2abs) + Disturbance + rescale(suitabilitypc1V2abs) + 
                                          Heating + rescale(awmpd) + rescale(OptGermRt) +
                                          (1|Species) + (1|Plot),
                                        family=poisson(link = "log"), 
                                        data=haueseryr2earlycount)

yr2early_count_GammaPhy_diff <- glmmTMB(surviveyr2early_num ~ scale(GammaPhy_diff) + suitabilitypc2V2abs + Disturbance + suitabilitypc1V2abs + 
                                          Heating + scale(awmpd) + OptGermRt +
                                          (1|Species) + (1|Plot),
                                        family=poisson(link = "log"), 
                                        data=haueseryr2earlycount)
summary(yr2early_count_GammaPhy_diff)
Anova(yr2early_count_GammaPhy_diff, type="III")
confint(yr2early_count_GammaPhy_diff)
r.squaredGLMM(yr2early_count_GammaPhy_diff)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = yr2early_count_GammaPhy_diff, plot=F)
plot(simulationOutput)
check_overdispersion(yr2early_count_GammaPhy_diff)
check_singularity(yr2early_count_GammaPhy_diff)

#Raw Median AlphaPhy

yr2early_count_Med_AlphaPhy <- glmmTMB(surviveyr2early_num ~ Med_AlphaPhy + suitabilitypc2V2abs + Disturbance + suitabilitypc1V2abs + 
                                         + Heating + awmpd + OptGermRt +
                                         (1|Species) + (1|Plot),
                                       family=poisson(link = "log"), 
                                       data=haueseryr2earlycount)
summary(yr2early_count_Med_AlphaPhy)
Anova(yr2early_count_Med_AlphaPhy, type="III")
confint(yr2early_count_Med_AlphaPhy)
r.squaredGLMM(yr2early_count_Med_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = yr2early_count_Med_AlphaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(yr2early_count_Med_AlphaPhy)
check_singularity(yr2early_count_Med_AlphaPhy)

#Median AlphaPhy difference 

yr2early_count_Med_AlphaPhy_diff <- glmmTMB(surviveyr2early_num ~ rescale(Med_AlphaPhy_diff) + rescale(suitabilitypc2V2abs) + Disturbance + rescale(suitabilitypc1V2abs) + 
                                              Heating + rescale(awmpd) + rescale(OptGermRt) +
                                              (1|Species) + (1|Plot),
                                            family=poisson(link = "log"), 
                                            data=haueseryr2earlycount)

yr2early_count_Med_AlphaPhy_diff <- glmmTMB(surviveyr2early_num ~ Med_AlphaPhy_diff + suitabilitypc2V2abs + Disturbance + suitabilitypc1V2abs + 
                                              Heating + scale(awmpd) + OptGermRt +
                                              (1|Species) + (1|Plot),
                                            family=poisson(link = "log"), 
                                            data=haueseryr2earlycount)
summary(yr2early_count_Med_AlphaPhy_diff)
Anova(yr2early_count_Med_AlphaPhy_diff, type="III")
confint(yr2early_count_Med_AlphaPhy_diff)
r.squaredGLMM(yr2early_count_Med_AlphaPhy_diff)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = yr2early_count_Med_AlphaPhy_diff, plot=F)
plot(simulationOutput)
check_overdispersion(yr2early_count_Med_AlphaPhy_diff)
check_singularity(yr2early_count_Med_AlphaPhy_diff)

####################################################################

#LATE SECOND YEAR SURVIVAL COUNT

haueseryr2latecount <- filter(haueserdata, surviveyr2late_num>0)

#Raw Mean AlphaPhy #no interactions

yr2late_count_Mean_AlphaPhy <- glmmTMB(surviveyr2late_num ~ Mean_AlphaPhy + suitabilitypc2V2abs + Disturbance + suitabilitypc1V2abs + 
                                         + Heating + awmpd + OptGermRt +
                                         (1|Species) + (1|Plot),
                                       family=poisson(link = "log"), 
                                       data=haueseryr2latecount)
summary(yr2late_count_Mean_AlphaPhy)
Anova(yr2late_count_Mean_AlphaPhy, type="III")
confint(yr2late_count_Mean_AlphaPhy)
r.squaredGLMM(yr2late_count_Mean_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = yr2late_count_Mean_AlphaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(yr2late_count_Mean_AlphaPhy)
check_singularity(yr2late_count_Mean_AlphaPhy)

#Mean AlphaPhy difference 

yr2late_count_Mean_AlphaPhy_diff <- glmmTMB(surviveyr2late_num ~ rescale(Mean_AlphaPhy_diff) + rescale(suitabilitypc2V2abs) + Disturbance + rescale(suitabilitypc1V2abs) + 
                                              Heating + rescale(awmpd) + rescale(OptGermRt) +
                                              (1|Species) + (1|Plot),
                                            family=poisson(link = "log"), 
                                            data=haueseryr2latecount)

yr2late_count_Mean_AlphaPhy_diff <- glmmTMB(surviveyr2late_num ~ Mean_AlphaPhy_diff + suitabilitypc2V2abs + Disturbance + suitabilitypc1V2abs + 
                                              Heating + scale(awmpd) + OptGermRt +
                                              (1|Species) + (1|Plot),
                                            family=poisson(link = "log"), 
                                            data=haueseryr2latecount)
summary(yr2late_count_Mean_AlphaPhy_diff)
Anova(yr2late_count_Mean_AlphaPhy_diff, type="III")
confint(yr2late_count_Mean_AlphaPhy_diff)
r.squaredGLMM(yr2late_count_Mean_AlphaPhy_diff)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = yr2late_count_Mean_AlphaPhy_diff, plot=F)
plot(simulationOutput)
check_overdispersion(yr2late_count_Mean_AlphaPhy_diff)
check_singularity(yr2late_count_Mean_AlphaPhy_diff)

#Testing whether native status has an effect

yr2late_count_Mean_AlphaPhy_diffN <- glmmTMB(surviveyr2late_num ~ Mean_AlphaPhy_diff + suitabilitypc2V2abs + Disturbance + suitabilitypc1V2abs + 
                                               Heating + scale(awmpd) + OptGermRt + Status +
                                               (1|Species) + (1|Plot),
                                             family=poisson(link = "log"), 
                                             data=haueseryr2latecount)
summary(yr2late_count_Mean_AlphaPhy_diffN)

AICtab(yr2late_count_Mean_AlphaPhy_diffN, yr2late_count_Mean_AlphaPhy_diff)

#Test effect of biogeographic origin

haueseryr2latecountBG <- filter(haueserdataBG, surviveyr2late_num>0)

yr2late_count_Mean_AlphaPhy_diffBG <- glmmTMB(surviveyr2late_num ~ Mean_AlphaPhy_diff + suitabilitypc2V2abs + Disturbance + suitabilitypc1V2abs + 
                                                Heating + scale(awmpd) + OptGermRt +
                                                (1|Species) + (1|Plot) + (1|BG_ID),
                                              family=poisson(link = "log"), 
                                              data=haueseryr2latecountBG)
summary(yr2late_count_Mean_AlphaPhy_diffBG)

AICtab(yr2late_count_Mean_AlphaPhy_diffBG, yr2late_count_Mean_AlphaPhy_diff)

#Raw Max AlphaPhy 

yr2late_count_Max_AlphaPhy <- glmmTMB(surviveyr2late_num ~ Max_AlphaPhy + Disturbance + suitabilitypc1V2abs + 
                                        suitabilitypc2V2abs + Heating + awmpd + OptGermRt +
                                        (1|Species) + (1|Plot),
                                      family=poisson(link = "log"), 
                                      data=haueseryr2latecount)
summary(yr2late_count_Max_AlphaPhy)
Anova(yr2late_count_Max_AlphaPhy, type="III")
confint(yr2late_count_Max_AlphaPhy)
r.squaredGLMM(yr2late_count_Max_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = yr2late_count_Max_AlphaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(yr2late_count_Max_AlphaPhy)
check_singularity(yr2late_count_Max_AlphaPhy)

#Max AlphaPhy difference 

yr2late_count_Max_AlphaPhy_diff <- glmmTMB(surviveyr2late_num ~ rescale(Max_AlphaPhy_diff) + Disturbance + rescale(suitabilitypc1V2abs) + 
                                             rescale(suitabilitypc2V2abs) + Heating + rescale(awmpd) + rescale(OptGermRt) +
                                             (1|Species) + (1|Plot),
                                           family=poisson(link = "log"), 
                                           data=haueseryr2latecount)

yr2late_count_Max_AlphaPhy_diff <- glmmTMB(surviveyr2late_num ~ Max_AlphaPhy_diff + Disturbance + suitabilitypc1V2abs + 
                                             suitabilitypc2V2abs + Heating + scale(awmpd) + OptGermRt +
                                             (1|Species) + (1|Plot),
                                           family=poisson(link = "log"), 
                                           data=haueseryr2latecount)
summary(yr2late_count_Max_AlphaPhy_diff)
Anova(yr2late_count_Max_AlphaPhy_diff, type="III")
confint(yr2late_count_Max_AlphaPhy_diff)
r.squaredGLMM(yr2late_count_Max_AlphaPhy_diff)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = yr2late_count_Max_AlphaPhy_diff, plot=F)
plot(simulationOutput)
check_overdispersion(yr2late_count_Max_AlphaPhy_diff)
check_singularity(yr2late_count_Max_AlphaPhy_diff)

##

#Raw GammaPhy

yr2late_count_GammaPhy <- glmmTMB(surviveyr2late_num ~ scale(GammaPhy)*suitabilitypc1V2abs + Disturbance + suitabilitypc2V2abs + 
                                    Heating + awmpd + OptGermRt +
                                    (1|Species) + (1|Plot),
                                  family=poisson(link = "log"), 
                                  data=haueseryr2latecount)
summary(yr2late_count_GammaPhy)
Anova(yr2late_count_GammaPhy, type="III")
confint(yr2late_count_GammaPhy)
r.squaredGLMM(yr2late_count_GammaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = yr2late_count_GammaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(yr2late_count_GammaPhy)
check_singularity(yr2late_count_GammaPhy)
check_model(yr2late_count_GammaPhy)

#GammaPhy difference

yr2late_count_GammaPhy_diff <- glmmTMB(surviveyr2late_num ~ scale(GammaPhy_diff) + rescale(suitabilitypc2V2abs) + Disturbance + rescale(suitabilitypc1V2abs) + 
                                         Heating + rescale(wmpd) + rescale(OptGermRt) +
                                         (1|Species) + (1|Plot),
                                       family=poisson(link = "log"), 
                                       data=haueseryr2latecount)

yr2late_count_GammaPhy_diff <- glmmTMB(surviveyr2late_num ~ scale(GammaPhy_diff) + suitabilitypc2V2abs + Disturbance + suitabilitypc1V2abs + 
                                         Heating + scale(awmpd) + OptGermRt +
                                         (1|Species) + (1|Plot),
                                       family=poisson(link = "log"), 
                                       data=haueseryr2latecount)
summary(yr2late_count_GammaPhy_diff)
Anova(yr2late_count_GammaPhy_diff, type="III")
confint(yr2late_count_GammaPhy_diff)
r.squaredGLMM(yr2late_count_GammaPhy_diff)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = yr2late_count_GammaPhy_diff, plot=F)
plot(simulationOutput)
check_overdispersion(yr2late_count_GammaPhy_diff)
check_singularity(yr2late_count_GammaPhy_diff)

#Raw Med AlphaPhy 

yr2late_count_Med_AlphaPhy <- glmmTMB(surviveyr2late_num ~ Med_AlphaPhy + suitabilitypc1V2abs + Disturbance + suitabilitypc2V2abs + 
                                        + Heating + awmpd + OptGermRt +
                                        (1|Species) + (1|Plot),
                                      family=poisson(link = "log"), 
                                      data=haueseryr2latecount)
summary(yr2late_count_Med_AlphaPhy)
Anova(yr2late_count_Med_AlphaPhy, type="III")
confint(yr2late_count_Med_AlphaPhy)
r.squaredGLMM(yr2late_count_Med_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = yr2late_count_Med_AlphaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(yr2late_count_Med_AlphaPhy)
check_singularity(yr2late_count_Med_AlphaPhy)

#Median AlphaPhy difference

yr2late_count_Med_AlphaPhy_diff <- glmmTMB(surviveyr2late_num ~ rescale(Med_AlphaPhy_diff)*rescale(suitabilitypc2V2abs) + Disturbance + rescale(suitabilitypc1V2abs) + 
                                             Heating + rescale(awmpd) + rescale(OptGermRt) +
                                             (1|Species) + (1|Plot),
                                           family=poisson(link = "log"), 
                                           data=haueseryr2latecount)

yr2late_count_Med_AlphaPhy_diff <- glmmTMB(surviveyr2late_num ~ Med_AlphaPhy_diff*suitabilitypc2V2abs + Disturbance + suitabilitypc1V2abs + 
                                             Heating + scale(awmpd) + OptGermRt +
                                             (1|Species) + (1|Plot),
                                           family=poisson(link = "log"), 
                                           data=haueseryr2latecount)
summary(yr2late_count_Med_AlphaPhy_diff)
Anova(yr2late_count_Med_AlphaPhy_diff, type="III")
confint(yr2late_count_Med_AlphaPhy_diff)
r.squaredGLMM(yr2late_count_Med_AlphaPhy_diff)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = yr2late_count_Med_AlphaPhy_diff, plot=F)
plot(simulationOutput)
check_overdispersion(yr2late_count_Med_AlphaPhy_diff)
check_singularity(yr2late_count_Med_AlphaPhy_diff)

####################################################################

#NUMBER OF FLOWERS

haueserflowers <- filter(haueserdisturb, max_flowers>0)
#ONLY FLOWERS IN DISTURBED PLOTS

#Raw Mean AlphaPhy 

flowercount_Mean_AlphaPhy <- glmmTMB(log(max_flowers) ~ Mean_AlphaPhy + suitabilitypc1V2abs + suitabilitypc2V2abs + 
                                       + Heating + awmpd + OptGermRt +
                                       (1|Species) + (1|Plot),
                                     family=poisson(link = "log"), 
                                     data=haueserflowers)
summary(flowercount_Mean_AlphaPhy)
Anova(flowercount_Mean_AlphaPhy, type="III")
confint(flowercount_Mean_AlphaPhy)
r.squaredGLMM(flowercount_Mean_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = flowercount_Mean_AlphaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(flowercount_Mean_AlphaPhy)
check_singularity(flowercount_Mean_AlphaPhy)
check_model(flowercount_Mean_AlphaPhy)

#Mean AlphaPhy difference

flowercount_Mean_AlphaPhy_diff <- glmmTMB(log(max_flowers) ~ rescale(Mean_AlphaPhy_diff) + rescale(suitabilitypc2V2abs) + rescale(suitabilitypc1V2abs) + 
                                            Heating + rescale(awmpd) + rescale(OptGermRt) +
                                            (1|Species) + (1|Plot),
                                          family=poisson(link = "log"), 
                                          data=haueserflowers)

flowercount_Mean_AlphaPhy_diff <- glmmTMB(log(max_flowers) ~ Mean_AlphaPhy_diff + suitabilitypc2V2abs + suitabilitypc1V2abs + 
                                            Heating + scale(awmpd) + OptGermRt +
                                            (1|Species) + (1|Plot),
                                          family=poisson(link = "log"), 
                                          data=haueserflowers)
summary(flowercount_Mean_AlphaPhy_diff)
Anova(flowercount_Mean_AlphaPhy_diff, type="III")
confint(flowercount_Mean_AlphaPhy_diff)
r.squaredGLMM(flowercount_Mean_AlphaPhy_diff)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = flowercount_Mean_AlphaPhy_diff, plot=F)
plot(simulationOutput)
check_overdispersion(flowercount_Mean_AlphaPhy_diff)
check_singularity(flowercount_Mean_AlphaPhy_diff)
check_model(flowercount_Mean_AlphaPhy_diff)

#Raw Max AlphaPhy 

flowercount_Max_AlphaPhy <- glmmTMB(log(max_flowers) ~ Max_AlphaPhy + suitabilitypc2V2abs + 
                                      suitabilitypc1V2abs + Heating + awmpd + OptGermRt +
                                      (1|Species) + (1|Plot),
                                    family=poisson(link = "log"), 
                                    data=haueserflowers)
summary(flowercount_Max_AlphaPhy)
Anova(flowercount_Max_AlphaPhy, type="III")
confint(flowercount_Max_AlphaPhy)
r.squaredGLMM(flowercount_Max_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = flowercount_Max_AlphaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(flowercount_Max_AlphaPhy)
check_singularity(flowercount_Max_AlphaPhy)
check_model(flowercount_Max_AlphaPhy)

#Max AlphaPhy difference

flowercount_Max_AlphaPhy_diff <- glmmTMB(log(max_flowers) ~ rescale(Max_AlphaPhy_diff) + rescale(suitabilitypc2V2abs) + 
                                           rescale(suitabilitypc1V2abs) + Heating + rescale(awmpd) + rescale(OptGermRt) +
                                           (1|Species) + (1|Plot),
                                         family=poisson(link = "log"), 
                                         data=haueserflowers)

flowercount_Max_AlphaPhy_diff <- glmmTMB(log(max_flowers) ~ Max_AlphaPhy_diff + suitabilitypc2V2abs + 
                                           suitabilitypc1V2abs + Heating + scale(awmpd) + OptGermRt +
                                           (1|Species) + (1|Plot),
                                         family=poisson(link = "log"), 
                                         data=haueserflowers)
summary(flowercount_Max_AlphaPhy_diff)
Anova(flowercount_Max_AlphaPhy_diff, type="III")
confint(flowercount_Max_AlphaPhy_diff)
r.squaredGLMM(flowercount_Max_AlphaPhy_diff)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = flowercount_Max_AlphaPhy_diff, plot=F)
plot(simulationOutput)
check_overdispersion(flowercount_Max_AlphaPhy_diff)
check_singularity(flowercount_Max_AlphaPhy_diff)

##

#Raw GammaPhy 

flowercount_GammaPhy <- glmmTMB(log(max_flowers) ~ scale(GammaPhy) + suitabilitypc2V2abs + suitabilitypc1V2abs + 
                                  Heating + awmpd + OptGermRt +
                                  (1|Species) + (1|Plot),
                                family=poisson(link = "log"), 
                                data=haueserflowers)
summary(flowercount_GammaPhy)
Anova(flowercount_GammaPhy, type="III")
confint(flowercount_GammaPhy)
r.squaredGLMM(flowercount_GammaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = flowercount_GammaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(flowercount_GammaPhy)
check_singularity(flowercount_GammaPhy)
check_model(flowercount_GammaPhy)

#GammaPhy difference 

flowercount_GammaPhy_diff <- glmmTMB(log(max_flowers) ~ scale(GammaPhy_diff) + rescale(suitabilitypc2V2abs) + rescale(suitabilitypc1V2abs) + 
                                       Heating + rescale(awmpd) + rescale(OptGermRt) +
                                       (1|Species) + (1|Plot),
                                     family=poisson(link = "log"), 
                                     data=haueserflowers)

flowercount_GammaPhy_diff <- glmmTMB(log(max_flowers) ~ scale(GammaPhy_diff) + suitabilitypc2V2abs + suitabilitypc1V2abs + 
                                       Heating + scale(awmpd) + OptGermRt +
                                       (1|Species) + (1|Plot),
                                     family=poisson(link = "log"), 
                                     data=haueserflowers)
summary(flowercount_GammaPhy_diff)
Anova(flowercount_GammaPhy_diff, type="III")
confint(flowercount_GammaPhy_diff)
r.squaredGLMM(flowercount_GammaPhy_diff)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = flowercount_GammaPhy_diff, plot=F)
plot(simulationOutput)
check_overdispersion(flowercount_GammaPhy_diff)
check_singularity(flowercount_GammaPhy_diff)

#Raw Med AlphaPhy

flowercount_Med_AlphaPhy <- glmmTMB(log(max_flowers) ~ Med_AlphaPhy + suitabilitypc1V2abs + suitabilitypc2V2abs + 
                                      + Heating + awmpd + OptGermRt +
                                      (1|Species) + (1|Plot),
                                    family=poisson(link = "log"), 
                                    data=haueserflowers)
summary(flowercount_Med_AlphaPhy)
Anova(flowercount_Med_AlphaPhy, type="III")
confint(flowercount_Med_AlphaPhy)
r.squaredGLMM(flowercount_Med_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = flowercount_Med_AlphaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(flowercount_Med_AlphaPhy)
check_singularity(flowercount_Med_AlphaPhy)
check_model(flowercount_Med_AlphaPhy)

#Mean AlphaPhy difference 

flowercount_Med_AlphaPhy_diff <- glmmTMB(log(max_flowers) ~ rescale(Med_AlphaPhy_diff) + rescale(suitabilitypc2V2abs) + rescale(suitabilitypc1V2abs) + 
                                           Heating + rescale(awmpd) + rescale(OptGermRt) +
                                           (1|Species) + (1|Plot),
                                         family=poisson(link = "log"), 
                                         data=haueserflowers)

flowercount_Med_AlphaPhy_diff <- glmmTMB(log(max_flowers) ~ Med_AlphaPhy_diff + suitabilitypc2V2abs + suitabilitypc1V2abs + 
                                           Heating + scale(awmpd) + OptGermRt +
                                           (1|Species) + (1|Plot),
                                         family=poisson(link = "log"), 
                                         data=haueserflowers)
summary(flowercount_Med_AlphaPhy_diff)
Anova(flowercount_Med_AlphaPhy_diff, type="III")
confint(flowercount_Med_AlphaPhy_diff)
r.squaredGLMM(flowercount_Med_AlphaPhy_diff)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = flowercount_Med_AlphaPhy_diff, plot=F)
plot(simulationOutput)
check_overdispersion(flowercount_Med_AlphaPhy_diff)
check_singularity(flowercount_Med_AlphaPhy_diff)
check_model(flowercount_Med_AlphaPhy_diff)

#########################################################################
#########################################################################

#Plotting all the results

#ALL EFFECT SIZES - MEAN ALPHA PD - Figure 4

hauesereffectsizes <- read.csv("haueserresults_effectsize_revised.csv", header=T, stringsAsFactors = T)
hauesereffectsizes$Stage = factor(hauesereffectsizes$Stage, levels = c("Second year","Overwinter","First year","Colonisation"))
hauesereffectsizes$Variable = factor(hauesereffectsizes$Variable, levels = c("Alpha PD diff. x Climate (temp.)", "Alpha PD diff. x Disturbed", "Heating", "Phylogenetic distance", "Climate dissimilarity (temp.)",
                                                                             "Climate dissimilarity (precip.)", "Disturbed", "Mean Alpha PD difference"))

hauesereffectplot <- ggplot(hauesereffectsizes, aes(x=effect, y=Variable, color=Stage, shape=significant)) +
  geom_point(size=2.3, position=position_dodge(0.5)) +
  scale_shape_manual(values=c(2, 1, 16)) +
  geom_errorbar(data=hauesereffectsizes, aes(y=Variable, xmin=effect_low, xmax=effect_high), 
                width=0, size=1, position=position_dodge(0.5)) +
  facet_grid(. ~ Metric, scales="free") +
  scale_x_continuous(trans='pseudo_log') +
  scale_color_viridis(discrete = TRUE, guide = guide_legend(reverse = TRUE)) +
  theme_bw() + geom_vline(xintercept=0, linetype="solid", color="grey", size=0.9) +
  xlab("Standardised effect size ± 95% C.I.") + guides(shape = "none") + theme(legend.position = "top") +
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=10), 
        legend.text=element_text(size=8), 
        legend.title=element_text(size=10))
hauesereffectplot

ggsave(hauesereffectplot, 
       filename = "Fig_4.svg",
       height = 140, width = 180, units = "mm")

#PLOT ALL TIME POINTS AND RESPONSES FOR MEAN ALPHA PHY - Figure 5

#Fig. 5a

it_vr <- visreg(germ_YN_Mean_AlphaPhy_diff, "Mean_AlphaPhy_diff", ylab="Germination Success",
                xlab="Mean_AlphaPhy", scale="response", partial=T)
it_vr_df <- as.data.frame(cbind(it_vr$fit$Mean_AlphaPhy_diff, it_vr$fit$visregFit,
                                it_vr$fit$visregLwr, it_vr$fit$visregUpr))
colnames(it_vr_df) <- c("it", "fit", "lwr", "upr")
it_vr_pr <- as.data.frame(cbind(it_vr$res$Mean_AlphaPhy_diff,
                                it_vr$res$visregRes))
colnames(it_vr_pr) <- c("it", "resids")

fig5a <- ggplot() +
  geom_jitter(aes(y=resids, x= it), alpha = 0.8, data = it_vr_pr, width=0.1, height=0.0001, colour="darkgray") +
  geom_line(aes(x=it, y=fit), data=it_vr_df, size = 1.2, linetype="solid") +
  theme_bw() +
  ylab("Colonisation") +
  xlab("Difference in Mean Alpha Phylogenetic Diversity") + 
  geom_vline(xintercept = 0, linetype="dashed", colour="grey", size=1) +
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=10), 
        legend.text=element_text(size=8), 
        legend.title=element_text(size=10)) + 
  xlim(-1.9, 0.9) + xlab(NULL)
fig5a

#Fig 5b

it_vr <- visreg(germ_count_Mean_AlphaPhy_diff, "Mean_AlphaPhy_diff", ylab="Germination Success",
                xlab="Mean_AlphaPhy",  partial=T) #scale="response",
it_vr_df <- as.data.frame(cbind(it_vr$fit$Mean_AlphaPhy_diff, it_vr$fit$visregFit,
                                it_vr$fit$visregLwr, it_vr$fit$visregUpr))
colnames(it_vr_df) <- c("it", "fit", "lwr", "upr")
it_vr_pr <- as.data.frame(cbind(it_vr$res$Mean_AlphaPhy_diff,
                                it_vr$res$visregRes))
colnames(it_vr_pr) <- c("it", "resids")

fig5b <- ggplot() +
  geom_jitter(aes(y=resids, x= it), alpha = 0.8, data = it_vr_pr, width=0.1, height=0.0001, colour="darkgray") +
  geom_line(aes(x=it, y=fit), data=it_vr_df, size = 1.2, linetype="solid") +
  theme_bw() +
  ylab("Colonisation") +
  xlab("Difference in Mean Alpha Phylogenetic Diversity") + 
  geom_vline(xintercept = 0, linetype="dashed", colour="grey", size=1) +
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=10), 
        legend.text=element_text(size=8), 
        legend.title=element_text(size=10)) + 
  xlim(-1.9, 0.9) + xlab(NULL) + ylab(NULL)
fig5b

#Fig. 5c

it_vr <- visreg(surviveyr1_YN_Mean_AlphaPhy_diff, "Mean_AlphaPhy_diff", ylab="Germination Success",
                xlab="Mean_AlphaPhy", scale="response", partial=T)
it_vr_df <- as.data.frame(cbind(it_vr$fit$Mean_AlphaPhy_diff, it_vr$fit$visregFit,
                                it_vr$fit$visregLwr, it_vr$fit$visregUpr))
colnames(it_vr_df) <- c("it", "fit", "lwr", "upr")
it_vr_pr <- as.data.frame(cbind(it_vr$res$Mean_AlphaPhy_diff,
                                it_vr$res$visregRes))
colnames(it_vr_pr) <- c("it", "resids")

fig5c <- ggplot() +
  geom_jitter(aes(y=resids, x= it), alpha = 0.8, data = it_vr_pr, width=0.1, height=0.0001, colour="darkgray") +
  geom_line(aes(x=it, y=fit), data=it_vr_df, size = 1.2, linetype="solid") +
  theme_bw() +
  ylab("First year") +
  xlab("Difference in Mean Alpha Phylogenetic Diversity") + 
  geom_vline(xintercept = 0, linetype="dashed", colour="grey", size=1) +
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=10), 
        legend.text=element_text(size=8), 
        legend.title=element_text(size=10)) +
  xlim(-1.9, 0.9) + xlab(NULL)
fig5c

#Fig 5d

it_vr <- visreg(yr1_count_Mean_AlphaPhy_diff, "Mean_AlphaPhy_diff", ylab="Germination Success",
                xlab="Mean_AlphaPhy",  partial=T) #scale="response",
it_vr_df <- as.data.frame(cbind(it_vr$fit$Mean_AlphaPhy_diff, it_vr$fit$visregFit,
                                it_vr$fit$visregLwr, it_vr$fit$visregUpr))
colnames(it_vr_df) <- c("it", "fit", "lwr", "upr")
it_vr_pr <- as.data.frame(cbind(it_vr$res$Mean_AlphaPhy_diff,
                                it_vr$res$visregRes))
colnames(it_vr_pr) <- c("it", "resids")

fig5d <- ggplot() +
  geom_jitter(aes(y=resids, x= it), alpha = 0.8, data = it_vr_pr, width=0.1, height=0.0001, colour="darkgray") +
  geom_line(aes(x=it, y=fit), data=it_vr_df, size = 1.2, linetype="solid") +
  theme_bw() +
  ylab("First year") +
  xlab("Difference in Mean Alpha Phylogenetic Diversity") + 
  geom_vline(xintercept = 0, linetype="dashed", colour="grey", size=1) +
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=10), 
        legend.text=element_text(size=8), 
        legend.title=element_text(size=10)) + 
  xlim(-1.9, 0.9) + xlab(NULL) + ylab(NULL)
fig5d

#Fig 5e

yr2earlymeantempint <- interact_plot(surviveyr2early_YN_Mean_AlphaPhy_diff, pred="Mean_AlphaPhy_diff", modx = "suitabilitypc1V2abs",
                                     plot.points=TRUE, partial.residuals = TRUE, legend.main = "Climate dissimilarity (temp.)") +
  xlab("Difference in Mean Alpha Phylogenetic Diversity") + ylab("Likelihood of any plant in plot surviving overwinter") +
  scale_linetype_manual(values=c("solid", "solid", "solid")) +
  
  theme_bw() + guides(linetype = "none") + theme(legend.position = "top")  +
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=10), 
        legend.text=element_text(size=8), 
        legend.title=element_text(size=10))
yr2earlymeantempint

#Fig 5f

it_vr <- visreg(yr2early_count_Mean_AlphaPhy_diff, "Mean_AlphaPhy_diff", ylab="Germination Success",
                xlab="Mean_AlphaPhy",  partial=T) #scale="response",
it_vr_df <- as.data.frame(cbind(it_vr$fit$Mean_AlphaPhy_diff, it_vr$fit$visregFit,
                                it_vr$fit$visregLwr, it_vr$fit$visregUpr))
colnames(it_vr_df) <- c("it", "fit", "lwr", "upr")
it_vr_pr <- as.data.frame(cbind(it_vr$res$Mean_AlphaPhy_diff,
                                it_vr$res$visregRes))
colnames(it_vr_pr) <- c("it", "resids")

fig5f <- ggplot() +
  geom_jitter(aes(y=resids, x= it), alpha = 0.8, data = it_vr_pr, width=0.1, height=0.0001, colour="darkgray") +
  #geom_line(aes(x=it, y=fit), data=it_vr_df, size = 1.2, linetype="solid") +
  theme_bw() +
  ylab("Overwinter") +
  xlab("Difference in Mean Alpha Phylogenetic Diversity") + 
  geom_vline(xintercept = 0, linetype="dashed", colour="grey", size=1) +
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=10), 
        legend.text=element_text(size=8), 
        legend.title=element_text(size=10))  + 
  xlim(-1.9, 0.9) + xlab(NULL) + ylab(NULL)
fig5f

#Fig 5g

yr2YNmeaninteraction <- interact_plot(surviveyr2late_YN_Mean_AlphaPhy_diff, pred="Mean_AlphaPhy_diff", modx = "Disturbance",
                                      plot.points=TRUE, partial.residuals = TRUE) +
  xlab("Difference in Mean Alpha Phylogenetic Diversity") + ylab("Likelihood of any plant in plot surviving to end of second year") +
  scale_color_manual(values=c("#999999", "#E69F00")) +
  scale_linetype_manual(values=c("solid", "solid")) +
  theme_bw() + guides(linetype = "none") + theme(legend.position = "top")  +
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=10), 
        legend.text=element_text(size=8), 
        legend.title=element_text(size=10))
yr2YNmeaninteraction

#Fig 5h

it_vr <- visreg(yr2late_count_Mean_AlphaPhy_diff, "Mean_AlphaPhy_diff", ylab="Germination Success",
                xlab="Mean_AlphaPhy",  partial=T) #scale="response",
it_vr_df <- as.data.frame(cbind(it_vr$fit$Mean_AlphaPhy_diff, it_vr$fit$visregFit,
                                it_vr$fit$visregLwr, it_vr$fit$visregUpr))
colnames(it_vr_df) <- c("it", "fit", "lwr", "upr")
it_vr_pr <- as.data.frame(cbind(it_vr$res$Mean_AlphaPhy_diff,
                                it_vr$res$visregRes))
colnames(it_vr_pr) <- c("it", "resids")

fig5h <- ggplot() +
  geom_jitter(aes(y=resids, x= it), alpha = 0.8, data = it_vr_pr, width=0.1, height=0.0001, colour="darkgray") +
  #geom_line(aes(x=it, y=fit), data=it_vr_df, size = 1.2, linetype="solid") +
  theme_bw() +
  ylab("Second year") +
  xlab("Difference in Mean Alpha Phylogenetic Diversity") + 
  geom_vline(xintercept = 0, linetype="dashed", colour="grey", size=1) +
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=10), 
        legend.text=element_text(size=8), 
        legend.title=element_text(size=10))  + 
  xlim(-1.9, 0.9) + ylab(NULL)
fig5h

Fig5 <- (fig5a | fig5b)/
  (fig5c | fig5d)/
  (yr2earlymeantempint | fig5f)/
  (yr2YNmeaninteraction | fig5h)
Fig5
ggsave(Fig5, 
       filename = "Figure_5.svg",
       height = 240, width = 180, units = "mm")

#MAX ALPHA PHY - INTERACTIONS - Figure S13

germYNmaxinteraction <- interact_plot(germ_YN_Max_AlphaPhy_diff, pred="Max_AlphaPhy_diff", modx = "Disturbance",
                                      plot.points=TRUE, partial.residuals = TRUE) +
  xlab("Difference in Maximum Alpha Phylogenetic Diversity") + ylab("Likelihood of any plant in plot colonising") +
  scale_color_manual(values=c("#999999", "#E69F00")) +
  scale_linetype_manual(values=c("solid", "solid")) +
  theme_bw() + guides(linetype = "none") + theme(legend.position = "top")  +
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=10), 
        legend.text=element_text(size=8), 
        legend.title=element_text(size=10))
germYNmaxinteraction

germcountmaxinteraction <- interact_plot(germ_count_Max_AlphaPhy_diff, pred="Max_AlphaPhy_diff", modx = "suitabilitypc2V2abs",
                                         plot.points=TRUE, partial.residuals = TRUE, legend.main = "Climate dissimilarity (precip.)") +
  xlab("Difference in Maximum Alpha Phylogenetic Diversity") + ylab("Number of plants colonising") +
  scale_linetype_manual(values=c("solid", "solid", "solid")) +
  scale_y_continuous(trans='log10') +
  theme_bw() + guides(linetype = "none") + theme(legend.position = "top")  +
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=10), 
        legend.text=element_text(size=8), 
        legend.title=element_text(size=10))
germcountmaxinteraction

yr1YNmaxdistinteraction <- interact_plot(surviveyr1_YN_Max_AlphaPhy_diff, pred="Max_AlphaPhy_diff", modx = "Disturbance",
                                         plot.points=TRUE, partial.residuals = TRUE) +
  xlab("Difference in Maximum Alpha Phylogenetic Diversity") + ylab("Likelihood of any plant in plot surviving to end of first year") +
  scale_color_manual(values=c("#999999", "#E69F00")) +
  scale_linetype_manual(values=c("solid", "solid")) +
  theme_bw() + guides(linetype = "none") + theme(legend.position = "top")  +
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=10), 
        legend.text=element_text(size=8), 
        legend.title=element_text(size=10))
yr1YNmaxdistinteraction

yr1countmaxinteraction <- interact_plot(yr1_count_Max_AlphaPhy_diff, pred="Max_AlphaPhy_diff", modx = "suitabilitypc2V2abs",
                                        plot.points=TRUE, partial.residuals = TRUE, legend.main = "Climate dissimilarity (precip.)") +
  xlab("Difference in Maximum Alpha Phylogenetic Diversity") + ylab("Number of plants surviving to end of first year") +
  scale_linetype_manual(values=c("solid", "solid", "solid")) +
  scale_y_continuous(trans='log10') +
  theme_bw() + guides(linetype = "none") + theme(legend.position = "top")  +
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=10), 
        legend.text=element_text(size=8), 
        legend.title=element_text(size=10))
yr1countmaxinteraction

#Group plot for Max AlphaPhy 

maxalphahaueserfig <- ggarrange(germYNmaxinteraction, germcountmaxinteraction,
                                yr1YNmaxdistinteraction, yr1countmaxinteraction, nrow=2, ncol=2)
maxalphahaueserfig
ggsave(maxalphahaueserfig, 
       filename = "Fig_S13.svg",
       height = 210, width = 210, units = "mm")

#MEDIAN ALPHAPHY - INTERACTIONS - Figure S14

medYNdistinteraction <- interact_plot(germ_YN_Med_AlphaPhy_diff, pred="Med_AlphaPhy_diff", modx = "Disturbance",
                                      plot.points=TRUE, partial.residuals = TRUE) +
  xlab("Difference in Median Phylogenetic Diversity") + ylab("Likelihood of any plant in plot colonising") +
  scale_color_manual(values=c("#999999", "#E69F00")) +
  scale_linetype_manual(values=c("solid", "solid")) +
  theme_bw() + guides(linetype = "none") + theme(legend.position = "top")  +
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=10), 
        legend.text=element_text(size=8), 
        legend.title=element_text(size=10))
medYNdistinteraction

medyr1YNdistinteraction <- interact_plot(surviveyr1_YN_Med_AlphaPhy_diff, pred="Med_AlphaPhy_diff", modx = "Disturbance",
                                         plot.points=TRUE, partial.residuals = TRUE) +
  xlab("Difference in Median Phylogenetic Diversity") + ylab("Likelihood of any plant in plot surviving to end of first year") +
  scale_color_manual(values=c("#999999", "#E69F00")) +
  scale_linetype_manual(values=c("solid", "solid")) +
  theme_bw() + guides(linetype = "none") + theme(legend.position = "top")  +
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=10), 
        legend.text=element_text(size=8), 
        legend.title=element_text(size=10))
medyr1YNdistinteraction

yr2earlyYNmedtempinteraction <- interact_plot(surviveyr2early_YN_Med_AlphaPhy_diff, pred="Med_AlphaPhy_diff", modx = "suitabilitypc1V2abs",
                                              plot.points=TRUE, partial.residuals = TRUE, legend.main = "Climate dissimilarity (temp.)") +
  xlab("Difference in Median Phylogenetic Diversity") + ylab("Likelihood of any plant in plot surviving overwinter") +
  scale_y_continuous(trans='log10') +
  scale_linetype_manual(values=c("solid", "solid", "solid")) +
  theme_bw() + guides(linetype = "none") + theme(legend.position = "top")  +
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=10), 
        legend.text=element_text(size=8), 
        legend.title=element_text(size=10))
yr2earlyYNmedtempinteraction

medyr2lateYNdistinteraction <- interact_plot(surviveyr2late_YN_Med_AlphaPhy_diff, pred="Med_AlphaPhy_diff", modx = "Disturbance",
                                             plot.points=TRUE, partial.residuals = TRUE) +
  xlab("Difference in Median Phylogenetic Diversity") + ylab("Likelihood of any plant in plot surviving to end of second year") +
  scale_color_manual(values=c("#999999", "#E69F00")) +
  scale_linetype_manual(values=c("solid", "solid")) +
  theme_bw() + guides(linetype = "none") + theme(legend.position = "top")  +
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=10), 
        legend.text=element_text(size=8), 
        legend.title=element_text(size=10))
medyr2lateYNdistinteraction

yr2latecountmedtempinteraction <- interact_plot(yr2late_count_Med_AlphaPhy_diff, pred="Med_AlphaPhy_diff", modx = "suitabilitypc2V2abs",
                                                plot.points=TRUE, partial.residuals = TRUE, legend.main = "Climate dissimilarity (precip.)") +
  xlab("Difference in Median Phylogenetic Diversity") + ylab("Number of plants surviving to end of second year") +
  scale_y_continuous(trans='log10') +
  scale_linetype_manual(values=c("solid", "solid", "solid")) +
  theme_bw() + guides(linetype = "none") + theme(legend.position = "top")  +
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=10), 
        legend.text=element_text(size=8), 
        legend.title=element_text(size=10))
yr2latecountmedtempinteraction

medianhaueserfig <- ggarrange(medYNdistinteraction, medyr1YNdistinteraction, yr2earlyYNmedtempinteraction,
                              medyr2lateYNdistinteraction, yr2latecountmedtempinteraction, nrow=3, ncol=2)
medianhaueserfig
ggsave(medianhaueserfig, 
       filename = "Fig_S14.svg",
       height = 315, width = 210, units = "mm")

#GAMMA PHY - INTERACTIONS - Figure S15

germYNgammadistinteraction <- interact_plot(germ_YN_GammaPhy_diff, pred="GammaPhy_diff", modx = "Disturbance",
                                            plot.points=TRUE, partial.residuals = TRUE) +
  xlab("Difference in Gamma Phylogenetic Diversity") + ylab("Likelihood of any plant in plot colonising") +
  scale_color_manual(values=c("#999999", "#E69F00")) +
  scale_linetype_manual(values=c("solid", "solid")) +
  theme_bw() + guides(linetype = "none") + theme(legend.position = "top")  +
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=10), 
        legend.text=element_text(size=8), 
        legend.title=element_text(size=10))
germYNgammadistinteraction

overwintercountgammainteraction <- interact_plot(surviveyr2early_YN_GammaPhy_diff, pred="GammaPhy_diff", modx = "suitabilitypc2V2abs",
                                                 plot.points=TRUE, partial.residuals = TRUE, legend.main = "Climate dissimilarity (precip.)") +
  xlab("Difference in Gamma Phylogenetic Diversity") + ylab("Likelihood of any plant in plot surviving overwinter") +
  scale_linetype_manual(values=c("solid", "solid", "solid")) +
  scale_y_continuous(trans='log10') +
  theme_bw() + guides(linetype = "none") + theme(legend.position = "top")  +
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=10), 
        legend.text=element_text(size=8), 
        legend.title=element_text(size=10))
overwintercountgammainteraction

#Group plot for GammaPhy 

gammahaueserfig <- ggarrange(germYNgammadistinteraction, overwintercountgammainteraction, nrow=1, ncol=2)
gammahaueserfig
ggsave(gammahaueserfig, 
       filename = "Fig_S15.svg",
       height = 105, width = 210, units = "mm")

##################################################################################################

#QUESTION THREE: Is high home-range-PD correlated with traits that could explain increased 
#survival of sown species? 

##################################################################################################

#Read in file of trait data and link with species information

fulltraits <- read.csv("traits_complete.csv", header=T) 

fulltraits$hardiness[fulltraits$hardiness==""] <- NA
fulltraits$Life_span[fulltraits$Life_span==""] <- NA

speciesinfo <- read.csv("allspeciesall3exp_clean.csv", header=T)
fulldata <- left_join(fulltraits, speciesinfo, by="plant_name_id")

BG <- read.csv("biogeographic_ID.csv", header=T) %>% dplyr::select(plant_name_id, BG_ID)
BG$BG_ID <- as.factor(BG$BG_ID)

fulldataBG <- left_join(fulldata, BG, by="plant_name_id")

fulldataBGfiltered <- fulldataBG %>% #Limit to those biogeographic syndromes that are native to Europe or USA only
  filter(BG_ID %in% c("7", "11", "12", "13", "17", "19", "21", "28", "30", "36", "38", "42", "43", "48", "54",
                      "65", "67", "68", "69", "70", "73", "77", "84", "85", "92", "93", "104", "105", "106", 
                      "107", "108", "109", "110", "111", "112", "113", "114", "115", "116", "117", "118", "119",
                      "142", "143", "144", "145", "146", "147", "148", "149", "150", "151", "152", "153", "154",
                      "155", "156", "157", "158", "159", "160", "161", "5", "9", "20", "32", "37", "41", "44", 
                      "46", "49", "51", "53", "57", "58", "59", "66", "74", "75", "78", "79", "83", "87", "90",
                      "98", "100", "101", "103", "121", "122", "124", "127", "128", "130", "132", "134", "139", "141"))

#For all traits, do three tests per PD metric: use all available data (corresponds to 'p-value (all data)' in Table S12);
#limit data to species that are native to North America or Europe only (the fulldataBGfiltered dataset, corresponds to
#'p-value (EU/NA)' in Table S12; use all data but include biogeographic syndrome as a random variable (corresponds to 
#'p-value (BG RV)' in Table S12). 

#Test relationship between resistance to herbivory and phylogenetic diversity

summary(lm(herbivory_resistance ~ Mean_AlphaPhy, data=fulldataBG))
summary(lm(herbivory_resistance ~ Mean_AlphaPhy, data=fulldataBGfiltered))
summary(glmmTMB(herbivory_resistance ~ Mean_AlphaPhy + (1|BG_ID), family=gaussian(), data=fulldataBG))

summary(lm(herbivory_resistance ~ Max_AlphaPhy, data=fulldataBG))
summary(lm(herbivory_resistance ~ Max_AlphaPhy, data=fulldataBGfiltered))
summary(glmmTMB(herbivory_resistance ~ Max_AlphaPhy + (1|BG_ID), family=gaussian(), data=fulldataBG))

summary(lm(herbivory_resistance ~ Med_AlphaPhy, data=fulldataBG))
summary(lm(herbivory_resistance ~ Med_AlphaPhy, data=fulldataBGfiltered))
summary(glmmTMB(herbivory_resistance ~ Med_AlphaPhy + (1|BG_ID), family=gaussian(), data=fulldataBG))

summary(lm(herbivory_resistance ~ log(GammaPhy), data=fulldataBG))
summary(lm(herbivory_resistance ~ log(GammaPhy), data=fulldataBGfiltered))
summary(glmmTMB(herbivory_resistance ~ log(GammaPhy) + (1|BG_ID), family=gaussian(), data=fulldataBG))

#Test relationship between competition index and phylogenetic diversity

summary(lm(competition_index ~ Mean_AlphaPhy, data=fulldataBG))
summary(lm(competition_index ~ Mean_AlphaPhy, data=fulldataBGfiltered))
summary(glmmTMB(competition_index ~ Mean_AlphaPhy + (1|BG_ID), family=gaussian(), data=fulldataBG))

summary(lm(competition_index ~ Max_AlphaPhy, data=fulldataBG))
summary(lm(competition_index ~ Max_AlphaPhy, data=fulldataBGfiltered))
summary(glmmTMB(competition_index ~ Max_AlphaPhy + (1|BG_ID), family=gaussian(), data=fulldataBG))

summary(lm(competition_index ~ Med_AlphaPhy, data=fulldataBG))
summary(lm(competition_index ~ Med_AlphaPhy, data=fulldataBGfiltered))
summary(glmmTMB(competition_index ~ Med_AlphaPhy + (1|BG_ID), family=gaussian(), data=fulldataBG))

summary(lm(competition_index ~ log(GammaPhy), data=fulldataBG))
summary(lm(competition_index ~ log(GammaPhy), data=fulldataBGfiltered))
summary(glmmTMB(competition_index ~ log(GammaPhy) + (1|BG_ID), family=gaussian(), data=fulldataBG))

#Test relationship with hardiness index 

Anova(lm(Mean_AlphaPhy ~ hardiness, data=fulldataBG), type="III")
Anova(lm(Mean_AlphaPhy ~ hardiness, data=fulldataBGfiltered), type="III")
Anova(glmmTMB(Mean_AlphaPhy ~ hardiness + (1|BG_ID), family=gaussian(), data=fulldataBG), type="III")

Anova(lm(Max_AlphaPhy ~ hardiness, data=fulldataBG), type="III")
Anova(lm(Max_AlphaPhy ~ hardiness, data=fulldataBGfiltered), type="III")
Anova(glmmTMB(Max_AlphaPhy ~ hardiness + (1|BG_ID), family=gaussian(), data=fulldataBG), type="III")

Anova(lm(Med_AlphaPhy ~ hardiness, data=fulldataBG), type="III")
Anova(lm(Med_AlphaPhy ~ hardiness, data=fulldataBGfiltered), type="III")
Anova(glmmTMB(Med_AlphaPhy ~ hardiness + (1|BG_ID), family=gaussian(), data=fulldataBG), type="III")

Anova(lm(log(GammaPhy) ~ hardiness, data=fulldataBG), type="III")
Anova(lm(log(GammaPhy) ~ hardiness, data=fulldataBGfiltered), type="III")
Anova(glmmTMB(log(GammaPhy) ~ hardiness + (1|BG_ID), family=gaussian(), data=fulldataBG), type="III")

#Annual or perennial

Anova(lm(Mean_AlphaPhy ~ Life_span, data=fulldataBG), type="III")
Anova(lm(Mean_AlphaPhy ~ Life_span, data=fulldataBGfiltered), type="III")
Anova(glmmTMB(Mean_AlphaPhy ~ Life_span + (1|BG_ID), family=gaussian(), data=fulldataBG), type="III")

Anova(lm(Max_AlphaPhy ~ Life_span, data=fulldataBG), type="III")
Anova(lm(Max_AlphaPhy ~ Life_span, data=fulldataBGfiltered), type="III")
Anova(glmmTMB(Max_AlphaPhy ~ Life_span + (1|BG_ID), family=gaussian(), data=fulldataBG), type="III")

Anova(lm(Med_AlphaPhy ~ Life_span, data=fulldataBG), type="III")
Anova(lm(Med_AlphaPhy ~ Life_span, data=fulldataBGfiltered), type="III")
Anova(glmmTMB(Med_AlphaPhy ~ Life_span + (1|BG_ID), family=gaussian(), data=fulldataBG), type="III")

Anova(lm(log(GammaPhy) ~ Life_span, data=fulldataBG), type="III")
Anova(lm(log(GammaPhy) ~ Life_span, data=fulldataBGfiltered), type="III")
Anova(glmmTMB(log(GammaPhy) ~ Life_span + (1|BG_ID), family=gaussian(), data=fulldataBG), type="III")

#Test if phylogenetic diversity linked to seed weight

summary(lm(log(Seed_wt) ~ Mean_AlphaPhy, data=fulldataBG))
summary(lm(log(Seed_wt) ~ Mean_AlphaPhy, data=fulldataBGfiltered))
summary(glmmTMB(log(Seed_wt) ~ Mean_AlphaPhy + (1|BG_ID), family=gaussian(), data=fulldataBG))

meanseedweight <- ggplot(fulldataBG, aes(x=Mean_AlphaPhy, y=Seed_wt)) +
  geom_point() + geom_smooth(method="lm", colour="black") +
  xlab("Mean Alpha Phylogenetic Diversity") + ylab("Seed weight (mg)") +
  theme_bw() + scale_y_continuous(trans='log10') 
meanseedweight

summary(lm(log(Seed_wt) ~ Max_AlphaPhy, data=fulldataBG))
summary(lm(log(Seed_wt) ~ Max_AlphaPhy, data=fulldataBGfiltered))
summary(glmmTMB(log(Seed_wt) ~ Max_AlphaPhy + (1|BG_ID), family=gaussian(), data=fulldataBG))

summary(lm(log(Seed_wt) ~ Med_AlphaPhy, data=fulldataBG))
summary(lm(log(Seed_wt) ~ Med_AlphaPhy, data=fulldataBGfiltered))
summary(glmmTMB(log(Seed_wt) ~ Med_AlphaPhy + (1|BG_ID), family=gaussian(), data=fulldataBG))

summary(lm(log(Seed_wt) ~ log(GammaPhy), data=fulldataBG))
summary(lm(log(Seed_wt) ~ log(GammaPhy), data=fulldataBGfiltered))
summary(glmmTMB(log(Seed_wt) ~ log(GammaPhy) + (1|BG_ID), family=gaussian(), data=fulldataBG))

#Test if optimal germination rate linked to phylogenetic distance

summary(glm(OptGermRt ~ Mean_AlphaPhy, data=fulldataBG, family=binomial(link="logit")))
summary(glm(OptGermRt ~ Mean_AlphaPhy, data=fulldataBGfiltered, family=binomial(link="logit")))
summary(glmmTMB(OptGermRt ~ Mean_AlphaPhy + (1|BG_ID), data=fulldataBG, family=binomial(link="logit")))

summary(glm(OptGermRt ~ Max_AlphaPhy, data=fulldataBG, family=binomial(link="logit")))
summary(glm(OptGermRt ~ Max_AlphaPhy, data=fulldataBGfiltered, family=binomial(link="logit")))
summary(glmmTMB(OptGermRt ~ Max_AlphaPhy + (1|BG_ID), data=fulldataBG, family=binomial(link="logit")))

summary(glm(OptGermRt ~ Med_AlphaPhy, data=fulldataBG, family=binomial(link="logit")))
summary(glm(OptGermRt ~ Med_AlphaPhy, data=fulldataBGfiltered, family=binomial(link="logit")))
summary(glmmTMB(OptGermRt ~ Med_AlphaPhy + (1|BG_ID), data=fulldataBG, family=binomial(link="logit")))

summary(glm(OptGermRt ~ log(GammaPhy), data=fulldataBG, family=binomial(link="logit")))
summary(glm(OptGermRt ~ log(GammaPhy), data=fulldataBGfiltered, family=binomial(link="logit")))
summary(glmmTMB(OptGermRt ~ log(GammaPhy) + (1|BG_ID), data=fulldataBG, family=binomial(link="logit")))

#Test if phylogenetic diversity linked to vegetative height

summary(lm(log(vegheight) ~ Mean_AlphaPhy, data=fulldataBG))
summary(lm(log(vegheight) ~ Mean_AlphaPhy, data=fulldataBGfiltered))
summary(glmmTMB(log(vegheight) ~ Mean_AlphaPhy + (1|BG_ID), family=gaussian(), data=fulldataBG))

summary(lm(log(vegheight) ~ Max_AlphaPhy, data=fulldataBG))
summary(lm(log(vegheight) ~ Max_AlphaPhy, data=fulldataBGfiltered))
summary(glmmTMB(log(vegheight) ~ Max_AlphaPhy + (1|BG_ID), family=gaussian(), data=fulldataBG))

summary(lm(log(vegheight) ~ Med_AlphaPhy, data=fulldataBG))
summary(lm(log(vegheight) ~ Med_AlphaPhy, data=fulldataBGfiltered))
summary(glmmTMB(log(vegheight) ~ Med_AlphaPhy + (1|BG_ID), family=gaussian(), data=fulldataBG))

summary(lm(log(vegheight)~ log(GammaPhy), data=fulldataBG))
summary(lm(log(vegheight) ~ log(GammaPhy), data=fulldataBGfiltered))
summary(glmmTMB(log(vegheight) ~ log(GammaPhy) + (1|BG_ID), family=gaussian(), data=fulldataBG))

#Test if phylogenetic diversity linked to generative height

summary(lm(log(genheight) ~ Mean_AlphaPhy, data=fulldataBG))
summary(lm(log(genheight) ~ Mean_AlphaPhy, data=fulldataBGfiltered))
summary(glmmTMB(log(genheight) ~ Mean_AlphaPhy + (1|BG_ID), family=gaussian(), data=fulldataBG))

summary(lm(log(genheight) ~ Max_AlphaPhy, data=fulldataBG))
summary(lm(log(genheight) ~ Max_AlphaPhy, data=fulldataBGfiltered))
summary(glmmTMB(log(genheight) ~ Max_AlphaPhy + (1|BG_ID), family=gaussian(), data=fulldataBG))

summary(lm(log(genheight) ~ Med_AlphaPhy, data=fulldataBG))
summary(lm(log(genheight) ~ Med_AlphaPhy, data=fulldataBGfiltered))
summary(glmmTMB(log(genheight) ~ Med_AlphaPhy + (1|BG_ID), family=gaussian(), data=fulldataBG))

summary(lm(log(genheight)~ log(GammaPhy), data=fulldataBG))
summary(lm(log(genheight) ~ log(GammaPhy), data=fulldataBGfiltered))
summary(glmmTMB(log(genheight) ~ log(GammaPhy) + (1|BG_ID), family=gaussian(), data=fulldataBG))

#Test if phylogenetic diversity linked to leaf N

summary(lm(leafN ~ Mean_AlphaPhy, data=fulldataBG))
summary(lm(leafN ~ Mean_AlphaPhy, data=fulldataBGfiltered))
summary(glmmTMB(leafN ~ Mean_AlphaPhy + (1|BG_ID), family=gaussian(), data=fulldataBG))

meanleafN <- ggplot(fulldataBG, aes(x=Mean_AlphaPhy, y=leafN)) +
  geom_point() + geom_smooth(method="lm", colour="black") +
  xlab("Mean Alpha Phylogenetic Diversity") + ylab("Leaf N (mg/g)") +
  theme_bw() 
meanleafN

summary(lm(leafN ~ Max_AlphaPhy, data=fulldataBG))
summary(lm(leafN ~ Max_AlphaPhy, data=fulldataBGfiltered))
summary(glmmTMB(leafN ~ Max_AlphaPhy + (1|BG_ID), family=gaussian(), data=fulldataBG))

summary(lm(leafN ~ Med_AlphaPhy, data=fulldataBG))
summary(lm(leafN ~ Med_AlphaPhy, data=fulldataBGfiltered))
summary(glmmTMB(leafN ~ Med_AlphaPhy + (1|BG_ID), family=gaussian(), data=fulldataBG))

summary(lm(leafN ~ log(GammaPhy), data=fulldataBG))
summary(lm(leafN ~ log(GammaPhy), data=fulldataBGfiltered))
summary(glmmTMB(leafN ~ log(GammaPhy) + (1|BG_ID), family=gaussian(), data=fulldataBG))

#Test if phylogenetic diversity linked to LDMC

summary(lm(LDMC ~ Mean_AlphaPhy, data=fulldataBG))
summary(lm(LDMC ~ Mean_AlphaPhy, data=fulldataBGfiltered))
summary(glmmTMB(LDMC ~ Mean_AlphaPhy + (1|BG_ID), family=gaussian(), data=fulldataBG))

meanLDMC <- ggplot(fulldataBG, aes(x=Mean_AlphaPhy, y=LDMC)) +
  geom_point() + geom_smooth(method="lm", colour="black") +
  xlab("Mean Alpha Phylogenetic Diversity") + ylab("LDMC (g/g)") +
  theme_bw() 
meanLDMC

summary(lm(LDMC ~ Max_AlphaPhy, data=fulldataBG))
summary(lm(LDMC ~ Max_AlphaPhy, data=fulldataBGfiltered))
summary(glmmTMB(LDMC ~ Max_AlphaPhy + (1|BG_ID), family=gaussian(), data=fulldataBG))

summary(lm(LDMC ~ Med_AlphaPhy, data=fulldataBG))
summary(lm(LDMC ~ Med_AlphaPhy, data=fulldataBGfiltered))
summary(glmmTMB(LDMC ~ Med_AlphaPhy + (1|BG_ID), family=gaussian(), data=fulldataBG))

summary(lm(LDMC ~ log(GammaPhy), data=fulldataBG))
summary(lm(LDMC ~ log(GammaPhy), data=fulldataBGfiltered))
summary(glmmTMB(LDMC ~ log(GammaPhy) + (1|BG_ID), family=gaussian(), data=fulldataBG))

#Test if phylogenetic diversity linked to SLA

summary(lm(SLA ~ Mean_AlphaPhy, data=fulldataBG))
summary(lm(SLA ~ Mean_AlphaPhy, data=fulldataBGfiltered))
summary(glmmTMB(SLA ~ Mean_AlphaPhy + (1|BG_ID), family=gaussian(), data=fulldataBG))

summary(lm(SLA ~ Max_AlphaPhy, data=fulldataBG))
summary(lm(SLA ~ Max_AlphaPhy, data=fulldataBGfiltered))
summary(glmmTMB(SLA ~ Max_AlphaPhy + (1|BG_ID), family=gaussian(), data=fulldataBG))

summary(lm(SLA ~ Med_AlphaPhy, data=fulldataBG))
summary(lm(SLA ~ Med_AlphaPhy, data=fulldataBGfiltered))
summary(glmmTMB(SLA ~ Med_AlphaPhy + (1|BG_ID), family=gaussian(), data=fulldataBG))

summary(lm(SLA ~ log(GammaPhy), data=fulldataBG))
summary(lm(SLA ~ log(GammaPhy), data=fulldataBGfiltered))
summary(glmmTMB(SLA ~ log(GammaPhy) + (1|BG_ID), family=gaussian(), data=fulldataBG))

Figure_6 <- ggarrange(meanLDMC, meanseedweight, meanleafN,
                      nrow=2, ncol=2)
Figure_6
ggsave(Figure_6, 
       filename = "Figure_6.svg",
       height = 210, width = 210, units = "mm")

##################################################################################################

#SUPPLEMENTARY ANALYSES

#As per main text results, final selected models are presented

##################################################################################################

#Re-do Question 1 models WITHOUT HAUESER to see if trends still hold
#K+M models in Tables S2-S5

germinationv2 <- readRDS("germinationrevised.RDS") %>%
  mutate(suitabilitypc1V2abs = abs(suitabilitypc1V2)) %>%
  mutate(suitabilitypc2V2abs = abs(suitabilitypc2V2)) %>%
  filter(Study != "Haueser")

germinationv2$Disturbance <- relevel(germinationv2$Disturbance, ref="0")

surviveyr1v2 <- readRDS("surviveyr1revised.RDS") %>%
  mutate(suitabilitypc1V2abs = abs(suitabilitypc1V2)) %>%
  mutate(suitabilitypc2V2abs = abs(suitabilitypc2V2))  %>%
  filter(Study != "Haueser")

surviveyr1v2$Disturbance <- relevel(surviveyr1v2$Disturbance, ref="0")

surviveyr2earlyv2 <- readRDS("surviveyr2earlyrevised.RDS") %>%
  mutate(suitabilitypc1V2abs = abs(suitabilitypc1V2)) %>%
  mutate(suitabilitypc2V2abs = abs(suitabilitypc2V2))  %>%
  filter(Study != "Haueser")

surviveyr2earlyv2$Disturbance <- relevel(surviveyr2earlyv2$Disturbance, ref="0")

surviveyr2latev2 <- readRDS("surviveyr2laterevised.RDS") %>%
  mutate(suitabilitypc1V2abs = abs(suitabilitypc1V2)) %>%
  mutate(suitabilitypc2V2abs = abs(suitabilitypc2V2))  %>%
  filter(Study != "Haueser")

surviveyr2latev2$Disturbance <- relevel(surviveyr2latev2$Disturbance, ref="0")

#DO SURVIVAL Y/N FOR ALL DATASETS

#GERMINATION YES NO

#Mean alpha diversity 

germ_YN_Mean_AlphaPhy <- glmmTMB(germ_YN ~ Mean_AlphaPhy + Disturbance + propagule_pressure +
                                   suitabilitypc1V2abs + suitabilitypc2V2abs +
                                   (1|Family) + (1|POWO.name) + 
                                   (1|Site) + offset(log(Density)), 
                                 family=binomial(link="logit"), 
                                 data=germinationv2)
summary(germ_YN_Mean_AlphaPhy)
Anova(germ_YN_Mean_AlphaPhy, type="III")
confint(germ_YN_Mean_AlphaPhy)
r.squaredGLMM(germ_YN_Mean_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = germ_YN_Mean_AlphaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(germ_YN_Mean_AlphaPhy)
check_singularity(germ_YN_Mean_AlphaPhy)

#Maximum alpha diversity

germ_YN_Max_AlphaPhy <- glmmTMB(germ_YN ~ Max_AlphaPhy + Disturbance + propagule_pressure +
                                  suitabilitypc1V2abs + suitabilitypc2V2abs +
                                  (1|Family) + (1|POWO.name) +
                                  (1|Site) + offset(log(Density)), 
                                family=binomial(link="logit"), 
                                data=germinationv2)
summary(germ_YN_Max_AlphaPhy)
Anova(germ_YN_Max_AlphaPhy, type="III")
confint(germ_YN_Max_AlphaPhy)
r.squaredGLMM(germ_YN_Max_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = germ_YN_Max_AlphaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(germ_YN_Max_AlphaPhy)
check_singularity(germ_YN_Max_AlphaPhy)

#Gamma diversity

germ_YN_GammaPhy <- glmmTMB(germ_YN ~ log(GammaPhy)*suitabilitypc1V2abs +  + Disturbance + propagule_pressure +
                              suitabilitypc2V2abs +
                              (1|Family) + (1|POWO.name) +
                              (1|Site) + offset(log(Density)), 
                            family=binomial(link="logit"), 
                            data=germinationv2)
summary(germ_YN_GammaPhy)
Anova(germ_YN_GammaPhy, type="III")
confint(germ_YN_GammaPhy)
r.squaredGLMM(germ_YN_GammaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = germ_YN_GammaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(germ_YN_GammaPhy)
check_singularity(germ_YN_GammaPhy)

#Median alpha diversity 

germ_YN_MedPhy <- glmmTMB(germ_YN ~ Med_AlphaPhy*Disturbance + propagule_pressure +
                            suitabilitypc1V2abs + suitabilitypc2V2abs +
                            (1|Family) + (1|POWO.name) +
                            (1|Site) + offset(log(Density)), 
                          family=binomial(link="logit"), 
                          data=germination)
summary(germ_YN_MedPhy)
Anova(germ_YN_MedPhy, type="III")
confint(germ_YN_MedPhy)
r.squaredGLMM(germ_YN_MedPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = germ_YN_MedPhy, plot=F)
plot(simulationOutput)
check_overdispersion(germ_YN_MedPhy)
check_singularity(germ_YN_MedPhy)

#FIRST YEAR SURVIVAL YES NO

#Mean alpha diversity 

surviveyr1_YN_Mean_AlphaPhy <- glmmTMB(surviveyr1_YN ~ Mean_AlphaPhy*suitabilitypc2V2abs+Disturbance + propagule_pressure +
                                         + suitabilitypc1V2abs + Herbivory +
                                         (1|Family) + (1|POWO.name) + 
                                         (1|Site) + offset(log(Density)), 
                                       family=binomial(link="logit"), 
                                       data=surviveyr1v2)
summary(surviveyr1_YN_Mean_AlphaPhy)
Anova(surviveyr1_YN_Mean_AlphaPhy, type="III")
confint(surviveyr1_YN_Mean_AlphaPhy)
r.squaredGLMM(surviveyr1_YN_Mean_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr1_YN_Mean_AlphaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr1_YN_Mean_AlphaPhy)
check_singularity(surviveyr1_YN_Mean_AlphaPhy)

#Maximum alpha diversity

surviveyr1_YN_Max_AlphaPhy <- glmmTMB(surviveyr1_YN ~ Max_AlphaPhy*Disturbance + propagule_pressure +
                                        suitabilitypc1V2abs + suitabilitypc2V2abs + Herbivory +
                                        (1|Family) + (1|POWO.name) + 
                                        (1|Site) + offset(log(Density)), 
                                      family=binomial(link="logit"), 
                                      data=surviveyr1v2)
summary(surviveyr1_YN_Max_AlphaPhy)
Anova(surviveyr1_YN_Max_AlphaPhy, type="III")
confint(surviveyr1_YN_Max_AlphaPhy)
r.squaredGLMM(surviveyr1_YN_Max_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr1_YN_Max_AlphaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr1_YN_Max_AlphaPhy)
check_singularity(surviveyr1_YN_Max_AlphaPhy)

#Gamma diversity 

surviveyr1_YN_GammaPhy <- glmmTMB(surviveyr1_YN ~ log(GammaPhy) + Disturbance + propagule_pressure +
                                    suitabilitypc1V2abs + suitabilitypc2V2abs + Herbivory +
                                    (1|Family) + (1|POWO.name) + 
                                    (1|Site) + offset(log(Density)), 
                                  family=binomial(link="logit"), 
                                  data=surviveyr1v2)
summary(surviveyr1_YN_GammaPhy)
Anova(surviveyr1_YN_GammaPhy, type="III")
confint(surviveyr1_YN_GammaPhy)
r.squaredGLMM(surviveyr1_YN_GammaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr1_YN_GammaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr1_YN_GammaPhy)
check_singularity(surviveyr1_YN_GammaPhy)

#Median alpha diversity 

surviveyr1_YN_MedPhy <- glmmTMB(surviveyr1_YN ~ Med_AlphaPhy*suitabilitypc2V2abs + Disturbance + propagule_pressure +
                                  suitabilitypc1V2abs + Herbivory +
                                  (1|Family) + (1|POWO.name) + 
                                  (1|Site) + offset(log(Density)), 
                                family=binomial(link="logit"), 
                                data=surviveyr1v2)
summary(surviveyr1_YN_MedPhy)
Anova(surviveyr1_YN_MedPhy, type="III")
confint(surviveyr1_YN_MedPhy)
r.squaredGLMM(surviveyr1_YN_MedPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr1_YN_MedPhy, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr1_YN_MedPhy)
check_singularity(surviveyr1_YN_MedPhy)

#SECOND YEAR EARLY SURVIVAL YES NO

#Mean alpha diversity

surviveyr2early_YN_Mean_AlphaPhy <- glmmTMB(surviveyr2early_YN ~ Mean_AlphaPhy + Disturbance + propagule_pressure +
                                              suitabilitypc1V2abs + suitabilitypc2V2abs + Herbivory +
                                              (1|Family) + (1|POWO.name) + 
                                              (1|Site) + offset(log(Density)), 
                                            family=binomial(link="logit"), 
                                            data=surviveyr2earlyv2)
summary(surviveyr2early_YN_Mean_AlphaPhy)
Anova(surviveyr2early_YN_Mean_AlphaPhy, type="III")
confint(surviveyr2early_YN_Mean_AlphaPhy)
r.squaredGLMM(surviveyr2early_YN_Mean_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2early_YN_Mean_AlphaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr2early_YN_Mean_AlphaPhy)
check_singularity(surviveyr2early_YN_Mean_AlphaPhy)

# Maximum alpha diversity

surviveyr2early_YN_Max_AlphaPhy <- glmmTMB(surviveyr2early_YN ~ Max_AlphaPhy*suitabilitypc2V2abs +  Disturbance + propagule_pressure +
                                             suitabilitypc1V2abs + Herbivory +
                                             (1|Family) + (1|POWO.name) + 
                                             (1|Site) + offset(log(Density)), 
                                           family=binomial(link="logit"), 
                                           data=surviveyr2earlyv2)
summary(surviveyr2early_YN_Max_AlphaPhy)
Anova(surviveyr2early_YN_Max_AlphaPhy, type="III")
confint(surviveyr2early_YN_Max_AlphaPhy)
r.squaredGLMM(surviveyr2early_YN_Max_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2early_YN_Max_AlphaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr2early_YN_Max_AlphaPhy)
check_singularity(surviveyr2early_YN_Max_AlphaPhy)

#Gamma diversity 

surviveyr2early_YN_GammaPhy <- glmmTMB(surviveyr2early_YN ~ log(GammaPhy) + suitabilitypc2V2abs + Disturbance + propagule_pressure +
                                         suitabilitypc1V2abs + Herbivory +
                                         (1|Family) + (1|POWO.name) + 
                                         (1|Site) + offset(log(Density)), 
                                       family=binomial(link="logit"), 
                                       data=surviveyr2earlyv2)
summary(surviveyr2early_YN_GammaPhy)
Anova(surviveyr2early_YN_GammaPhy, type="III")
confint(surviveyr2early_YN_GammaPhy)
r.squaredGLMM(surviveyr2early_YN_GammaPhy)
#interaction with PC2 0.051

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2early_YN_GammaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr2early_YN_GammaPhy)
check_singularity(surviveyr2early_YN_GammaPhy)

#Median alpha diversity 

surviveyr2early_YN_MedPhy <- glmmTMB(surviveyr2early_YN ~ Med_AlphaPhy*suitabilitypc2V2abs + Disturbance + propagule_pressure +
                                       suitabilitypc1V2abs + Herbivory +
                                       (1|Family) + (1|POWO.name) + 
                                       (1|Site) + offset(log(Density)), 
                                     family=binomial(link="logit"), 
                                     data=surviveyr2earlyv2)
summary(surviveyr2early_YN_MedPhy)
Anova(surviveyr2early_YN_MedPhy, type="III")
confint(surviveyr2early_YN_MedPhy)
r.squaredGLMM(surviveyr2early_YN_MedPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2early_YN_MedPhy, plot=F)
plot(simulationOutput)
testOutliers(simulationOutput)
check_overdispersion(surviveyr2early_YN_MedPhy)
check_singularity(surviveyr2early_YN_MedPhy)

check_model(surviveyr2early_YN_MedPhy)

#SECOND YEAR LATE SURVIVAL YES NO

#Mean alpha diversity 

surviveyr2late_YN_Mean_AlphaPhy <- glmmTMB(surviveyr2late_YN ~ Mean_AlphaPhy + Disturbance + propagule_pressure +
                                             suitabilitypc1V2abs + suitabilitypc2V2abs + 
                                             (1|Family) + (1|POWO.name) + 
                                             (1|Site) + offset(log(Density)), 
                                           family=binomial(link="logit"), 
                                           data=surviveyr2latev2)
summary(surviveyr2late_YN_Mean_AlphaPhy)
Anova(surviveyr2late_YN_Mean_AlphaPhy, type="III")
confint(surviveyr2late_YN_Mean_AlphaPhy)
r.squaredGLMM(surviveyr2late_YN_Mean_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2late_YN_Mean_AlphaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr2late_YN_Mean_AlphaPhy)
check_singularity(surviveyr2late_YN_Mean_AlphaPhy)

#Maximum alpha diversity

surviveyr2late_YN_Max_AlphaPhy <- glmmTMB(surviveyr2late_YN ~ Max_AlphaPhy + Disturbance + propagule_pressure +
                                            suitabilitypc1V2abs + suitabilitypc2V2abs + 
                                            (1|Family) + (1|POWO.name) + 
                                            (1|Site) + offset(log(Density)), 
                                          family=binomial(link="logit"), 
                                          data=surviveyr2latev2)
summary(surviveyr2late_YN_Max_AlphaPhy)
Anova(surviveyr2late_YN_Max_AlphaPhy, type="III")
confint(surviveyr2late_YN_Max_AlphaPhy)
r.squaredGLMM(surviveyr2late_YN_Max_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2late_YN_Max_AlphaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr2late_YN_Max_AlphaPhy)
check_singularity(surviveyr2late_YN_Max_AlphaPhy)

#Gamma diversity 

surviveyr2late_YN_GammaPhy <- glmmTMB(surviveyr2late_YN ~ GammaPhy + Disturbance + propagule_pressure +
                                        suitabilitypc1V2abs + suitabilitypc2V2abs + 
                                        (1|POWO.name) + 
                                        (1|Site) + offset(log(Density)), 
                                      family=binomial(link="logit"), 
                                      data=surviveyr2latev2)
summary(surviveyr2late_YN_GammaPhy)
Anova(surviveyr2late_YN_GammaPhy, type="III")
confint(surviveyr2late_YN_GammaPhy)
r.squaredGLMM(surviveyr2late_YN_GammaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2late_YN_GammaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr2late_YN_GammaPhy)
check_singularity(surviveyr2late_YN_GammaPhy)

#Median alpha diversity 

surviveyr2late_YN_MedPhy <- glmmTMB(surviveyr2late_YN ~ Med_AlphaPhy + Disturbance + propagule_pressure +
                                      suitabilitypc1V2abs + suitabilitypc2V2abs + 
                                      (1|Family) + (1|POWO.name) + 
                                      (1|Site) + offset(log(Density)), 
                                    family=binomial(link="logit"), 
                                    data=surviveyr2latev2)
summary(surviveyr2late_YN_MedPhy)
Anova(surviveyr2late_YN_MedPhy, type="III")
confint(surviveyr2late_YN_MedPhy)
r.squaredGLMM(surviveyr2late_YN_MedPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2late_YN_MedPhy, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr2late_YN_MedPhy)
check_singularity(surviveyr2late_YN_MedPhy)

#######################################################################################

#DO NUMBER SURVIVING FOR ALL DATASETS

#GERMINATION COUNT

germinationcountv2 <- filter(germinationv2, germ_num>0)

#Mean alpha diversity

germ_count_Mean_AlphaPhy <- glmmTMB(germ_num ~ Mean_AlphaPhy*Disturbance + propagule_pressure +
                                      suitabilitypc1V2abs + suitabilitypc2V2abs + Density +
                                      (1|Family) + (1|POWO.name) + (1|Site),
                                    family=nbinom2(link = "log"), 
                                    data=germinationcountv2)
summary(germ_count_Mean_AlphaPhy)
Anova(germ_count_Mean_AlphaPhy, type="III")
confint(germ_count_Mean_AlphaPhy)
r.squaredGLMM(germ_count_Mean_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = germ_count_Mean_AlphaPhy, plot=F)
plot(simulationOutput)

check_overdispersion(germ_count_Mean_AlphaPhy)
check_singularity(germ_count_Mean_AlphaPhy)
check_model(germ_count_Mean_AlphaPhy)

#Maximum alpha diversity

germ_count_Max_AlphaPhy <- glmmTMB(germ_num ~ Max_AlphaPhy*Disturbance + propagule_pressure +
                                     suitabilitypc1V2abs + suitabilitypc2V2abs + Density +
                                     (1|Family) + (1|POWO.name) + (1|Site),
                                   family=nbinom2(link = "log"), 
                                   data=germinationcountv2)
summary(germ_count_Max_AlphaPhy)
Anova(germ_count_Max_AlphaPhy, type="III")
confint(germ_count_Max_AlphaPhy)
r.squaredGLMM(germ_count_Max_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = germ_count_Max_AlphaPhy, plot=F)
plot(simulationOutput)

check_overdispersion(germ_count_Max_AlphaPhy)
check_singularity(germ_count_Max_AlphaPhy)

#Gamma diversity 

germ_count_GammaPhy <- glmmTMB(germ_num ~ log(GammaPhy) + Disturbance + propagule_pressure +
                                 suitabilitypc1V2abs + suitabilitypc2V2abs + Density +
                                 (1|Family) + (1|POWO.name) + (1|Site),
                               family=nbinom2(link = "log"), 
                               data=germinationcountv2)
summary(germ_count_GammaPhy)
Anova(germ_count_GammaPhy, type="III")
confint(germ_count_GammaPhy)
r.squaredGLMM(germ_count_GammaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = germ_count_GammaPhy, plot=F)
plot(simulationOutput)

check_overdispersion(germ_count_GammaPhy)
check_singularity(germ_count_GammaPhy)

#Median Alpha diversity

germ_count_MedPhy <- glmmTMB(germ_num ~ Med_AlphaPhy*Disturbance + propagule_pressure +
                               suitabilitypc2V2abs + suitabilitypc1V2abs + Density +
                               (1|Family) + (1|POWO.name) + (1|Site),
                             family=nbinom2(link = "log"), 
                             data=germinationcountv2)
summary(germ_count_MedPhy)
Anova(germ_count_MedPhy, type="III")
confint(germ_count_MedPhy)
r.squaredGLMM(germ_count_MedPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = germ_count_MedPhy, plot=F)
plot(simulationOutput)

check_overdispersion(germ_count_MedPhy)
check_singularity(germ_count_MedPhy)

#FIRST YEAR SURVIVAL COUNT

surviveyr1countv2 <- filter(surviveyr1v2, surviveyr1_num>0)

#Mean alpha diversity 

surviveyr1_count_Mean_AlphaPhy <- glmmTMB(surviveyr1_num ~ Mean_AlphaPhy*Disturbance + propagule_pressure +
                                            suitabilitypc1V2abs +  Density + Herbivory +
                                            (1|POWO.name) + (1|Site),
                                          family=nbinom2(link = "log"), 
                                          data=surviveyr1countv2)
summary(surviveyr1_count_Mean_AlphaPhy)
Anova(surviveyr1_count_Mean_AlphaPhy, type="III")
confint(surviveyr1_count_Mean_AlphaPhy)
r.squaredGLMM(surviveyr1_count_Mean_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr1_count_Mean_AlphaPhy, plot=F)
plot(simulationOutput)

check_overdispersion(surviveyr1_count_Mean_AlphaPhy)
check_singularity(surviveyr1_count_Mean_AlphaPhy)

#Maximum alpha diversity

surviveyr1_count_Max_AlphaPhy <- glmmTMB(surviveyr1_num ~ Max_AlphaPhy*Disturbance + propagule_pressure +
                                           suitabilitypc1V2abs + suitabilitypc2V2abs + Density + Herbivory +
                                           (1|POWO.name) + (1|Site),
                                         family=nbinom2(link = "log"), 
                                         data=surviveyr1countv2)
summary(surviveyr1_count_Max_AlphaPhy)
Anova(surviveyr1_count_Max_AlphaPhy, type="III")
confint(surviveyr1_count_Max_AlphaPhy)
r.squaredGLMM(surviveyr1_count_Max_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr1_count_Max_AlphaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr1_count_Max_AlphaPhy)
check_singularity(surviveyr1_count_Max_AlphaPhy)

#Gamma diversity

surviveyr1_count_GammaPhy <- glmmTMB(surviveyr1_num ~ log(GammaPhy) + Disturbance + propagule_pressure +
                                       suitabilitypc1V2abs + suitabilitypc2V2abs + Density + Herbivory +
                                       (1|POWO.name) + (1|Site),
                                     family=nbinom2(link = "log"), 
                                     data=surviveyr1countv2)
summary(surviveyr1_count_GammaPhy)
Anova(surviveyr1_count_GammaPhy, type="III")
confint(surviveyr1_count_GammaPhy)
r.squaredGLMM(surviveyr1_count_GammaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr1_count_GammaPhy, plot=F)
plot(simulationOutput)

check_overdispersion(surviveyr1_count_GammaPhy)
check_singularity(surviveyr1_count_GammaPhy)

#Median alpha diversity

surviveyr1_count_MedPhy <- glmmTMB(surviveyr1_num ~ Med_AlphaPhy*Disturbance + propagule_pressure +
                                     suitabilitypc1V2abs + suitabilitypc2V2abs + Density + Herbivory +
                                     (1|POWO.name) + (1|Site),
                                   family=nbinom2(link = "log"), 
                                   data=surviveyr1countv2)
summary(surviveyr1_count_MedPhy)
Anova(surviveyr1_count_MedPhy, type="III")
confint(surviveyr1_count_MedPhy)
r.squaredGLMM(surviveyr1_count_MedPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr1_count_MedPhy, plot=F)
plot(simulationOutput)

check_overdispersion(surviveyr1_count_MedPhy)
check_singularity(surviveyr1_count_MedPhy)

#EARLY SECOND YEAR SURVIVAL COUNT

surviveyr2earlycountv2 <- filter(surviveyr2earlyv2, surviveyr2early_num>0)

#Mean alpha diversity 

surviveyr2early_count_Mean_AlphaPhy <- glmmTMB(surviveyr2early_num ~ Mean_AlphaPhy + Disturbance + propagule_pressure +
                                                 suitabilitypc1V2abs + suitabilitypc2V2abs + Density + Herbivory +
                                                 (1|POWO.name) + (1|Site),
                                               family=nbinom2(link = "log"), 
                                               data=surviveyr2earlycountv2)
summary(surviveyr2early_count_Mean_AlphaPhy)
Anova(surviveyr2early_count_Mean_AlphaPhy, type="III")
confint(surviveyr2early_count_Mean_AlphaPhy)
r.squaredGLMM(surviveyr2early_count_Mean_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2early_count_Mean_AlphaPhy, plot=F)
plot(simulationOutput)

check_overdispersion(surviveyr2early_count_Mean_AlphaPhy)
check_singularity(surviveyr2early_count_Mean_AlphaPhy)

#Max alpha diversity 

surviveyr2early_count_Max_AlphaPhy <- glmmTMB(surviveyr2early_num ~ Max_AlphaPhy + Disturbance + propagule_pressure +
                                                suitabilitypc1V2abs + suitabilitypc2V2abs + Density + Herbivory +
                                                (1|POWO.name) + (1|Site),
                                              family=nbinom2(link = "log"), 
                                              data=surviveyr2earlycountv2)
summary(surviveyr2early_count_Max_AlphaPhy)
Anova(surviveyr2early_count_Max_AlphaPhy, type="III")
confint(surviveyr2early_count_Max_AlphaPhy)
r.squaredGLMM(surviveyr2early_count_Max_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2early_count_Max_AlphaPhy, plot=F)
plot(simulationOutput)

check_overdispersion(surviveyr2early_count_Max_AlphaPhy)
check_singularity(surviveyr2early_count_Max_AlphaPhy)

#Gamma diversity 

surviveyr2early_count_GammaPhy <- glmmTMB(surviveyr2early_num ~ log(GammaPhy) + Disturbance + propagule_pressure +
                                            suitabilitypc1V2abs + suitabilitypc2V2abs + Density + Herbivory +
                                            (1|POWO.name) + (1|Site),
                                          family=nbinom2(link = "log"), 
                                          data=surviveyr2earlycountv2)
summary(surviveyr2early_count_GammaPhy)
Anova(surviveyr2early_count_GammaPhy, type="III")
confint(surviveyr2early_count_GammaPhy)
r.squaredGLMM(surviveyr2early_count_GammaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2early_count_GammaPhy, plot=F)
plot(simulationOutput)

check_overdispersion(surviveyr2early_count_GammaPhy)
check_singularity(surviveyr2early_count_GammaPhy)

#Median alpha diversity 

surviveyr2early_count_MedPhy <- glmmTMB(surviveyr2early_num ~ Med_AlphaPhy + Disturbance + propagule_pressure +
                                          suitabilitypc1V2abs + suitabilitypc2V2abs + Density + Herbivory +
                                          (1|POWO.name) + (1|Site),
                                        family=nbinom2(link = "log"), 
                                        data=surviveyr2earlycountv2)
summary(surviveyr2early_count_MedPhy)
Anova(surviveyr2early_count_MedPhy, type="III")
confint(surviveyr2early_count_MedPhy)
r.squaredGLMM(surviveyr2early_count_MedPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2early_count_MedPhy, plot=F)
plot(simulationOutput)

check_overdispersion(surviveyr2early_count_MedPhy)
check_singularity(surviveyr2early_count_MedPhy)

#LATE SECOND YEAR SURVIVAL COUNT

surviveyr2latecountv2 <- filter(surviveyr2latev2, surviveyr2late_num>0)

#Mean alpha diversity 

surviveyr2late_count_Mean_AlphaPhy <- glmmTMB(surviveyr2late_num ~ Mean_AlphaPhy*suitabilitypc2V2abs + Disturbance +
                                                suitabilitypc1V2abs + suitabilitypc2V2abs +
                                                (1|POWO.name),
                                              family=nbinom2(link = "log"), 
                                              data=surviveyr2latecountv2)
summary(surviveyr2late_count_Mean_AlphaPhy)
Anova(surviveyr2late_count_Mean_AlphaPhy, type="III")
confint(surviveyr2late_count_Mean_AlphaPhy)
r.squaredGLMM(surviveyr2late_count_Mean_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2late_count_Mean_AlphaPhy, plot=F)
plot(simulationOutput)

check_overdispersion(surviveyr2late_count_Mean_AlphaPhy)
check_singularity(surviveyr2late_count_Mean_AlphaPhy)

#Max alpha diversity 

surviveyr2late_count_Max_AlphaPhy <- glmmTMB(surviveyr2late_num ~ Max_AlphaPhy + Disturbance + 
                                               suitabilitypc1V2abs + suitabilitypc2V2abs +
                                               (1|POWO.name),
                                             family=nbinom2(link = "log"), 
                                             data=surviveyr2latecountv2)
summary(surviveyr2late_count_Max_AlphaPhy)
Anova(surviveyr2late_count_Max_AlphaPhy, type="III")
confint(surviveyr2late_count_Max_AlphaPhy)
r.squaredGLMM(surviveyr2late_count_Max_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2late_count_Max_AlphaPhy, plot=F)
plot(simulationOutput)

check_overdispersion(surviveyr2late_count_Max_AlphaPhy)
check_singularity(surviveyr2late_count_Max_AlphaPhy)

#Gamma diversity 

surviveyr2late_count_GammaPhy <- glmmTMB(surviveyr2late_num ~ log(GammaPhy) + Disturbance + 
                                           suitabilitypc1V2abs + suitabilitypc2V2abs + 
                                           (1|POWO.name),
                                         family=nbinom2(link = "log"), 
                                         data=surviveyr2latecountv2)
summary(surviveyr2late_count_GammaPhy)
Anova(surviveyr2late_count_GammaPhy, type="III")
confint(surviveyr2late_count_GammaPhy)
r.squaredGLMM(surviveyr2late_count_GammaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2late_count_GammaPhy, plot=F)
plot(simulationOutput)

check_overdispersion(surviveyr2late_count_GammaPhy)
check_singularity(surviveyr2late_count_GammaPhy)

#median alpha diversity 

surviveyr2late_count_MedPhy <- glmmTMB(surviveyr2late_num ~ Med_AlphaPhy + suitabilitypc1V2abs + Disturbance +
                                         suitabilitypc2V2abs + 
                                         (1|POWO.name),
                                       family=nbinom2(link = "log"), 
                                       data=surviveyr2latecountv2)
summary(surviveyr2late_count_MedPhy)
Anova(surviveyr2late_count_MedPhy, type="III")
confint(surviveyr2late_count_MedPhy)
r.squaredGLMM(surviveyr2late_count_MedPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2late_count_MedPhy, plot=F)
plot(simulationOutput)

check_overdispersion(surviveyr2late_count_MedPhy)
check_singularity(surviveyr2late_count_MedPhy)

#################################################################

#Now doing it again looking at difference at a site level for
#Haueser and Kempel combined
#The K+H models in Tables S2-S5

#Do a version of difference between planted and resident diversity, just using overall average
#at site level using Kempel and Haueser (as do have site-level for Kempel)

#Haueser just average PD per plot

cover_H <- read.csv("plot_cover.csv", header=T, stringsAsFactors = T) %>%
  pivot_longer(cols=Achillea.millefolium:Veronica.persica, names_to="Species", values_to="cover", values_drop_na = T)
cover_H$Species=gsub("\\."," ",cover_H$Species)
cover_H$Species <- as.factor(cover_H$Species)

#Read in diversity metrics Haueser

comm_diversity_H <- read.csv("community_diversities.csv", header=T, stringsAsFactors = T)

#Join the two dataframes and calculate average of phylogenetic diversity Haueser

comm_diversity_H <- left_join(cover_H, comm_diversity_H, by="Species") %>%
  group_by(Plot) %>%
  dplyr::summarise(Max_AlphaPhy_mean = mean(Max_AlphaPhy), 
                   GammaPhy_mean = mean(GammaPhy),
                   Mean_AlphaPhy_mean = mean(Mean_AlphaPhy),
                   Med_AlphaPhy_mean = mean(Med_AlphaPhy)) %>%
  dplyr::rename(Site=Plot)

comm_diversity_H$Site <- plyr::revalue(comm_diversity_H$Site, c("1A"="01A", "1B"="01B", "2A"="02A",
                                                                "2B"="02B", "4A"="04A", "4B"="04B",
                                                                "5A"="05A", "5B"="05B", "6A"="06A",
                                                                "6B"="06B", "7A"="07A", "7B"="07B",
                                                                "8A"="08A", "8B"="08B", "9A"="09A",
                                                                "9B"="09B"))
#Kempel average PD per plot

cover_K <- read.csv("site_cover_Kempel.csv", header=T, stringsAsFactors = T) %>%
  pivot_longer(cols=Bueren:Heimiswil, names_to="Site", values_to="cover", values_drop_na = T)
cover_K$Site <- as.factor(cover_K$Site)

#Read in diversity metrics Kempel

comm_diversity_K <- read.csv("community_diversities_Kempel.csv", header=T, stringsAsFactors = T)

#Join the two dataframes and calculate average of phylogenetic diversity Kempel

comm_diversity_K <- left_join(cover_K, comm_diversity_K, by="Species") %>%
  group_by(Site) %>%
  dplyr::summarise(Max_AlphaPhy_mean = mean(Max_AlphaPhy), 
                   GammaPhy_mean = mean(GammaPhy),
                   Mean_AlphaPhy_mean = mean(Mean_AlphaPhy), 
                   Med_AlphaPhy_mean = mean(Med_AlphaPhy))

#Combine and bind to the datasets

comm_diversity_total <- bind_rows(comm_diversity_H, comm_diversity_K)

germinationNMrevised <- germination %>%
  filter(Study != "Muller") %>%
  droplevels()
germinationNMrevised <- left_join(germinationNMrevised, comm_diversity_total, by="Site") %>%
  mutate(Max_AlphaPhy_diff = Max_AlphaPhy - Max_AlphaPhy_mean,
         Mean_AlphaPhy_diff = Mean_AlphaPhy - Mean_AlphaPhy_mean,
         GammaPhy_diff = GammaPhy - GammaPhy_mean,
         Med_AlphaPhy_diff = Med_AlphaPhy - Med_AlphaPhy_mean) %>%
  mutate(suitabilitypc1V2abs = abs(suitabilitypc1V2)) %>%
  mutate(suitabilitypc2V2abs = abs(suitabilitypc2V2))

surviveyr1NMrevised <- surviveyr1  %>%
  filter(Study != "Muller") %>%
  droplevels()
surviveyr1NMrevised <- left_join(surviveyr1NMrevised, comm_diversity_total, by="Site") %>%
  mutate(Max_AlphaPhy_diff = Max_AlphaPhy - Max_AlphaPhy_mean,
         Mean_AlphaPhy_diff = Mean_AlphaPhy - Mean_AlphaPhy_mean,
         GammaPhy_diff = GammaPhy - GammaPhy_mean,
         Med_AlphaPhy_diff = Med_AlphaPhy - Med_AlphaPhy_mean) %>%
  mutate(suitabilitypc1V2abs = abs(suitabilitypc1V2)) %>%
  mutate(suitabilitypc2V2abs = abs(suitabilitypc2V2))

surviveyr2earlyNMrevised <- surviveyr2early  %>%
  filter(Study != "Muller") %>%
  droplevels()
surviveyr2earlyNMrevised <- left_join(surviveyr2earlyNMrevised, comm_diversity_total, by="Site") %>%
  mutate(Max_AlphaPhy_diff = Max_AlphaPhy - Max_AlphaPhy_mean,
         Mean_AlphaPhy_diff = Mean_AlphaPhy - Mean_AlphaPhy_mean,
         GammaPhy_diff = GammaPhy - GammaPhy_mean,
         Med_AlphaPhy_diff = Med_AlphaPhy - Med_AlphaPhy_mean) %>%
  mutate(suitabilitypc1V2abs = abs(suitabilitypc1V2)) %>%
  mutate(suitabilitypc2V2abs = abs(suitabilitypc2V2))

#surviveyr2late is already just Kempel and Haueser
surviveyr2lateNMrevised <- left_join(surviveyr2late, comm_diversity_total, by="Site") %>%
  mutate(Max_AlphaPhy_diff = Max_AlphaPhy - Max_AlphaPhy_mean,
         Mean_AlphaPhy_diff = Mean_AlphaPhy - Mean_AlphaPhy_mean,
         GammaPhy_diff = GammaPhy - GammaPhy_mean,
         Med_AlphaPhy_diff = Med_AlphaPhy - Med_AlphaPhy_mean) %>%
  mutate(suitabilitypc1V2abs = abs(suitabilitypc1V2)) %>%
  mutate(suitabilitypc2V2abs = abs(suitabilitypc2V2))

#DO SURVIVAL Y/N FOR ALL DATASETS

#GERMINATION YES NO

#Mean alpha diversity 

germ_YN_Mean_AlphaPhy <- glmmTMB(germ_YN ~ Mean_AlphaPhy_diff + Disturbance + propagule_pressure +
                                   suitabilitypc1V2abs + suitabilitypc2V2abs +
                                   (1|Family) + (1|POWO.name) + 
                                   (1|Site) + offset(log(Density)), 
                                 family=binomial(link="logit"), 
                                 data=germinationNMrevised)
summary(germ_YN_Mean_AlphaPhy)
Anova(germ_YN_Mean_AlphaPhy, type="III")
confint(germ_YN_Mean_AlphaPhy)
r.squaredGLMM(germ_YN_Mean_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = germ_YN_Mean_AlphaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(germ_YN_Mean_AlphaPhy)
check_singularity(germ_YN_Mean_AlphaPhy)

#Median alpha diversity 

germ_YN_Med_AlphaPhy <- glmmTMB(germ_YN ~ Med_AlphaPhy_diff + Disturbance + propagule_pressure +
                                  suitabilitypc1V2abs + suitabilitypc2V2abs +
                                  (1|Family) + (1|POWO.name) + 
                                  (1|Site) + offset(log(Density)), 
                                family=binomial(link="logit"), 
                                data=germinationNMrevised)
summary(germ_YN_Med_AlphaPhy)
Anova(germ_YN_Med_AlphaPhy, type="III")
confint(germ_YN_Med_AlphaPhy)
r.squaredGLMM(germ_YN_Med_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = germ_YN_Med_AlphaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(germ_YN_Med_AlphaPhy)
check_singularity(germ_YN_Med_AlphaPhy)

#Maximum alpha diversity

germ_YN_Max_AlphaPhy <- glmmTMB(germ_YN ~ Max_AlphaPhy_diff + Disturbance + propagule_pressure +
                                  suitabilitypc1V2abs + suitabilitypc2V2abs +
                                  (1|Family) + (1|POWO.name) +
                                  (1|Site) + offset(log(Density)), 
                                family=binomial(link="logit"), 
                                data=germinationNMrevised)
summary(germ_YN_Max_AlphaPhy)
Anova(germ_YN_Max_AlphaPhy, type="III")
confint(germ_YN_Max_AlphaPhy)
r.squaredGLMM(germ_YN_Max_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = germ_YN_Max_AlphaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(germ_YN_Max_AlphaPhy)
check_singularity(germ_YN_Max_AlphaPhy)

#Gamma diversity

germ_YN_GammaPhy <- glmmTMB(germ_YN ~ scale(GammaPhy_diff) + Disturbance + propagule_pressure +
                              suitabilitypc1V2abs + suitabilitypc2V2abs +
                              (1|Family) + (1|POWO.name) +
                              (1|Site) + offset(log(Density)), 
                            family=binomial(link="logit"), 
                            data=germinationNMrevised)
summary(germ_YN_GammaPhy)
Anova(germ_YN_GammaPhy, type="III")
confint(germ_YN_GammaPhy)
r.squaredGLMM(germ_YN_GammaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = germ_YN_GammaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(germ_YN_GammaPhy)
check_singularity(germ_YN_GammaPhy)

#FIRST YEAR SURVIVAL YES NO

#Mean alpha diversity 

surviveyr1_YN_Mean_AlphaPhy <- glmmTMB(surviveyr1_YN ~ Mean_AlphaPhy_diff + suitabilitypc2V2abs+Disturbance + propagule_pressure +
                                         + suitabilitypc1V2abs +
                                         (1|Family) + (1|POWO.name) + 
                                         (1|Site) + offset(log(Density)), 
                                       family=binomial(link="logit"), 
                                       data=surviveyr1NMrevised)
summary(surviveyr1_YN_Mean_AlphaPhy)
Anova(surviveyr1_YN_Mean_AlphaPhy, type="III")
confint(surviveyr1_YN_Mean_AlphaPhy)
r.squaredGLMM(surviveyr1_YN_Mean_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr1_YN_Mean_AlphaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr1_YN_Mean_AlphaPhy)
check_singularity(surviveyr1_YN_Mean_AlphaPhy)

#Median alpha diversity 

surviveyr1_YN_Med_AlphaPhy <- glmmTMB(surviveyr1_YN ~ Med_AlphaPhy_diff + suitabilitypc2V2abs+Disturbance + propagule_pressure +
                                        + suitabilitypc1V2abs +
                                        (1|Family) + (1|POWO.name) + 
                                        (1|Site) + offset(log(Density)), 
                                      family=binomial(link="logit"), 
                                      data=surviveyr1NMrevised)
summary(surviveyr1_YN_Med_AlphaPhy)
Anova(surviveyr1_YN_Med_AlphaPhy, type="III")
confint(surviveyr1_YN_Med_AlphaPhy)
r.squaredGLMM(surviveyr1_YN_Med_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr1_YN_Med_AlphaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr1_YN_Med_AlphaPhy)
check_singularity(surviveyr1_YN_Med_AlphaPhy)

#Maximum alpha diversity

surviveyr1_YN_Max_AlphaPhy <- glmmTMB(surviveyr1_YN ~ Max_AlphaPhy_diff + Disturbance + propagule_pressure +
                                        suitabilitypc1V2abs + suitabilitypc2V2abs + 
                                        (1|Family) + (1|POWO.name) + 
                                        (1|Site) + offset(log(Density)), 
                                      family=binomial(link="logit"), 
                                      data=surviveyr1NMrevised)
summary(surviveyr1_YN_Max_AlphaPhy)
Anova(surviveyr1_YN_Max_AlphaPhy, type="III")
confint(surviveyr1_YN_Max_AlphaPhy)
r.squaredGLMM(surviveyr1_YN_Max_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr1_YN_Max_AlphaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr1_YN_Max_AlphaPhy)
check_singularity(surviveyr1_YN_Max_AlphaPhy)

#Gamma diversity 

surviveyr1_YN_GammaPhy <- glmmTMB(surviveyr1_YN ~ scale(GammaPhy_diff)*suitabilitypc1V2abs + Disturbance + propagule_pressure +
                                    suitabilitypc2V2abs +
                                    (1|POWO.name) + 
                                    (1|Site) + offset(log(Density)), 
                                  family=binomial(link="logit"), 
                                  data=surviveyr1NMrevised)
summary(surviveyr1_YN_GammaPhy)
Anova(surviveyr1_YN_GammaPhy, type="III")
confint(surviveyr1_YN_GammaPhy)
r.squaredGLMM(surviveyr1_YN_GammaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr1_YN_GammaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr1_YN_GammaPhy)
check_singularity(surviveyr1_YN_GammaPhy)

#SECOND YEAR EARLY SURVIVAL YES NO

#Mean alpha diversity

surviveyr2early_YN_Mean_AlphaPhy <- glmmTMB(surviveyr2early_YN ~ Mean_AlphaPhy_diff*suitabilitypc2V2abs + Disturbance + propagule_pressure +
                                              suitabilitypc1V2abs +
                                              (1|Family) + (1|POWO.name) + 
                                              (1|Site) + offset(log(Density)), 
                                            family=binomial(link="logit"), 
                                            data=surviveyr2earlyNMrevised)
summary(surviveyr2early_YN_Mean_AlphaPhy)
Anova(surviveyr2early_YN_Mean_AlphaPhy, type="III")
confint(surviveyr2early_YN_Mean_AlphaPhy)
r.squaredGLMM(surviveyr2early_YN_Mean_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2early_YN_Mean_AlphaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr2early_YN_Mean_AlphaPhy)
check_singularity(surviveyr2early_YN_Mean_AlphaPhy)

#Median alpha diversity 

surviveyr2early_YN_Med_AlphaPhy <- glmmTMB(surviveyr2early_YN ~ Med_AlphaPhy_diff*suitabilitypc1V2abs + Disturbance + propagule_pressure +
                                             suitabilitypc2V2abs +
                                             (1|Family) + (1|POWO.name) + 
                                             (1|Site) + offset(log(Density)), 
                                           family=binomial(link="logit"), 
                                           data=surviveyr2earlyNMrevised)
summary(surviveyr2early_YN_Med_AlphaPhy)
Anova(surviveyr2early_YN_Med_AlphaPhy, type="III")
confint(surviveyr2early_YN_Med_AlphaPhy)
r.squaredGLMM(surviveyr2early_YN_Med_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2early_YN_Med_AlphaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr2early_YN_Med_AlphaPhy)
check_singularity(surviveyr2early_YN_Med_AlphaPhy)

# Maximum alpha diversity

surviveyr2early_YN_Max_AlphaPhy <- glmmTMB(surviveyr2early_YN ~ Max_AlphaPhy_diff*suitabilitypc2V2abs + Disturbance + propagule_pressure +
                                             suitabilitypc1V2abs + 
                                             (1|Family) + (1|POWO.name) + 
                                             (1|Site) + offset(log(Density)), 
                                           family=binomial(link="logit"), 
                                           data=surviveyr2earlyNMrevised)
summary(surviveyr2early_YN_Max_AlphaPhy)
Anova(surviveyr2early_YN_Max_AlphaPhy, type="III")
confint(surviveyr2early_YN_Max_AlphaPhy)
r.squaredGLMM(surviveyr2early_YN_Max_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2early_YN_Max_AlphaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr2early_YN_Max_AlphaPhy)
check_singularity(surviveyr2early_YN_Max_AlphaPhy)

#Gamma diversity 

surviveyr2early_YN_GammaPhy <- glmmTMB(surviveyr2early_YN ~ scale(GammaPhy_diff) + suitabilitypc2V2abs + Disturbance + propagule_pressure +
                                         suitabilitypc1V2abs +
                                         (1|POWO.name) + 
                                         offset(log(Density)), 
                                       family=binomial(link="logit"), 
                                       data=surviveyr2earlyNMrevised)
summary(surviveyr2early_YN_GammaPhy)
Anova(surviveyr2early_YN_GammaPhy, type="III")
confint(surviveyr2early_YN_GammaPhy)
r.squaredGLMM(surviveyr2early_YN_GammaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2early_YN_GammaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr2early_YN_GammaPhy)
check_singularity(surviveyr2early_YN_GammaPhy)

#SECOND YEAR LATE SURVIVAL YES NO

#Mean alpha diversity 

surviveyr2late_YN_Mean_AlphaPhy <- glmmTMB(surviveyr2late_YN ~ Mean_AlphaPhy_diff*suitabilitypc1V2abs + Disturbance + propagule_pressure +
                                             suitabilitypc2V2abs + 
                                             (1|Family) + (1|POWO.name) + 
                                             (1|Site) + offset(log(Density)), 
                                           family=binomial(link="logit"), 
                                           data=surviveyr2lateNMrevised)
summary(surviveyr2late_YN_Mean_AlphaPhy)
Anova(surviveyr2late_YN_Mean_AlphaPhy, type="III")
confint(surviveyr2late_YN_Mean_AlphaPhy)
r.squaredGLMM(surviveyr2late_YN_Mean_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2late_YN_Mean_AlphaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr2late_YN_Mean_AlphaPhy)
check_singularity(surviveyr2late_YN_Mean_AlphaPhy)

#Median alpha diversity (#checked and no interaction between diversity and disturbance or climate suitability)

surviveyr2late_YN_Med_AlphaPhy <- glmmTMB(surviveyr2late_YN ~ Med_AlphaPhy_diff*suitabilitypc1V2abs + Disturbance + propagule_pressure +
                                            suitabilitypc2V2abs + 
                                            (1|Family) + (1|POWO.name) + 
                                            (1|Site) + offset(log(Density)), 
                                          family=binomial(link="logit"), 
                                          data=surviveyr2lateNMrevised)
summary(surviveyr2late_YN_Med_AlphaPhy)
Anova(surviveyr2late_YN_Med_AlphaPhy, type="III")
confint(surviveyr2late_YN_Med_AlphaPhy)
r.squaredGLMM(surviveyr2late_YN_Med_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2late_YN_Med_AlphaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr2late_YN_Med_AlphaPhy)
check_singularity(surviveyr2late_YN_Med_AlphaPhy)

#Maximum alpha diversity

surviveyr2late_YN_Max_AlphaPhy <- glmmTMB(surviveyr2late_YN ~ Max_AlphaPhy_diff + Disturbance + propagule_pressure +
                                            suitabilitypc1V2abs + suitabilitypc2V2abs + 
                                            (1|Family) + (1|POWO.name) + 
                                            (1|Site) + offset(log(Density)), 
                                          family=binomial(link="logit"), 
                                          data=surviveyr2lateNMrevised)
summary(surviveyr2late_YN_Max_AlphaPhy)
Anova(surviveyr2late_YN_Max_AlphaPhy, type="III")
confint(surviveyr2late_YN_Max_AlphaPhy)
r.squaredGLMM(surviveyr2late_YN_Max_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2late_YN_Max_AlphaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr2late_YN_Max_AlphaPhy)
check_singularity(surviveyr2late_YN_Max_AlphaPhy)

#Gamma diversity 

surviveyr2late_YN_GammaPhy <- glmmTMB(surviveyr2late_YN ~ scale(GammaPhy_diff) + Disturbance + propagule_pressure +
                                        suitabilitypc1V2abs + suitabilitypc2V2abs + 
                                        (1|POWO.name) + 
                                        (1|Site) + offset(log(Density)), 
                                      family=binomial(link="logit"), 
                                      data=surviveyr2lateNMrevised)
summary(surviveyr2late_YN_GammaPhy)
Anova(surviveyr2late_YN_GammaPhy, type="III")
confint(surviveyr2late_YN_GammaPhy)
r.squaredGLMM(surviveyr2late_YN_GammaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2late_YN_GammaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr2late_YN_GammaPhy)
check_singularity(surviveyr2late_YN_GammaPhy)

#######################################################################################

#DO NUMBER SURVIVING FOR ALL DATASETS

#GERMINATION COUNT

germinationcountNMrevised <- filter(germinationNMrevised, germ_num>0)

#Mean alpha diversity 

germ_count_Mean_AlphaPhy <- glmmTMB(germ_num ~ Mean_AlphaPhy_diff*Disturbance + propagule_pressure +
                                      suitabilitypc1V2abs + suitabilitypc2V2abs + Density +
                                      (1|Family) + (1|POWO.name) + (1|Site),
                                    family=nbinom2(link = "log"), 
                                    data=germinationcountNMrevised)
summary(germ_count_Mean_AlphaPhy)
Anova(germ_count_Mean_AlphaPhy, type="III")
confint(germ_count_Mean_AlphaPhy)
r.squaredGLMM(germ_count_Mean_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = germ_count_Mean_AlphaPhy, plot=F)
plot(simulationOutput)

check_overdispersion(germ_count_Mean_AlphaPhy)
check_singularity(germ_count_Mean_AlphaPhy)
check_model(germ_count_Mean_AlphaPhy)

#Median alpha diversity 

germ_count_Med_AlphaPhy <- glmmTMB(germ_num ~ Med_AlphaPhy_diff*suitabilitypc1V2abs + Med_AlphaPhy_diff*Disturbance + propagule_pressure +
                                     suitabilitypc2V2abs + Density +
                                     (1|Family) + (1|POWO.name) + (1|Site),
                                   family=nbinom2(link = "log"), 
                                   data=germinationcountNMrevised)
summary(germ_count_Med_AlphaPhy)
Anova(germ_count_Med_AlphaPhy, type="III")
confint(germ_count_Med_AlphaPhy)
r.squaredGLMM(germ_count_Med_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = germ_count_Med_AlphaPhy, plot=F)
plot(simulationOutput)

check_overdispersion(germ_count_Med_AlphaPhy)
check_singularity(germ_count_Med_AlphaPhy)
check_model(germ_count_Med_AlphaPhy)

#Maximum alpha diversity

germ_count_Max_AlphaPhy <- glmmTMB(germ_num ~ Max_AlphaPhy_diff*suitabilitypc2V2abs + Max_AlphaPhy_diff*Disturbance + propagule_pressure +
                                     suitabilitypc1V2abs + Density +
                                     (1|Family) + (1|POWO.name) + (1|Site),
                                   family=nbinom2(link = "log"), 
                                   data=germinationcountNMrevised)
summary(germ_count_Max_AlphaPhy)
Anova(germ_count_Max_AlphaPhy, type="III")
confint(germ_count_Max_AlphaPhy)
r.squaredGLMM(germ_count_Max_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = germ_count_Max_AlphaPhy, plot=F)
plot(simulationOutput)

check_overdispersion(germ_count_Max_AlphaPhy)
check_singularity(germ_count_Max_AlphaPhy)

#Gamma diversity 

germ_count_GammaPhy <- glmmTMB(germ_num ~ scale(GammaPhy_diff) + Disturbance + propagule_pressure +
                                 suitabilitypc1V2abs + suitabilitypc2V2abs + Density +
                                 (1|Family) + (1|POWO.name) + (1|Site),
                               family=nbinom2(link = "log"), 
                               data=germinationcountNMrevised)
summary(germ_count_GammaPhy)
Anova(germ_count_GammaPhy, type="III")
confint(germ_count_GammaPhy)
r.squaredGLMM(germ_count_GammaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = germ_count_GammaPhy, plot=F)
plot(simulationOutput)

check_overdispersion(germ_count_GammaPhy)
check_singularity(germ_count_GammaPhy)

#FIRST YEAR SURVIVAL COUNT

surviveyr1countNMrevised <- filter(surviveyr1NMrevised, surviveyr1_num>0)

#Mean alpha diversity 

surviveyr1_count_Mean_AlphaPhy <- glmmTMB(surviveyr1_num ~ Mean_AlphaPhy_diff*Disturbance + propagule_pressure +
                                            suitabilitypc1V2abs + suitabilitypc2V2abs + Density +
                                            (1|POWO.name) + (1|Site),
                                          family=nbinom2(link = "log"), 
                                          data=surviveyr1countNMrevised)
summary(surviveyr1_count_Mean_AlphaPhy)
Anova(surviveyr1_count_Mean_AlphaPhy, type="III")
confint(surviveyr1_count_Mean_AlphaPhy)
r.squaredGLMM(surviveyr1_count_Mean_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr1_count_Mean_AlphaPhy, plot=F)
plot(simulationOutput)

check_overdispersion(surviveyr1_count_Mean_AlphaPhy)
check_singularity(surviveyr1_count_Mean_AlphaPhy)

#Median alpha diversity 

surviveyr1_count_Med_AlphaPhy <- glmmTMB(surviveyr1_num ~ Med_AlphaPhy_diff*Disturbance + propagule_pressure +
                                           suitabilitypc1V2abs + suitabilitypc2V2abs + Density +
                                           (1|POWO.name) + (1|Site),
                                         family=nbinom2(link = "log"), 
                                         data=surviveyr1countNMrevised)
summary(surviveyr1_count_Med_AlphaPhy)
Anova(surviveyr1_count_Med_AlphaPhy, type="III")
confint(surviveyr1_count_Med_AlphaPhy)
r.squaredGLMM(surviveyr1_count_Med_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr1_count_Med_AlphaPhy, plot=F)
plot(simulationOutput)

check_overdispersion(surviveyr1_count_Med_AlphaPhy)
check_singularity(surviveyr1_count_Med_AlphaPhy)

#Maximum alpha diversity

surviveyr1_count_Max_AlphaPhy <- glmmTMB(surviveyr1_num ~ Max_AlphaPhy_diff*Disturbance + propagule_pressure +
                                           suitabilitypc1V2abs + suitabilitypc2V2abs + Density +
                                           (1|POWO.name) + (1|Site),
                                         family=nbinom2(link = "log"), 
                                         data=surviveyr1countNMrevised)
summary(surviveyr1_count_Max_AlphaPhy)
Anova(surviveyr1_count_Max_AlphaPhy, type="III")
confint(surviveyr1_count_Max_AlphaPhy)
r.squaredGLMM(surviveyr1_count_Max_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr1_count_Max_AlphaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr1_count_Max_AlphaPhy)
check_singularity(surviveyr1_count_Max_AlphaPhy)

#Gamma diversity

surviveyr1_count_GammaPhy <- glmmTMB(surviveyr1_num ~ scale(GammaPhy_diff) + Disturbance + propagule_pressure +
                                       suitabilitypc1V2abs + suitabilitypc2V2abs + Density +
                                       (1|POWO.name) + (1|Site),
                                     family=nbinom2(link = "log"), 
                                     data=surviveyr1countNMrevised)
summary(surviveyr1_count_GammaPhy)
Anova(surviveyr1_count_GammaPhy, type="III")
confint(surviveyr1_count_GammaPhy)
r.squaredGLMM(surviveyr1_count_GammaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr1_count_GammaPhy, plot=F)
plot(simulationOutput)

check_overdispersion(surviveyr1_count_GammaPhy)
check_singularity(surviveyr1_count_GammaPhy)

#EARLY SECOND YEAR SURVIVAL COUNT

surviveyr2earlycountNMrevised <- filter(surviveyr2earlyNMrevised, surviveyr2early_num>0)

#Mean alpha diversity 

surviveyr2early_count_Mean_AlphaPhy <- glmmTMB(surviveyr2early_num ~ Mean_AlphaPhy_diff + Disturbance + propagule_pressure +
                                                 suitabilitypc1V2abs + suitabilitypc2V2abs + Density +
                                                 (1|POWO.name) + (1|Site),
                                               family=nbinom2(link = "log"), 
                                               data=surviveyr2earlycountNMrevised)
summary(surviveyr2early_count_Mean_AlphaPhy)
Anova(surviveyr2early_count_Mean_AlphaPhy, type="III")
confint(surviveyr2early_count_Mean_AlphaPhy)
r.squaredGLMM(surviveyr2early_count_Mean_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2early_count_Mean_AlphaPhy, plot=F)
plot(simulationOutput)

check_overdispersion(surviveyr2early_count_Mean_AlphaPhy)
check_singularity(surviveyr2early_count_Mean_AlphaPhy)

#Median alpha diversity 

surviveyr2early_count_Med_AlphaPhy <- glmmTMB(surviveyr2early_num ~ Med_AlphaPhy_diff + suitabilitypc1V2abs + Disturbance + propagule_pressure +
                                                suitabilitypc2V2abs + Density +
                                                (1|POWO.name) + (1|Site),
                                              family=nbinom2(link = "log"), 
                                              data=surviveyr2earlycountNMrevised)
summary(surviveyr2early_count_Med_AlphaPhy)
Anova(surviveyr2early_count_Med_AlphaPhy, type="III")
confint(surviveyr2early_count_Med_AlphaPhy)
r.squaredGLMM(surviveyr2early_count_Med_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2early_count_Med_AlphaPhy, plot=F)
plot(simulationOutput)

check_overdispersion(surviveyr2early_count_Med_AlphaPhy)
check_singularity(surviveyr2early_count_Med_AlphaPhy)

#Max alpha diversity 

surviveyr2early_count_Max_AlphaPhy <- glmmTMB(surviveyr2early_num ~ Max_AlphaPhy_diff + Disturbance + propagule_pressure +
                                                suitabilitypc1V2abs + suitabilitypc2V2abs + Density +
                                                (1|POWO.name) + (1|Site),
                                              family=nbinom2(link = "log"), 
                                              data=surviveyr2earlycountNMrevised)
summary(surviveyr2early_count_Max_AlphaPhy)
Anova(surviveyr2early_count_Max_AlphaPhy, type="III")
confint(surviveyr2early_count_Max_AlphaPhy)
r.squaredGLMM(surviveyr2early_count_Max_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2early_count_Max_AlphaPhy, plot=F)
plot(simulationOutput)

check_overdispersion(surviveyr2early_count_Max_AlphaPhy)
check_singularity(surviveyr2early_count_Max_AlphaPhy)

#Gamma diversity 

surviveyr2early_count_GammaPhy <- glmmTMB(surviveyr2early_num ~ scale(GammaPhy_diff) + Disturbance + propagule_pressure +
                                            suitabilitypc1V2abs + suitabilitypc2V2abs + Density + 
                                            (1|POWO.name),
                                          family=nbinom2(link = "log"), 
                                          data=surviveyr2earlycountNMrevised)
summary(surviveyr2early_count_GammaPhy)
Anova(surviveyr2early_count_GammaPhy, type="III")
confint(surviveyr2early_count_GammaPhy)
r.squaredGLMM(surviveyr2early_count_GammaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2early_count_GammaPhy, plot=F)
plot(simulationOutput)

check_overdispersion(surviveyr2early_count_GammaPhy)
check_singularity(surviveyr2early_count_GammaPhy)

#LATE SECOND YEAR SURVIVAL COUNT

surviveyr2latecountNMrevised <- filter(surviveyr2lateNMrevised, surviveyr2late_num>0)

#Mean alpha diversity 

surviveyr2late_count_Mean_AlphaPhy <- glmmTMB(surviveyr2late_num ~ Mean_AlphaPhy_diff + suitabilitypc2V2abs + Disturbance +
                                                suitabilitypc1V2abs + suitabilitypc2V2abs + Density + propagule_pressure +
                                                (1|POWO.name),
                                              family=nbinom2(link = "log"), 
                                              data=surviveyr2latecountNMrevised)
summary(surviveyr2late_count_Mean_AlphaPhy)
Anova(surviveyr2late_count_Mean_AlphaPhy, type="III")
confint(surviveyr2late_count_Mean_AlphaPhy)
r.squaredGLMM(surviveyr2late_count_Mean_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2late_count_Mean_AlphaPhy, plot=F)
plot(simulationOutput)

check_overdispersion(surviveyr2late_count_Mean_AlphaPhy)
check_singularity(surviveyr2late_count_Mean_AlphaPhy)

#Median alpha diversity 

surviveyr2late_count_Med_AlphaPhy <- glmmTMB(surviveyr2late_num ~ Med_AlphaPhy_diff*suitabilitypc2V2abs + Disturbance +
                                               suitabilitypc1V2abs + suitabilitypc2V2abs + Density + propagule_pressure +
                                               (1|POWO.name),
                                             family=nbinom2(link = "log"), 
                                             data=surviveyr2latecountNMrevised)
summary(surviveyr2late_count_Med_AlphaPhy)
Anova(surviveyr2late_count_Med_AlphaPhy, type="III")
confint(surviveyr2late_count_Med_AlphaPhy)
r.squaredGLMM(surviveyr2late_count_Med_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2late_count_Med_AlphaPhy, plot=F)
plot(simulationOutput)

check_overdispersion(surviveyr2late_count_Med_AlphaPhy)
check_singularity(surviveyr2late_count_Med_AlphaPhy)

#Max alpha diversity 

surviveyr2late_count_Max_AlphaPhy <- glmmTMB(surviveyr2late_num ~ Max_AlphaPhy_diff + Disturbance + 
                                               suitabilitypc1V2abs + suitabilitypc2V2abs + Density + propagule_pressure +
                                               (1|POWO.name),
                                             family=nbinom2(link = "log"), 
                                             data=surviveyr2latecountNMrevised)
summary(surviveyr2late_count_Max_AlphaPhy)
Anova(surviveyr2late_count_Max_AlphaPhy, type="III")
confint(surviveyr2late_count_Max_AlphaPhy)
r.squaredGLMM(surviveyr2late_count_Max_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2late_count_Max_AlphaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr2late_count_Max_AlphaPhy)
check_singularity(surviveyr2late_count_Max_AlphaPhy)

#Gamma diversity 

surviveyr2late_count_GammaPhy <- glmmTMB(surviveyr2late_num ~ scale(GammaPhy_diff) + Disturbance + 
                                           suitabilitypc1V2abs + suitabilitypc2V2abs +  Density + propagule_pressure +
                                           (1|POWO.name),
                                         family=nbinom2(link = "log"), 
                                         data=surviveyr2latecountNMrevised)
summary(surviveyr2late_count_GammaPhy)
Anova(surviveyr2late_count_GammaPhy, type="III")
confint(surviveyr2late_count_GammaPhy)
r.squaredGLMM(surviveyr2late_count_GammaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2late_count_GammaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr2late_count_GammaPhy)
check_singularity(surviveyr2late_count_GammaPhy)

#Plot the climate interactions for the K+H difference models

yr2earlyYNmeanprecipinteraction <- interact_plot(surviveyr2early_YN_Mean_AlphaPhy, pred="Mean_AlphaPhy_diff", modx = "suitabilitypc2V2abs",
                                                 plot.points=TRUE, partial.residuals = TRUE, legend.main = "Climate dissimilarity (precip.)") +
  scale_y_continuous(trans='log10') +
  xlab("Difference in Mean Alpha Phylogenetic Diversity") + ylab("Likelihood of any plants in plot surviving overwinter") +
  scale_linetype_manual(values=c("solid", "solid", "solid")) +
  theme_bw() + guides(linetype = "none") + theme(legend.position = "top")  +
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=10), 
        legend.text=element_text(size=8), 
        legend.title=element_text(size=10))
yr2earlyYNmeanprecipinteraction

yr2earlyYNmaxprecipinteraction <- interact_plot(surviveyr2early_YN_Max_AlphaPhy, pred="Max_AlphaPhy_diff", modx = "suitabilitypc2V2abs",
                                                plot.points=TRUE, partial.residuals = TRUE, legend.main = "Climate dissimilarity (precip.)") +
  scale_y_continuous(trans='log10') +
  xlab("Difference in Maximum Alpha Phylogenetic Diversity") + ylab("Likelihood of any plants in plot surviving overwinter") +
  scale_linetype_manual(values=c("solid", "solid", "solid")) +
  theme_bw() + guides(linetype = "none") + theme(legend.position = "top")  +
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=10), 
        legend.text=element_text(size=8), 
        legend.title=element_text(size=10))
yr2earlyYNmaxprecipinteraction

yr2earlyYNmedtempinteraction <- interact_plot(surviveyr2early_YN_Med_AlphaPhy, pred="Med_AlphaPhy_diff", modx = "suitabilitypc1V2abs",
                                              plot.points=TRUE, partial.residuals = TRUE, legend.main = "Climate dissimilarity (temp.)") +
  scale_y_continuous(trans='log10') +
  xlab("Difference in Median Alpha Phylogenetic Diversity") + ylab("Likelihood of any plants in plot surviving overwinter") +
  scale_linetype_manual(values=c("solid", "solid", "solid")) +
  theme_bw() + guides(linetype = "none") + theme(legend.position = "top")  +
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=10), 
        legend.text=element_text(size=8), 
        legend.title=element_text(size=10))
yr2earlyYNmedtempinteraction

KHfig <- ggarrange(yr2earlyYNmeanprecipinteraction, yr2earlyYNmaxprecipinteraction, 
                   yr2earlyYNmedtempinteraction, nrow=2, ncol=2)
KHfig
ggsave(KHfig, 
       filename = "Fig_S16.svg",
       height = 210, width = 210, units = "mm")

####################################################################################################

#DO IT AGAIN BUT WITH SINGLE CLIMATE VARIABLES 
#The 'SC' models in Tables S2-S9

####################################################################################################

#Doing it for Question 1

#Read in all four data files

germinationv2 <- readRDS("germinationv2revised.RDS") %>%
  mutate(suitabilitypc1V3abs = abs(suitabilitypc1V3)) %>%
  mutate(suitabilitypc2V3abs = abs(suitabilitypc2V3))
str(germinationv2)

germinationv2$Disturbance <- relevel(germinationv2$Disturbance, ref="0")

surviveyr1v2 <- readRDS("surviveyr1v2revised.RDS") %>%
  mutate(suitabilitypc1V3abs = abs(suitabilitypc1V3)) %>%
  mutate(suitabilitypc2V3abs = abs(suitabilitypc2V3))
str(surviveyr1v2)

surviveyr1v2$Disturbance <- relevel(surviveyr1v2$Disturbance, ref="0")

surviveyr2earlyv2 <- readRDS("surviveyr2earlyv2revised.RDS") %>%
  mutate(suitabilitypc1V3abs = abs(suitabilitypc1V3)) %>%
  mutate(suitabilitypc2V3abs = abs(suitabilitypc2V3))
str(surviveyr2earlyv2)

surviveyr2earlyv2$Disturbance <- relevel(surviveyr2earlyv2$Disturbance, ref="0")

surviveyr2latev2 <- readRDS("surviveyr2latev2revised.RDS") %>%
  mutate(suitabilitypc1V3abs = abs(suitabilitypc1V3)) %>%
  mutate(suitabilitypc2V3abs = abs(suitabilitypc2V3))
str(surviveyr2latev2)

surviveyr2latev2$Disturbance <- relevel(surviveyr2latev2$Disturbance, ref="0")

hist(germinationv2$suitabilitypc1V3abs)
hist(germinationv2$suitabilitypc2V3abs) #need to log these to help meet assumptions

################################################################################

#DO SURVIVAL Y/N FOR ALL DATASETS

#GERMINATION YES NO

#Mean alpha diversity

germ_YN_Mean_AlphaPhy <- glmmTMB(germ_YN ~ Mean_AlphaPhy + Disturbance + propagule_pressure +
                                   log(suitabilitypc1V3abs) + log(suitabilitypc2V3abs) +
                                   (1|Family) + (1|POWO.name) + 
                                   (1|Site) + offset(log(Density)), 
                                 family=binomial(link="logit"), 
                                 data=germinationv2)
summary(germ_YN_Mean_AlphaPhy)
Anova(germ_YN_Mean_AlphaPhy, type="III")
confint(germ_YN_Mean_AlphaPhy)
r.squaredGLMM(germ_YN_Mean_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = germ_YN_Mean_AlphaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(germ_YN_Mean_AlphaPhy)
check_singularity(germ_YN_Mean_AlphaPhy)

#Maximum alpha diversity

germ_YN_Max_AlphaPhy <- glmmTMB(germ_YN ~ Max_AlphaPhy + Disturbance + propagule_pressure +
                                  log(suitabilitypc1V3abs) + log(suitabilitypc2V3abs) +
                                  (1|Family) + (1|POWO.name) +
                                  (1|Site) + offset(log(Density)), 
                                family=binomial(link="logit"), 
                                data=germinationv2)
summary(germ_YN_Max_AlphaPhy)
Anova(germ_YN_Max_AlphaPhy, type="III")
confint(germ_YN_Max_AlphaPhy)
r.squaredGLMM(germ_YN_Max_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = germ_YN_Max_AlphaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(germ_YN_Max_AlphaPhy)
check_singularity(germ_YN_Max_AlphaPhy)

#Gamma diversity

germ_YN_GammaPhy <- glmmTMB(germ_YN ~ log(GammaPhy) + Disturbance + propagule_pressure +
                              log(suitabilitypc1V3abs) + log(suitabilitypc2V3abs) +
                              (1|Family) + (1|POWO.name) +
                              (1|Site) + offset(log(Density)), 
                            family=binomial(link="logit"), 
                            data=germinationv2)
summary(germ_YN_GammaPhy)
Anova(germ_YN_GammaPhy, type="III")
confint(germ_YN_GammaPhy)
r.squaredGLMM(germ_YN_GammaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = germ_YN_GammaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(germ_YN_GammaPhy)
check_singularity(germ_YN_GammaPhy)

#Median alpha diversity 

germ_YN_MedPhy <- glmmTMB(germ_YN ~ Med_AlphaPhy*Disturbance + propagule_pressure +
                            log(suitabilitypc1V3abs) + log(suitabilitypc2V3abs) +
                            (1|Family) + (1|POWO.name) +
                            (1|Site) + offset(log(Density)), 
                          family=binomial(link="logit"), 
                          data=germinationv2)
summary(germ_YN_MedPhy)
Anova(germ_YN_MedPhy, type="III")
confint(germ_YN_MedPhy)
r.squaredGLMM(germ_YN_MedPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = germ_YN_MedPhy, plot=F)
plot(simulationOutput)
check_overdispersion(germ_YN_MedPhy)
check_singularity(germ_YN_MedPhy)

#FIRST YEAR SURVIVAL YES NO

#Mean alpha diversity 

surviveyr1_YN_Mean_AlphaPhy <- glmmTMB(surviveyr1_YN ~ Mean_AlphaPhy + Disturbance + propagule_pressure +
                                         log(suitabilitypc1V3abs) + log(suitabilitypc2V3abs) + Herbivory +
                                         (1|Family) + (1|POWO.name) + 
                                         (1|Site) + offset(log(Density)), 
                                       family=binomial(link="logit"), 
                                       data=surviveyr1v2)
summary(surviveyr1_YN_Mean_AlphaPhy)
Anova(surviveyr1_YN_Mean_AlphaPhy, type="III")
confint(surviveyr1_YN_Mean_AlphaPhy)
r.squaredGLMM(surviveyr1_YN_Mean_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr1_YN_Mean_AlphaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr1_YN_Mean_AlphaPhy)
check_singularity(surviveyr1_YN_Mean_AlphaPhy)

#Maximum alpha diversity

surviveyr1_YN_Max_AlphaPhy <- glmmTMB(surviveyr1_YN ~ Max_AlphaPhy*Disturbance + propagule_pressure +
                                        log(suitabilitypc1V3abs) + log(suitabilitypc2V3abs) + Herbivory +
                                        (1|Family) + (1|POWO.name) + 
                                        (1|Site) + offset(log(Density)), 
                                      family=binomial(link="logit"), 
                                      data=surviveyr1v2)
summary(surviveyr1_YN_Max_AlphaPhy)
Anova(surviveyr1_YN_Max_AlphaPhy, type="III")
confint(surviveyr1_YN_Max_AlphaPhy)
r.squaredGLMM(surviveyr1_YN_Max_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr1_YN_Max_AlphaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr1_YN_Max_AlphaPhy)
check_singularity(surviveyr1_YN_Max_AlphaPhy)

#Gamma diversity 

surviveyr1_YN_GammaPhy <- glmmTMB(surviveyr1_YN ~ log(GammaPhy) + log(suitabilitypc1V3abs) + Disturbance + propagule_pressure +
                                    log(suitabilitypc2V3abs) + Herbivory +
                                    (1|Family) + (1|POWO.name) + 
                                    (1|Site) + offset(log(Density)), 
                                  family=binomial(link="logit"), 
                                  data=surviveyr1v2)
summary(surviveyr1_YN_GammaPhy)
Anova(surviveyr1_YN_GammaPhy, type="III")
confint(surviveyr1_YN_GammaPhy)
r.squaredGLMM(surviveyr1_YN_GammaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr1_YN_GammaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr1_YN_GammaPhy)
check_singularity(surviveyr1_YN_GammaPhy)

#Median alpha diversity 

surviveyr1_YN_MedPhy <- glmmTMB(surviveyr1_YN ~ Med_AlphaPhy + Disturbance + propagule_pressure +
                                  log(suitabilitypc1V3abs) + log(suitabilitypc2V3abs) + Herbivory +
                                  (1|Family) + (1|POWO.name) + 
                                  (1|Site) + offset(log(Density)), 
                                family=binomial(link="logit"), 
                                data=surviveyr1v2)
summary(surviveyr1_YN_MedPhy)
Anova(surviveyr1_YN_MedPhy, type="III")
confint(surviveyr1_YN_MedPhy)
r.squaredGLMM(surviveyr1_YN_MedPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr1_YN_MedPhy, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr1_YN_MedPhy)
check_singularity(surviveyr1_YN_MedPhy)

#SECOND YEAR EARLY SURVIVAL YES NO

#Mean alpha diversity

surviveyr2early_YN_Mean_AlphaPhy <- glmmTMB(surviveyr2early_YN ~ Mean_AlphaPhy + Disturbance + propagule_pressure +
                                              log(suitabilitypc1V3abs) + log(suitabilitypc2V3abs) + Herbivory +
                                              (1|Family) + (1|POWO.name) + 
                                              (1|Site) + offset(log(Density)), 
                                            family=binomial(link="logit"), 
                                            data=surviveyr2earlyv2)
summary(surviveyr2early_YN_Mean_AlphaPhy)
Anova(surviveyr2early_YN_Mean_AlphaPhy, type="III")
confint(surviveyr2early_YN_Mean_AlphaPhy)
r.squaredGLMM(surviveyr2early_YN_Mean_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2early_YN_Mean_AlphaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr2early_YN_Mean_AlphaPhy)
check_singularity(surviveyr2early_YN_Mean_AlphaPhy)

# Maximum alpha diversity

surviveyr2early_YN_Max_AlphaPhy <- glmmTMB(surviveyr2early_YN ~ Max_AlphaPhy*log(suitabilitypc1V3abs) + Disturbance + propagule_pressure +
                                             log(suitabilitypc2V3abs) + Herbivory +
                                             (1|Family) + (1|POWO.name) + 
                                             (1|Site) + offset(log(Density)), 
                                           family=binomial(link="logit"), 
                                           data=surviveyr2earlyv2)
summary(surviveyr2early_YN_Max_AlphaPhy)
Anova(surviveyr2early_YN_Max_AlphaPhy, type="III")
confint(surviveyr2early_YN_Max_AlphaPhy)
r.squaredGLMM(surviveyr2early_YN_Max_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2early_YN_Max_AlphaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr2early_YN_Max_AlphaPhy)
check_singularity(surviveyr2early_YN_Max_AlphaPhy)

#Gamma diversity 

surviveyr2early_YN_GammaPhy <- glmmTMB(surviveyr2early_YN ~ log(GammaPhy) + Disturbance + propagule_pressure +
                                         log(suitabilitypc1V3abs) + log(suitabilitypc2V3abs) + Herbivory +
                                         (1|Family) + (1|POWO.name) + 
                                         (1|Site) + offset(log(Density)), 
                                       family=binomial(link="logit"), 
                                       data=surviveyr2earlyv2)
summary(surviveyr2early_YN_GammaPhy)
Anova(surviveyr2early_YN_GammaPhy, type="III")
confint(surviveyr2early_YN_GammaPhy)
r.squaredGLMM(surviveyr2early_YN_GammaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2early_YN_GammaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr2early_YN_GammaPhy)
check_singularity(surviveyr2early_YN_GammaPhy)

#Median alpha diversity 

surviveyr2early_YN_MedPhy <- glmmTMB(surviveyr2early_YN ~ Med_AlphaPhy + Disturbance + propagule_pressure +
                                       log(suitabilitypc1V3abs) + log(suitabilitypc2V3abs) + Herbivory +
                                       (1|Family) + (1|POWO.name) + 
                                       (1|Site) + offset(log(Density)), 
                                     family=binomial(link="logit"), 
                                     data=surviveyr2earlyv2)
summary(surviveyr2early_YN_MedPhy)
Anova(surviveyr2early_YN_MedPhy, type="III")
confint(surviveyr2early_YN_MedPhy)
r.squaredGLMM(surviveyr2early_YN_MedPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2early_YN_MedPhy, plot=F)
plot(simulationOutput)
testOutliers(simulationOutput)
check_overdispersion(surviveyr2early_YN_MedPhy)
check_singularity(surviveyr2early_YN_MedPhy)

check_model(surviveyr2early_YN_MedPhy)

#SECOND YEAR LATE SURVIVAL YES NO

#Mean alpha diversity 

surviveyr2late_YN_Mean_AlphaPhy <- glmmTMB(surviveyr2late_YN ~ Mean_AlphaPhy*Disturbance + log(suitabilitypc1V3abs) + propagule_pressure +
                                             log(suitabilitypc2V3abs) + 
                                             (1|Family) + (1|POWO.name) + 
                                             (1|Site) + offset(log(Density)), 
                                           family=binomial(link="logit"), 
                                           data=surviveyr2latev2)
summary(surviveyr2late_YN_Mean_AlphaPhy)
Anova(surviveyr2late_YN_Mean_AlphaPhy, type="III")
confint(surviveyr2late_YN_Mean_AlphaPhy)
r.squaredGLMM(surviveyr2late_YN_Mean_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2late_YN_Mean_AlphaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr2late_YN_Mean_AlphaPhy)
check_singularity(surviveyr2late_YN_Mean_AlphaPhy)

#Maximum alpha diversity

surviveyr2late_YN_Max_AlphaPhy <- glmmTMB(surviveyr2late_YN ~ Max_AlphaPhy*Disturbance + propagule_pressure +
                                            log(suitabilitypc1V3abs) + log(suitabilitypc2V3abs) + 
                                            (1|Family) + (1|POWO.name) + 
                                            (1|Site) + offset(log(Density)), 
                                          family=binomial(link="logit"), 
                                          data=surviveyr2latev2)
summary(surviveyr2late_YN_Max_AlphaPhy)
Anova(surviveyr2late_YN_Max_AlphaPhy, type="III")
confint(surviveyr2late_YN_Max_AlphaPhy)
r.squaredGLMM(surviveyr2late_YN_Max_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2late_YN_Max_AlphaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr2late_YN_Max_AlphaPhy)
check_singularity(surviveyr2late_YN_Max_AlphaPhy)

#Gamma diversity 

surviveyr2late_YN_GammaPhy <- glmmTMB(surviveyr2late_YN ~ log(suitabilitypc1V3abs) + log(GammaPhy)*Disturbance + propagule_pressure +
                                        log(suitabilitypc2V3abs) + 
                                        (1|Family) + (1|POWO.name) + 
                                        (1|Site) + offset(log(Density)), 
                                      family=binomial(link="logit"), 
                                      data=surviveyr2latev2)
summary(surviveyr2late_YN_GammaPhy)
Anova(surviveyr2late_YN_GammaPhy, type="III")
confint(surviveyr2late_YN_GammaPhy)
r.squaredGLMM(surviveyr2late_YN_GammaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2late_YN_GammaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr2late_YN_GammaPhy)
check_singularity(surviveyr2late_YN_GammaPhy)

#Median alpha diversity 

surviveyr2late_YN_MedPhy <- glmmTMB(surviveyr2late_YN ~ Med_AlphaPhy + Disturbance + propagule_pressure +
                                      log(suitabilitypc1V3abs) + log(suitabilitypc2V3abs) + 
                                      (1|Family) + (1|POWO.name) + 
                                      (1|Site) + offset(log(Density)), 
                                    family=binomial(link="logit"), 
                                    data=surviveyr2latev2)
summary(surviveyr2late_YN_MedPhy)
Anova(surviveyr2late_YN_MedPhy, type="III")
confint(surviveyr2late_YN_MedPhy)
r.squaredGLMM(surviveyr2late_YN_MedPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2late_YN_MedPhy, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr2late_YN_MedPhy)
check_singularity(surviveyr2late_YN_MedPhy)

#######################################################################################

#DO NUMBER SURVIVING FOR ALL DATASETS

#GERMINATION COUNT

germinationcountv2 <- filter(germinationv2, germ_num>0)

#Mean alpha diversity

germ_count_Mean_AlphaPhy <- glmmTMB(germ_num ~ Mean_AlphaPhy*Disturbance + propagule_pressure +
                                      log(suitabilitypc1V3abs) + log(suitabilitypc2V3abs) + Density +
                                      (1|Family) + (1|POWO.name) + (1|Site),
                                    family=nbinom2(link = "log"), 
                                    data=germinationcountv2)
summary(germ_count_Mean_AlphaPhy)
Anova(germ_count_Mean_AlphaPhy, type="III")
confint(germ_count_Mean_AlphaPhy)
r.squaredGLMM(germ_count_Mean_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = germ_count_Mean_AlphaPhy, plot=F)
plot(simulationOutput)

check_overdispersion(germ_count_Mean_AlphaPhy)
check_singularity(germ_count_Mean_AlphaPhy)
check_model(germ_count_Mean_AlphaPhy)

#Maximum alpha diversity

germ_count_Max_AlphaPhy <- glmmTMB(germ_num ~ Max_AlphaPhy*Disturbance + propagule_pressure +
                                     log(suitabilitypc1V3abs) + log(suitabilitypc2V3abs) + Density +
                                     (1|Family) + (1|POWO.name) + (1|Site),
                                   family=nbinom2(link = "log"), 
                                   data=germinationcountv2)
summary(germ_count_Max_AlphaPhy)
Anova(germ_count_Max_AlphaPhy, type="III")
confint(germ_count_Max_AlphaPhy)
r.squaredGLMM(germ_count_Max_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = germ_count_Max_AlphaPhy, plot=F)
plot(simulationOutput)

check_overdispersion(germ_count_Max_AlphaPhy)
check_singularity(germ_count_Max_AlphaPhy)

#Gamma diversity 

germ_count_GammaPhy <- glmmTMB(germ_num ~ log(GammaPhy) + Disturbance + propagule_pressure +
                                 log(suitabilitypc1V3abs) + log(suitabilitypc2V3abs) + Density +
                                 (1|Family) + (1|POWO.name) + (1|Site),
                               family=nbinom2(link = "log"), 
                               data=germinationcountv2)
summary(germ_count_GammaPhy)
Anova(germ_count_GammaPhy, type="III")
confint(germ_count_GammaPhy)
r.squaredGLMM(germ_count_GammaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = germ_count_GammaPhy, plot=F)
plot(simulationOutput)

check_overdispersion(germ_count_GammaPhy)
check_singularity(germ_count_GammaPhy)

#Median Alpha diversity

germ_count_MedPhy <- glmmTMB(germ_num ~ Med_AlphaPhy*Disturbance + log(suitabilitypc1V3abs) + propagule_pressure +
                               log(suitabilitypc2V3abs) + Density +
                               (1|Family) + (1|POWO.name) + (1|Site),
                             family=nbinom2(link = "log"), 
                             data=germinationcountv2)
summary(germ_count_MedPhy)
Anova(germ_count_MedPhy, type="III")
confint(germ_count_MedPhy)
r.squaredGLMM(germ_count_MedPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = germ_count_MedPhy, plot=F)
plot(simulationOutput)

check_overdispersion(germ_count_MedPhy)
check_singularity(germ_count_MedPhy)

#FIRST YEAR SURVIVAL COUNT

surviveyr1countv2 <- filter(surviveyr1v2, surviveyr1_num>0)

#Mean alpha diversity 

surviveyr1_count_Mean_AlphaPhy <- glmmTMB(surviveyr1_num ~ Mean_AlphaPhy*Disturbance + propagule_pressure +
                                            log(suitabilitypc1V3abs) + log(suitabilitypc2V3abs) + Density + Herbivory +
                                            (1|POWO.name) + (1|Site),
                                          family=nbinom2(link = "log"), 
                                          data=surviveyr1countv2)
summary(surviveyr1_count_Mean_AlphaPhy)
Anova(surviveyr1_count_Mean_AlphaPhy, type="III")
confint(surviveyr1_count_Mean_AlphaPhy)
r.squaredGLMM(surviveyr1_count_Mean_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr1_count_Mean_AlphaPhy, plot=F)
plot(simulationOutput)

check_overdispersion(surviveyr1_count_Mean_AlphaPhy)
check_singularity(surviveyr1_count_Mean_AlphaPhy)

#Maximum alpha diversity

surviveyr1_count_Max_AlphaPhy <- glmmTMB(surviveyr1_num ~ Max_AlphaPhy*Disturbance + propagule_pressure +
                                           log(suitabilitypc1V3abs) + log(suitabilitypc2V3abs) + Density + Herbivory +
                                           (1|POWO.name) + (1|Site),
                                         family=nbinom2(link = "log"), 
                                         data=surviveyr1countv2)
summary(surviveyr1_count_Max_AlphaPhy)
Anova(surviveyr1_count_Max_AlphaPhy, type="III")
confint(surviveyr1_count_Max_AlphaPhy)
r.squaredGLMM(surviveyr1_count_Max_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr1_count_Max_AlphaPhy, plot=F)
plot(simulationOutput)

check_overdispersion(surviveyr1_count_Max_AlphaPhy)
check_singularity(surviveyr1_count_Max_AlphaPhy)

#Gamma diversity

surviveyr1_count_GammaPhy <- glmmTMB(surviveyr1_num ~ log(GammaPhy) + Disturbance + propagule_pressure +
                                       log(suitabilitypc1V3abs) + log(suitabilitypc2V3abs) + Density + Herbivory +
                                       (1|POWO.name) + (1|Site),
                                     family=nbinom2(link = "log"), 
                                     data=surviveyr1countv2)
summary(surviveyr1_count_GammaPhy)
Anova(surviveyr1_count_GammaPhy, type="III")
confint(surviveyr1_count_GammaPhy)
r.squaredGLMM(surviveyr1_count_GammaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr1_count_GammaPhy, plot=F)
plot(simulationOutput)

check_overdispersion(surviveyr1_count_GammaPhy)
check_singularity(surviveyr1_count_GammaPhy)

#Median alpha diversity

surviveyr1_count_MedPhy <- glmmTMB(surviveyr1_num ~ Med_AlphaPhy*Disturbance + propagule_pressure +
                                     log(suitabilitypc1V3abs) + log(suitabilitypc2V3abs) + Density + Herbivory +
                                     (1|POWO.name) + (1|Site),
                                   family=nbinom2(link = "log"), 
                                   data=surviveyr1countv2)
summary(surviveyr1_count_MedPhy)
Anova(surviveyr1_count_MedPhy, type="III")
confint(surviveyr1_count_MedPhy)
r.squaredGLMM(surviveyr1_count_MedPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr1_count_MedPhy, plot=F)
plot(simulationOutput)

check_overdispersion(surviveyr1_count_MedPhy)
check_singularity(surviveyr1_count_MedPhy)

#EARLY SECOND YEAR SURVIVAL COUNT

surviveyr2earlycountv2 <- filter(surviveyr2earlyv2, surviveyr2early_num>0)

#Mean alpha diversity 

surviveyr2early_count_Mean_AlphaPhy <- glmmTMB(surviveyr2early_num ~ Mean_AlphaPhy + Disturbance + propagule_pressure +
                                                 log(suitabilitypc1V3abs) + log(suitabilitypc2V3abs) + Density + Herbivory +
                                                 (1|POWO.name) + (1|Site),
                                               family=nbinom2(link = "log"), 
                                               data=surviveyr2earlycountv2)
summary(surviveyr2early_count_Mean_AlphaPhy)
Anova(surviveyr2early_count_Mean_AlphaPhy, type="III")
confint(surviveyr2early_count_Mean_AlphaPhy)
r.squaredGLMM(surviveyr2early_count_Mean_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2early_count_Mean_AlphaPhy, plot=F)
plot(simulationOutput)

check_overdispersion(surviveyr2early_count_Mean_AlphaPhy)
check_singularity(surviveyr2early_count_Mean_AlphaPhy)

#Max alpha diversity 

surviveyr2early_count_Max_AlphaPhy <- glmmTMB(surviveyr2early_num ~ Max_AlphaPhy + Disturbance + propagule_pressure +
                                                log(suitabilitypc1V3abs) + log(suitabilitypc2V3abs) + Density + Herbivory +
                                                (1|POWO.name) + (1|Site),
                                              family=nbinom2(link = "log"), 
                                              data=surviveyr2earlycountv2)
summary(surviveyr2early_count_Max_AlphaPhy)
Anova(surviveyr2early_count_Max_AlphaPhy, type="III")
confint(surviveyr2early_count_Max_AlphaPhy)
r.squaredGLMM(surviveyr2early_count_Max_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2early_count_Max_AlphaPhy, plot=F)
plot(simulationOutput)

check_overdispersion(surviveyr2early_count_Max_AlphaPhy)
check_singularity(surviveyr2early_count_Max_AlphaPhy)

#Gamma diversity 

surviveyr2early_count_GammaPhy <- glmmTMB(surviveyr2early_num ~ log(GammaPhy) + Disturbance + propagule_pressure +
                                            log(suitabilitypc1V3abs) + log(suitabilitypc2V3abs) + Density + Herbivory +
                                            (1|POWO.name) + (1|Site),
                                          family=nbinom2(link = "log"), 
                                          data=surviveyr2earlycountv2)
summary(surviveyr2early_count_GammaPhy)
Anova(surviveyr2early_count_GammaPhy, type="III")
confint(surviveyr2early_count_GammaPhy)
r.squaredGLMM(surviveyr2early_count_GammaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2early_count_GammaPhy, plot=F)
plot(simulationOutput)

check_overdispersion(surviveyr2early_count_GammaPhy)
check_singularity(surviveyr2early_count_GammaPhy)

#Median alpha diversity 

surviveyr2early_count_MedPhy <- glmmTMB(surviveyr2early_num ~ Med_AlphaPhy + Disturbance + propagule_pressure +
                                          log(suitabilitypc1V3abs) + log(suitabilitypc2V3abs) + Density + Herbivory +
                                          (1|POWO.name) + (1|Site),
                                        family=nbinom2(link = "log"), 
                                        data=surviveyr2earlycountv2)
summary(surviveyr2early_count_MedPhy)
Anova(surviveyr2early_count_MedPhy, type="III")
confint(surviveyr2early_count_MedPhy)
r.squaredGLMM(surviveyr2early_count_MedPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2early_count_MedPhy, plot=F)
plot(simulationOutput)

check_overdispersion(surviveyr2early_count_MedPhy)
check_singularity(surviveyr2early_count_MedPhy)

#LATE SECOND YEAR SURVIVAL COUNT

surviveyr2latecountv2 <- filter(surviveyr2latev2, surviveyr2late_num>0)

#Mean alpha diversity 

surviveyr2late_count_Mean_AlphaPhy <- glmmTMB(surviveyr2late_num ~ Mean_AlphaPhy + Disturbance + propagule_pressure +
                                                log(suitabilitypc1V3abs) + log(suitabilitypc2V3abs) + Density +
                                                (1|POWO.name),
                                              family=nbinom2(link = "log"), 
                                              data=surviveyr2latecountv2)
summary(surviveyr2late_count_Mean_AlphaPhy)
Anova(surviveyr2late_count_Mean_AlphaPhy, type="III")
confint(surviveyr2late_count_Mean_AlphaPhy)
r.squaredGLMM(surviveyr2late_count_Mean_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2late_count_Mean_AlphaPhy, plot=F)
plot(simulationOutput)

check_overdispersion(surviveyr2late_count_Mean_AlphaPhy)
check_singularity(surviveyr2late_count_Mean_AlphaPhy)

#Max alpha diversity 

surviveyr2late_count_Max_AlphaPhy <- glmmTMB(surviveyr2late_num ~ Max_AlphaPhy + Disturbance + propagule_pressure +
                                               log(suitabilitypc1V3abs) + log(suitabilitypc2V3abs) + Density +
                                               (1|POWO.name),
                                             family=nbinom2(link = "log"), 
                                             data=surviveyr2latecountv2)
summary(surviveyr2late_count_Max_AlphaPhy)
Anova(surviveyr2late_count_Max_AlphaPhy, type="III")
confint(surviveyr2late_count_Max_AlphaPhy)
r.squaredGLMM(surviveyr2late_count_Max_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2late_count_Max_AlphaPhy, plot=F)
plot(simulationOutput)

check_overdispersion(surviveyr2late_count_Max_AlphaPhy)
check_singularity(surviveyr2late_count_Max_AlphaPhy)

#Gamma diversity 

surviveyr2late_count_GammaPhy <- glmmTMB(surviveyr2late_num ~ log(GammaPhy) + Disturbance + propagule_pressure +
                                           log(suitabilitypc1V3abs) + log(suitabilitypc2V3abs) + Density +
                                           (1|POWO.name),
                                         family=nbinom2(link = "log"), 
                                         data=surviveyr2latecountv2)
summary(surviveyr2late_count_GammaPhy)
Anova(surviveyr2late_count_GammaPhy, type="III")
confint(surviveyr2late_count_GammaPhy)
r.squaredGLMM(surviveyr2late_count_GammaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2late_count_GammaPhy, plot=F)
plot(simulationOutput)

check_overdispersion(surviveyr2late_count_GammaPhy)
check_singularity(surviveyr2late_count_GammaPhy)

#median alpha diversity 

surviveyr2late_count_MedPhy <- glmmTMB(surviveyr2late_num ~ Med_AlphaPhy*log(suitabilitypc2V3abs) + Disturbance + propagule_pressure +
                                         log(suitabilitypc1V3abs) + Density +
                                         (1|POWO.name),
                                       family=nbinom2(link = "log"), 
                                       data=surviveyr2latecountv2)
summary(surviveyr2late_count_MedPhy)
Anova(surviveyr2late_count_MedPhy, type="III")
confint(surviveyr2late_count_MedPhy)
r.squaredGLMM(surviveyr2late_count_MedPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2late_count_MedPhy, plot=F)
plot(simulationOutput)

check_overdispersion(surviveyr2late_count_MedPhy)
check_singularity(surviveyr2late_count_MedPhy)

#Plot the disturbance interactions for the single climate variable responses

yr2latemeanSCinteraction <- interact_plot(surviveyr2late_YN_Mean_AlphaPhy, pred="Mean_AlphaPhy", modx = "Disturbance",
                                          plot.points=TRUE, partial.residuals = TRUE) +
  scale_y_continuous(trans='log10') +
  xlab("Mean Alpha Phylogenetic Diversity") + ylab("Likelihood of any plant in plot surviving to end of second year") +
  scale_color_manual(values=c("#999999", "#E69F00")) +
  scale_linetype_manual(values=c("solid", "solid")) +
  theme_bw() + guides(linetype = "none") + theme(legend.position = "top")  +
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=10), 
        legend.text=element_text(size=8), 
        legend.title=element_text(size=10))
yr2latemeanSCinteraction

yr2latemaxSCinteraction <- interact_plot(surviveyr2late_YN_Max_AlphaPhy, pred="Max_AlphaPhy", modx = "Disturbance",
                                         plot.points=TRUE, partial.residuals = TRUE) +
  scale_y_continuous(trans='log10') +
  xlab("Maximum Alpha Phylogenetic Diversity") + ylab("Likelihood of any plant in plot surviving to end of second year") +
  scale_color_manual(values=c("#999999", "#E69F00")) +
  scale_linetype_manual(values=c("solid", "solid")) +
  theme_bw() + guides(linetype = "none") + theme(legend.position = "top")  +
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=10), 
        legend.text=element_text(size=8), 
        legend.title=element_text(size=10))
yr2latemaxSCinteraction

yr2lategammaSCinteraction <- interact_plot(surviveyr2late_YN_GammaPhy, pred="GammaPhy", modx = "Disturbance",
                                           plot.points=TRUE, partial.residuals = TRUE) +
  scale_y_continuous(trans='log10') +
  xlab("Gamma Phylogenetic Diversity") + ylab("Likelihood of any plant in plot surviving to end of second year") +
  scale_color_manual(values=c("#999999", "#E69F00")) +
  scale_linetype_manual(values=c("solid", "solid")) +
  theme_bw() + guides(linetype = "none") + theme(legend.position = "top")  +
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=10), 
        legend.text=element_text(size=8), 
        legend.title=element_text(size=10))
yr2lategammaSCinteraction


SCfig <- ggarrange(yr2latemeanSCinteraction, yr2latemaxSCinteraction, 
                   yr2lategammaSCinteraction, nrow=2, ncol=2)
SCfig
ggsave(SCfig, 
       filename = "Fig_S12.svg",
       height = 210, width = 210, units = "mm")

#######################################################################################

#Doing it for Question 2

haueserdatav2 <- readRDS("haueserdatav2.RDS") %>%
  mutate(suitabilitypc1V3abs = abs(suitabilitypc1V3)) %>%
  mutate(suitabilitypc2V3abs = abs(suitabilitypc2V3))
#this contains data for germination, 1st year survival, 2nd year survival
#when analysing, will need to filter when analysing 1st and 2nd year data to only allow plants that germinated
haueserdatav2$awmpd <- as.numeric(haueserdatav2$awmpd)
haueserdatav2$Species <- as.factor(haueserdatav2$Species)
haueserdatav2$Plot <- as.factor(haueserdatav2$Plot)
str(haueserdatav2)

haueserdatav2$Disturbance <- relevel(haueserdatav2$Disturbance, ref="0")

hist(haueserdatav2$suitabilitypc1V3abs)
hist(haueserdatav2$suitabilitypc2V3abs) #need to log these to help meet assumptions

####################################################################

#BINARY DATA

####################################################################

#GERMINATION YES/NO 

#Difference in Mean AlphaPhy 

germ_YN_Mean_AlphaPhy_diff <- glmmTMB(germ_YN ~ Mean_AlphaPhy_diff + Disturbance + log(suitabilitypc1V3abs) +
                                        log(suitabilitypc2V3abs) + Heating + scale(awmpd) + OptGermRt +
                                        (1|Species) + (1|Plot), 
                                      family=binomial(link="logit"), 
                                      data=haueserdatav2)
summary(germ_YN_Mean_AlphaPhy_diff)
Anova(germ_YN_Mean_AlphaPhy_diff, type="III")
confint(germ_YN_Mean_AlphaPhy_diff)
r.squaredGLMM(germ_YN_Mean_AlphaPhy_diff)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = germ_YN_Mean_AlphaPhy_diff, plot=F)
plot(simulationOutput)
check_overdispersion(germ_YN_Mean_AlphaPhy_diff)
check_singularity(germ_YN_Mean_AlphaPhy_diff)

#Difference in Max AlphaPhy 

germ_YN_Max_AlphaPhy_diff <- glmmTMB(germ_YN ~ Max_AlphaPhy_diff*Disturbance + log(suitabilitypc1V3abs) +
                                       log(suitabilitypc2V3abs) + Heating + awmpd + OptGermRt +
                                       (1|Species) + (1|Plot), 
                                     family=binomial(link="logit"), 
                                     data=haueserdatav2)
summary(germ_YN_Max_AlphaPhy_diff)
Anova(germ_YN_Max_AlphaPhy_diff, type="III")
confint(germ_YN_Max_AlphaPhy_diff)
r.squaredGLMM(germ_YN_Max_AlphaPhy_diff)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = germ_YN_Max_AlphaPhy_diff, plot=F)
plot(simulationOutput)
check_overdispersion(germ_YN_Max_AlphaPhy_diff)
check_singularity(germ_YN_Max_AlphaPhy_diff)

#Difference in GammaPhy 

germ_YN_GammaPhy_diff <- glmmTMB(germ_YN ~ log(suitabilitypc2V3abs) + scale(GammaPhy_diff)*Disturbance + 
                                   log(suitabilitypc1V3abs) + Heating + awmpd + OptGermRt +
                                   (1|Species) + (1|Plot), 
                                 family=binomial(link="logit"), 
                                 data=haueserdatav2)
summary(germ_YN_GammaPhy_diff)
Anova(germ_YN_GammaPhy_diff, type="III")
confint(germ_YN_GammaPhy_diff)
r.squaredGLMM(germ_YN_GammaPhy_diff)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = germ_YN_GammaPhy_diff, plot=F)
plot(simulationOutput)
check_overdispersion(germ_YN_GammaPhy_diff)
check_singularity(germ_YN_GammaPhy_diff)

#Difference in Median AlphaPhy 

germ_YN_Med_AlphaPhy_diff <- glmmTMB(germ_YN ~ log(suitabilitypc2V3abs) + Med_AlphaPhy_diff*Disturbance + 
                                       log(suitabilitypc1V3abs) + Heating + awmpd + OptGermRt +
                                       (1|Species) + (1|Plot), 
                                     family=binomial(link="logit"), 
                                     data=haueserdatav2)
summary(germ_YN_Med_AlphaPhy_diff)
Anova(germ_YN_Med_AlphaPhy_diff, type="III")
confint(germ_YN_Med_AlphaPhy_diff)
r.squaredGLMM(germ_YN_Med_AlphaPhy_diff)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = germ_YN_Med_AlphaPhy_diff, plot=F)
plot(simulationOutput)
check_overdispersion(germ_YN_Med_AlphaPhy_diff)
check_singularity(germ_YN_Med_AlphaPhy_diff)

##########################################################

#FIRST YEAR SURVIVAL

haueserdatafirstyearv2 <- haueserdatav2 %>% filter(germ_YN>0 | surviveyr1_YN>0)

#Difference in Mean AlphaPhy 

surviveyr1_YN_Mean_AlphaPhy_diff <- glmmTMB(surviveyr1_YN ~ Mean_AlphaPhy_diff + log(suitabilitypc1V3abs) + Disturbance + 
                                              log(suitabilitypc2V3abs) + Heating + scale(awmpd) + OptGermRt +
                                              (1|Species) + (1|Plot), 
                                            family=binomial(link="logit"), 
                                            data=haueserdatafirstyearv2)
summary(surviveyr1_YN_Mean_AlphaPhy_diff)
Anova(surviveyr1_YN_Mean_AlphaPhy_diff, type="III")
confint(surviveyr1_YN_Mean_AlphaPhy_diff)
r.squaredGLMM(surviveyr1_YN_Mean_AlphaPhy_diff)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr1_YN_Mean_AlphaPhy_diff, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr1_YN_Mean_AlphaPhy_diff)
check_singularity(surviveyr1_YN_Mean_AlphaPhy_diff)
#slight issue with singularity, but this likely due to the plot random variable, which is still 
#important to have in there

#Difference in Max AlphaPhy 

surviveyr1_YN_Max_AlphaPhy_diff <- glmmTMB(surviveyr1_YN ~ Max_AlphaPhy_diff*Disturbance + log(suitabilitypc2V3abs) +    
                                             log(suitabilitypc1V3abs) + Heating + scale(awmpd) + OptGermRt +
                                             (1|Species) + (1|Plot), 
                                           family=binomial(link="logit"), 
                                           data=haueserdatafirstyearv2)
summary(surviveyr1_YN_Max_AlphaPhy_diff)
Anova(surviveyr1_YN_Max_AlphaPhy_diff, type="III")
confint(surviveyr1_YN_Max_AlphaPhy_diff)
r.squaredGLMM(surviveyr1_YN_Max_AlphaPhy_diff)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr1_YN_Max_AlphaPhy_diff, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr1_YN_Max_AlphaPhy_diff)
check_singularity(surviveyr1_YN_Max_AlphaPhy_diff)
#slight issue with singularity, but this likely due to the plot random variable, which is still 
#important to have in there

#Difference in GammaPhy 

surviveyr1_YN_GammaPhy_diff <- glmmTMB(surviveyr1_YN ~ scale(GammaPhy_diff) + log(suitabilitypc2V3abs) + Disturbance +   
                                         log(suitabilitypc1V3abs) + Heating + scale(awmpd) + OptGermRt +
                                         (1|Species) + (1|Plot), 
                                       family=binomial(link="logit"), 
                                       data=haueserdatafirstyearv2)
summary(surviveyr1_YN_GammaPhy_diff)
Anova(surviveyr1_YN_GammaPhy_diff, type="III")
confint(surviveyr1_YN_GammaPhy_diff)
r.squaredGLMM(surviveyr1_YN_GammaPhy_diff)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr1_YN_GammaPhy_diff, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr1_YN_GammaPhy_diff)
check_singularity(surviveyr1_YN_GammaPhy_diff)
#slight issue with singularity, but this likely due to the plot random variable, which is still 
#important to have in there

#Difference in Median AlphaPhy 

surviveyr1_YN_Med_AlphaPhy_diff <- glmmTMB(surviveyr1_YN ~ Med_AlphaPhy_diff*Disturbance + log(suitabilitypc1V3abs) +    
                                             log(suitabilitypc2V3abs) + Heating + awmpd + OptGermRt +
                                             (1|Species) + (1|Plot), 
                                           family=binomial(link="logit"), 
                                           data=haueserdatafirstyearv2)
summary(surviveyr1_YN_Med_AlphaPhy_diff)
Anova(surviveyr1_YN_Med_AlphaPhy_diff, type="III")
confint(surviveyr1_YN_Med_AlphaPhy_diff)
r.squaredGLMM(surviveyr1_YN_Med_AlphaPhy_diff)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr1_YN_Med_AlphaPhy_diff, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr1_YN_Med_AlphaPhy_diff)
check_singularity(surviveyr1_YN_Med_AlphaPhy_diff)
#slight issue with singularity, but this likely due to the plot random variable, which is still 
#important to have in there

##########################################################

#OVERWINTER SURVIVAL

haueserdataoverwinterv2 <- haueserdatav2 %>% filter(surviveyr1_YN>0 | surviveyr2early_YN > 0)

#Difference in Mean AlphaPhy 

surviveyr2early_YN_Mean_AlphaPhy_diff <- glmmTMB(surviveyr2early_YN ~ Mean_AlphaPhy_diff + log(suitabilitypc2V3abs) + Disturbance + 
                                                   log(suitabilitypc1V3abs) + Heating + scale(awmpd) + OptGermRt +
                                                   (1|Species) + (1|Plot), 
                                                 family=binomial(link="logit"), 
                                                 data=haueserdataoverwinterv2)
summary(surviveyr2early_YN_Mean_AlphaPhy_diff)
Anova(surviveyr2early_YN_Mean_AlphaPhy_diff, type="III")
confint(surviveyr2early_YN_Mean_AlphaPhy_diff)
r.squaredGLMM(surviveyr2early_YN_Mean_AlphaPhy_diff)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2early_YN_Mean_AlphaPhy_diff, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr2early_YN_Mean_AlphaPhy_diff)
check_singularity(surviveyr2early_YN_Mean_AlphaPhy_diff)
#slight issue with singularity, but this likely due to the plot random variable, which is still 
#important to have in there

#Difference in Max AlphaPhy 

surviveyr2early_YN_Max_AlphaPhy_diff <- glmmTMB(surviveyr2early_YN ~ Max_AlphaPhy_diff + log(suitabilitypc1V3abs) + Disturbance + 
                                                  log(suitabilitypc2V3abs) + Heating + scale(awmpd) + OptGermRt +
                                                  (1|Species) + (1|Plot), 
                                                family=binomial(link="logit"), 
                                                data=haueserdataoverwinterv2)
summary(surviveyr2early_YN_Max_AlphaPhy_diff)
Anova(surviveyr2early_YN_Max_AlphaPhy_diff, type="III")
confint(surviveyr2early_YN_Max_AlphaPhy_diff)
r.squaredGLMM(surviveyr2early_YN_Max_AlphaPhy_diff)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2early_YN_Max_AlphaPhy_diff, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr2early_YN_Max_AlphaPhy_diff)
check_singularity(surviveyr2early_YN_Max_AlphaPhy_diff)
#slight issue with singularity, but this likely due to the plot random variable, which is still 
#important to have in there

#Difference in GammaPhy #no interactions

surviveyr2early_YN_GammaPhy_diff <- glmmTMB(surviveyr2early_YN ~ scale(GammaPhy_diff) + log(suitabilitypc2V3abs) +  Disturbance + 
                                              log(suitabilitypc1V3abs) + Heating + scale(awmpd) + OptGermRt +
                                              (1|Species) + (1|Plot), 
                                            family=binomial(link="logit"), 
                                            data=haueserdataoverwinterv2)
summary(surviveyr2early_YN_GammaPhy_diff)
Anova(surviveyr2early_YN_GammaPhy_diff, type="III")
confint(surviveyr2early_YN_GammaPhy_diff)
r.squaredGLMM(surviveyr2early_YN_GammaPhy_diff)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2early_YN_GammaPhy_diff, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr2early_YN_GammaPhy_diff)
check_singularity(surviveyr2early_YN_GammaPhy_diff)
#slight issue with singularity, but this likely due to the plot random variable, which is still 
#important to have in there

#Difference in Median AlphaPhy 

surviveyr2early_YN_Med_AlphaPhy_diff <- glmmTMB(surviveyr2early_YN ~ Med_AlphaPhy_diff + Disturbance + log(suitabilitypc1V3abs) + 
                                                  log(suitabilitypc2V3abs) + Heating + scale(awmpd) + OptGermRt +
                                                  (1|Species) + (1|Plot), 
                                                family=binomial(link="logit"), 
                                                data=haueserdataoverwinterv2)
summary(surviveyr2early_YN_Med_AlphaPhy_diff)
Anova(surviveyr2early_YN_Med_AlphaPhy_diff, type="III")
confint(surviveyr2early_YN_Med_AlphaPhy_diff)
r.squaredGLMM(surviveyr2early_YN_Med_AlphaPhy_diff)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2early_YN_Med_AlphaPhy_diff, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr2early_YN_Med_AlphaPhy_diff)
check_singularity(surviveyr2early_YN_Med_AlphaPhy_diff)
#slight issue with singularity, but this likely due to the plot random variable, which is still 
#important to have in there

##########################################################

#LATE SECOND YEAR SURVIVAL

haueserdatayr2latev2 <- haueserdatav2 %>% filter(surviveyr2early_YN > 0 | surviveyr2late_YN>0)

#Difference in Mean AlphaPhy 

surviveyr2late_YN_Mean_AlphaPhy_diff <- glmmTMB(surviveyr2late_YN ~ Mean_AlphaPhy_diff*Disturbance + log(suitabilitypc1V3abs) + 
                                                  log(suitabilitypc2V3abs) + Heating + scale(awmpd) + OptGermRt +
                                                  (1|Species) + (1|Plot), 
                                                family=binomial(link="logit"), 
                                                data=haueserdatayr2latev2)
summary(surviveyr2late_YN_Mean_AlphaPhy_diff)
Anova(surviveyr2late_YN_Mean_AlphaPhy_diff, type="III")
confint(surviveyr2late_YN_Mean_AlphaPhy_diff)
r.squaredGLMM(surviveyr2late_YN_Mean_AlphaPhy_diff)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2late_YN_Mean_AlphaPhy_diff, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr2late_YN_Mean_AlphaPhy_diff)
check_singularity(surviveyr2late_YN_Mean_AlphaPhy_diff)

#Difference in Max AlphaPhy

surviveyr2late_YN_Max_AlphaPhy_diff <- glmmTMB(surviveyr2late_YN ~ Max_AlphaPhy_diff + Disturbance + log(suitabilitypc1V3abs) +  
                                                 log(suitabilitypc2V3abs) + Heating + scale(awmpd) + OptGermRt +
                                                 (1|Species) + (1|Plot), 
                                               family=binomial(link="logit"), 
                                               data=haueserdatayr2latev2)
summary(surviveyr2late_YN_Max_AlphaPhy_diff)
Anova(surviveyr2late_YN_Max_AlphaPhy_diff, type="III")
confint(surviveyr2late_YN_Max_AlphaPhy_diff)
r.squaredGLMM(surviveyr2late_YN_Max_AlphaPhy_diff)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2late_YN_Max_AlphaPhy_diff, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr2late_YN_Max_AlphaPhy_diff)
check_singularity(surviveyr2late_YN_Max_AlphaPhy_diff)

#Difference in GammaPhy

surviveyr2late_YN_GammaPhy_diff <- glmmTMB(surviveyr2late_YN ~ scale(GammaPhy_diff) +  Disturbance + log(suitabilitypc1V3abs) + 
                                             log(suitabilitypc2V3abs) + Heating + scale(awmpd) + OptGermRt +
                                             (1|Species) + (1|Plot), 
                                           family=binomial(link="logit"), 
                                           data=haueserdatayr2latev2)
summary(surviveyr2late_YN_GammaPhy_diff)
Anova(surviveyr2late_YN_GammaPhy_diff, type="III")
confint(surviveyr2late_YN_GammaPhy_diff)
r.squaredGLMM(surviveyr2late_YN_GammaPhy_diff)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2late_YN_GammaPhy_diff, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr2late_YN_GammaPhy_diff)
check_singularity(surviveyr2late_YN_GammaPhy_diff)

#Difference in Median AlphaPhy 

surviveyr2late_YN_Med_AlphaPhy_diff <- glmmTMB(surviveyr2late_YN ~ Med_AlphaPhy_diff*Disturbance + log(suitabilitypc1V3abs) + 
                                                 log(suitabilitypc2V3abs) + Heating + scale(awmpd) + OptGermRt +
                                                 (1|Species) + (1|Plot), 
                                               family=binomial(link="logit"), 
                                               data=haueserdatayr2latev2)
summary(surviveyr2late_YN_Med_AlphaPhy_diff)
Anova(surviveyr2late_YN_Med_AlphaPhy_diff, type="III")
confint(surviveyr2late_YN_Med_AlphaPhy_diff)
r.squaredGLMM(surviveyr2late_YN_Med_AlphaPhy_diff)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2late_YN_Med_AlphaPhy_diff, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr2late_YN_Med_AlphaPhy_diff)
check_singularity(surviveyr2late_YN_Med_AlphaPhy_diff)

##########################################################

#FLOWERING

#there was NO flowering in undisturbed plots
haueserdisturbv2 <- filter(haueserdatav2, Disturbance=="1")

#Difference in Mean AlphaPhy 

flowers_YN_Mean_AlphaPhy_diff <- glmmTMB(flowers_YN ~ Mean_AlphaPhy_diff + log(suitabilitypc2V3abs) +
                                           log(suitabilitypc1V3abs) + Heating + scale(awmpd) + OptGermRt +
                                           (1|Species) + (1|Plot), 
                                         family=binomial(link="logit"), 
                                         data=haueserdisturbv2)
summary(flowers_YN_Mean_AlphaPhy_diff)
Anova(flowers_YN_Mean_AlphaPhy_diff, type="III")
confint(flowers_YN_Mean_AlphaPhy_diff)
r.squaredGLMM(flowers_YN_Mean_AlphaPhy_diff)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = flowers_YN_Mean_AlphaPhy_diff, plot=F)
plot(simulationOutput)
check_overdispersion(flowers_YN_Mean_AlphaPhy_diff)
check_singularity(flowers_YN_Mean_AlphaPhy_diff)

#Difference in Max AlphaPhy 

flowers_YN_Max_AlphaPhy_diff <- glmmTMB(flowers_YN ~ Max_AlphaPhy_diff + log(suitabilitypc2V3abs) +  
                                          log(suitabilitypc1V3abs) + Heating + scale(awmpd) + OptGermRt +
                                          (1|Species) + (1|Plot), 
                                        family=binomial(link="logit"), 
                                        data=haueserdisturbv2)
summary(flowers_YN_Max_AlphaPhy_diff)
Anova(flowers_YN_Max_AlphaPhy_diff, type="III")
confint(flowers_YN_Max_AlphaPhy_diff)
r.squaredGLMM(flowers_YN_Max_AlphaPhy_diff)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = flowers_YN_Max_AlphaPhy_diff, plot=F)
plot(simulationOutput)
check_overdispersion(flowers_YN_Max_AlphaPhy_diff)
check_singularity(flowers_YN_Max_AlphaPhy_diff)

#Difference in GammaPhy

flowers_YN_GammaPhy_diff <- glmmTMB(flowers_YN ~ scale(GammaPhy_diff) + log(suitabilitypc2V3abs) + 
                                      log(suitabilitypc1V3abs) + Heating + scale(awmpd) + OptGermRt +
                                      (1|Species) + (1|Plot), 
                                    family=binomial(link="logit"), 
                                    data=haueserdisturbv2)
summary(flowers_YN_GammaPhy_diff)
Anova(flowers_YN_GammaPhy_diff, type="III")
confint(flowers_YN_GammaPhy_diff)
r.squaredGLMM(flowers_YN_GammaPhy_diff)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = flowers_YN_GammaPhy_diff, plot=F)
plot(simulationOutput)
check_overdispersion(flowers_YN_GammaPhy_diff)
check_singularity(flowers_YN_GammaPhy_diff)

#Difference in Med AlphaPhy 

flowers_YN_Med_AlphaPhy_diff <- glmmTMB(flowers_YN ~ Med_AlphaPhy_diff + log(suitabilitypc1V3abs) +  
                                          log(suitabilitypc2V3abs) + Heating + scale(awmpd) + OptGermRt +
                                          (1|Species) + (1|Plot), 
                                        family=binomial(link="logit"), 
                                        data=haueserdisturbv2)
summary(flowers_YN_Med_AlphaPhy_diff)
Anova(flowers_YN_Med_AlphaPhy_diff, type="III")
confint(flowers_YN_Med_AlphaPhy_diff)
r.squaredGLMM(flowers_YN_Med_AlphaPhy_diff)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = flowers_YN_Med_AlphaPhy_diff, plot=F)
plot(simulationOutput)
check_overdispersion(flowers_YN_Med_AlphaPhy_diff)
check_singularity(flowers_YN_Med_AlphaPhy_diff)

####################################################################

#COUNT DATA

####################################################################

#GERMINATION COUNT

hauesergermcountv2 <- filter(haueserdatav2, germ_num>0)

#Mean AlphaPhy difference 

germ_count_Mean_AlphaPhy_diff <- glmmTMB(germ_num ~ Mean_AlphaPhy_diff + Disturbance + log(suitabilitypc1V3abs) + 
                                           log(suitabilitypc2V3abs) + Heating + scale(awmpd) + OptGermRt +
                                           (1|Species) + (1|Plot),
                                         family=poisson(link = "log"), 
                                         data=hauesergermcountv2)
summary(germ_count_Mean_AlphaPhy_diff)
Anova(germ_count_Mean_AlphaPhy_diff, type="III")
confint(germ_count_Mean_AlphaPhy_diff)
r.squaredGLMM(germ_count_Mean_AlphaPhy_diff)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = germ_count_Mean_AlphaPhy_diff, plot=F)
plot(simulationOutput)
check_overdispersion(germ_count_Mean_AlphaPhy_diff)
check_singularity(germ_count_Mean_AlphaPhy_diff)
check_model(germ_count_Mean_AlphaPhy_diff)

#Max AlphaPhy difference

germ_count_Max_AlphaPhy_diff <- glmmTMB(germ_num ~ Max_AlphaPhy_diff + log(suitabilitypc2V3abs) + 
                                          Disturbance + log(suitabilitypc1V3abs) + 
                                          Heating + scale(awmpd) + OptGermRt +
                                          (1|Species) + (1|Plot),
                                        family=poisson(link = "log"), 
                                        data=hauesergermcountv2)
summary(germ_count_Max_AlphaPhy_diff)
Anova(germ_count_Max_AlphaPhy_diff, type="III")
confint(germ_count_Max_AlphaPhy_diff)
r.squaredGLMM(germ_count_Max_AlphaPhy_diff)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = germ_count_Max_AlphaPhy_diff, plot=F)
plot(simulationOutput)
check_overdispersion(germ_count_Max_AlphaPhy_diff)
check_singularity(germ_count_Max_AlphaPhy_diff)

#GammaPhy difference 

germ_count_GammaPhy_diff <- glmmTMB(germ_num ~ scale(GammaPhy_diff) + log(suitabilitypc1V3abs) + Disturbance +  
                                      log(suitabilitypc2V3abs) + Heating + scale(awmpd) + OptGermRt +
                                      (1|Species) + (1|Plot),
                                    family=poisson(link = "log"), 
                                    data=hauesergermcountv2)
summary(germ_count_GammaPhy_diff)
Anova(germ_count_GammaPhy_diff, type="III")
confint(germ_count_GammaPhy_diff)
r.squaredGLMM(germ_count_GammaPhy_diff)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = germ_count_GammaPhy_diff, plot=F)
plot(simulationOutput)
check_overdispersion(germ_count_GammaPhy_diff)
check_singularity(germ_count_GammaPhy_diff)

#Med AlphaPhy difference

germ_count_Med_AlphaPhy_diff <- glmmTMB(germ_num ~ Med_AlphaPhy_diff + log(suitabilitypc1V3abs) + Disturbance +  
                                          log(suitabilitypc2V3abs) + Heating + scale(awmpd) + OptGermRt +
                                          (1|Species) + (1|Plot),
                                        family=poisson(link = "log"), 
                                        data=hauesergermcountv2)
summary(germ_count_Med_AlphaPhy_diff)
Anova(germ_count_Med_AlphaPhy_diff, type="III")
confint(germ_count_Med_AlphaPhy_diff)
r.squaredGLMM(germ_count_Med_AlphaPhy_diff)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = germ_count_Med_AlphaPhy_diff, plot=F)
plot(simulationOutput)
check_overdispersion(germ_count_Med_AlphaPhy_diff)
check_singularity(germ_count_Med_AlphaPhy_diff)
check_model(germ_count_Med_AlphaPhy_diff)

####################################################################

#FIRST YEAR SURVIVAL COUNT

haueseryr1countv2 <- filter(haueserdatav2, surviveyr1_num>0)

#Mean AlphaPhy difference 

yr1_count_Mean_AlphaPhy_diff <- glmmTMB(surviveyr1_num ~ Mean_AlphaPhy_diff + Disturbance + log(suitabilitypc1V3abs) + 
                                          log(suitabilitypc2V3abs) + Heating + scale(awmpd) + OptGermRt +
                                          (1|Species) + (1|Plot),
                                        family=poisson(link = "log"), 
                                        data=haueseryr1countv2)
summary(yr1_count_Mean_AlphaPhy_diff)
Anova(yr1_count_Mean_AlphaPhy_diff, type="III")
confint(yr1_count_Mean_AlphaPhy_diff)
r.squaredGLMM(yr1_count_Mean_AlphaPhy_diff)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = yr1_count_Mean_AlphaPhy_diff, plot=F)
plot(simulationOutput)
check_overdispersion(yr1_count_Mean_AlphaPhy_diff)
check_singularity(yr1_count_Mean_AlphaPhy_diff)

#Max AlphaPhy difference 

yr1_count_Max_AlphaPhy_diff <- glmmTMB(surviveyr1_num ~ Max_AlphaPhy_diff*log(suitabilitypc2V3abs) + Disturbance +  
                                         log(suitabilitypc1V3abs) + Heating + scale(awmpd) + OptGermRt +
                                         (1|Species) + (1|Plot),
                                       family=poisson(link = "log"), 
                                       data=haueseryr1countv2)
summary(yr1_count_Max_AlphaPhy_diff)
Anova(yr1_count_Max_AlphaPhy_diff, type="III")
confint(yr1_count_Max_AlphaPhy_diff)
r.squaredGLMM(yr1_count_Max_AlphaPhy_diff)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = yr1_count_Max_AlphaPhy_diff, plot=F)
plot(simulationOutput)
check_overdispersion(yr1_count_Max_AlphaPhy_diff)
check_singularity(yr1_count_Max_AlphaPhy_diff)

#GammaPhy difference 

yr1_count_GammaPhy_diff <- glmmTMB(surviveyr1_num ~ scale(GammaPhy_diff)*log(suitabilitypc1V3abs) + Disturbance +  
                                     log(suitabilitypc2V3abs) + Heating + scale(awmpd) + OptGermRt +
                                     (1|Species) + (1|Plot),
                                   family=poisson(link = "log"), 
                                   data=haueseryr1countv2)
summary(yr1_count_GammaPhy_diff)
Anova(yr1_count_GammaPhy_diff, type="III")
confint(yr1_count_GammaPhy_diff)
r.squaredGLMM(yr1_count_GammaPhy_diff)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = yr1_count_GammaPhy_diff, plot=F)
plot(simulationOutput)
check_overdispersion(yr1_count_GammaPhy_diff)
check_singularity(yr1_count_GammaPhy_diff)

#Med AlphaPhy difference 

yr1_count_Med_AlphaPhy_diff <- glmmTMB(surviveyr1_num ~ Med_AlphaPhy_diff + log(suitabilitypc1V3abs) +  Disturbance + 
                                         log(suitabilitypc2V3abs) + Heating + scale(awmpd) + OptGermRt +
                                         (1|Species) + (1|Plot),
                                       family=poisson(link = "log"), 
                                       data=haueseryr1countv2)
summary(yr1_count_Med_AlphaPhy_diff)
Anova(yr1_count_Med_AlphaPhy_diff, type="III")
confint(yr1_count_Med_AlphaPhy_diff)
r.squaredGLMM(yr1_count_Med_AlphaPhy_diff)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = yr1_count_Med_AlphaPhy_diff, plot=F)
plot(simulationOutput)
check_overdispersion(yr1_count_Med_AlphaPhy_diff)
check_singularity(yr1_count_Med_AlphaPhy_diff)

####################################################################

#OVERWINTER SURVIVAL COUNT

haueseryr2earlycountv2 <- filter(haueserdatav2, surviveyr2early_num>0)

#Mean AlphaPhy difference 

yr2early_count_Mean_AlphaPhy_diff <- glmmTMB(surviveyr2early_num ~ Mean_AlphaPhy_diff + log(suitabilitypc2V3abs) +
                                               Disturbance + log(suitabilitypc1V3abs) + 
                                               Heating + scale(awmpd) + OptGermRt +
                                               (1|Species) + (1|Plot),
                                             family=poisson(link = "log"), 
                                             data=haueseryr2earlycountv2)
summary(yr2early_count_Mean_AlphaPhy_diff)
Anova(yr2early_count_Mean_AlphaPhy_diff, type="III")
confint(yr2early_count_Mean_AlphaPhy_diff)
r.squaredGLMM(yr2early_count_Mean_AlphaPhy_diff)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = yr2early_count_Mean_AlphaPhy_diff, plot=F)
plot(simulationOutput)
check_overdispersion(yr2early_count_Mean_AlphaPhy_diff)
check_singularity(yr2early_count_Mean_AlphaPhy_diff)

#Max AlphaPhy difference 

yr2early_count_Max_AlphaPhy_diff <- glmmTMB(surviveyr2early_num ~ Max_AlphaPhy_diff + Disturbance + log(suitabilitypc1V3abs) + 
                                              log(suitabilitypc2V3abs) + Heating + scale(awmpd) + OptGermRt +
                                              (1|Species) + (1|Plot),
                                            family=poisson(link = "log"), 
                                            data=haueseryr2earlycountv2)
summary(yr2early_count_Max_AlphaPhy_diff)
Anova(yr2early_count_Max_AlphaPhy_diff, type="III")
confint(yr2early_count_Max_AlphaPhy_diff)
r.squaredGLMM(yr2early_count_Max_AlphaPhy_diff)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = yr2early_count_Max_AlphaPhy_diff, plot=F)
plot(simulationOutput)
check_overdispersion(yr2early_count_Max_AlphaPhy_diff)
check_singularity(yr2early_count_Max_AlphaPhy_diff)

#GammaPhy difference 

yr2early_count_GammaPhy_diff <- glmmTMB(surviveyr2early_num ~ scale(GammaPhy_diff) + log(suitabilitypc2V3abs)
                                        + Disturbance + log(suitabilitypc1V3abs) + 
                                          Heating + scale(awmpd) + OptGermRt +
                                          (1|Species) + (1|Plot),
                                        family=poisson(link = "log"), 
                                        data=haueseryr2earlycountv2)
summary(yr2early_count_GammaPhy_diff)
Anova(yr2early_count_GammaPhy_diff, type="III")
confint(yr2early_count_GammaPhy_diff)
r.squaredGLMM(yr2early_count_GammaPhy_diff)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = yr2early_count_GammaPhy_diff, plot=F)
plot(simulationOutput)
check_overdispersion(yr2early_count_GammaPhy_diff)
check_singularity(yr2early_count_GammaPhy_diff)

#Median AlphaPhy difference 

yr2early_count_Med_AlphaPhy_diff <- glmmTMB(surviveyr2early_num ~ Med_AlphaPhy_diff + log(suitabilitypc2V3abs) + 
                                              Disturbance + log(suitabilitypc1V3abs) + 
                                              Heating + scale(awmpd) + OptGermRt +
                                              (1|Species) + (1|Plot),
                                            family=poisson(link = "log"), 
                                            data=haueseryr2earlycountv2)
summary(yr2early_count_Med_AlphaPhy_diff)
Anova(yr2early_count_Med_AlphaPhy_diff, type="III")
confint(yr2early_count_Med_AlphaPhy_diff)
r.squaredGLMM(yr2early_count_Med_AlphaPhy_diff)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = yr2early_count_Med_AlphaPhy_diff, plot=F)
plot(simulationOutput)
check_overdispersion(yr2early_count_Med_AlphaPhy_diff)
check_singularity(yr2early_count_Med_AlphaPhy_diff)

####################################################################

#LATE SECOND YEAR SURVIVAL COUNT

haueseryr2latecountv2 <- filter(haueserdatav2, surviveyr2late_num>0)

#Mean AlphaPhy difference

yr2late_count_Mean_AlphaPhy_diff <- glmmTMB(surviveyr2late_num ~ Mean_AlphaPhy_diff + log(suitabilitypc2V3abs) +
                                              Disturbance + log(suitabilitypc1V3abs) + 
                                              Heating + scale(awmpd) + OptGermRt +
                                              (1|Species) + (1|Plot),
                                            family=poisson(link = "log"), 
                                            data=haueseryr2latecountv2)
summary(yr2late_count_Mean_AlphaPhy_diff)
Anova(yr2late_count_Mean_AlphaPhy_diff, type="III")
confint(yr2late_count_Mean_AlphaPhy_diff)
r.squaredGLMM(yr2late_count_Mean_AlphaPhy_diff)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = yr2late_count_Mean_AlphaPhy_diff, plot=F)
plot(simulationOutput)
check_overdispersion(yr2late_count_Mean_AlphaPhy_diff)
check_singularity(yr2late_count_Mean_AlphaPhy_diff)

#Max AlphaPhy difference

yr2late_count_Max_AlphaPhy_diff <- glmmTMB(surviveyr2late_num ~ Max_AlphaPhy_diff + Disturbance + log(suitabilitypc1V3abs) + 
                                             log(suitabilitypc2V3abs) + Heating + scale(awmpd) + OptGermRt +
                                             (1|Species) + (1|Plot),
                                           family=poisson(link = "log"), 
                                           data=haueseryr2latecountv2)
summary(yr2late_count_Max_AlphaPhy_diff)
Anova(yr2late_count_Max_AlphaPhy_diff, type="III")
confint(yr2late_count_Max_AlphaPhy_diff)
r.squaredGLMM(yr2late_count_Max_AlphaPhy_diff)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = yr2late_count_Max_AlphaPhy_diff, plot=F)
plot(simulationOutput)
check_overdispersion(yr2late_count_Max_AlphaPhy_diff)
check_singularity(yr2late_count_Max_AlphaPhy_diff)

#GammaPhy difference 

yr2late_count_GammaPhy_diff <- glmmTMB(surviveyr2late_num ~ scale(GammaPhy_diff)*log(suitabilitypc1V3abs) +
                                         Disturbance + log(suitabilitypc2V3abs) + 
                                         Heating + scale(awmpd) + OptGermRt +
                                         (1|Species) + (1|Plot),
                                       family=poisson(link = "log"), 
                                       data=haueseryr2latecountv2)
summary(yr2late_count_GammaPhy_diff)
Anova(yr2late_count_GammaPhy_diff, type="III")
confint(yr2late_count_GammaPhy_diff)
r.squaredGLMM(yr2late_count_GammaPhy_diff)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = yr2late_count_GammaPhy_diff, plot=F)
plot(simulationOutput)
check_overdispersion(yr2late_count_GammaPhy_diff)
check_singularity(yr2late_count_GammaPhy_diff)

#Median AlphaPhy difference 

yr2late_count_Med_AlphaPhy_diff <- glmmTMB(surviveyr2late_num ~ Med_AlphaPhy_diff*log(suitabilitypc2V3abs) + 
                                             Disturbance + log(suitabilitypc1V3abs) + 
                                             Heating + scale(awmpd) + OptGermRt +
                                             (1|Species) + (1|Plot),
                                           family=poisson(link = "log"), 
                                           data=haueseryr2latecountv2)
summary(yr2late_count_Med_AlphaPhy_diff)
Anova(yr2late_count_Med_AlphaPhy_diff, type="III")
confint(yr2late_count_Med_AlphaPhy_diff)
r.squaredGLMM(yr2late_count_Med_AlphaPhy_diff)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = yr2late_count_Med_AlphaPhy_diff, plot=F)
plot(simulationOutput)
check_overdispersion(yr2late_count_Med_AlphaPhy_diff)
check_singularity(yr2late_count_Med_AlphaPhy_diff)

####################################################################

#NUMBER OF FLOWERS

haueserflowersv2 <- filter(haueserdisturbv2, max_flowers>0)
#ONLY FLOWERS IN DISTURBED PLOTS

#Mean AlphaPhy difference 

flowercount_Mean_AlphaPhy_diff <- glmmTMB(log(max_flowers) ~ Mean_AlphaPhy_diff + log(suitabilitypc2V3abs) + log(suitabilitypc1V3abs) + 
                                            Heating + scale(awmpd) + OptGermRt +
                                            (1|Species) + (1|Plot),
                                          family=poisson(link = "log"), 
                                          data=haueserflowersv2)
summary(flowercount_Mean_AlphaPhy_diff)
Anova(flowercount_Mean_AlphaPhy_diff, type="III")
confint(flowercount_Mean_AlphaPhy_diff)
r.squaredGLMM(flowercount_Mean_AlphaPhy_diff)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = flowercount_Mean_AlphaPhy_diff, plot=F)
plot(simulationOutput)
check_overdispersion(flowercount_Mean_AlphaPhy_diff)
check_singularity(flowercount_Mean_AlphaPhy_diff)
check_model(flowercount_Mean_AlphaPhy_diff)

#Max AlphaPhy difference

flowercount_Max_AlphaPhy_diff <- glmmTMB(log(max_flowers) ~ Max_AlphaPhy_diff + log(suitabilitypc2V3abs) + 
                                           log(suitabilitypc1V3abs) + Heating + scale(awmpd) + OptGermRt +
                                           (1|Species) + (1|Plot),
                                         family=poisson(link = "log"), 
                                         data=haueserflowersv2)
summary(flowercount_Max_AlphaPhy_diff)
Anova(flowercount_Max_AlphaPhy_diff, type="III")
confint(flowercount_Max_AlphaPhy_diff)
r.squaredGLMM(flowercount_Max_AlphaPhy_diff)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = flowercount_Max_AlphaPhy_diff, plot=F)
plot(simulationOutput)
check_overdispersion(flowercount_Max_AlphaPhy_diff)
check_singularity(flowercount_Max_AlphaPhy_diff)

#GammaPhy difference 

flowercount_GammaPhy_diff <- glmmTMB(log(max_flowers) ~ scale(GammaPhy_diff) + log(suitabilitypc2V3abs) + 
                                       log(suitabilitypc1V3abs) + 
                                       Heating + scale(awmpd) + OptGermRt +
                                       (1|Species) + (1|Plot),
                                     family=poisson(link = "log"), 
                                     data=haueserflowersv2)
summary(flowercount_GammaPhy_diff)
Anova(flowercount_GammaPhy_diff, type="III")
confint(flowercount_GammaPhy_diff)
r.squaredGLMM(flowercount_GammaPhy_diff)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = flowercount_GammaPhy_diff, plot=F)
plot(simulationOutput)
check_overdispersion(flowercount_GammaPhy_diff)
check_singularity(flowercount_GammaPhy_diff)

#Med AlphaPhy difference #no interaction

flowercount_Med_AlphaPhy_diff <- glmmTMB(log(max_flowers) ~ Med_AlphaPhy_diff + log(suitabilitypc2V3abs) +
                                           log(suitabilitypc1V3abs) + 
                                           Heating + scale(awmpd) + OptGermRt +
                                           (1|Species) + (1|Plot),
                                         family=poisson(link = "log"), 
                                         data=haueserflowersv2)
summary(flowercount_Med_AlphaPhy_diff)
Anova(flowercount_Med_AlphaPhy_diff, type="III")
confint(flowercount_Med_AlphaPhy_diff)
r.squaredGLMM(flowercount_Med_AlphaPhy_diff)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = flowercount_Med_AlphaPhy_diff, plot=F)
plot(simulationOutput)
check_overdispersion(flowercount_Med_AlphaPhy_diff)
check_singularity(flowercount_Med_AlphaPhy_diff)
check_model(flowercount_Med_AlphaPhy_diff)

####################################################################################################

#STRUCTURAL EQUATION MODELS

####################################################################################################

#All three studies

#Read in all four data files (output of 'data cleaning all 3 experiments.R')

germination <- readRDS("germinationrevised.RDS") %>%
  mutate(suitabilitypc1V2abs = abs(suitabilitypc1V2)) %>%
  mutate(suitabilitypc2V2abs = abs(suitabilitypc2V2))

germination2trait <- left_join(germination, fulltraits, by="plant_name_id") %>%
  drop_na(vegheight, SLA)
germination3trait <- germination2trait %>% drop_na(Seed_wt)

surviveyr1 <- readRDS("surviveyr1revised.RDS") %>%
  mutate(suitabilitypc1V2abs = abs(suitabilitypc1V2)) %>%
  mutate(suitabilitypc2V2abs = abs(suitabilitypc2V2))

surviveyr12trait <- left_join(surviveyr1, fulltraits, by="plant_name_id") %>%
  drop_na(vegheight, SLA)
surviveyr13trait <- surviveyr12trait %>% drop_na(Seed_wt)

surviveyr2early <- readRDS("surviveyr2earlyrevised.RDS") %>%
  mutate(suitabilitypc1V2abs = abs(suitabilitypc1V2)) %>%
  mutate(suitabilitypc2V2abs = abs(suitabilitypc2V2))

surviveyr2early2trait <- left_join(surviveyr2early, fulltraits, by="plant_name_id") %>%
  drop_na(vegheight, SLA)
surviveyr2early3trait <- surviveyr2early2trait %>% drop_na(Seed_wt)

surviveyr2late <- readRDS("surviveyr2laterevised.RDS") %>%
  mutate(suitabilitypc1V2abs = abs(suitabilitypc1V2)) %>%
  mutate(suitabilitypc2V2abs = abs(suitabilitypc2V2))

surviveyr2late2trait <- left_join(surviveyr2late, fulltraits, by="plant_name_id") %>%
  drop_na(vegheight, SLA)
surviveyr2late3trait <- surviveyr2late2trait %>% drop_na(Seed_wt)

#Note, models were fit here using lme4::glmer() not glmmTMB::glmmTMB() as the piecewiseSEM version I
#have can't take glmmTMB models; this didn't change the results of the original models

#Below we provide exemplar code for the structure of the analyses. Analyses were carried out by running the same 
#model structure but removing Mean_AlphaPhy from the main models, and then comparing the AIC of the overall 
#SEMs with and without this removal. This was done for both the '2trait' and '3trait' datasets (see Supplementary Methods).
#Here we show the procedure for SEMs using all three traits (vegetative height, SLA, seed weight) i.e. using the '3trait' 
#datasets - for the two-trait results the same procedure can be followed but using the '2trait' datasets and excluding 
#the 'seedmod' and Seed_wt from the models in each case

#All results are in Table S11.

########################

#Germination SEM

#Germination Y/N

#Original model:
##germ_YN_Mean_AlphaPhy <- glmmTMB(germ_YN ~ Mean_AlphaPhy + Disturbance + propagule_pressure +
##                                   suitabilitypc1V2abs + suitabilitypc2V2abs +
##                                   (1|Family) + (1|POWO.name) +
##                                   (1|Site) + offset(log(Density)), 
##                                 family=binomial(link="logit"), 
##                                 data=germination)

SLAmod <- lm(SLA ~ Mean_AlphaPhy, data=germination3trait)
plot(SLAmod)
vegheightmod <- lm(log(vegheight) ~ Mean_AlphaPhy, data=germination3trait)
plot(vegheightmod)
seedmod <- lm(log(Seed_wt) ~ Mean_AlphaPhy, data=germinationcount3trait)

germ_YN_Mean_AlphaPhy <- glmer(germ_YN ~ Mean_AlphaPhy + Disturbance + propagule_pressure + #i.e. run this with and without Mean_AlphaPhy and compare AIC from the summary() call
                                 suitabilitypc1V2abs + suitabilitypc2V2abs + SLA + log(vegheight) +
                                 (1|Family) + (1|POWO.name) + log(Seed_wt) +
                                 (1|Site) + offset(log(Density)), 
                               family=binomial(link="logit"), 
                               data=germination3trait)
summary(germ_YN_Mean_AlphaPhy)

mod <- psem(SLAmod, vegheightmod, seedmod, germ_YN_Mean_AlphaPhy)
summary(mod)

#Germination count

germinationcount2trait <- filter(germination2trait, germ_num>0)
germinationcount3trait <- filter(germination3trait, germ_num>0)

SLAmod <- lm(SLA ~ Mean_AlphaPhy, data=germinationcount3trait)
plot(SLAmod)
vegheightmod <- lm(log(vegheight) ~ Mean_AlphaPhy, data=germinationcount3trait)
plot(vegheightmod)
seedmod <- lm(log(Seed_wt) ~ Mean_AlphaPhy, data=germinationcount3trait)

germ_count_Mean_AlphaPhy <- glmer(germ_num ~ Mean_AlphaPhy + Disturbance + propagule_pressure + #i.e. run this with and without Mean_AlphaPhy and compare AIC from the summary() call
                                    suitabilitypc1V2abs + suitabilitypc2V2abs + Density +
                                    SLA + log(vegheight) + log(Seed_wt) +
                                    (1|Family) + (1|POWO.name) + (1|Study) + (1|Site),
                                  family=poisson(link="log"), 
                                  data=germinationcount3trait)
summary(germ_count_Mean_AlphaPhy)

simulationOutput <- simulateResiduals(fittedModel = germ_count_Mean_AlphaPhy, plot=F)
plot(simulationOutput)

mod <- psem(SLAmod, vegheightmod, seedmod, germ_count_Mean_AlphaPhy)
summary(mod)

#Survive yr1 Y/N SEM

SLAmod <- lm(SLA ~ Mean_AlphaPhy, data=surviveyr13trait)
plot(SLAmod)
vegheightmod <- lm(log(vegheight) ~ Mean_AlphaPhy, data=surviveyr13trait)
plot(vegheightmod)
seedmod <- lm(log(Seed_wt) ~ Mean_AlphaPhy, data=surviveyr13trait)

surviveyr1_YN_Mean_AlphaPhy <- glmer(surviveyr1_YN ~ Mean_AlphaPhy + Disturbance + propagule_pressure + #i.e. run this with and without Mean_AlphaPhy and compare AIC from the summary() call
                                       suitabilitypc1V2abs + suitabilitypc2V2abs + Herbivory + 
                                       SLA + log(vegheight) + log(Seed_wt) +
                                       (1|Family) + (1|POWO.name) + 
                                       (1|Site) + offset(log(Density)), 
                                     family=binomial(link="logit"), 
                                     data=surviveyr13trait)
summary(surviveyr1_YN_Mean_AlphaPhy)

mod <- psem(SLAmod, vegheightmod, seedmod, surviveyr1_YN_Mean_AlphaPhy)
summary(mod)

#Survive year 1 count SEM

surviveyr1count3trait <- filter(surviveyr13trait, surviveyr1_num>0)

SLAmod <- lm(SLA ~ Mean_AlphaPhy, data=surviveyr1count3trait)
plot(SLAmod)
vegheightmod <- lm(log(vegheight) ~ Mean_AlphaPhy, data=surviveyr1count3trait)
plot(vegheightmod)
seedmod <- lm(log(Seed_wt) ~ Mean_AlphaPhy, data=surviveyr1count3trait)

surviveyr1_count_Mean_AlphaPhy <- glmer(surviveyr1_num ~ Mean_AlphaPhy + Disturbance + propagule_pressure + #i.e. run this with and without Mean_AlphaPhy and compare AIC from the summary() call
                                          suitabilitypc1V2abs + suitabilitypc2V2abs + Density + Herbivory +
                                          SLA + log(vegheight) + log(Seed_wt) +
                                          (1|Family) + (1|POWO.name) + (1|Site) + (1|Study),
                                        family=poisson(link="log"), 
                                        data=surviveyr1count3trait)
summary(surviveyr1_count_Mean_AlphaPhy)

mod <- psem(SLAmod, vegheightmod, seedmod, surviveyr1_count_Mean_AlphaPhy)
summary(mod)

##################################################################################

#Haueser et al. only

haueserdata <- readRDS("haueserdata.RDS") %>%
  mutate(suitabilitypc1V2abs = abs(suitabilitypc1V2)) %>%
  mutate(suitabilitypc2V2abs = abs(suitabilitypc2V2))
haueserdata$awmpd <- as.numeric(haueserdata$awmpd)
haueserdata$Species <- as.factor(haueserdata$Species)
haueserdata$Plot <- as.factor(haueserdata$Plot)

haueserdata2trait <- left_join(haueserdata, fulltraits, by="plant_name_id") %>%
  drop_na(vegheight, SLA)
#Don't need to filter further as all have seed weight

str(haueserdata)

#Germination yes/no SEMs

SLAmod <- lm(SLA ~ Mean_AlphaPhy_diff, data=haueserdata2trait)
plot(SLAmod)
vegheightmod <- lm(log(vegheight) ~ Mean_AlphaPhy_diff, data=haueserdata2trait)
plot(vegheightmod)
seedmod <- lm(log(Seed_wt.x) ~ Mean_AlphaPhy_diff, data=haueserdata2trait)

germ_YN_Mean_AlphaPhy <- glmer(germ_YN ~ Mean_AlphaPhy_diff + Disturbance + suitabilitypc1V2abs + #i.e. run this with and without Mean_AlphaPhy_diff and compare AIC from the summary() call
                                 suitabilitypc2V2abs + Heating + scale(awmpd) + OptGermRt.x +
                                 SLA + log(vegheight) + log(Seed_wt.x) +
                                 (1|Species.x) + (1|Plot), 
                               family=binomial(link="logit"), 
                               data=haueserdata2trait)
summary(germ_YN_Mean_AlphaPhy)

mod <- psem(SLAmod, vegheightmod, seedmod, germ_YN_Mean_AlphaPhy)
summary(mod)

#Germination count SEMs

hauesergermcount2trait <- filter(haueserdata2trait, germ_num>0)

SLAmod <- lm(SLA ~ Mean_AlphaPhy_diff, data=hauesergermcount2trait)
plot(SLAmod)
vegheightmod <- lm(log(vegheight) ~ Mean_AlphaPhy_diff, data=hauesergermcount2trait)
plot(vegheightmod)
seedmod <- lm(log(Seed_wt.x) ~ Mean_AlphaPhy_diff, data=hauesergermcount2trait)

germ_count_Mean_AlphaPhy_diff <- glmer(germ_num ~ Mean_AlphaPhy_diff + Disturbance + suitabilitypc1V2abs + #i.e. run this with and without Mean_AlphaPhy_diff and compare AIC from the summary() call
                                         suitabilitypc2V2abs + Heating + scale(awmpd) + OptGermRt.x +
                                         SLA + log(vegheight) + log(Seed_wt.x) +
                                         (1|Species.x) + (1|Plot),
                                       family=poisson(link = "log"), 
                                       data=hauesergermcount2trait)
summary(germ_count_Mean_AlphaPhy_diff)

mod <- psem(SLAmod, vegheightmod, seedmod, germ_count_Mean_AlphaPhy_diff)
summary(mod)

#Survive yr 1 SEMS

haueserdatafirstyear2trait <- haueserdata2trait %>% filter(germ_YN>0 | surviveyr1_YN>0)

SLAmod <- lm(SLA ~ Mean_AlphaPhy_diff, data=haueserdatafirstyear2trait)
plot(SLAmod)
vegheightmod <- lm(log(vegheight) ~ Mean_AlphaPhy_diff, data=haueserdatafirstyear2trait)
plot(vegheightmod)
seedmod <- lm(log(Seed_wt.x) ~ Mean_AlphaPhy_diff, data=haueserdatafirstyear2trait)

surviveyr1_YN_Mean_AlphaPhy_diff <- glmer(surviveyr1_YN ~ Mean_AlphaPhy_diff + suitabilitypc1V2abs + Disturbance + #i.e. run this with and without Mean_AlphaPhy_diff and compare AIC from the summary() call
                                            suitabilitypc2V2abs + Heating + scale(awmpd) + OptGermRt.x +
                                            SLA + log(vegheight) + log(Seed_wt.x) +
                                            (1|Species.x) + (1|Plot), 
                                          family=binomial(link="logit"), 
                                          data=haueserdatafirstyear2trait)
summary(surviveyr1_YN_Mean_AlphaPhy_diff)

mod <- psem(SLAmod, vegheightmod, seedmod, surviveyr1_YN_Mean_AlphaPhy_diff)
summary(mod)

#Survive yr 1 count 

haueseryr1count2trait <- filter(haueserdata2trait, surviveyr1_num>0)

SLAmod <- lm(SLA ~ Mean_AlphaPhy_diff, data=haueseryr1count2trait)
plot(SLAmod)
vegheightmod <- lm(log(vegheight) ~ Mean_AlphaPhy_diff, data=haueseryr1count2trait)
plot(vegheightmod)
seedmod <- lm(log(Seed_wt.x) ~ Mean_AlphaPhy_diff, data=haueseryr1count2trait)

yr1_count_Mean_AlphaPhy_diff <- glmer(surviveyr1_num ~ Mean_AlphaPhy_diff + Disturbance + suitabilitypc1V2abs + #i.e. run this with and without Mean_AlphaPhy_diff and compare AIC from the summary() call
                                        suitabilitypc2V2abs + Heating + scale(awmpd) + OptGermRt.x +
                                        SLA + log(vegheight) + log(Seed_wt.x) +
                                        (1|Species.x) + (1|Plot),
                                      family=poisson(link = "log"), 
                                      data=haueseryr1count2trait)
summary(yr1_count_Mean_AlphaPhy_diff)

mod <- psem(SLAmod, vegheightmod, seedmod, yr1_count_Mean_AlphaPhy_diff)
summary(mod)

#Oh my golly gosh, I can't believe you made it this far! Here's a joke for you.

#Where do bad desserts get taken?

#Into custardy




