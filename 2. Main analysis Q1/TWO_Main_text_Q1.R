#Code for 'Plants that evolved under high phylogenetic diversity have higher invasion success, particularly in undisturbed communities'
#Joshua Brian, Mark van Kleunen, Wayne Dawson, Anne Kempel, Weihan Zhao, Jane Catford

#Code prepared by Joshua Brian; joshua.brian@kcl.ac.uk OR jshbrian@gmail.com

##PART TWO: ANALYSIS FOR QUESTION 1 

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
library(arm) #version 1.13-1
library(patchwork) #version 1.3.2.9000

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
  xlab("Standardised effect size ? 95% C.I.") + guides(shape = "none") + theme(legend.position = "top") +
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