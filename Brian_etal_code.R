#Code for 'Plants from regions of high phylogenetic diversity have higher invasion success, particularly in undisturbed communities'
#Joshua Brian, Mark van Kleunen, Wayne Dawson, Anne Kempel, Jane Catford

#Code prepared by Joshua Brian; joshua.brian@kcl.ac.uk OR jshbrian@gmail.com

library(tidyverse)
library(glmmTMB)
library(DHARMa)
library(car)
library(MuMIn)
library(bbmle)
library(visreg)
library(bbmle)
library(lme4)
library(performance)
library(interactions)
library(viridis)
library(ggpubr)
library(cowplot)
library(rtry)

#This code provides all details required to perform the statistical analyses and obtain the 
#results described in the text. It also provides code to produce all results figures 
#(Figs 1-5 and S7-S13).
#Models for all four response metrics are included (mean alpha PD, max alpha PD, gamma PD and 
#median alpha PD) in that order. Results for mean alpha PD correspond to those reported and 
#discussed in the main text, results for the other three metrics are reported in the
#Supporting Information.

#In all cases, models described below are the final selected models 
#(e.g. only interactions found to be significant were maintained in
#models)

#Main analyses for all three research questions are presented first. Supplementary analyses
#as described in the Supporting Information then follow. See the Supporting Information for 
#full details. Note the 'Species origin' and 'Hnd' supplementary analyses are included in Main text
#Question 1 and 2 analyses.

#Full contents:
#Main text Question 1: Line 49
#Main text Question 2: Line 1079
#Main text Question 3: Line 3022
#Supplementary analysis K+M models: Line 3371
#Supplementary analysis K+H models: Line 4037
#Supplementary analysis SC models: Line 4775

##################################################################################################

#QUESTION ONE: Do sown species with higher home-range-PD have higher rates of colonisation and 
#survival than species with low home-range-PD?

##################################################################################################

#Read in all four data files 

germination <- readRDS("germination.RDS") %>%
  mutate(suitabilitypc1V2abs = abs(suitabilitypc1V2)) %>%
  mutate(suitabilitypc2V2abs = abs(suitabilitypc2V2))
str(germination)

surviveyr1 <- readRDS("surviveyr1.RDS") %>%
  mutate(suitabilitypc1V2abs = abs(suitabilitypc1V2)) %>%
  mutate(suitabilitypc2V2abs = abs(suitabilitypc2V2))
str(surviveyr1)

surviveyr2early <- readRDS("surviveyr2early.RDS") %>%
  mutate(suitabilitypc1V2abs = abs(suitabilitypc1V2)) %>%
  mutate(suitabilitypc2V2abs = abs(suitabilitypc2V2))
str(surviveyr2early)

surviveyr2late <- readRDS("surviveyr2late.RDS") %>%
  mutate(suitabilitypc1V2abs = abs(suitabilitypc1V2)) %>%
  mutate(suitabilitypc2V2abs = abs(suitabilitypc2V2))
str(surviveyr2late)

################################################################################

#FIRST, BINARY RESPONSE VARIABLE FOR ALL DATASETS (did any plants in plot survive)

#COLONISATION

#Mean alpha diversity 

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

#Testing whether worth adding native/exotic status

germ_YN_Mean_AlphaPhyN <- glmmTMB(germ_YN ~ Mean_AlphaPhy + Disturbance + propagule_pressure +
                                    suitabilitypc1V2abs + suitabilitypc2V2abs + Status +
                                    (1|Family) + (1|POWO.name) +
                                    (1|Site) + offset(log(Density)), 
                                  family=binomial(link="logit"), 
                                  data=germination)

AICtab(germ_YN_Mean_AlphaPhyN, germ_YN_Mean_AlphaPhy)

#Maximum alpha diversity

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

#FIRST GROWING SEASON

#Mean alpha diversity 

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

#Maximum alpha diversity

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

#OVERWINTER

#Mean alpha diversity

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

# Maximum alpha diversity

surviveyr2early_YN_Max_AlphaPhy <- glmmTMB(surviveyr2early_YN ~ Max_AlphaPhy + Disturbance + propagule_pressure +
                                             suitabilitypc1V2abs + suitabilitypc2V2abs + Herbivory +
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

surviveyr2early_YN_MedPhy <- glmmTMB(surviveyr2early_YN ~ Med_AlphaPhy + Disturbance + propagule_pressure +
                                       suitabilitypc1V2abs + suitabilitypc2V2abs + Herbivory +
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

#SECOND GROWING SEASON

#Mean alpha diversity

surviveyr2late_YN_Mean_AlphaPhy <- glmmTMB(surviveyr2late_YN ~ Mean_AlphaPhy*suitabilitypc1V2abs + Disturbance + propagule_pressure +
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

#Maximum alpha diversity

surviveyr2late_YN_Max_AlphaPhy <- glmmTMB(surviveyr2late_YN ~ Max_AlphaPhy*Disturbance + propagule_pressure +
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

surviveyr2late_YN_GammaPhy <- glmmTMB(surviveyr2late_YN ~ log(GammaPhy)*suitabilitypc1V2abs + log(GammaPhy)*Disturbance + propagule_pressure +
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

surviveyr2late_YN_MedPhy <- glmmTMB(surviveyr2late_YN ~ Med_AlphaPhy + Disturbance + propagule_pressure +
                                      suitabilitypc1V2abs + suitabilitypc2V2abs + 
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

#SECOND, COUNT RESPONSE VARIABLE FOR ALL DATASETS (how many plants in a plot survived?)

#COLONISATION

germinationcount <- filter(germination, germ_num>0)

#Mean alpha diversity 

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

#Maximum alpha diversity

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

#FIRST GROWING SEASON

surviveyr1count <- filter(surviveyr1, surviveyr1_num>0)

#Mean alpha diversity 

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

#Maximum alpha diversity

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

#OVERWINTER 

surviveyr2earlycount <- filter(surviveyr2early, surviveyr2early_num>0)

#Mean alpha diversity 

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

#Max alpha diversity 

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

#SECOND GROWING SEASON

surviveyr2latecount <- filter(surviveyr2late, surviveyr2late_num>0)

#Mean alpha diversity 

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

#Max alpha diversity 

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

#median alpha diversity 

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

#Plot all effect sizes for Mean_AlphaPhy

effectsizes <- read.csv("all3results_effectsize_v2.csv", header=T, stringsAsFactors = T)
effectsizes$Stage = factor(effectsizes$Stage, levels = c("Second growing season","Overwinter","First growing season","Colonisation"))
effectsizes$Variable = factor(effectsizes$Variable, levels = c("Alpha PD x Climate (temp.)", "Alpha PD x Undisturbed","Herbivory", "Climate dissimilarity (temp.)",
                                                               "Climate dissimilarity (precip.)", "Undisturbed","Mean Alpha PD"))

Figure_1 <- ggplot(effectsizes, aes(x=effect, y=Variable, color=Stage, shape=significant)) +
  geom_point(size=2.3, position=position_dodge(0.5)) +
  scale_shape_manual(values=c(2, 1, 16)) +
  geom_errorbar(data=effectsizes, aes(y=Variable, xmin=effect_low, xmax=effect_high), 
                width=0, size=1, position=position_dodge(0.5)) +
  facet_grid(. ~ Metric, scales="free") +
  scale_color_viridis(discrete = TRUE, guide = guide_legend(reverse = TRUE)) +
  theme_bw() + geom_vline(xintercept=0, linetype="solid", color="grey", size=0.9) +
  xlab("Effect size ± 95% C.I.") + guides(shape = "none") + theme(legend.position = "top") +
  theme(axis.text=element_text(size=7),
        axis.title=element_text(size=9), 
        legend.text=element_text(size=7), 
        legend.title=element_text(size=9))
Figure_1

ggsave(Figure_1, 
       filename = "Figure_1.svg",
       height = 130, width = 180, units = "mm")

#MEAN ALPHA PHY - INTERACTIONS

germmeanalphacountinteraction <- interact_plot(germ_count_Mean_AlphaPhy, pred="Mean_AlphaPhy",
                                               modx = "Disturbance",
                                               plot.points=TRUE, partial.residuals = TRUE) +
  scale_y_continuous(trans='log10') + 
  xlab("Mean Alpha Phylogenetic Diversity") + ylab("Number of plants colonising") +
  scale_color_manual(values=c("#999999", "#E69F00")) +
  scale_linetype_manual(values=c("solid", "solid")) +
  theme_bw() + guides(linetype = "none") + theme(legend.position = "top")  +
  theme(axis.text=element_text(size=7),
        axis.title=element_text(size=9), 
        legend.text=element_text(size=7), 
        legend.title=element_text(size=9))
germmeanalphacountinteraction

yr1meanalphacountinteraction <- interact_plot(surviveyr1_count_Mean_AlphaPhy, pred="Mean_AlphaPhy", modx = "Disturbance",
                                              plot.points=TRUE, partial.residuals = TRUE) +
  scale_y_continuous(trans='log10') + 
  xlab("Mean Alpha Phylogenetic Diversity") + ylab("Number of plants surviving to end of first year") +
  scale_color_manual(values=c("#999999", "#E69F00")) +
  scale_linetype_manual(values=c("solid", "solid")) +
  theme_bw() + guides(linetype = "none") + theme(legend.position = "top")  +
  theme(axis.text=element_text(size=7),
        axis.title=element_text(size=9), 
        legend.text=element_text(size=7), 
        legend.title=element_text(size=9))
yr1meanalphacountinteraction

yr2YNtempinteraction <- interact_plot(surviveyr2late_YN_Mean_AlphaPhy, pred="Mean_AlphaPhy", modx = "suitabilitypc1V2abs",
                                      plot.points=TRUE, partial.residuals = TRUE, legend.main = "Climate dissimilarity (temp.)") +
  scale_y_continuous(trans='log10') +
  xlab("Mean Alpha Phylogenetic Diversity") + ylab("Likelihood of any plant in plot surviving to end of second year") +
  scale_linetype_manual(values=c("solid", "solid", "solid")) +
  theme_bw() + guides(linetype = "none") + theme(legend.position = "top")  +
  theme(axis.text=element_text(size=7),
        axis.title=element_text(size=9), 
        legend.text=element_text(size=7), 
        legend.title=element_text(size=9))
yr2YNtempinteraction

Figure_2 <- ggarrange(germmeanalphacountinteraction, yr1meanalphacountinteraction, yr2YNtempinteraction, nrow=2, ncol=2)
Figure_2
ggsave(Figure_2, 
       filename = "Figure_2.svg",
       height = 210, width = 210, units = "mm")

#MAX ALPHA PHY - INTERACTIONS

germmaxalphacountinteraction <- interact_plot(germ_count_Max_AlphaPhy, pred="Max_AlphaPhy", modx = "Disturbance",
                                              plot.points=TRUE, partial.residuals = TRUE) +
  scale_y_continuous(trans='log10') + 
  xlab("Maximum Alpha Phylogenetic Diversity") + ylab("Number of plants colonising") +
  scale_color_manual(values=c("#999999", "#E69F00")) +
  scale_linetype_manual(values=c("solid", "solid")) +
  theme_bw() + guides(linetype = "none") + theme(legend.position = "top")  +
  theme(axis.text=element_text(size=7),
        axis.title=element_text(size=9), 
        legend.text=element_text(size=7), 
        legend.title=element_text(size=9))
germmaxalphacountinteraction

yr1YNmaxinteraction <- interact_plot(surviveyr1_YN_Max_AlphaPhy, pred="Max_AlphaPhy", modx = "Disturbance",
                                     plot.points=TRUE, partial.residuals = TRUE) +
  scale_y_continuous(trans='log10') + 
  xlab("Maximum Alpha Phylogenetic Diversity") + ylab("Likelihood of any plant in plot surviving to end of first year") +
  scale_color_manual(values=c("#999999", "#E69F00")) +
  scale_linetype_manual(values=c("solid", "solid")) +
  theme_bw() + guides(linetype = "none") + theme(legend.position = "top")  +
  theme(axis.text=element_text(size=7),
        axis.title=element_text(size=9), 
        legend.text=element_text(size=7), 
        legend.title=element_text(size=9))
yr1YNmaxinteraction

yr1maxalphacountinteraction <- interact_plot(surviveyr1_count_Max_AlphaPhy, pred="Max_AlphaPhy", modx = "Disturbance",
                                             plot.points=TRUE, partial.residuals = TRUE) +
  scale_y_continuous(trans='log10') + 
  xlab("Maximum Alpha Phylogenetic Diversity") + ylab("Number of plants surviving to end of first year") +
  scale_color_manual(values=c("#999999", "#E69F00")) +
  scale_linetype_manual(values=c("solid", "solid")) +
  theme_bw() + guides(linetype = "none") + theme(legend.position = "top")  +
  theme(axis.text=element_text(size=7),
        axis.title=element_text(size=9), 
        legend.text=element_text(size=7), 
        legend.title=element_text(size=9))
yr1maxalphacountinteraction

yr2YNmaxinteraction <- interact_plot(surviveyr2late_YN_Max_AlphaPhy, pred="Max_AlphaPhy", modx = "Disturbance",
                                     plot.points=TRUE, partial.residuals = TRUE) +
  scale_y_continuous(trans='log10') +
  xlab("Maximum Alpha Phylogenetic Diversity") + ylab("Likelihood of any plant in plot surviving to end of second year") +
  scale_color_manual(values=c("#999999", "#E69F00")) +
  scale_linetype_manual(values=c("solid", "solid")) +
  theme_bw() + guides(linetype = "none") + theme(legend.position = "top")  +
  theme(axis.text=element_text(size=7),
        axis.title=element_text(size=9), 
        legend.text=element_text(size=7), 
        legend.title=element_text(size=9))
yr2YNmaxinteraction

Figure_S7 <- ggarrange(germmaxalphacountinteraction, yr1YNmaxinteraction, 
                         yr1maxalphacountinteraction, yr2YNmaxinteraction, nrow=2, ncol=2)
Figure_S7
ggsave(Figure_S7, 
       filename = "Figure_S7.svg",
       height = 210, width = 210, units = "mm")

#MEDIAN ALPHA PHY - INTERACTIONS

germYNmedinteraction <- interact_plot(germ_YN_MedPhy, pred="Med_AlphaPhy", modx = "Disturbance",
                                      plot.points=TRUE, partial.residuals = TRUE) +
  scale_y_continuous(trans='log10') +
  xlab("Median Alpha Phylogenetic Diversity") + ylab("Likelihood of any plant in plot colonising") +
  scale_color_manual(values=c("#999999", "#E69F00")) +
  scale_linetype_manual(values=c("solid", "solid")) +
  theme_bw() + guides(linetype = "none") + theme(legend.position = "top")  +
  theme(axis.text=element_text(size=7),
        axis.title=element_text(size=9), 
        legend.text=element_text(size=7), 
        legend.title=element_text(size=9))
germYNmedinteraction

germmedalphacountinteraction <- interact_plot(germ_count_MedPhy, pred="Med_AlphaPhy", modx = "Disturbance",
                                              plot.points=TRUE, partial.residuals = TRUE) +
  scale_y_continuous(trans='log10') + 
  xlab("Median Alpha Phylogenetic Diversity") + ylab("Number of plants colonising") +
  scale_color_manual(values=c("#999999", "#E69F00")) +
  scale_linetype_manual(values=c("solid", "solid")) +
  theme_bw() + guides(linetype = "none") + theme(legend.position = "top")  +
  theme(axis.text=element_text(size=7),
        axis.title=element_text(size=9), 
        legend.text=element_text(size=7), 
        legend.title=element_text(size=9))
germmedalphacountinteraction

germcountmedtempinteraction <- interact_plot(germ_count_MedPhy, pred="Med_AlphaPhy", modx = "suitabilitypc1V2abs",
                                             plot.points=TRUE, partial.residuals = TRUE, legend.main = "Climate dissimilarity (temp.)") +
  scale_y_continuous(trans='log10') +
  xlab("Median Alpha Phylogenetic Diversity") + ylab("Number of plants colonising") +
  scale_linetype_manual(values=c("solid", "solid", "solid")) +
  theme_bw() + guides(linetype = "none") + theme(legend.position = "top")  +
  theme(axis.text=element_text(size=7),
        axis.title=element_text(size=9), 
        legend.text=element_text(size=7), 
        legend.title=element_text(size=9))
germcountmedtempinteraction

yr1medalphacountinteraction <- interact_plot(surviveyr1_count_MedPhy, pred="Med_AlphaPhy", modx = "Disturbance",
                                             plot.points=TRUE, partial.residuals = TRUE) +
  scale_y_continuous(trans='log10') + 
  xlab("Median Alpha Phylogenetic Diversity") + ylab("Number of plants surviving to end of first year") +
  scale_color_manual(values=c("#999999", "#E69F00")) +
  scale_linetype_manual(values=c("solid", "solid")) +
  theme_bw() + guides(linetype = "none") + theme(legend.position = "top")  +
  theme(axis.text=element_text(size=7),
        axis.title=element_text(size=9), 
        legend.text=element_text(size=7), 
        legend.title=element_text(size=9))
yr1medalphacountinteraction

yr2latecountmedtempinteraction <- interact_plot(surviveyr2late_count_MedPhy, pred="Med_AlphaPhy", modx = "suitabilitypc2V2abs",
                                                plot.points=TRUE, partial.residuals = TRUE, legend.main = "Climate dissimilarity (precip.)") +
  scale_y_continuous(trans='log10') +
  xlab("Median Alpha Phylogenetic Diversity") + ylab("Number of plants surviving to end of second year") +
  scale_linetype_manual(values=c("solid", "solid", "solid")) +
  theme_bw() + guides(linetype = "none") + theme(legend.position = "top")  +
  theme(axis.text=element_text(size=7),
        axis.title=element_text(size=9), 
        legend.text=element_text(size=7), 
        legend.title=element_text(size=9))
yr2latecountmedtempinteraction

Figure_S8 <- ggarrange(germYNmedinteraction, germmedalphacountinteraction,
                            germcountmedtempinteraction, yr1medalphacountinteraction,
                            yr2latecountmedtempinteraction, nrow=3, ncol=2)
Figure_S8
ggsave(Figure_S8, 
       filename = "Figure_S8.svg",
       height = 315, width = 210, units = "mm")

#GAMMAPHY - INTERACTIONS

yr1YNgammatempinteraction <- interact_plot(surviveyr1_YN_GammaPhy, pred="GammaPhy", modx = "suitabilitypc1V2abs",
                                           plot.points=TRUE, partial.residuals = TRUE, legend.main = "Climate dissimilarity (temp.)") +
  scale_y_continuous(trans='log10') +
  xlab("Gamma Phylogenetic Diversity") + ylab("Likelihood of any plant in plot surviving to end of first year") +
  scale_linetype_manual(values=c("solid", "solid", "solid")) +
  theme_bw() + guides(linetype = "none") + theme(legend.position = "top")  +
  theme(axis.text=element_text(size=7),
        axis.title=element_text(size=9), 
        legend.text=element_text(size=7), 
        legend.title=element_text(size=9))
yr1YNgammatempinteraction

yr2earlygammacountinteraction <- interact_plot(surviveyr2early_count_GammaPhy, pred="GammaPhy", modx = "Disturbance",
                                               plot.points=TRUE, partial.residuals = TRUE) +
  scale_y_continuous(trans='log10') + 
  xlab("Gamma Phylogenetic Diversity") + ylab("Number of plants surviving overwinter") +
  scale_color_manual(values=c("#999999", "#E69F00")) +
  scale_linetype_manual(values=c("solid", "solid")) +
  theme_bw() + guides(linetype = "none") + theme(legend.position = "top")  +
  theme(axis.text=element_text(size=7),
        axis.title=element_text(size=9), 
        legend.text=element_text(size=7), 
        legend.title=element_text(size=9))
yr2earlygammacountinteraction

yr2YNgammainteraction <- interact_plot(surviveyr2late_YN_GammaPhy, pred="GammaPhy", modx = "Disturbance",
                                       plot.points=TRUE, partial.residuals = TRUE) +
  scale_y_continuous(trans='log10') +
  xlab("Gamma Phylogenetic Diversity") + ylab("Likelihood of any plant in plot surviving to end of second year") +
  scale_color_manual(values=c("#999999", "#E69F00")) +
  scale_linetype_manual(values=c("solid", "solid")) +
  theme_bw() + guides(linetype = "none") + theme(legend.position = "top")  +
  theme(axis.text=element_text(size=7),
        axis.title=element_text(size=9), 
        legend.text=element_text(size=7), 
        legend.title=element_text(size=9))
yr2YNgammainteraction

yr2YNgammatempinteraction <- interact_plot(surviveyr2late_YN_GammaPhy, pred="GammaPhy", modx = "suitabilitypc1V2abs",
                                           plot.points=TRUE, partial.residuals = TRUE, legend.main = "Climate dissimilarity (temp.)") +
  scale_y_continuous(trans='log10') +
  xlab("Gamma Phylogenetic Diversity") + ylab("Likelihood of any plant in plot surviving to end of second year") +
  scale_linetype_manual(values=c("solid", "solid", "solid")) +
  theme_bw() + guides(linetype = "none") + theme(legend.position = "top")  +
  theme(axis.text=element_text(size=7),
        axis.title=element_text(size=9), 
        legend.text=element_text(size=7), 
        legend.title=element_text(size=9))
yr2YNgammatempinteraction

Figure_S9 <- ggarrange(yr1YNgammatempinteraction, yr2earlygammacountinteraction, 
                      yr2YNgammainteraction, yr2YNgammatempinteraction, nrow=2, ncol=2)
Figure_S9
ggsave(Figure_S9, 
       filename = "Figure_S9.svg",
       height = 210, width = 210, units = "mm")

##################################################################################################

#QUESTION TWO: Do sown species with higher relative home-range-PD have higher rates of 
#colonisation and survival than species with lower relative home-range-PD?

##################################################################################################

haueserdata <- readRDS("haueserdata.RDS") %>%
  mutate(suitabilitypc1V2abs = abs(suitabilitypc1V2)) %>%
  mutate(suitabilitypc2V2abs = abs(suitabilitypc2V2))
#this contains data for germination, 1st year survival, 2nd year survival
#when analysing, will need to filter when analysing 1st and 2nd year data to only allow plants that germinated
haueserdata$awmpd <- as.numeric(haueserdata$awmpd)
haueserdata$Species <- as.factor(haueserdata$Species)
haueserdata$Plot <- as.factor(haueserdata$Plot)
str(haueserdata)

#For these models, try both absolute values of PD ('raw', as per Q1), and relative PD (the 
#"diff' models). The relative PD diff models are those presented in the final paper,
#and the 'raw' models are those presented in the Supporting Information ('Hnd' models 
#in Tables S4-S8)

###################################################################################

#FIRST, BINARY RESPONSE VARIABLE FOR ALL DATASETS (did any plants in plot survive)

#COLONISATION

#Raw Mean AlphaPhy 

germ_YN_Mean_AlphaPhy <- glmmTMB(germ_YN ~ Mean_AlphaPhy +  Disturbance + suitabilitypc1V2abs + 
                                   suitabilitypc2V2abs + Heating + awmpd + OptGermRt +
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

#FIRST GROWING SEASON

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

##

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

#OVERWINTER 

haueserdataoverwinter <- haueserdata %>% filter(germ_YN>0 | surviveyr1_YN>0 | surviveyr2early_YN > 0)

#Raw Mean AlphaPhy 

surviveyr2early_YN_Mean_AlphaPhy <- glmmTMB(surviveyr2early_YN ~ Mean_AlphaPhy +  Disturbance + suitabilitypc1V2abs +
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

surviveyr2early_YN_Mean_AlphaPhy_diff <- glmmTMB(surviveyr2early_YN ~ Mean_AlphaPhy_diff*suitabilitypc2V2abs + Disturbance + 
                                                   suitabilitypc1V2abs + Heating + scale(awmpd) + OptGermRt +
                                                   (1|Species) + (1|Plot), 
                                                 family=binomial(link="logit"), 
                                                 data=haueserdataoverwinter)
summary(surviveyr2early_YN_Mean_AlphaPhy_diff)
Anova(surviveyr2early_YN_Mean_AlphaPhy_diff, type="III")
confint(surviveyr2early_YN_Mean_AlphaPhy_diff)
r.squaredGLMM(surviveyr2early_YN_Mean_AlphaPhy_diff)

###

#Testing whether native status has an effect

surviveyr2early_YN_Mean_AlphaPhy_diffN <- glmmTMB(surviveyr2early_YN ~ Mean_AlphaPhy_diff*suitabilitypc2V2abs + Disturbance + 
                                                    suitabilitypc1V2abs + Heating + scale(awmpd) + OptGermRt + Status +
                                                    (1|Species) + (1|Plot), 
                                                  family=binomial(link="logit"), 
                                                  data=haueserdataoverwinter)
summary(surviveyr2early_YN_Mean_AlphaPhy_diffN)

AICtab(surviveyr2early_YN_Mean_AlphaPhy_diffN, surviveyr2early_YN_Mean_AlphaPhy_diff)

#Segue: checking climate suitability

surviveyr2early_YN_Mean_AlphaPhy_diffsq <- glmmTMB(surviveyr2early_YN ~ Mean_AlphaPhy_diff + I(suitabilitypc2V2^2) + Disturbance + 
                                                     suitabilitypc1V2 + Heating + scale(awmpd) + OptGermRt +
                                                     (1|Species) + (1|Plot), 
                                                   family=binomial(link="logit"), 
                                                   data=haueserdataoverwinter)
summary(surviveyr2early_YN_Mean_AlphaPhy_diffsq)

AICtab(surviveyr2early_YN_Mean_AlphaPhy_diff, surviveyr2early_YN_Mean_AlphaPhy_diffsq)

Figure_S13a <- visreg(surviveyr2early_YN_Mean_AlphaPhy_diffsq, 
                     "suitabilitypc2V2", gg=TRUE, xlab="Untransformed precipitation dissimilarity",
                     ylab="Scaled likelihood of overwinter survival", line=list(col="black"),
                     fill=list(fill=rgb(red=0.8, green=0.8, blue=0.8, alpha=0.6)),
                     points=list(size=2, pch=16, col="black")) + theme_bw() + 
  theme(axis.text=element_text(size=7), axis.title=element_text(size=9),
        legend.text=element_text(size=7), legend.title=element_text(size=9))
Figure_S13a

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

surviveyr2early_YN_GammaPhy <- glmmTMB(surviveyr2early_YN ~ scale(GammaPhy) +  Disturbance + suitabilitypc1V2abs +
                                         suitabilitypc2V2abs + Heating + awmpd + OptGermRt +
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

surviveyr2early_YN_Med_AlphaPhy <- glmmTMB(surviveyr2early_YN ~ Med_AlphaPhy +  Disturbance + suitabilitypc1V2abs +
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

surviveyr2early_YN_Med_AlphaPhy_diff <- glmmTMB(surviveyr2early_YN ~ Med_AlphaPhy_diff + Disturbance + suitabilitypc1V2abs + 
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

#SECOND GROWING SEASON

haueserdatayr2late <- haueserdata %>% filter(germ_YN>0 | surviveyr1_YN>0 | surviveyr2early_YN > 0 | surviveyr2late_YN>0)

#Raw Mean AlphaPhy 

surviveyr2late_YN_Mean_AlphaPhy <- glmmTMB(surviveyr2late_YN ~ Mean_AlphaPhy*Disturbance + suitabilitypc1V2abs +
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

#Segue: checking climate suitability

surviveyr2late_YN_Mean_AlphaPhy_diffsq <- glmmTMB(surviveyr2late_YN ~ Mean_AlphaPhy_diff + I(suitabilitypc2V2^2) + Disturbance + 
                                                    suitabilitypc1V2 + Heating + scale(awmpd) + OptGermRt +
                                                    (1|Species) + (1|Plot), 
                                                  family=binomial(link="logit"), 
                                                  data=haueserdatayr2late)
summary(surviveyr2late_YN_Mean_AlphaPhy_diffsq)

AICtab(surviveyr2late_YN_Mean_AlphaPhy_diff, surviveyr2late_YN_Mean_AlphaPhy_diffsq)

Figure_S13b <- visreg(surviveyr2late_YN_Mean_AlphaPhy_diffsq, 
                         "suitabilitypc2V2", gg=TRUE, xlab="Untransformed precipitation dissimilarity",
                         ylab="Scaled likelihood of survival to end of second growing season", line=list(col="black"),
                         fill=list(fill=rgb(red=0.8, green=0.8, blue=0.8, alpha=0.6)),
                         points=list(size=2, pch=16, col="black")) + theme_bw() +
  theme(axis.text=element_text(size=7),
        axis.title=element_text(size=9), 
        legend.text=element_text(size=7), 
        legend.title=element_text(size=9))
Figure_S13b

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

##

#Raw GammaPhy 

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

#Difference in GammaPhy #no interactions

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

#Raw Median AlphaPhy #no interactions 

surviveyr2late_YN_Med_AlphaPhy <- glmmTMB(surviveyr2late_YN ~ Med_AlphaPhy + Disturbance + suitabilitypc1V2abs +
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

surviveyr2late_YN_Med_AlphaPhy_diff <- glmmTMB(surviveyr2late_YN ~ Med_AlphaPhy_diff +  Disturbance + suitabilitypc1V2abs + 
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

##

#Raw Max AlphaPhy 

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

##

#Raw GammaPhy  

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

###################################################################################

#COUNT RESPONSE VARIABLE FOR ALL DATASETS (number of plants in plot that survived)

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

##

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

#Mean AlphaPhy difference #interaction with precipitation

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

##

#Raw GammaPhy #no interactions

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

#Med AlphaPhy difference #marginal interaction with temp

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

#FIRST GROWING SEASON

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

##

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

#Max AlphaPhy difference #no interaction

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

##

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

#GammaPhy difference #no interactions

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

#OVERWINTER 

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

#Mean AlphaPhy difference #no interaction

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

##

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

#Max AlphaPhy difference #no interaction

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

##

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

#SECOND GROWING SEASON

haueseryr2latecount <- filter(haueserdata, surviveyr2late_num>0)

#Raw Mean AlphaPhy 

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

##

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

##

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

#ALL EFFECT SIZES - MEAN ALPHA

hauesereffectsizes <- read.csv("haueserresults_effectsizev2.csv", header=T, stringsAsFactors = T)
hauesereffectsizes$Stage = factor(hauesereffectsizes$Stage, levels = c("Second growing season","Overwinter","First growing season","Colonisation"))
hauesereffectsizes$Variable = factor(hauesereffectsizes$Variable, levels = c("Alpha PD diff. x Climate (precip.)", "Alpha PD diff. x Undisturbed", "Heating", "Phylogenetic distance", "Climate dissimilarity (temp.)",
                                                                             "Climate dissimilarity (precip.)", "Undisturbed", "Mean Alpha PD difference"))

Figure_3 <- ggplot(hauesereffectsizes, aes(x=effect, y=Variable, color=Stage, shape=significant)) +
  geom_point(size=2.3, position=position_dodge(0.5)) +
  scale_shape_manual(values=c(2, 1, 16)) +
  geom_errorbar(data=hauesereffectsizes, aes(y=Variable, xmin=effect_low, xmax=effect_high), 
                width=0, size=1, position=position_dodge(0.5)) +
  facet_grid(. ~ Metric, scales="free") +
  scale_x_continuous(trans='pseudo_log') +
  scale_color_viridis(discrete = TRUE, guide = guide_legend(reverse = TRUE)) +
  theme_bw() + geom_vline(xintercept=0, linetype="solid", color="grey", size=0.9) +
  xlab("Effect size ± 95% C.I.") + guides(shape = "none") + theme(legend.position = "top") +
  theme(axis.text=element_text(size=7),
        axis.title=element_text(size=9), 
        legend.text=element_text(size=7), 
        legend.title=element_text(size=9))
Figure_3

ggsave(Figure_3, 
       filename = "Figure_3.svg",
       height = 140, width = 180, units = "mm")

#MEAN ALPHA PHY - INTERACTIONS

yr2earlymeantempint <- interact_plot(surviveyr2early_YN_Mean_AlphaPhy_diff, pred="Mean_AlphaPhy_diff", modx = "suitabilitypc2V2abs",
                                     plot.points=TRUE, partial.residuals = TRUE, legend.main = "Climate dissimilarity (precip.)") +
  xlab("Difference in Mean Alpha Phylogenetic Diversity") + ylab("Likelihood of any plant in plot surviving overwinter") +
  scale_linetype_manual(values=c("solid", "solid", "solid")) +
  scale_y_continuous(trans='log10') +
  theme_bw() + guides(linetype = "none") + theme(legend.position = "top")  +
  theme(axis.text=element_text(size=7),
        axis.title=element_text(size=9), 
        legend.text=element_text(size=7), 
        legend.title=element_text(size=9))
yr2earlymeantempint

yr2YNmeaninteraction <- interact_plot(surviveyr2late_YN_Mean_AlphaPhy_diff, pred="Mean_AlphaPhy_diff", modx = "Disturbance",
                                      plot.points=TRUE, partial.residuals = TRUE) +
  xlab("Difference in Mean Alpha Phylogenetic Diversity") + ylab("Likelihood of any plant in plot surviving to end of second year") +
  scale_color_manual(values=c("#999999", "#E69F00")) +
  scale_linetype_manual(values=c("solid", "solid")) +
  theme_bw() + guides(linetype = "none") + theme(legend.position = "top")  +
  theme(axis.text=element_text(size=7),
        axis.title=element_text(size=9), 
        legend.text=element_text(size=7), 
        legend.title=element_text(size=9))
yr2YNmeaninteraction

Figure_4 <- ggarrange(yr2earlymeantempint, yr2YNmeaninteraction,
                                   nrow=1, ncol=2)
Figure_4
ggsave(Figure_4, 
       filename = "Figure_4.svg",
       height = 105, width = 210, units = "mm")

#MAX ALPHA PHY - INTERACTIONS

germYNmaxinteraction <- interact_plot(germ_YN_Max_AlphaPhy_diff, pred="Max_AlphaPhy_diff", modx = "Disturbance",
                                      plot.points=TRUE, partial.residuals = TRUE) +
  xlab("Difference in Maximum Alpha Phylogenetic Diversity") + ylab("Likelihood of any plant in plot colonising") +
  scale_color_manual(values=c("#999999", "#E69F00")) +
  scale_linetype_manual(values=c("solid", "solid")) +
  theme_bw() + guides(linetype = "none") + theme(legend.position = "top")  +
  theme(axis.text=element_text(size=7),
        axis.title=element_text(size=9), 
        legend.text=element_text(size=7), 
        legend.title=element_text(size=9))
germYNmaxinteraction

germcountmaxinteraction <- interact_plot(germ_count_Max_AlphaPhy_diff, pred="Max_AlphaPhy_diff", modx = "suitabilitypc2V2abs",
                                         plot.points=TRUE, partial.residuals = TRUE, legend.main = "Climate dissimilarity (precip.)") +
  xlab("Difference in Maximum Alpha Phylogenetic Diversity") + ylab("Number of plants colonising") +
  scale_linetype_manual(values=c("solid", "solid", "solid")) +
  scale_y_continuous(trans='log10') +
  theme_bw() + guides(linetype = "none") + theme(legend.position = "top")  +
  theme(axis.text=element_text(size=7),
        axis.title=element_text(size=9), 
        legend.text=element_text(size=7), 
        legend.title=element_text(size=9))
germcountmaxinteraction

yr1YNmaxdistinteraction <- interact_plot(surviveyr1_YN_Max_AlphaPhy_diff, pred="Max_AlphaPhy_diff", modx = "Disturbance",
                                         plot.points=TRUE, partial.residuals = TRUE) +
  xlab("Difference in Maximum Alpha Phylogenetic Diversity") + ylab("Likelihood of any plant in plot surviving to end of first year") +
  scale_color_manual(values=c("#999999", "#E69F00")) +
  scale_linetype_manual(values=c("solid", "solid")) +
  theme_bw() + guides(linetype = "none") + theme(legend.position = "top")  +
  theme(axis.text=element_text(size=7),
        axis.title=element_text(size=9), 
        legend.text=element_text(size=7), 
        legend.title=element_text(size=9))
yr1YNmaxdistinteraction

yr1countmaxinteraction <- interact_plot(yr1_count_Max_AlphaPhy_diff, pred="Max_AlphaPhy_diff", modx = "suitabilitypc2V2abs",
                                        plot.points=TRUE, partial.residuals = TRUE, legend.main = "Climate dissimilarity (precip.)") +
  xlab("Difference in Maximum Alpha Phylogenetic Diversity") + ylab("Number of plants surviving to end of first year") +
  scale_linetype_manual(values=c("solid", "solid", "solid")) +
  scale_y_continuous(trans='log10') +
  theme_bw() + guides(linetype = "none") + theme(legend.position = "top")  +
  theme(axis.text=element_text(size=7),
        axis.title=element_text(size=9), 
        legend.text=element_text(size=7), 
        legend.title=element_text(size=9))
yr1countmaxinteraction

#Group plot for Max AlphaPhy 

Figure_S10 <- ggarrange(germYNmaxinteraction, germcountmaxinteraction,
                                yr1YNmaxdistinteraction, yr1countmaxinteraction, nrow=2, ncol=2)
Figure_S10
ggsave(Figure_S10, 
       filename = "maxalphahaueserfigv2.svg",
       height = 210, width = 210, units = "mm")

#GAMMA PHY - INTERACTIONS

germYNgammadistinteraction <- interact_plot(germ_YN_GammaPhy_diff, pred="GammaPhy_diff", modx = "Disturbance",
                                            plot.points=TRUE, partial.residuals = TRUE) +
  xlab("Difference in Gamma Phylogenetic Diversity") + ylab("Likelihood of any plant in plot colonising") +
  scale_color_manual(values=c("#999999", "#E69F00")) +
  scale_linetype_manual(values=c("solid", "solid")) +
  theme_bw() + guides(linetype = "none") + theme(legend.position = "top")  +
  theme(axis.text=element_text(size=7),
        axis.title=element_text(size=9), 
        legend.text=element_text(size=7), 
        legend.title=element_text(size=9))
germYNgammadistinteraction

overwintercountgammainteraction <- interact_plot(surviveyr2early_YN_GammaPhy_diff, pred="GammaPhy_diff", modx = "suitabilitypc2V2abs",
                                                 plot.points=TRUE, partial.residuals = TRUE, legend.main = "Climate dissimilarity (precip.)") +
  xlab("Difference in Gamma Phylogenetic Diversity") + ylab("Likelihood of any plant in plot surviving overwinter") +
  scale_linetype_manual(values=c("solid", "solid", "solid")) +
  scale_y_continuous(trans='log10') +
  theme_bw() + guides(linetype = "none") + theme(legend.position = "top")  +
  theme(axis.text=element_text(size=7),
        axis.title=element_text(size=9), 
        legend.text=element_text(size=7), 
        legend.title=element_text(size=9))
overwintercountgammainteraction

#Group plot for GammaPhy 

Figure_S12 <- ggarrange(germYNgammadistinteraction, overwintercountgammainteraction, nrow=1, ncol=2)
Figure_S12
ggsave(Figure_S12, 
       filename = "Figure_S12.svg",
       height = 105, width = 210, units = "mm")

#MEDIAN ALPHAPHY - INTERACTIONS

medYNdistinteraction <- interact_plot(germ_YN_Med_AlphaPhy_diff, pred="Med_AlphaPhy_diff", modx = "Disturbance",
                                      plot.points=TRUE, partial.residuals = TRUE) +
  xlab("Difference in Median Phylogenetic Diversity") + ylab("Likelihood of any plant in plot colonising") +
  scale_color_manual(values=c("#999999", "#E69F00")) +
  scale_linetype_manual(values=c("solid", "solid")) +
  theme_bw() + guides(linetype = "none") + theme(legend.position = "top")  +
  theme(axis.text=element_text(size=7),
        axis.title=element_text(size=9), 
        legend.text=element_text(size=7), 
        legend.title=element_text(size=9))
medYNdistinteraction

medyr1YNdistinteraction <- interact_plot(surviveyr1_YN_Med_AlphaPhy_diff, pred="Med_AlphaPhy_diff", modx = "Disturbance",
                                         plot.points=TRUE, partial.residuals = TRUE) +
  xlab("Difference in Median Phylogenetic Diversity") + ylab("Likelihood of any plant in plot surviving to end of first year") +
  scale_color_manual(values=c("#999999", "#E69F00")) +
  scale_linetype_manual(values=c("solid", "solid")) +
  theme_bw() + guides(linetype = "none") + theme(legend.position = "top")  +
  theme(axis.text=element_text(size=7),
        axis.title=element_text(size=9), 
        legend.text=element_text(size=7), 
        legend.title=element_text(size=9))
medyr1YNdistinteraction

yr2latecountmedtempinteraction <- interact_plot(yr2late_count_Med_AlphaPhy_diff, pred="Med_AlphaPhy_diff", modx = "suitabilitypc2V2abs",
                                                plot.points=TRUE, partial.residuals = TRUE, legend.main = "Climate dissimilarity (precip.)") +
  xlab("Difference in Median Phylogenetic Diversity") + ylab("Number of plants surviving to end of second year") +
  scale_y_continuous(trans='log10') +
  scale_linetype_manual(values=c("solid", "solid", "solid")) +
  theme_bw() + guides(linetype = "none") + theme(legend.position = "top")  +
  theme(axis.text=element_text(size=7),
        axis.title=element_text(size=9), 
        legend.text=element_text(size=7), 
        legend.title=element_text(size=9))
yr2latecountmedtempinteraction

Figure_S11 <- ggarrange(medYNdistinteraction, medyr1YNdistinteraction, 
                              yr2latecountmedtempinteraction, nrow=2, ncol=2)
Figure_S11
ggsave(Figure_S11, 
       filename = "Figure_S11.svg",
       height = 210, width = 210, units = "mm")

##################################################################################################

#QUESTION Three: Does high home-range-PD act on specific traits to improve plant performance?

##################################################################################################

#Read in data from Kempel et al. 2013 and test the traits just available in Kempel

kempeltraits <- read.csv("Kempel.csv", header=T, stringsAsFactors = T) %>%
  select(plant_name_id, seed1000, Competitionindex_lnR, Germ.green.mean, herb.all) %>%
  mutate(Seed_wt = seed1000/1000) %>%
  unique() %>%
  rename(OptGermRt=Germ.green.mean)
kempeltraits$Competitionindex_lnR  <- as.numeric(as.character(kempeltraits$Competitionindex_lnR))
kempeltraits$herb.all  <- as.numeric(as.character(kempeltraits$herb.all))
kempeltraits$experiment <- rep("kempel", times=93)

climatekempel <- readRDS("climatedfv2.RDS") %>% filter(experiment_country=="switzerland")

kempeltraits <- left_join(kempeltraits, climatekempel, by="plant_name_id")

#Test relationship between resistance to herbivory and phylogenetic diversity

summary(lm(herb.all ~ Mean_AlphaPhy, data=kempeltraits))

summary(lm(herb.all ~ Max_AlphaPhy, data=kempeltraits))

summary(lm(herb.all ~ Med_AlphaPhy, data=kempeltraits))

summary(lm(herb.all ~ log(GammaPhy), data=kempeltraits))

#Test relationship between competition index and phylogenetic diversity

summary(lm(Competitionindex_lnR ~ Mean_AlphaPhy, data=kempeltraits))

summary(lm(Competitionindex_lnR ~ Max_AlphaPhy, data=kempeltraits))

summary(lm(Competitionindex_lnR ~ Med_AlphaPhy, data=kempeltraits))

summary(lm(Competitionindex_lnR ~ log(GammaPhy), data=kempeltraits))

#Link with data from Haueser et al. 2017 data and analyse traits available in both

hauesertraits <- readRDS("haueserdata.RDS")

#First just do Haueser traits

haueseranalysis <- hauesertraits %>%
  select(plant_name_id, Seed_wt, Status, OptGermRt, Max_AlphaPhy, GammaPhy, Mean_AlphaPhy, Med_AlphaPhy, HZ, Life_span) %>%
  unique()

#Hardiness index 

mod1 <- lm(Mean_AlphaPhy ~ HZ, data=haueseranalysis)
Anova(mod1, type="III")

mod2 <- lm(Max_AlphaPhy ~ HZ, data=haueseranalysis)
Anova(mod2, type="III")

mod3 <- lm(Med_AlphaPhy ~ HZ, data=haueseranalysis)
Anova(mod3, type="III")

mod4 <- lm(GammaPhy ~ HZ, data=haueseranalysis)
Anova(mod4, type="III")

#Annual or perennial

mod5 <- lm(Mean_AlphaPhy ~ Life_span, data=haueseranalysis)
Anova(mod5, type="III")

mod6 <- lm(Max_AlphaPhy ~ Life_span, data=haueseranalysis)
Anova(mod6, type="III")

mod7 <- lm(Med_AlphaPhy ~ Life_span, data=haueseranalysis)
Anova(mod7, type="III")

mod8 <- lm(GammaPhy ~ Life_span, data=haueseranalysis)
Anova(mod8, type="III")

#Now join Kempel with Haueser

haueserjoining <- hauesertraits %>%
  select(plant_name_id, Seed_wt, Status, OptGermRt, Max_AlphaPhy, GammaPhy, Mean_AlphaPhy, Med_AlphaPhy) %>%
  unique()
kempeljoining <- kempeltraits %>% 
  select(plant_name_id, Seed_wt, Status, OptGermRt, Max_AlphaPhy, GammaPhy, Mean_AlphaPhy, Med_AlphaPhy)

grouptraits <- bind_rows(haueserjoining, kempeljoining)

#Test if phylogenetic diversity linked to seed weight

summary(lm(log(Seed_wt) ~ Mean_AlphaPhy, data=grouptraits))

meanseedweight <- ggplot(grouptraits, aes(x=Mean_AlphaPhy, y=Seed_wt)) +
  geom_point() + geom_smooth(method="lm", colour="black") +
  xlab("Mean Alpha Phylogenetic Diversity") + ylab("Seed weight (mg)") +
  theme_bw() + scale_y_continuous(trans='log10') 
meanseedweight

summary(lm(log(Seed_wt) ~ Max_AlphaPhy, data=grouptraits))

summary(lm(log(Seed_wt) ~ Med_AlphaPhy, data=grouptraits))

summary(lm(log(Seed_wt) ~ log(GammaPhy), data=grouptraits))

#Test if optimal germination rate linked to phylogenetic distance

summary(glm(OptGermRt ~ Mean_AlphaPhy, data=grouptraits, family=binomial(link="logit")))

summary(glm(OptGermRt ~ Max_AlphaPhy, data=grouptraits, family=binomial(link="logit")))

summary(glm(OptGermRt ~ Med_AlphaPhy, data=grouptraits, family=binomial(link="logit")))

summary(glm(OptGermRt ~ log(GammaPhy), data=grouptraits, family=binomial(link="logit")))

#Doing height using TRY
#The following code shows the processing from the raw TRY data - the processed file is provided
#and the analysis continues at line 3209

##heights <- rtry_import("42106.txt") %>%
##  filter(TraitID=="3106"|TraitID=="3107")

##trynames <- unique(heights$AccSpeciesName)
##naughtynames <- unique(heights$SpeciesName)

##experimentdata <- read.csv("allspeciesall3exp_clean.csv") %>%
##  select(POWO.name, Status, Max_AlphaPhy, GammaPhy, Med_AlphaPhy, Mean_AlphaPhy) %>%
##  unique()

##experimentnames <- unique(experimentdata$POWO.name)

##specieswithdata <- experimentnames[experimentnames %in% trynames]
##specieswithoutdata <- experimentnames[!experimentnames %in% trynames]

#check if the ones missing are in the TRY original names

##specieswithoutdata[specieswithoutdata %in% naughtynames]

##somespecies <- heights %>%
##  filter(SpeciesName %in% specieswithoutdata)

##species1 <- filter(heights, AccSpeciesName=="Pseudofumaria lutea")
##species2 <- filter(heights, AccSpeciesName=="Digitalis lanata")
##species3 <- filter(heights, AccSpeciesName=="Jacobaea maritima")
##species6 <- filter(heights, AccSpeciesName=="Lychnis chalcedonica")

##mostspecies <- heights %>%
##  filter(AccSpeciesName %in% specieswithdata)

##fullheightdata <- bind_rows(somespecies, mostspecies, species1, species2, 
##                            species3, species6)

#Rename the necessary species in experimentdata to match AccSpeciesName

##experimentdata$POWO.name <- fct_recode(experimentdata$POWO.name,
##                                       "Erigeron andicola" = "Aster bellidiastrum", 
##                                       "ATROPA BELLA-DONNA" = "Atropa bella-donna",
##                                       "Galeopsis ladanum" = "Galeopsis angustifolia", 
##                                       "LYCHNIS FLOS-CUCULI" = "Silene flos-cuculi", 
##                                       "Pseudofumaria lutea" = "Pseudo-fumaria lutea",
##                                       "Digitalis lanata"="Digitalis lanata subsp. trojana",
##                                       "Jacobaea maritima"="Jacobaea maritima subsp. bicolor",
##                                       "Lychnis chalcedonica"="Silene chalcedonica")

#Calculate average height data for each species

#First do vegetative height (trait 3106)

##vegheight <- fullheightdata %>%
##  filter(TraitID==3106) %>% 
##  group_by(AccSpeciesName) %>%
##  summarise(vegheight=mean(StdValue, na.rm=TRUE)) %>%
##  filter(AccSpeciesName!="Lithotoma axillaris") %>%
##  rename(POWO.name = AccSpeciesName)
#no data for Lithotoma axillaris

#Now do for generative height (trait 3107)

##genheight <- fullheightdata %>%
##  filter(TraitID==3107) %>% 
##  group_by(AccSpeciesName) %>%
##  summarise(genheight=mean(StdValue, na.rm=TRUE)) %>%
##  rename(POWO.name = AccSpeciesName)

##combinedheights <- left_join(experimentdata, vegheight, by="POWO.name")
##combinedheights <- left_join(combinedheights, genheight, by="POWO.name")

##write.csv(combinedheights, "heighttraits.csv")

combinedheights <- read.csv("heighttraits.csv", header=T)

#Test if phylogenetic diversity linked to vegetative height

summary(lm(log(vegheight) ~ Mean_AlphaPhy, data=combinedheights))

summary(lm(log(vegheight) ~ Max_AlphaPhy, data=combinedheights))

summary(lm(log(vegheight) ~ Med_AlphaPhy, data=combinedheights))

summary(lm(log(vegheight) ~ log(GammaPhy), data=combinedheights))

#Test if phylogenetic diversity linked to generative height

summary(lm(log(genheight) ~ Mean_AlphaPhy, data=combinedheights))

summary(lm(log(genheight) ~ Max_AlphaPhy, data=combinedheights))

summary(lm(log(genheight) ~ Med_AlphaPhy, data=combinedheights))

summary(lm(log(genheight) ~ log(GammaPhy), data=combinedheights))

#Now do other TRY traits - leaf N (14), LDMC (47) and SLA (3117)
#The following code shows the processing from the raw TRY data - the processed file is provided
#and the analysis continues at line 3312

##manytraits <- rtry_import("42495.txt") %>%
##  filter(TraitID=="14"|TraitID=="47"|TraitID=="3117")

##trynames <- unique(manytraits$AccSpeciesName)
##naughtynames <- unique(manytraits$SpeciesName)

##experimentdata <- read.csv("allspeciesall3exp_clean.csv") %>%
##  select(POWO.name, Status, Max_AlphaPhy, GammaPhy, Med_AlphaPhy, Mean_AlphaPhy) %>%
##  unique()

##experimentnames <- unique(experimentdata$POWO.name)

##specieswithdata <- experimentnames[experimentnames %in% trynames]
##specieswithoutdata <- experimentnames[!experimentnames %in% trynames]

#check if the ones missing are in the TRY original names

##specieswithoutdata[specieswithoutdata %in% naughtynames]

##somespecies <- manytraits %>%
##  filter(SpeciesName %in% specieswithoutdata)

##species1 <- filter(manytraits, AccSpeciesName=="Pseudofumaria lutea")
##species2 <- filter(manytraits, AccSpeciesName=="Digitalis lanata")
##species3 <- filter(manytraits, AccSpeciesName=="Jacobaea maritima")
##species4 <- filter(manytraits, AccSpeciesName=="Lychnis chalcedonica")

##mostspecies <- manytraits %>%
##  filter(AccSpeciesName %in% specieswithdata)

##fulltraitdata <- bind_rows(somespecies, mostspecies, species1, species2, 
##                           species3, species4)

#Rename the necessary species in experimentdata to match AccSpeciesName

##experimentdata$POWO.name <- fct_recode(experimentdata$POWO.name,
##                                       "Erigeron andicola" = "Aster bellidiastrum", 
##                                       "ATROPA BELLA-DONNA" = "Atropa bella-donna",
##                                       "Galeopsis ladanum" = "Galeopsis angustifolia", 
##                                       "LYCHNIS FLOS-CUCULI" = "Silene flos-cuculi", 
##                                       "Pseudofumaria lutea" = "Pseudo-fumaria lutea",
##                                       "Digitalis lanata"="Digitalis lanata subsp. trojana",
##                                       "Jacobaea maritima"="Jacobaea maritima subsp. bicolor",
##                                       "Lychnis chalcedonica"="Silene chalcedonica")

#Calculate average trait data for each species

#First do leaf N (trait 14)

##leafN <- fulltraitdata %>%
##  filter(TraitID==14) %>% 
##  group_by(AccSpeciesName) %>%
##  summarise(leafN=mean(StdValue, na.rm=TRUE)) %>%
##  rename(POWO.name = AccSpeciesName)

#Now do for LDMC (trait 47)

##LDMC <- fulltraitdata %>%
##  filter(TraitID==47) %>% 
##  group_by(AccSpeciesName) %>%
##  summarise(LDMC=mean(StdValue, na.rm=TRUE)) %>%
##  rename(POWO.name = AccSpeciesName)

#Now do for SLA (trait 3117)

##SLA <- fulltraitdata %>%
##  filter(TraitID==3117) %>% 
##  group_by(AccSpeciesName) %>%
##  summarise(SLA=mean(StdValue, na.rm=TRUE)) %>%
##  rename(POWO.name = AccSpeciesName)

##combinedtraits <- left_join(experimentdata, leafN, by="POWO.name")
##combinedtraits <- left_join(combinedtraits, LDMC, by="POWO.name")
##combinedtraits <- left_join(combinedtraits, SLA, by="POWO.name")

##write.csv(combinedtraits, "leaftraits.csv")

combinedtraits <- read.csv("leaftraits.csv", header=T)

#Test if phylogenetic diversity linked to leaf N

meanleafN <- ggplot(combinedtraits, aes(x=Mean_AlphaPhy, y=leafN)) +
  geom_point() + geom_smooth(method="lm", colour="black") +
  xlab("Mean Alpha Phylogenetic Diversity") + ylab("Leaf N (mg/g)") +
  theme_bw() 
meanleafN

summary(lm(leafN ~ Mean_AlphaPhy, data=combinedtraits))

summary(lm(leafN ~ Max_AlphaPhy, data=combinedtraits))

summary(lm(leafN ~ Med_AlphaPhy, data=combinedtraits))

summary(lm(leafN ~ log(GammaPhy), data=combinedtraits))

#Test if phylogenetic diversity linked to LDMC

meanLDMC <- ggplot(combinedtraits, aes(x=Mean_AlphaPhy, y=LDMC)) +
  geom_point() + geom_smooth(method="lm", colour="black") +
  xlab("Mean Alpha Phylogenetic Diversity") + ylab("LDMC (g/g)") +
  theme_bw() 
meanLDMC

summary(lm(LDMC ~ Mean_AlphaPhy, data=combinedtraits))

summary(lm(LDMC ~ Max_AlphaPhy, data=combinedtraits))

summary(lm(LDMC ~ Med_AlphaPhy, data=combinedtraits))

summary(lm(LDMC ~ log(GammaPhy), data=combinedtraits))

#Test if phylogenetic diversity linked to SLA

summary(lm(SLA ~ Mean_AlphaPhy, data=combinedtraits))

summary(lm(SLA ~ Max_AlphaPhy, data=combinedtraits))

summary(lm(SLA ~ Med_AlphaPhy, data=combinedtraits))

summary(lm(SLA ~ log(GammaPhy), data=combinedtraits))

Figure_5 <- ggarrange(meanseedweight, meanLDMC, meanleafN,
                      nrow=2, ncol=2)
Figure_5
ggsave(Figure_5, 
       filename = "Figure_5.svg",
       height = 210, width = 210, units = "mm")

##################################################################################################

#SUPPLEMENTARY ANALYSES

#As per main text results, final selected models are presented

##################################################################################################

#Re-do Question 1 models WITHOUT HAUESER to see if trends still hold
#K+M models in Tables S1-S4

germinationv2 <- readRDS("germination.RDS") %>%
  mutate(suitabilitypc1V2abs = abs(suitabilitypc1V2)) %>%
  mutate(suitabilitypc2V2abs = abs(suitabilitypc2V2)) %>%
  filter(Study != "Haueser")

surviveyr1v2 <- readRDS("surviveyr1.RDS") %>%
  mutate(suitabilitypc1V2abs = abs(suitabilitypc1V2)) %>%
  mutate(suitabilitypc2V2abs = abs(suitabilitypc2V2))  %>%
  filter(Study != "Haueser")

surviveyr2earlyv2 <- readRDS("surviveyr2early.RDS") %>%
  mutate(suitabilitypc1V2abs = abs(suitabilitypc1V2)) %>%
  mutate(suitabilitypc2V2abs = abs(suitabilitypc2V2))  %>%
  filter(Study != "Haueser")

surviveyr2latev2 <- readRDS("surviveyr2late.RDS") %>%
  mutate(suitabilitypc1V2abs = abs(suitabilitypc1V2)) %>%
  mutate(suitabilitypc2V2abs = abs(suitabilitypc2V2))  %>%
  filter(Study != "Haueser")

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

surviveyr2early_YN_Max_AlphaPhy <- glmmTMB(surviveyr2early_YN ~ Max_AlphaPhy + Disturbance + propagule_pressure +
                                             suitabilitypc1V2abs + suitabilitypc2V2abs + Herbivory +
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

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2early_YN_GammaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr2early_YN_GammaPhy)
check_singularity(surviveyr2early_YN_GammaPhy)

#Median alpha diversity 

surviveyr2early_YN_MedPhy <- glmmTMB(surviveyr2early_YN ~ Med_AlphaPhy + Disturbance + propagule_pressure +
                                       suitabilitypc1V2abs + suitabilitypc2V2abs + Herbivory +
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

surviveyr2late_YN_Max_AlphaPhy <- glmmTMB(surviveyr2late_YN ~ Max_AlphaPhy*Disturbance + propagule_pressure +
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

#Now doing it agaib looking at difference at a site level for
#Haueser and Kempel combined
#The K+H models in Tables S1-S4

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

germinationNM <- germination %>%
  filter(Study != "Muller") %>%
  droplevels()
germinationNM <- left_join(germinationNM, comm_diversity_total, by="Site") %>%
  mutate(Max_AlphaPhy_diff = Max_AlphaPhy - Max_AlphaPhy_mean,
         Mean_AlphaPhy_diff = Mean_AlphaPhy - Mean_AlphaPhy_mean,
         GammaPhy_diff = GammaPhy - GammaPhy_mean,
         Med_AlphaPhy_diff = Med_AlphaPhy - Med_AlphaPhy_mean) %>%
  mutate(suitabilitypc1V2abs = abs(suitabilitypc1V2)) %>%
  mutate(suitabilitypc2V2abs = abs(suitabilitypc2V2))

surviveyr1NM <- surviveyr1  %>%
  filter(Study != "Muller") %>%
  droplevels()
surviveyr1NM <- left_join(surviveyr1NM, comm_diversity_total, by="Site") %>%
  mutate(Max_AlphaPhy_diff = Max_AlphaPhy - Max_AlphaPhy_mean,
         Mean_AlphaPhy_diff = Mean_AlphaPhy - Mean_AlphaPhy_mean,
         GammaPhy_diff = GammaPhy - GammaPhy_mean,
         Med_AlphaPhy_diff = Med_AlphaPhy - Med_AlphaPhy_mean) %>%
  mutate(suitabilitypc1V2abs = abs(suitabilitypc1V2)) %>%
  mutate(suitabilitypc2V2abs = abs(suitabilitypc2V2))

surviveyr2earlyNM <- surviveyr2early  %>%
  filter(Study != "Muller") %>%
  droplevels()
surviveyr2earlyNM <- left_join(surviveyr2earlyNM, comm_diversity_total, by="Site") %>%
  mutate(Max_AlphaPhy_diff = Max_AlphaPhy - Max_AlphaPhy_mean,
         Mean_AlphaPhy_diff = Mean_AlphaPhy - Mean_AlphaPhy_mean,
         GammaPhy_diff = GammaPhy - GammaPhy_mean,
         Med_AlphaPhy_diff = Med_AlphaPhy - Med_AlphaPhy_mean) %>%
  mutate(suitabilitypc1V2abs = abs(suitabilitypc1V2)) %>%
  mutate(suitabilitypc2V2abs = abs(suitabilitypc2V2))

#surviveyr2late is already just Kempel and Haueser
surviveyr2lateNM <- left_join(surviveyr2late, comm_diversity_total, by="Site") %>%
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
                                 data=germinationNM)
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
                                data=germinationNM)
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
                                data=germinationNM)
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
                            data=germinationNM)
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
                                       data=surviveyr1NM)
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
                                      data=surviveyr1NM)
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
                                      data=surviveyr1NM)
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
                                  data=surviveyr1NM)
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

surviveyr2early_YN_Mean_AlphaPhy <- glmmTMB(surviveyr2early_YN ~ Mean_AlphaPhy_diff + Disturbance + propagule_pressure +
                                              suitabilitypc1V2abs + suitabilitypc2V2abs +
                                              (1|Family) + (1|POWO.name) + 
                                              (1|Site) + offset(log(Density)), 
                                            family=binomial(link="logit"), 
                                            data=surviveyr2earlyNM)
summary(surviveyr2early_YN_Mean_AlphaPhy)
Anova(surviveyr2early_YN_Mean_AlphaPhy, type="III")
confint(surviveyr2early_YN_Mean_AlphaPhy)
r.squaredGLMM(surviveyr2early_YN_Mean_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2early_YN_Mean_AlphaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr2early_YN_Mean_AlphaPhy)
check_singularity(surviveyr2early_YN_Mean_AlphaPhy)

#Median alpha diversity (#checked and no interaction between diversity and disturbance or climate suitability)

surviveyr2early_YN_Med_AlphaPhy <- glmmTMB(surviveyr2early_YN ~ Med_AlphaPhy_diff + Disturbance + propagule_pressure +
                                             suitabilitypc1V2abs + suitabilitypc2V2abs +
                                             (1|Family) + (1|POWO.name) + 
                                             (1|Site) + offset(log(Density)), 
                                           family=binomial(link="logit"), 
                                           data=surviveyr2earlyNM)
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

surviveyr2early_YN_Max_AlphaPhy <- glmmTMB(surviveyr2early_YN ~ Max_AlphaPhy_diff + Disturbance + propagule_pressure +
                                             suitabilitypc1V2abs + suitabilitypc2V2abs + 
                                             (1|Family) + (1|POWO.name) + 
                                             (1|Site) + offset(log(Density)), 
                                           family=binomial(link="logit"), 
                                           data=surviveyr2earlyNM)
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
                                       data=surviveyr2earlyNM)
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
                                           data=surviveyr2lateNM)
summary(surviveyr2late_YN_Mean_AlphaPhy)
Anova(surviveyr2late_YN_Mean_AlphaPhy, type="III")
confint(surviveyr2late_YN_Mean_AlphaPhy)
r.squaredGLMM(surviveyr2late_YN_Mean_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2late_YN_Mean_AlphaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr2late_YN_Mean_AlphaPhy)
check_singularity(surviveyr2late_YN_Mean_AlphaPhy)

#Median alpha diversity 

surviveyr2late_YN_Med_AlphaPhy <- glmmTMB(surviveyr2late_YN ~ Med_AlphaPhy_diff + Disturbance + propagule_pressure +
                                            suitabilitypc1V2abs + suitabilitypc2V2abs + 
                                            (1|Family) + (1|POWO.name) + 
                                            (1|Site) + offset(log(Density)), 
                                          family=binomial(link="logit"), 
                                          data=surviveyr2lateNM)
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

surviveyr2late_YN_Max_AlphaPhy <- glmmTMB(surviveyr2late_YN ~ Max_AlphaPhy_diff*Disturbance + propagule_pressure +
                                            suitabilitypc1V2abs + suitabilitypc2V2abs + 
                                            (1|Family) + (1|POWO.name) + 
                                            (1|Site) + offset(log(Density)), 
                                          family=binomial(link="logit"), 
                                          data=surviveyr2lateNM)
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

surviveyr2late_YN_GammaPhy <- glmmTMB(surviveyr2late_YN ~ scale(GammaPhy_diff)*Disturbance + propagule_pressure +
                                        suitabilitypc1V2abs + suitabilitypc2V2abs + 
                                        (1|POWO.name) + 
                                        (1|Site) + offset(log(Density)), 
                                      family=binomial(link="logit"), 
                                      data=surviveyr2lateNM)
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

germinationcountNM <- filter(germinationNM, germ_num>0)

#Mean alpha diversity

germ_count_Mean_AlphaPhy <- glmmTMB(germ_num ~ Mean_AlphaPhy_diff*Disturbance + propagule_pressure +
                                      suitabilitypc1V2abs + suitabilitypc2V2abs + Density +
                                      (1|Family) + (1|POWO.name) + (1|Site),
                                    family=nbinom2(link = "log"), 
                                    data=germinationcountNM)
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
                                   data=germinationcountNM)
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
                                   data=germinationcountNM)
summary(germ_count_Max_AlphaPhy)
Anova(germ_count_Max_AlphaPhy, type="III")
confint(germ_count_Max_AlphaPhy)
r.squaredGLMM(germ_count_Max_AlphaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = germ_count_Max_AlphaPhy, plot=F)
plot(simulationOutput)

check_overdispersion(germ_count_Max_AlphaPhy)
check_singularity(germ_count_Max_AlphaPhy)

#Gamma diversity (note no interaction unlike other two)

germ_count_GammaPhy <- glmmTMB(germ_num ~ scale(GammaPhy_diff) + Disturbance + propagule_pressure +
                                 suitabilitypc1V2abs + suitabilitypc2V2abs + Density +
                                 (1|Family) + (1|POWO.name) + (1|Site),
                               family=nbinom2(link = "log"), 
                               data=germinationcountNM)
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

surviveyr1countNM <- filter(surviveyr1NM, surviveyr1_num>0)

#Mean alpha diversity 

surviveyr1_count_Mean_AlphaPhy <- glmmTMB(surviveyr1_num ~ Mean_AlphaPhy_diff*Disturbance + propagule_pressure +
                                            suitabilitypc1V2abs + suitabilitypc2V2abs + Density +
                                            (1|POWO.name) + (1|Site),
                                          family=nbinom2(link = "log"), 
                                          data=surviveyr1countNM)
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
                                         data=surviveyr1countNM)
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
                                         data=surviveyr1countNM)
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
                                     data=surviveyr1countNM)
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

surviveyr2earlycountNM <- filter(surviveyr2earlyNM, surviveyr2early_num>0)

#Mean alpha diversity 

surviveyr2early_count_Mean_AlphaPhy <- glmmTMB(surviveyr2early_num ~ Mean_AlphaPhy_diff + Disturbance + propagule_pressure +
                                                 suitabilitypc1V2abs + suitabilitypc2V2abs + Density +
                                                 (1|POWO.name) + (1|Site),
                                               family=nbinom2(link = "log"), 
                                               data=surviveyr2earlycountNM)
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
                                              data=surviveyr2earlycountNM)
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
                                              data=surviveyr2earlycountNM)
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
                                          data=surviveyr2earlycountNM)
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

surviveyr2latecountNM <- filter(surviveyr2lateNM, surviveyr2late_num>0)

#Mean alpha diversity 

surviveyr2late_count_Mean_AlphaPhy <- glmmTMB(surviveyr2late_num ~ Mean_AlphaPhy_diff + suitabilitypc2V2abs + Disturbance +
                                                suitabilitypc1V2abs + suitabilitypc2V2abs + Density + propagule_pressure +
                                                (1|POWO.name),
                                              family=nbinom2(link = "log"), 
                                              data=surviveyr2latecountNM)
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
                                             data=surviveyr2latecountNM)
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
                                             data=surviveyr2latecountNM)
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
                                         data=surviveyr2latecountNM)
summary(surviveyr2late_count_GammaPhy)
Anova(surviveyr2late_count_GammaPhy, type="III")
confint(surviveyr2late_count_GammaPhy)
r.squaredGLMM(surviveyr2late_count_GammaPhy)

#diagnostics 
simulationOutput <- simulateResiduals(fittedModel = surviveyr2late_count_GammaPhy, plot=F)
plot(simulationOutput)
check_overdispersion(surviveyr2late_count_GammaPhy)
check_singularity(surviveyr2late_count_GammaPhy)

####################################################################################################

#DO IT AGAIN BUT WITH SINGLE CLIMATE VARIABLES 
#The 'SC' models in Tables S1-S8

####################################################################################################

#Doing it for Question 1

#Read in all four data files (output of 'data cleaning all 3 experiments.R')

germinationv2 <- readRDS("germinationv2.RDS") %>%
  mutate(suitabilitypc1V3abs = abs(suitabilitypc1V3)) %>%
  mutate(suitabilitypc2V3abs = abs(suitabilitypc2V3))
str(germinationv2)

surviveyr1v2 <- readRDS("surviveyr1v2.RDS") %>%
  mutate(suitabilitypc1V3abs = abs(suitabilitypc1V3)) %>%
  mutate(suitabilitypc2V3abs = abs(suitabilitypc2V3))
str(surviveyr1v2)

surviveyr2earlyv2 <- readRDS("surviveyr2earlyv2.RDS") %>%
  mutate(suitabilitypc1V3abs = abs(suitabilitypc1V3)) %>%
  mutate(suitabilitypc2V3abs = abs(suitabilitypc2V3))
str(surviveyr2earlyv2)

surviveyr2latev2 <- readRDS("surviveyr2latev2.RDS") %>%
  mutate(suitabilitypc1V3abs = abs(suitabilitypc1V3)) %>%
  mutate(suitabilitypc2V3abs = abs(suitabilitypc2V3))
str(surviveyr2latev2)

hist(log(germinationv2$suitabilitypc1V3abs))
hist(log(germinationv2$suitabilitypc2V3abs))

################################################################################

#NEW DATASET WITH NEW CLIMATE AND MEDIAN ALPHA

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

surviveyr2early_YN_Max_AlphaPhy <- glmmTMB(surviveyr2early_YN ~ Max_AlphaPhy + Disturbance + propagule_pressure +
                                             log(suitabilitypc1V3abs) + log(suitabilitypc2V3abs) + Herbivory +
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

haueserdataoverwinterv2 <- haueserdatav2 %>% filter(germ_YN>0 | surviveyr1_YN>0 | surviveyr2early_YN > 0)

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

#Difference in Median AlphaPhy #no interactions

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

haueserdatayr2latev2 <- haueserdatav2 %>% filter(germ_YN>0 | surviveyr1_YN>0 | surviveyr2early_YN > 0 | surviveyr2late_YN>0)

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

#Difference in GammaPhy #no interactions

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

surviveyr2late_YN_Med_AlphaPhy_diff <- glmmTMB(surviveyr2late_YN ~ Med_AlphaPhy_diff +  Disturbance + log(suitabilitypc1V3abs) + 
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

#Difference in GammaPhy #no interactions

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

#Med AlphaPhy difference 

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
