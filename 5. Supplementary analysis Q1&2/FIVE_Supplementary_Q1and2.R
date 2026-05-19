#Code for 'Plants that evolved under high phylogenetic diversity have higher invasion success, particularly in undisturbed communities'
#Joshua Brian, Mark van Kleunen, Wayne Dawson, Anne Kempel, Weihan Zhao, Jane Catford

#Code prepared by Joshua Brian; joshua.brian@kcl.ac.uk OR jshbrian@gmail.com

##PART FIVE: SUPPLEMENTARY ANALYSIS FOR QUESTIONS 1 AND 2 

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

#SUPPLEMENTARY ANALYSES FOR QUESTIONS 1 AND 2

#As per main text results, final selected models are presented

##################################################################################################

#Re-do Question 1 models WITHOUT HAUESER to see if trends still hold
#K+M models in Tables S2-S5

germination <- readRDS("germinationrevised.RDS") %>%
  mutate(suitabilitypc1V2abs = abs(suitabilitypc1V2)) %>%
  mutate(suitabilitypc2V2abs = abs(suitabilitypc2V2)) 
germinationv2 <- germination %>%
  filter(Study != "Haueser")
germination$Disturbance <- relevel(germination$Disturbance, ref="0")
germinationv2$Disturbance <- relevel(germinationv2$Disturbance, ref="0")

surviveyr1 <- readRDS("surviveyr1revised.RDS") %>%
  mutate(suitabilitypc1V2abs = abs(suitabilitypc1V2)) %>%
  mutate(suitabilitypc2V2abs = abs(suitabilitypc2V2))
surviveyr1v2 <- surviveyr1 %>%
  filter(Study != "Haueser")
surviveyr1$Disturbance <- relevel(surviveyr1$Disturbance, ref="0")
surviveyr1v2$Disturbance <- relevel(surviveyr1v2$Disturbance, ref="0")

surviveyr2early <- readRDS("surviveyr2earlyrevised.RDS") %>%
  mutate(suitabilitypc1V2abs = abs(suitabilitypc1V2)) %>%
  mutate(suitabilitypc2V2abs = abs(suitabilitypc2V2))
surviveyr2earlyv2 <- surviveyr2early %>%
  filter(Study != "Haueser")
surviveyr2early$Disturbance <- relevel(surviveyr2early$Disturbance, ref="0")
surviveyr2earlyv2$Disturbance <- relevel(surviveyr2earlyv2$Disturbance, ref="0")

surviveyr2late <- readRDS("surviveyr2laterevised.RDS") %>%
  mutate(suitabilitypc1V2abs = abs(suitabilitypc1V2)) %>%
  mutate(suitabilitypc2V2abs = abs(suitabilitypc2V2)) 
surviveyr2latev2 <- surviveyr2late %>%
  filter(Study != "Haueser")
surviveyr2late$Disturbance <- relevel(surviveyr2late$Disturbance, ref="0")
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