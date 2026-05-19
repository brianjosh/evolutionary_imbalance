#Code for 'Plants that evolved under high phylogenetic diversity have higher invasion success, particularly in undisturbed communities'
#Joshua Brian, Mark van Kleunen, Wayne Dawson, Anne Kempel, Weihan Zhao, Jane Catford

#Code prepared by Joshua Brian; joshua.brian@kcl.ac.uk OR jshbrian@gmail.com

##PART THREE: ANALYSIS FOR QUESTION 2 

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
  xlab("Standardised effect size ? 95% C.I.") + guides(shape = "none") + theme(legend.position = "top") +
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