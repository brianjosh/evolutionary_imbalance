#Code for 'Plants that evolved under high phylogenetic diversity have higher invasion success, particularly in undisturbed communities'
#Joshua Brian, Mark van Kleunen, Wayne Dawson, Anne Kempel, Weihan Zhao, Jane Catford

#Code prepared by Joshua Brian; joshua.brian@kcl.ac.uk OR jshbrian@gmail.com

##PART SIX: SUPPLEMENTARY ANALYSIS FOR QUESTION 3

#This code was run using R v.4.1.1 and the following packages:
library(tidyverse) #version 1.3.2
library(glmmTMB) #version 1.1.4
library(car) #version 3.1-1
library(MuMIn) #version 1.46.0
library(bbmle) #version 1.0.25
library(lme4) #version 1.1-31
library(ggpubr) #version 0.4.0
library(cowplot) #version 1.1.1
library(piecewiseSEM) #version 2.1.2

####################################################################################################

##SUPPLEMENTARY ANALYSES FOR QUESTION 3: STRUCTURAL EQUATION MODELS

####################################################################################################

#All three studies

#Read in trait data

fulltraits <- read.csv("traits_complete.csv", header=T) 

#Read in all four data files 

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
