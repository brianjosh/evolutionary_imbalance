#Code for 'Plants that evolved under high phylogenetic diversity have higher invasion success, particularly in undisturbed communities'
#Joshua Brian, Mark van Kleunen, Wayne Dawson, Anne Kempel, Weihan Zhao, Jane Catford

#Code prepared by Joshua Brian; joshua.brian@kcl.ac.uk OR jshbrian@gmail.com

##PART FOUR: ANALYSIS FOR QUESTION 3 

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