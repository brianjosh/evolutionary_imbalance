5. Supplementary analyses for Q1 and Q2 

This folder describes the supplementary analyses for questions 1 and 2.

It contains the following files:

FIVE_Supplementary_Q1and2.R - code file for running analyses

plot_cover.csv - datafile containing the percent cover of each resident species in all plots in the Haueser experiment

site_cover_Kempel.csv - datafile containing the presence or absence of resident vegetation at all sites in the Kempel experiment

community_diversities.csv - datafile containing the PD metrics for all resident species in the Haueser experiment

community_diversities_Kempel.csv - datafile containing the PD metrics for all resident species in the Kempel experiment

germinationrevised.RDS - datafile describing whether sown species successfully germinated (~6 weeks after sowing)

surviveyr1revised.RDS - datafile describing whether sown species successfully survived to the end of the first year

surviveyr2earlyrevised.RDS - datafile describing whether sown species successfully survived overwinter

surviveyr2earlyrevised.RDS - datafile describing whether sown species successfully survived overwinter

haueserdata.RDS - datafile describing success of sown species through time

germinationv2revised.RDS - datafile describing whether sown species successfully germinated (~6 weeks after sowing)

surviveyr1v2revised.RDS - datafile describing whether sown species successfully survived to the end of the first year

surviveyr2v2earlyrevised.RDS - datafile describing whether sown species successfully survived overwinter

surviveyr2v2earlyrevised.RDS - datafile describing whether sown species successfully survived overwinter

haueserdatav2.RDS - datafile describing success of sown species through time

These latter ten datafiles contain the full dataset from all three experiments (Kempel, Muller, Haueser). 
Note that the 'v2' files are identical in all aspects EXCEPT they contain 'suitabilitypc1v3abs' and 'suitabilitypc2v3abs', rather than the original versions of these climate dissimilarity metrics. These climate dissimilarity metrics are drawn directly from bio6 and bio17 respectively from the WorldClim database, rather than from principal component analysis of all climate variables.

The datasets collectively have the following relevant columns:
- 'plant_name_id' = The identifying number of each species in the World Checklist of Vascular Plants database
- 'Mean_AlphaPhy' = The mean alpha phylogenetic diversity of a species’ home range 
- 'Mean_AlphaPhy_diff' = The difference between the mean alpha phylogenetic diversity of a species’ home range and the abundance-weighted value for species in the resident plot
- 'Max_AlphaPhy' = The maximum alpha phylogenetic diversity of a species’ home range
- 'Max_AlphaPhy_diff' = The difference between the maximum alpha phylogenetic diversity of a species’ home range and the abundance-weighted value for species in the resident plot
- 'Med_AlphaPhy' = The median alpha phylogenetic diversity of a species’ home range
- 'Med_AlphaPhy_diff' = The difference between the median alpha phylogenetic diversity of a species’ home range and the abundance-weighted value for species in the resident plot
- 'GammaPhy' = The gamma phylogenetic diversity of a species’ home range
- 'GammaPhy_diff' = The difference between the gamma phylogenetic diversity of a species’ home range and the abundance-weighted value for species in the resident plot
- 'germ_YN' = Whether any sown plants in a plot colonised the plot successfully (1=yes, 0=no)
- 'surviveyr1_YN' = Whether any sown plants in a plot survived to the end of the first growing season of the experiment (1=yes, 0=no)
- 'surviveyr2early_YN' = Whether any sown plants in a plot survived overwinter (1=yes, 0=no)
- 'surviveyr2late_YN' = Whether any sown plants in a plot survived to the end of the second growing season of the experiment (1=yes, 0=no)
- 'flowers_YN' = Whether a plant produced any flowers (1=yes, 0=no)
- 'germ_num' = The number of sown plants that colonised a plot
- 'surviveyr1_num' = The number of sown plants that survived in a plot to the end of the first growing season
- 'surviveyr2early_num' = The number of sown plants that survived in a plot overwinter
- 'surviveyr2late_num' = The number of sown plants that survived in a plot to the end of the second growing season
- 'max_flowers' = The maximum number of flowers produced by a plant
- 'Disturbance' = Whether a plot was disturbed (1=yes, 0=no)
- 'propagule_pressure' = The seeding rate of sown species (seeds cm^-2)
- 'suitabilitypc1v2abs' = The absolute value of temperature dissimilarity between the home range and experimental climate of all sown species in the experiment (based on PC axis 1 of climate variability)
- 'suitabilitypc2v2abs' = The absolute value of precipitation dissimilarity between the home range and experimental climate of all sown species in the experiment (based on PC axis 2 of climate variability)
- 'Herbivory' = Whether or not a plot was exposed to insect herbivores (open=exposed to herbivores, closed=protected from herbivores)
- 'Status' = Whether the sown species was native or exotic to the region
- 'awMPD' = Abundance-weighted mean phylogenetic distance from the sown species to the resident species in the plot
- 'Heating' = Whether or not the plot was heated 
- 'Density' = The total number of seeds sown in a plot
- 'Family' = The family (Linnean taxonomy) a given species belongs to
- 'POWO.name' = The species name as listed in the Plants of the World database
- 'Site' = The site an experiment took place at within a given study
- 'Study' = The study from which the data was drawn (Haueser, Kempel, Muller)