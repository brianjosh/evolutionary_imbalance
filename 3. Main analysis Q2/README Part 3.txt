3. Main analysis - Q2

This folder describes the analyses for question 2 in the main text: Does high relative indigenous-PD increase the likelihood of successful colonisation and survival of sown species? 

It contains the following files:

THREE_Main_text_Q2.R - code file for running analyses

biogeographic_ID.csv - datafile giving the 'biogeographic syndrome' of each species (combination of continents each species is native to). See Fristoe et al. 2023 (Nature Ecol. Evol. 7(10): 1633-1644) for details

haueserresults_effectsize_revised.csv - cleaned results file enabling efficient production of Fig. 4 in main text

haueserdata.RDS - datafile describing success of sown species through time

This latter datafile contain the full dataset from the Haueser experiment, which recorded resident plot diversity and thus allowed us to answer Q2 (effect of the relative relationship between sown species and resident species diversity). It has the following relevant columns:
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
- 'Heating' = Whether or not the plot was heated 
- 'Status' = Whether the sown species was native or exotic to the region
- 'Density' = The total number of seeds sown in a plot
- 'awMPD' = Abundance-weighted mean phylogenetic distance from the sown species to the resident species in the plot
- 'Family' = The family (Linnean taxonomy) a given species belongs to
- 'POWO.name' = The species name as listed in the Plants of the World database