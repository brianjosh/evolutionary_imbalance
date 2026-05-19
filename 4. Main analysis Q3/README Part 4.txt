4. Main analysis - Q3

This folder describes the analyses for question 3 in the main text: Is high indigenous-PD correlated with traits that could explain increased survival of sown species?  

It contains the following files:

FOUR_Main_text_Q3.R - code file for running analyses

biogeographic_ID.csv - datafile giving the 'biogeographic syndrome' of each species (combination of continents each species is native to). See Fristoe et al. 2023 (Nature Ecol. Evol. 7(10): 1633-1644) for details

allspeciesall3exp_clean.csv - datafile listing all experimental species, their WCVP ID (plant_name_id), and their four phylogenetic indigenous range metrics ('Mean_AlphaPhy', 'Max_AlphaPhy', 'Med_AlphaPhy', 'GammaPhy')

traits_complete.csv - datafile containing trait data for all species. It has the following columns:
- 'plant_name_id' = The identifying number of each species in the World Checklist of Vascular Plants database
- 'OptGermRt' = The optimal germination rate of a species, determined in a pot experiment
- 'Seed_wt' = Weight of a single seed (mg)
- 'competition_index' = The competitive ability of a species when tested in a pairwise experiment against a heterospecific 
- 'herbivory_resistance' = Change in biomass of the generalist herbivore Spodoptera littoralis when feeding on the plant
- 'Life_span' = The life history strategy of a species (annual, perennial)
- 'hardiness' = The hardiness zone of a species, describing how tolerant it is to winter conditions (1, 2, 3, 4, 5, with 5 being the most tolerant)
- 'vegheight' = The vegetative height of a species (m)
- 'genheight' = The generative (i.e. reproductive) height of a species (m)
- 'leafN' = The proportional nitrogen in a leaf of a given species (mg g^-1)
- 'SLA' = The specific leaf area of a given species (mm^2 mg^-1)
- 'LDMC' = The dry matter content of a leaf of a given species (g g^-1)




