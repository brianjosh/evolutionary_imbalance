#Code for 'Plants that evolved under high phylogenetic diversity have higher invasion success, particularly in undisturbed communities'
#Joshua Brian, Mark van Kleunen, Wayne Dawson, Anne Kempel, Weihan Zhao, Jane Catford

#Code prepared by Joshua Brian; joshua.brian@kcl.ac.uk OR jshbrian@gmail.com

#PART ONE: CALCULATING PD VALUES AND VISUALISING GLOBAL PD

library(tidyverse) #version 1.3.2
library(funrar) #version 1.5.0
library(picante) #version 1.8.2
library(ape) #version 5.8-1
library(sf) #version 1.0-12
library(tibble) #version 3.3.1

######################################################

#Calculate maximum alpha PD and median alpha PD

######################################################

load("wcvp.v11.tree.wcvp.id.sub.RData") #full phylogenetic tree for all accepted vascular seed plants
load("SpecReg2.RData") #lists the native range TDWG3 regions for all vascular seed plants, that are also in the tree

#Calculate PD for each TDWG3 region

mat <- SpecReg2
rownames(mat) <- mat[[1]]
mat <- mat[ , -1]

matv2 <- mat %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("variable")
rownames(matv2) <- matv2[[1]]
matv2 <- matv2[ , -1]

#Calculate Faith's PD
pd_res <- pd(matv2, wcvp.v11.tree.wcvp.id.sub, include.root = FALSE) #THIS CODE TAKES A LOT OF TIME

pd_df <- data.frame(
  Code = rownames(pd_res),
  PD = pd_res$PD,
  richness = pd_res$SR
)

#Read in a file with the areas for each TDWG3 region, join with PD for each region

tdwg3_area <- read.csv("tdwg3_areas.csv", header=T, stringsAsFactors = T)
pd_tdwg3_area <- inner_join(tdwg3_area, pd_df, by="Code")

#Regress the PD of each region (y-axis) against region size (x-axis), extract residuals (to get area-corrected PD)

mod <- lm(log(PD) ~ log(Area), data = pd_tdwg3_area)
pd_tdwg3_area$ResidualPD <- residuals(mod)
pd_tdwg3_area <- pd_tdwg3_area[, c("Code", "ResidualPD")]

#For each species in the dataset, calculate maximum alpha PD and median alpha PD 

max_med_df <- data.frame(plant_name_id = as.numeric(), 
                                   Max_AlphaPD = as.numeric(), Med_AlphaPD = as.numeric())

for(i in 1:length(SpecReg2)){
  
  #Get plant name ID for each species
  
  plant_name_id <- SpecReg2[i, ]$plant_name_id
  
  #Get TDWG native range for that species
  
  regions <- filter(SpecReg2, plant_name_id==SpecReg2[i, ]$plant_name_id) %>%
    select(where(~is.numeric(.x) && any(.x != 0))) %>%
    colnames() #'regions' now has the vector of TDWG3 regions for that species
  
  #Filter the area-corrected PD values for that species
  
  pd_tdwg3_area_species <- pd_tdwg3_area %>%
    filter(Code %in% regions)
  
  #Calculate maximum alpha PD and median alpha PD
  
  Max_AlphaPD <- max(pd_tdwg3_area_species$RedisualPD)   
  Med_AlphaPD <- median(pd_tdwg3_area_species$RedisualPD)
  
  values <- data.frame(plant_name_id, Max_AlphaPD, Med_AlphaPD)
  
  max_med_df <- bind_rows(max_med_df, values)
  
}

######################################################

#Calculate mean alpha PD and Gamma PD

######################################################

#Function to compute PD from a species set

calc_pd_species <- function(species, tree) {
  if (length(species) == 0) return(NA)
  sub_tree <- drop.tip(tree, setdiff(tree$tip.label, species))
  sum(sub_tree$edge.length)
}

#Calculate Gamma PD (total PD across all regions where a species occurs)

species_PD_union <- apply(mat, 1, function(x) { #THIS CODE TAKES A LOT OF TIME
  
  # Regions where focal species occurs
  regions <- colnames(mat)[x == 1]
  
  if (length(regions) == 0) return(NA)
  
  # Species occurring in any of those regions 
  sp_union <- rownames(mat)[rowSums(mat[, regions, drop = FALSE]) > 0]
  
  # Calculate PD for this union
  calc_pd_species(sp_union, wcvp.v11.tree.wcvp.id.sub)
})


gamma_df <- data.frame(
  plant_name_id = rownames(mat),
  Gamma_PD = species_PD_union
)

#Calculate Mean Alpha PD

#First, calculate total area of species' native range

species_areas_df <- data.frame(plant_name_id = as.numeric(), 
                         region_area = as.numeric())

for(i in 1:length(SpecReg2)){ #THIS CODE TAKES A LOT OF TIME
  
  #Get plant name ID for each species
  
  plant_name_id <- SpecReg2[i, ]$plant_name_id
  
  #Get the names of TDWG3 regions within the range of that species
  
  species_regions <- filter(SpecReg2, plant_name_id==SpecReg2[i, ]$plant_name_id) %>%
    select(where(~is.numeric(.x) && any(.x != 0))) %>%
    colnames()
  
  #Get the areas of those TDWG3 regions
  
  tdwg3_area_species <- tdwg3_area %>%
    filter(Code %in% species_regions)
  
  #Calculate total area for each species
  
  region_area <- sum(tdwg3_area_species$Area)
  
  values <- data.frame(plant_name_id, region_area)
  
  species_areas_df <- bind_rows(species_areas_df, values)
  
}

#Join with Gamma PD

gamma_mean_df <- inner_join(gamma_df, species_areas_df, by="plant_name_id")

##Regress the PD of each species' native range (y-axis) against range size (x-axis), extract residuals (to get area-corrected PD across
#whole range i.e. Mean Alpha PD)

mod <- lm(log(Gamma_PD) ~ log(region_area), data = gamma_mean_df)
gamma_mean_df$Mean_AlphaPD <- residuals(mod)

######################################################

#Make final dataframe

######################################################

PD_dataframe <- inner_join(gamma_mean_df, max_med_df, by="plant_name_id")

######################################################

#METHODS FIGURE: Mean Alpha PD and number of experimental species by region

#######################################################

#Read in the shapefile of regions
sf_tdwg3_ea <- st_read("level3.shp") |> #note the .shx, .dbf and .prj files also need to be in the working directory
  st_make_valid() |>
  st_wrap_dateline() |>
  st_transform("+proj=eqearth")
plot(sf_tdwg3_ea["geometry"])
st_crs(sf_tdwg3_ea)

#Read in phylogenetic diversity by region - cleaned version of pd_tdwg3_area above, with addition of number of 
#experimental species native to each TDWG3 region

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



