1. PD Calculations

This folder describes how to calculate PD metrics, and provides a visualisation of global PD trends (Fig. 1 main text).

It contains the following files:

ONE_PD_calculations.R - code file for running analyses

SpecReg2.Rdata - file containing WCVP_ID names for all accepted seed plants and which TDWG3 regions of the world they are native to

wcvp.v11.tree.wcvp.id.sub.RData - phylogenetic supertree containing all species in SpecReg2.RData

tdwg3_areas.csv - file containing the area (km^2) of every TDWG3 region

TDWG3_PD.csv - file containing the Mean Alpha PD for each TDWG3 region

level3.shp/shx/prj/dbf - shapefile and associated files containing TDWG3 regions

NOTE: This code file takes a long time to run, and is provided here for reproducibility. The PD values for all focal species are provided in subsequent data files (Parts 2-6 of this repository) - it is NOT necessary
to run this file in order to carry out the analyses described in the other parts, or to produce Fig. 1, for which a cleaned file is provided.
