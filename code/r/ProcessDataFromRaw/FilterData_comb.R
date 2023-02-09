########################################################################################
## Filter Data for Publication for Comb file ###########################################
################################################################################
# The data use policy of the Nutrient Network is to publish the minimum dataset required to produce the analyses, from the derived data. 
# This code is used to filter to the data to meet their requirements, though we also provide the full code to show all steps from our data processing from the raw data for reproducibility.
# The raw data for unmanipulated plots that were not included in the analyses, because they did not meet the inclusion criteria, 
# are available under restricted access for which permission can be obtained by contacting the Nutrient Network at https://nutnet.org. 

# Written by Laura Dee Feb 9 2023

#Close graphics and clear local memory
rm(list = ls())

# Define project directory
setwd("~/Documents/GitHub/NutNetCausalinf/")  

# load packages; version numbers are noted for each package used.
library(dplyr) 
library(plyr) # 1.8.6
library(data.table) # v 1.13.6
library(broom)  # v 0.7.4
library(tidyverse)  # v 1.3.0 
library(RColorBrewer) #1.1-2


comb <- fread("./data/processed/NutNetControlPlotData_v201804.csv",na.strings='NA')
comb$V1 = NULL
head(comb)
colnames(comb)

######################################################################################################################
### Biomass update and fix for marc and comp sites ######################################################################
######################################################################################################################
#load in new data  for comp.pt 2015 data and marc.ar 2011 and 2012 (see read me doc about this update: NutNet Biomass Issue log_August12_2022.docx)
updates <- fread("./data/marc-comp-comb-by.csv", na.strings = 'NA')
setnames(updates, old=c("live_mass","dead_mass","total_mass"), new=c("u_live_mass","u_dead_mass","u_total_mass"))
updates[,update:=T]

# merge with existing data
comb = merge(comb, updates, by=c("site_code","year","plot"), all.x=T)

# overwrite old biomass info with updated info for relevant updated recs
comb[update==T, `:=`(live_mass=u_live_mass,
                     dead_mass=u_dead_mass,
                     total_mass=u_total_mass)]
# remove temporary columns
comb[,`:=`(u_live_mass=NULL,
           u_dead_mass=NULL,
           u_total_mass=NULL,
           update=NULL)]

# write.csv(comb, "./data/NutNetControlPlotData_biomassfixAug22.csv")


combfiltered = comb[, c("site_code", "year", "plot", "newplotid", "block", "region", "country", "habitat",  "continent",   "country",  "region","managed" ,   "burned", "grazed",  "anthropogenic",                    
               "habitat", "elevation" , "TEMP_VAR_v2", "MIN_TEMP_v2", "MAX_TEMP_v2", "TEMP_WET_Q_v2", "TEMP_DRY_Q_v2", "TEMP_WARM_Q_v2", "TEMP_COLD_Q_v2",
               "pct_C" , "pct_N", "ppm_P", "ppm_K",  "ppm_Na",  "ppm_Mg", "ppm_S", "ppm_Na", "ppm_Zn",  "ppm_Mn", "ppm_Fe","ppm_Cu", "ppm_B",  "proportion_par", "avg.trt.neigh.rich.within.block",
               "pH", "PercentSand", "PercentSilt", "PercentClay", "total_cover", "rich", "even", "simpson", "shan", "has.5.yrs.data" , "is.PretreatmentYr" , "is.control",
               "newblockid", "laggedrich" ,  "laggedlive_mass", "live_mass", "site_year_rich", "site_richness", "site_native_richness" , "site_introduced_richness", "changeEvenness",
               "site_live_mass.yr",  "ave_site_live_mass", "initial_site_rich", "laggedrich","changerich", "changelive_mass", "lagged_changelive_mass","lagged_changerich",
               "rich_increase", "rich_decrease","rich_nochange", "total_cover" ) ]
 

write.csv(comb, "./data/NutNetControlPlotData_derived.csv")

