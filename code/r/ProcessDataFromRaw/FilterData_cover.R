########################################################################################
## Filter Data for Publication for Cover file ###########################################
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

#read file with extra columns that need to be filtered per the Nutrient Network data sharing guidelines:
cover <- fread("data/NutNetCoverData_ProcessedFinal.csv",na.strings='NA')

cover = cover[trt == "Control",]

### Confirm that only control plots are in the data
table(cover$trt)
### Confirm the # of years 
table(cover$year)

###############################################################################
#### Filter Data by years   ################################################
################################################################################ 

# in cover - filter to control plots with 5 years of data to match the main analyses
#Get maximum number of treatment years per site
cover[,max.trt.yr:=max(year_trt), by=.(site_code)]
#Find min value of treatment years; 0 = pre-treatment data
cover[,min.trt.yr:=min(year_trt), by=.(site_code)]

#filter to those sites and years
cover[, has.5.yrs.data:=(min.trt.yr == 0 & max.trt.yr >= 4)|(min.trt.yr==1 & max.trt.yr>=5)] #this keeps control plots without a pre-treatment survey
cover = cover[has.5.yrs.data==T,]

## the filter missed two sites that only have 4 years of data - but over the span of 5 years 
# Azi is missing live_mass data in 2012. Barta.us is missing data for all variables 2010, so only has 4 years of data.
# remove these 2 sites from the data with the filter that AT LEAST 5 years of data is needed for the site to be included:
cover = cover[site_code != "azi.cn",]
cover = cover[site_code != "barta.us",]

#################################################################################################
## Grab the summary columns from cover to merge  #######################################################
################################################################################################

#to subset columns and also remove duplicate rows from the cover file so that there is one observation per plot and year 
# and the data isn't artificially replicated 
coversummaries = unique(cover[, .(site_code, year,  site_name,  plot,  year_trt , trt, totplotcover.yr.live, sr_NA,
                                  sr_INT, sr_NAT, sr_UNK, sr_INT.site, 
                                  sr_non.nat_rare,  sr_nat_rare, sr_non.rare_non.nat, sr_non.rare_nat, sr_nat_dom, sr_non.nat_dom, 
                                  sr_nat_unk_rare, ## 2. Including the unknown spp origin all as native: ####
                                  sr_non.nat_unk_rare, # 3.Including them all as non-native: 
                                  sr_non.rare_nat_unk, ## 2. Include the unknown spp origin all as native: ####
                                  sr_non.rare_non.nat_unk , # 3.Including them all as non-native
                                  sr_non_rare_spp.Freq, sr_non.rare_nat.Freq, sr_non.rare_non.nat.Freq,
                                  sr_rare_non.nat.Freq, sr_rare_nat.Freq,
                                  sr_rare_unk_nat.Freq , sr_non.rare_nat_unk.Freq, ## 2. Including the unknown spp origin all as native: ####
                                  sr_non.rare_non.nat_unk.Freq, sr_rare_non.nat_unk.Freq,   # 3.Including them all as non-native
                                  sr_non.rare_nat2, sr_non.rare_non.nat2 , sr_nat_rare2, sr_non.nat_rare2 , #include variables for cut-off 2 as sensitivity test for main model for Figure 5
                                  sr_non.rare_nat3, sr_non.rare_non.nat3 , sr_nat_rare3, sr_non.nat_rare3 #include variables for cut-off 2 as sensitvity test for main model for Figure 5
)])
nrow(coversummaries)

write.csv(coversummaries, "./data/NutNetCover_derived.csv")
