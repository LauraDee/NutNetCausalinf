########################################################################
## NutNet Main Text Models ###########################################
## #######################################################################
##########################################################################
rm(list = ls())

# Define project directory
#  cdir <- "C:/GitHub/NutNetCausalinf/"
#  setwd(cdir)
setwd("~/Documents/GitHub/NutNetCausalinf/")  

######################################################
## BELOW THIS POINT, code should just run ##

# load packages; version numbers are noted for each package used.
library(dplyr) 
require(ggplot2) # 3.3.3
library(plyr) # 1.8.6
library(data.table) # v 1.13.6
library(sandwich) #3.0-0
library(foreign) # 0.8-80
library(car)  #v  3.0-10
library(fixest)  # v 0.8.2
library(lme4)  # 1.1-26
library(texreg) # 1.37. 5  
library(broom)  # v 0.7.4
library(tidyverse)  # v 1.3.0 
library(RColorBrewer) #1.1-2
library(cowplot) # 1.1.1
library(corrplot)  # 0.84
library(gridExtra)


###  *** Need to run *** 
### purpose built functions 
source("./code/r/useful_functions.R")

##########################################
## Analysis using small complete data ####
combonly <- TRUE  # combonly -> finalprocess_and_datachecks
comb <- fread("./data/processed/NutNetControlPlotData_v201804.csv",na.strings='NA')

###### Note on data version ######
# these analyses used the 'comb-by-plot-clim-soil-diversity-09-Apr-2018.csv' with updates to the biomass data.
# Updated with notice of biomass issue August 2022. 
# Updated biomass data for comp.pt 2015 data and marc.ar 2011 and 2012 live mass data from Peter W. 
#replace data in NutNetControlPlotData_v201804.csv for comp.pt 2015 data and marc.ar 2011 and 2012 with data in the following file

#load in new data  for comp.pt 2015 data and marc.ar 2011 and 2012 (see read me doc about this update)
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

write.csv(comb, "./data/NutNetControlPlotData_biomassfixAug22.csv")

source("./code/r/finalprocess_and_datachecks.R") ## Produces Table S1

source("./code/r/analysis_main.R") ## Produces Figures 2A, 2B, 3, and Tables S2, S3, and Figure S4 
source("./code/r/analysis_sm5.R") ## Produces Tables S4, S5, and S6
source("./code/r/FiguresS1_S2_S3.R") #Produces Figures S1, S2 and S3
###########################################################################################################################
## Analysis using large cover dataset from Nutrient Network, processed from version 'full-cover-09-April-2018.csv'    ####
############################################################################################################################
## *cover data is in separate location (not on Github) due its size!!!!!*

rm(list=setdiff(ls(),c("cdir","ihs","tidy")))
setwd("~/Documents/GitHub/NutNetCausalinf/")  
source("./code/r/useful_functions.R")
combonly <- FALSE
#comb <- fread("./data/processed/NutNetControlPlotData_v201804.csv",na.strings='NA')
comb <- fread("./data/NutNetControlPlotData_biomassfixAug22.csv",na.strings='NA')
source("./code/r/finalprocess_and_datachecks.R") # Doesn't produce Table S1 this time

### Load cover Data  - data processed using the code "ProcessNutNet_coverData_FINAL_public.R" which created all variables of SR counts and groupings 
cover <- fread("./data/NutNetCoverData_ProcessedFinal.csv")

source("./code/r/finalprocess_coverdata.R") # **need to run this line to prep the data for running the models** 
#the finalprocess_coverdata.R also produces Figure S10.

source("./code/r/analysis_fig5_smsection8.R") ## Produces Figure 5, and tables and supplemental results for section S8 of the supplemental materials

#code to create other supplemental figures
#code to create rank abundance curves for each site
source("./code/r/rank_abundance_curves.R") # need to run the finalprocess_coverdata.R to produce. 
source("./code/r/analysis_fig5_smsection8.R") #need to run the finalprocess_coverdata.R to produce.
source("./code/r/FiguresS1_S2_S3.R")
