########################################################################
## NutNet Main Text Models ###########################################
## #######################################################################
##########################################################################
rm(list = ls())

# Define project directory
# cdir <- "C:/GitHub/NutNetCausalinf/"
# setwd(cdir)

setwd("~/Documents/GitHub/NutNetCausalinf/")  

######################################################
## BELOW THIS POINT, code should just run ##

# load packages; version numbers are noted for each package used.
require(ggplot2) # 3.3.3
library(plyr) # 1.8.6
library(data.table) # v 1.13.6
library(AER) # v 1.2-9 ### Do we use this anymore??? 
library(sandwich) #3.0-0
library(foreign) # 0.8-80
library(car)  #v  3.0-10
library(fixest)  # v 0.8.2
library(lme4)  # 1.1-26
library(texreg) # 1.37. 5  #### Do we use this anymore??? 
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
source("./code/r/finalprocess_and_datachecks.R") ## Produces Table S1

source("./code/r/analysis_main.R") ## Produces Figures 2A, 2B, 3, and Tables S2, S3, and Figure S4 
source("./code/r/analysis_sm5.R") ## Produces Tables S4, S5, and S6

###########################################################################################################################
## Analysis using large cover dataset from Nutrient Network, processed from version 'full-cover-09-April-2018.csv'    ####
############################################################################################################################
## *cover data is in separate location due its size!!!!!*

rm(list=setdiff(ls(),c("cdir","ihs","tidy")))
combonly <- FALSE
comb <- fread("./data/processed/NutNetControlPlotData_v201804.csv",na.strings='NA')
source("./code/r/finalprocess_and_datachecks.R") # Doesn't produce Table S1 this time

### Load cover Data  - data processed using the code "ProcessNutNet_coverData_FINAL_pulic.R" which created all variables of SR counts and groupings 
cover <- fread("./data/NutNetCoverData_ProcessedFinal.csv")

source("./code/r/finalprocess_coverdata.R") # **need to run this line to prep the data for running the models** 
#the finalprocess_coverdata.R also produces Figure S10.

source("./code/r/analysis_fig5_smsection7.R") ## Produces Figure 5, and tables and supplemental results for section S9 of the supplemental materials

#code to create other supplemental figures
source("./code/r/analysis_fig5_smsection7.R")
source("./code/r/SM_Figures.R")

#code to create rank abudance curves for each site 

#experimental analysis code, see the folder: 

