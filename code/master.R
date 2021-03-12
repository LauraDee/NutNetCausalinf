########################################################################
## NutNet Main Text Models ###########################################
## Laura Dee - in LFE ###############################################
##########################################################################
# updated Dec 27, 2019
#updated May 11, 2020 to include Simpsons div models 
#cleaned up to focus on main results for Figures 2 & 3 on Feb 5 2021 

#notes on using lfe versus previous implementation of models:
#it does, especially if you're in a "FE nested within clusters" setting
#which the code I'd given you before did not account for
# since i was unaware of it at the time
# basically when doing degrees of freedom corrections, lfe won't count fixed effects nested
# in clusters b/c that would sort of be double counting

#plotting coefficient estimates from felm objects:
# https://raw.githack.com/uo-ec607/lectures/master/08-regression/08-regression.html#high_dimensional_fes_and_(multiway)_clustering

#plotting resurces: http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
# theme(legend.position="top") 
#http://www.sthda.com/english/wiki/ggplot2-legend-easy-steps-to-change-the-position-and-the-appearance-of-a-graph-legend-in-r-software

rm(list = ls())

# Define project directory
cdir <- "C:/GitHub/NutNetCausalinf/"

###########################
## BELOW THIS POINT, code should just run ##

setwd(cdir)
setwd("~/Documents/GitHub/NutNetCausalinf/")
# load packages; version numbers are noted for each package used.
require(ggplot2) # 3.3.3
library(plyr) # 1.8.6
library(data.table) # v 1.13.6
library(AER) # v 1.2-9 ### Do we use this anymore??? 
library(sandwich) #3.0-0
library(foreign) # 0.8-80
library(car)  #v  3.0-10
library(lfe)  # 2.8-6
library(fixest)  # v 0.8.2
library(lme4)  # 1.1-26
library(texreg) # 1.37. 5  #### Do we use this anymore??? 
library(broom)  # v 0.7.4
library(tidyverse)  # v 1.3.0 
library(RColorBrewer) #1.1-2
library(cowplot) # 1.1.1
library(corrplot)  # 0.84

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


##########################################
## Analysis using large coverage data ####
## cover data is in separate location!!!!!

rm(list=setdiff(ls(),c("cdir","ihs","tidy")))
combonly <- FALSE
comb <- fread("./data/processed/NutNetControlPlotData_v201804.csv",na.strings='NA')

cover <- fread("~/Dropbox/IV in ecology/NutNet/NutNetCoverData_ProcessedAug2019.csv")  

source("./code/r/finalprocess_and_datachecks.R") # Doesn't produce Table S1 this time
source("./code/r/finalprocess_coverdata.R") 

source("./code/r/analysis_fig5.R") ## Produces Figure 5
source("./code/r/analysis_sm7.R") ## In progress <_> LAURA TO WORK ON.




