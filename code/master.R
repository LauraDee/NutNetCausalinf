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

require(ggplot2)
library(plyr)
library(data.table)
library(AER)
library(sandwich)
library(foreign)
library(car)
library(lfe)
library(fixest)
library(lme4)
library(texreg)
library(broom)
library(tidyverse)
library(RColorBrewer)
library(cowplot)

### purpose built functions ###
source("./code/r/useful_functions.R")

##########################################
## Analysis using small complete data ####
combonly <- TRUE  # combonly -> finalprocess_and_datachecks
comb <- fread("./data/processed/NutNetControlPlotData_v201804.csv",na.strings='NA')
source("./code/r/finalprocess_and_datachecks.R") ## Produces Table S1

source("./code/r/analysis_main.R") ## Produces only Figures 2A, 2B, 3, and Table S2 (maybe others???)

##########################################
## Analysis using large coverage data ####
## cover data is in separate location!!!!!
rm(list=setdiff(ls(),c("cdir","ihs","tidy")))
combonly <- FALSE
comb <- fread("./data/processed/NutNetControlPlotData_v201804.csv",na.strings='NA')

cover <- fread("C:/Dropbox/Dropbox/Research_Active/IV in ecology/NutNet/NutNetCoverData_ProcessedAug2019.csv")  

source("./code/r/finalprocess_and_datachecks.R") # Doesn't produce Table S1 this time
source("./code/r/finalprocess_coverdata.R") # Doesn't produce Table S1 this time

source("./code/r/analysis_fig5.R") ## Produces only Figures 2A, 2B, 3, and Table S2 (maybe others???)
