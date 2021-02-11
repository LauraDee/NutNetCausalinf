########################################################################################################
## NutNet Supplemental Abalayses - SM Section 5 - Laura Dee   ###########################################
########################################################################################################

## Code for analyses presented in the Supplemental Meterials for Section 5:
# 5b: Checking assumptions about functional forms for the effect of biodiversity on productivity
# 5c:
# 5d: 

#Close graphics and clear local memory
#graphics.off()
rm(list=ls())

#load libraries
require(ggplot2)
library(plyr)
library(data.table)
library(AER)
library(sandwich)
library(foreign)
library(car)
library(lfe)
library(texreg)
library(broom)
library(tidyverse)
library(RColorBrewer)
library(cowplot)

# When log(0) need to use inverse hyperbolic sine transformation (REF NEEDED)
#https://en.wikipedia.org/wiki/Inverse_hyperbolic_functions#Inverse_hyperbolic_sine
ihs = function(x) {
  return(log(x + sqrt(x^2+1)))
}
# v.s. log(x+1) <- is defined for a negative x. 

##Load processed Data, processed from version 'comb-by-plot-clim-soil-diversity-09-Apr-2018.csv'
setwd("~/Dropbox/IV in ecology/NutNet")
comb <- fread("NutNetControlPlotDataToUseApril2018.csv",na.strings='NA')

length(comb$live_mass) #74 NAs for live mass
summary(comb$live_mass)
length(comb$rich)
summary(comb$rich) # 2 NAs for richness

#################################################################
### Process data to prep to use it in the models ###############
###############################################################
comb$site <- comb$site_code

# Filter data to records with non-NA live_mass and non-NA richness.
# Models can generally handle NAs by implicitly dropping, but clustered SE
# function currently does not.
comb = comb[!is.na(live_mass) & !is.na(rich)]

#count for rich and live_mass
length(comb$rich)
# 1252
length(comb$live_mass)
# 1252
# Determine the number of Obs.
nrow(comb) #1252

## Confirm that only control plots are in the data
table(comb$trt)
## Confirm the # of years 
table(comb$year)

## Check to make sure that site and year are factors/characters
class(comb$year)
class(comb$site_code) # 43
class(comb$newplotid) #172

## 
length(unique(comb$newplotid)) #172
comP = comb[site_code == "comp.pt ",]

#see if there are some plots with only one year - one control observation:
comb[,.N, by=c("newplotid")][N==1,]  # there are 21 that need to be removed.

#Create a new version of comb filtering out singletons- plots with only one year that is one control observation:
# flag the singleton
comb[,singleton:=(.N==1), by=c("newplotid")]
comb[singleton==T, .(newplotid)]  # these are the singletons

#remove the singletons, select rows that are NOT singletons 
comb = comb[singleton == F, ]

# Determine the number of Obs. removing the singletons and obs with NA or rich or live_mass
nrow(comb) #1231 
length(unique(comb$newplotid)) #151

##############################################################################
### Make dummy variables for the panel regression  analyses #################
##############################################################################

# make year a character, to be a dummy variable: 
comb$year <- as.character(comb$year)

#make a factor that is site by year
comb[, site.by.yeardummy := paste(site_code, year, sep = "_")]

############################################################################################################################
#### RUN MAIN DESIGN (as in Figure 2) Panel FE Plus: Plot-FE and site*year effects  ############################################
############################################################################################################################
######################################################################################################################
### Main Models (Models 1). Log-Log and fixed effects/dummies only. ################################################################
##################################################################################################################
##### Models with productivity as the Y variable as log live mass ##########
#A.  Log-log and fixed effects/dummies only.
ModPFE <- felm(log(live_mass) ~ log(rich)  | newplotid + site.by.yeardummy | 0 | newplotid, data = comb)
summary(ModPFE, robust = TRUE)

##################################################################################################
###### SM S5b: FUNCTIONAL FORM ASSUMPTION CHECK  ###############################################
####################################################################################################
### Check Different Functional form assumptions, though the log-log for the effect of richness on productivity 
# is supported by theory and past work. e.g., reviewed in Cardinale, B.J., et al. (2011). The functional role of producer diversity in ecosystems. Am. J. Bot., 98, 572–592.
#  Levels means untransformmed variables.

## Run these models for results in SM Table S3.

#B. Log-Level
ModPFE.loglevel <- felm(log(live_mass) ~ rich | newplotid + site.by.yeardummy | 0 | newplotid, data = comb, exactDOF='rM')
summary(ModPFE.loglevel , robust = TRUE)

#C. Levels #B.
ModPFE.levels <- felm(live_mass ~ rich | newplotid + site.by.yeardummy| 0 | newplotid, data = comb, exactDOF='rM')
summary(ModPFE.levels, robust = TRUE)

#print results
#output models into a single table
screenreg(list(ModPFE, ModPFE.levels , ModPFE.loglevel), 
          custom.model.names=c("Log-Log", "Levels", "Log-Level" ))

#D - Quadratic form 
ModPFE.quad <- felm(live_mass ~ rich + I(rich^2) | newplotid + site.by.yeardummy | 0 | newplotid, data = comb, exactDOF='rM')
summary(ModPFE.quad, robust = TRUE)
# With the quadratic specification, we find that the estimated effect or richness on productivity 
# only turns positive in plots over 31 species, which represents only 14 observations and 1.14% of the data 

#E. Cubic 
ModPFE.cubic <- felm(live_mass ~ rich + I(rich^2) +I(rich^3)   | newplotid + site.by.yeardummy | 0 | newplotid, data = comb, exactDOF='rM')
summary(ModPFE.cubic, robust = TRUE)
# The results imply that the cubic specification does not fit well (and thus arent reported in Table S3)
# Coefficients for level and squared richness are nearly identical to values in quadratic specification (ModPFE.quad)
# but very imprecisely estimated (p>0.4 for the squared term)
# and cubic coefficient is nearly zero (.001078, p=0.94).

#print all functional form results into a single table
screenreg(list(ModPFE, ModPFE.levels, ModPFE.loglevel, ModPFE.quad, ModPFE.cubic), 
          custom.model.names=c("Log-Log (Main)", "Levels", "Log-Level", "Quadratic", "Cubic"))


#######################################################################################################################################
### Supplemental Analyses: Test robustness to moderators: SM sections S5c and S5d ################################################################################################
#######################################################################################################################################

##### Analyses for Section S5c: Moderating Effect of site-level species richness ####
## See Table S6: code to reproduce Table S6 is below.

#average site richness as moderator
ModModSite1 <- felm(log(live_mass) ~ log(rich) + log(rich):site_richness | newplotid + site.by.yeardummy | 0 | newplotid, data = comb, exactDOF='rM')
summary(ModModSite1 , robust = TRUE)

#ave site INT (introduced) richness  as moderator 
ModModSite2 <- felm(log(live_mass) ~ log(rich) + log(rich):site_introduced_richness | newplotid + site.by.yeardummy | 0 | newplotid, data = comb, exactDOF='rM')
summary(ModModSite2 , robust = TRUE)

#yearly site  richness as moderator
ModModSite3 <- felm(log(live_mass) ~ log(rich) + log(rich):site_year_rich | newplotid + site.by.yeardummy| 0 | newplotid, data = comb, exactDOF='rM')
summary(ModModSite3 , robust = TRUE)

#ave site NATIVE richness  as moderator
ModModSite4 <- felm(log(live_mass) ~ log(rich) + log(rich):site_native_richness | newplotid + site.by.yeardummy| 0 | newplotid, data = comb, exactDOF='rM')
summary(ModModSite4, robust = TRUE, cluster = TRUE)

#print results - simple table for main results 
screenreg(list(ModModSite1, ModModSite2, ModModSite3, ModModSite4),     # object with results 
          custom.model.names= c("Mean Site SR" , "Site Introduced SR", "Site Yearly SR", "Site Native SR"))


##** I CUT THIS FROM THE SM SO MAYBE WE SHOULD CUT?????????***
###########################################################################
## With total cover: ########################################################
###########################################################################
ModModCover1 <- felm(log(live_mass) ~ log(rich) + log(rich):total_cover | newplotid + site.by.yeardummy | 0 | newplotid, data = comb, exactDOF='rM')
summary(ModModCover1, robust = TRUE, cluster = TRUE)

ModModCover2 <- felm(log(live_mass) ~ log(rich) + total_cover | newplotid + site.by.yeardummy, data = comb, exactDOF='rM')
summary(ModModCover2, robust = TRUE, cluster = TRUE)

screenreg(list( ModModCover2,  ModModCover1 ),   # object with results 
          custom.model.names= c("Total Cover", "Cover interacted with Richness " ))

# **** CUT in between????? *******

#######################################################################################################################################
### Supplemental Analyses:  S5d Moderating Effect of site-level productivity  ################################################################################################
#######################################################################################################################################

## See Table S5 & 6: code to reproduce results

###########################################################
### Compute average site-level productivity (live mass) ###
###########################################################
comb[, ave_site_livemass := ave(live_mass, na.rm = T), by = .(site_code)]
hist(comb$ave_site_livemass)

comb[, ave_site_livemass.peryr := ave(live_mass, na.rm = T), by = .(site_code, year)]
hist(comb$ave_site_livemass.peryr)

##Test for heterogeneous treatment effects moderated by site-level productivity levels - analyses for Table S5
#A.  Log-log and fixed effects/dummies only.
ModPFE.prod <- felm(log(live_mass) ~ log(rich) + log(rich):ave_site_livemass  | newplotid + site.by.yeardummy, data = comb, exactDOF='rM')
summary(ModPFE.prod, robust = TRUE, cluster = TRUE)

#with the average site-level live mass per year (likely the one we want!)
ModPFE.prod2 <- felm(log(live_mass) ~ log(rich) + log(rich):ave_site_livemass.peryr  | newplotid + site.by.yeardummy, data = comb, exactDOF='rM')
summary(ModPFE.prod2, robust = TRUE, cluster = TRUE)

#output productivity moderator results into a single table - Table S5
screenreg(list(ModPFE.prod2, ModPFE.prod ),       # object with results from clx
          custom.model.names=c("Average Prod per site &  year", "Average Prod per site"))



#### Analyses for Table S6 using cut-offs akin to Wang et al. 2019 Nature Communications

## Using the productivity groups from Wang et al 2019 Nature Communications
# Hi Laura,
# Here is Yongfan’s response:
# The 151 grids in HerbDivNet data were divided into three equal groups,
# depending on their mean productivity: low, medium, and high productivity (with 50 to 51 grids each).
# The cutoff values of primary productivity of each group:
#   
# Low (51 grids): 30.18-238.73 (g/m^2)
# Medium (50 grids): 239.67-409.69 (g/m^2)
# High (50 grids): 414.29-1382.42 (g/m^2)
# 
# Best,
# Michel

# do the groups based on the overall average since Wang et al was a cross-section so that is more comparable
summary(comb$ave_site_livemass.peryr)
summary(comb$ave_site_livemass)
#check that it worked length(unique(comb$ave_site_livemass))  == 43 

## Do a different cut off Low, Medium, High groups of site livemass, based on Wang et al cut-offs
comb[, ProdGroup_WangCutoffs:=cut(ave_site_livemass, breaks=c(30.18,239.67,414.20,1609), labels=c("Low","Medium","High")), by = .(site_code) ]
head(comb)
table(comb$ProdGroup_WangCutoffs)
summary(comb$ProdGroup_WangCutoffs)
plot(comb$ProdGroup_WangCutoffs, main = "Productivity Groups (Average Across Years)")

## the Wang et al. paper uses a single year of data. Thus, we also need to create these cut-offs for each site and year in our data:
comb[, ProdGroup := cut(ave_site_livemass.peryr, breaks=c(30.18,239.67,414.20,1609), labels=c("Low","Medium","High"))]
head(comb)
table(comb$ProdGroup)
plot(comb$ProdGroup, main = "Productivity Groups Per Year")

## For results in Table S6 -- column 1. 
ModPFE.prodgroup.peryr <- felm(log(live_mass) ~ log(rich) + log(rich):ProdGroup  | newplotid + site.by.yeardummy, data = comb, exactDOF='rM')
summary(ModPFE.prodgroup.peryr , robust = TRUE, cluster = TRUE)

## For results in Table S6 -- column 2. 
ModPFE.prodgroup.wang <- felm(log(live_mass) ~ log(rich) + log(rich):ProdGroup_WangCutoffs  | newplotid + site.by.yeardummy, data = comb, exactDOF='rM')
summary(ModPFE.prodgroup.wang , robust = TRUE, cluster = TRUE)

#output productivity moderator results into a single table: Table S6
screenreg(list(ModPFE.prodgroup.peryr , ModPFE.prodgroup.wang ),       # object with results from clx
          custom.model.names=c("Average Prod per site &  year", "Average Prod per site"))

