########################################################################
## Plot Rare species richness  ###########################################
## Laura Dee  - Update Feb 18 2021 ###############################################
##########################################################################


#Close graphics and clear local memory
#graphics.off()
rm(list=ls())

#load libraries
require(ggplot2)
library(plyr)
library(data.table)
library(sandwich)
library(foreign)
library(car)
library(lfe)
library(texreg)
library(broom)
library(tidyverse)
library(cowplot)

setwd("~/Dropbox/IV in ecology/NutNet")

## Run functions
# When log(0) need to use inverse hyperbolic sine transformation (Bellemare & Wichman 2020 Oxford Bulletin of Economics and Statistics)
#https://en.wikipedia.org/wiki/Inverse_hyperbolic_functions#Inverse_hyperbolic_sine
ihs = function(x) {
  return(log(x + sqrt(x^2+1)))
} # v.s. log(x+1) <- is defined for a negative x. 


### Load Data 
##Load Processed Cover Data, processed from version 'full-cover-09-April-2018.csv'
#code for data processing is NutNet_coverData_Sept2018.R
cover <- fread("NutNetCoverData_ProcessedAug2019.csv")  

##Load processed Data, processed from version 'comb-by-plot-clim-soil-diversity-28-Apr-2017.csv'
comb <- fread("NutNetControlPlotDataToUseApril2018.csv",na.strings='NA')
comb$site <- comb$site_code

length(comb$live_mass) #74 NAs for live mass
summary(comb$live_mass)
length(comb$rich)
summary(comb$rich) # 2 NAs for richness

################################################################# ############### ###############
### Process comb data to prep to use it in the models and harmonize with the main analyses ###############
############################################################### ############### ############### ###############
comb$site <- comb$site_code

# Filter data to records with non-NA live_mass and non-NA richness.
# Models can generally handle NAs by implicitly dropping, but clustered SE
# function currently does not.
comb = comb[!is.na(live_mass) & !is.na(rich)]

#count for rich and live_mass
length(comb$rich) # 1252
length(comb$live_mass) # 1252
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

####################################################################################################
## Filter cover data to only control plots & prep for merge and models ###############################
##################################################################################################
#comb = comb[is.control==T & has.5.yrs.data==T,]  #this should be done already in the saved file from "NutNetControlPlotDataToUseApril2018.csv"

### Confirm that only control plots are in the data
table(comb$trt)
# ## Confirm the # of years 
table(comb$year)

cover = cover[trt == "Control",]

### Confirm that only control plots are in the data
table(cover$trt)
### Confirm the # of years 
# table(cover$year)

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

# create a variable for one of the analyes for Fig 4
cover[, non_rare_spp.DI2 := sr_non.rare_non.nat2 + sr_non.rare_nat2]

#test
# cover = cover[subplot== "A",]
#trouble shooting
# table(cover$subplot, cover$site_code)

##########################################################################################
### Process comb data to prep to use it in the models & merge #############################
############################################################################################
## make year a character in both datasets, to be a dummy variable: 
comb$year <- as.character(comb$year)
cover$year <- as.character(cover$year)

# same with plot
comb$plot <- as.character(comb$plot)
cover$plot <- as.character(cover$plot)
cover$V1 = NULL

#check that the DIgroups 1 & 2 look OK #
table(cover$DIgroup)
table(cover$DIgroup2)

#################################################################################################
## Grab the summary columns from cover to merge  #######################################################
################################################################################################

#these variables also need to be removed because non_rare_spp2 is a series of TRUE or FALSE per species so its adding extra rows:
# non_rare_spp2

#to subset columns and also remove duplicate rows from the cover file so that there is one observation per plot and year 
# and the data isn't artificially replicated 
coversummaries = unique(cover[, .(site_code, year,  site_name,  plot,  year_trt , trt, totplotcover.yr.live, LegumePercentcover.yr, cover_nat_dom, cover_nat_sub,
                                  sr_nat_sub, sr_non.nat_sub, cover_tot_non.rare, sr_INT, sr_NAT, sr_domspp, sr_rarespp, sr_subordspp, sr_non_rare_spp, 
                                  sr_non.nat_rare,  sr_nat_rare, sr_non.rare_non.nat, sr_non.rare_nat, sr_nat_dom, sr_non.nat_dom, relabund_sr_domspp,
                                  sr_non_rare_spp.RelA, 
                                  sr_non_rare_spp.Freq, sr_non.rare_nat.Freq, sr_non.rare_non.nat.Freq,
                                  #rare_spp.DI2,
                                  sr_rare_non.nat.Freq, sr_rare_nat.Freq,
                                  sr_domspp2, sr_rarespp2 , sr_subordspp2, non_rare_spp.DI2,
                                  relabund_sr_rarespp, relabund_sr_subordspp,  sr_nat_dom.Freq, sr_non.nat_dom.Freq,
                                  sr_nat_sub.Freq, sr_non.nat_sub.Freq, 
                                  sr_non.rare_nat.RelA, sr_non.rare_non.nat.RelA, sr_rare_non.nat.RelA, sr_rare_nat.RelA,
                                  # sr_non.rare_nat.Freq,  sr_rare_non.nat.Freq, sr_rare_nat.Freq,
                                  sr_non.rare_nat.Freq2, sr_non.rare_non.nat.Freq2, sr_rare_non.nat.Freq2, sr_rare_nat.Freq2,
                                  sr_non.rare_nat.RelA2, sr_non.rare_non.nat.RelA2, sr_rare_non.nat.RelA2, sr_rare_nat.RelA2, # non_rare_spp2,
                                  sr_non.rare_nat2, sr_non.rare_non.nat2, sr_non.nat_rare2, sr_nat_rare2, sr_Nfixer, 
                                  # sr_rarespp2, 
                                  #rare_spp.DI2,  
                                  sr_non.Nfixer, N_fixer_cover.yr,
                                  NonNative_cover.yr , Native_cover.yr , sr_annual, sr_peren, sr_null.lspan, sr_biennial ,
                                  sr_indeter, AnnualPercentcover.yr, PerenPercentcover.yr , sr_graminoid, sr_forbs, sr_woody,
                                  sr_legume, sr_bryophyte , sr_cactus, Dom_cover.yr, 
                                  GrassPercentcover.yr ,ForbPercentcover.yr, WoodyPercentcover.yr, 
                                  freq_sr_domspp, freq_sr_rarespp, freq_sr_subordspp
)])
#make sure number of rows isnt inflated 
nrow(coversummaries)

#################################################################################################
## Merge comb with Processed Cover Data #######################################################
################################################################################################
# Merge the cover data with the comb data. Processed the April 2018 versions of the datasets.
# merge to keep the data that is in the comb file (i.e., sites with at least 5 years of data for control plots)
mech.data = merge(comb, coversummaries, by=c("site_code","plot","year"), all.x=T)

#what I had been doing: mech.data = merge(comb, cover, by=c("site_code","plot","year"), all.x=T)
nrow(mech.data) #should be 1231
dim(mech.data)
nrow(comb)
nrow(cover)

#to write out this data file:
# write.csv(mech.data, "merged_NutNet_Oct2019.csv")

##########################################################################################
### Process merge data to prep to use it in the models  #################################
##########################################################################################
# make year a character, to be a dummy variable: 
mech.data$year <- as.character(mech.data$year)
## *Check to make sure that site and year are factors/characters* important*
class(mech.data$year)
class(mech.data$site_code)
class(mech.data$newplotid)

# Filter data to records with non-NA live_mass and non-NA richness.
# Models can generally handle NAs by implicitly dropping, but clustered SE function currently does not.
mech.data = mech.data[!is.na(live_mass) & !is.na(rich)]

## Confirm that only control plots are in the data
table(mech.data$trt.x)

## Confirm the # of years 
table(mech.data$year)

## Confirm its the same site list
list(unique(mech.data$site_code))

#########################################################################################################
### ***** RUN For models****** ################################################################################
###########################################################################
#make a factor that is site by year
mech.data[, site.by.yeardummy := paste(site_code, year, sep = "_")]

mech.data[, even_year0 := even[year_trt.x == "0" ], by = .(plot, site_code)]

#########################################################################################################
### Create rare species richness change variables  ################################################################################
###################################################################################################

################################################################################################
### Create a first-difference ########################################################################
#################################################################################################
#as for changes you can use the shift function in data.table; but things need to be sorted first
# if you want changes, rows need to be sorted by time. Sort your data table by whatever the time column is. then use the shift function.
#if you already computed laggedrich, you don't need the by part or the order part. 
# they are only needed for constructing the lags   #comb[, changerich := rich-laggedrich,]

mech.data[order(year), changerich := rich-shift(rich), by =.(plot, site_code)]
comb[order(year), changelive_mass := live_mass-shift(live_mass), by =.(plot, site_code)]

#non-rare native species change 
mech.data[order(year), change_nonrare.native := sr_non.rare_nat-shift(sr_non.rare_nat), by =.(plot, site_code)]
hist(mech.data$change_nonrare.native, main = "Change in native, non-rare species richness")
#rare native species change 
mech.data[order(year), change_rare.native := sr_nat_rare -shift(sr_nat_rare), by =.(plot, site_code)]
hist(mech.data$change_rare.native, main = "Change in native, rare species richness")
summary(mech.data$change_rare.native)

#non-rare non-native species change - sr_non.rare_non.nat
mech.data[order(year), change_nonrare.nonnative :=sr_non.rare_non.nat -shift(sr_non.rare_non.nat), by =.(plot, site_code)]
hist(mech.data$change_nonrare.nonnative, main = "Change in non-native, non-rare species richness")

#rare non-native species change  sr_non.nat_rare
mech.data[order(year), change_rare.nonnative := sr_non.nat_rare -shift(sr_non.nat_rare), by =.(plot, site_code)]
hist(mech.data$change_rare.nonnative,  main = "Change in non-native, rare species richness")


# # function currently does not.
mech.data = mech.data[!is.na(change_rare.native)]

## plot change in rare species by site 
FigSX <- ggplot(data = mech.data, aes(x = change_rare.native)) + geom_histogram()+ facet_wrap(~site) + theme_bw() +
  geom_vline(xintercept=c(0,0), color = "blue", linetype="dashed") +
  labs(x = "Plot-level change in rare native species richness year to year") +  theme_bw() +
  theme(axis.title.y= element_text(size=14)) + theme(axis.title.x= element_text(size=12)) +
  theme(axis.text.y = element_text(size = 14)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(size=14)) 
FigSX


FigSX <- ggplot(data = mech.data, aes(x = change_nonrare.nonnative)) + geom_histogram()+ facet_wrap(~site) + theme_bw() +
  geom_vline(xintercept=c(0,0), color = "blue", linetype="dashed") +
  labs(x = "Plot-level change in nonrare nonnative species richness year to year") +  theme_bw() +
  theme(axis.title.y= element_text(size=14)) + theme(axis.title.x= element_text(size=12)) +
  theme(axis.text.y = element_text(size = 14)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(size=14)) 
FigSX



## plot change in rare non native species by site 
FigSX <- ggplot(data = mech.data, aes(x = change_rare.nonnative)) + geom_histogram()+ facet_wrap(~site) + theme_bw() +
  geom_vline(xintercept=c(0,0), color = "blue", linetype="dashed") +
  labs(x = "Plot-level change in rare native speciesrichness year to year") +  theme_bw() +
  theme(axis.title.y= element_text(size=14)) + theme(axis.title.x= element_text(size=12)) +
  theme(axis.text.y = element_text(size = 14)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(size=14)) 
FigSX

# # function currently does not.
mech.data = mech.data[!is.na(changerich)]

#plot change richness vs change rare spp
FigSX <- ggplot(data = mech.data, aes(x = change_rare.nonnative, y = changerich)) + geom_histogram()+ facet_wrap(~site) + theme_bw() +
  geom_vline(xintercept=c(0,0), color = "blue", linetype="dashed") +
  labs(x = "Plot-level change in rare native species vs s richness year to year") +  theme_bw() +
  theme(axis.title.y= element_text(size=14)) + theme(axis.title.x= element_text(size=12)) +
  theme(axis.text.y = element_text(size = 14)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(size=14)) 
FigSX

plot(mech.data$changerich, mech.data$change_rare.native)
plot(mech.data$changerich, mech.data$change_rare.nonnative)
plot(mech.data$changerich, mech.data$change_nonrare.nonnative)
plot(mech.data$changerich, mech.data$change_nonrare.native)




## create richness bins 
mech.data[ ,RichGroup := cut(rich, breaks=c(1, 5,10,37), labels=c("Low","Medium","High"))]
plot(mech.data$RichGroup)

mech.data[ ,RichGroup2 := cut(rich, breaks=c(1, 10,37), labels=c("Low","High"))]
plot(mech.data$RichGroup2)

#create proportion of rare native species:
#mech.data[ , proportion_rare  := (sr_non.nat_rare + sr_nat_rare)]
mech.data[ , proportion_rare  := (sr_nat_rare)]
mech.data[ , proportion_rare  := (proportion_rare/rich)]
hist(mech.data$proportion_rare)
hist(mech.data$proportion_rare)

#create proportion of change in richness as change in rare native species:
#mech.data[ , proportion_rare  := (sr_non.nat_rare + sr_nat_rare)]

#****THIS DINDT WORK ****
mech.data[ , proportion_changerare  := (change_rare.native/changerich)]
hist(mech.data$proportion_changerare)


#proportion non native spp
mech.data[ , proportion_nonnative  := (sr_INT/rich)]
hist(mech.data$proportion_nonnative )

summary(mech.data$site_year_rich)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#8.0    24.0    36.0    37.7    51.0    88.0 


## plot proportion of rare species in richness group
FigSX <- ggplot(data = mech.data, aes(x = proportion_rare)) + geom_histogram()+ facet_wrap(~RichGroup2) + theme_bw() +
  geom_vline(xintercept=c(0,0), color = "blue", linetype="dashed") +
  labs(x = "proportion rare native species out of total number of species per plot and yr ") +  theme_bw() +
  theme(axis.title.y= element_text(size=14)) + theme(axis.title.x= element_text(size=12)) +
  theme(axis.text.y = element_text(size = 14)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(size=14)) 
FigSX

FigSX <- ggplot(data = mech.data, aes(x = proportion_nonnative)) + geom_histogram()+ facet_wrap(~RichGroup2) + theme_bw() +
  geom_vline(xintercept=c(0,0), color = "blue", linetype="dashed") +
 # labs(x = "Plot-level change in rare native species richness year to year") +  theme_bw() +
  theme(axis.title.y= element_text(size=14)) + theme(axis.title.x= element_text(size=12)) +
  theme(axis.text.y = element_text(size = 14)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(size=14)) 
FigSX


#proportion rare by site
FigSX <- ggplot(data = mech.data, aes(x = proportion_rare)) + geom_histogram()+ facet_wrap(~site_code) + theme_bw() +
  geom_vline(xintercept=c(0,0), color = "blue", linetype="dashed") +
  labs(x = "Plot-level change in rare native species richness year to year") +  theme_bw() +
  theme(axis.title.y= element_text(size=14)) + theme(axis.title.x= element_text(size=12)) +
  theme(axis.text.y = element_text(size = 14)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(size=14)) 
FigSX


ModrMod_YrSiteProd  <- feols(log(live_mass) ~ log(rich) + log(rich):RichGroup  | newplotid + site.by.yeardummy, mech.data)



###########
### A. Grouped based on the Dominance Indicator (DI) and cutoffs of:  breaks=c(0.0,0.2,0.8,1.0),
##########
# these are the variable names/groups:
# sr_non.rare_nat  + sr_non.rare_non.nat + sr_non.nat_rare +  sr_nat_rare
Mod4A.1 <- felm(log(live_mass) ~ ihs(sr_non.rare_nat) + ihs(sr_non.rare_non.nat)  + ihs(sr_non.nat_rare) +  ihs(sr_nat_rare) | newplotid + site.by.yeardummy | 0 | newplotid, data = mech.data)
summary(Mod4A.1, robust = TRUE)

## Hypoth

