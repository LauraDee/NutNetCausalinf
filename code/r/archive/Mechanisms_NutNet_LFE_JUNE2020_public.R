########################################################################
## NutNet Mechanisms Analyses ###########################################
## Laura Dee - in LFE ###############################################
##########################################################################
# updated Dec 28, 2019 to run models with different cut-offs and do some robustness checks of 
# main Panel FE+ model with cover variables and to plot results

#updated June 22 2020 to finalize analyses; and add a rare vs non rare analysis with DIgroup2. 

#notes on using lfe versus previous implementation of models:
#it does, especially if you're in a "FE nested within clusters" setting
#which the code I'd given you before did not account for
# since i was unaware of it at the time
# basically when doing degrees of freedom corrections, lfe won't count fixed effects nested
# in clusters b/c that would sort of be double counting

#plotting coefficient estimates from felm objects:
# https://raw.githack.com/uo-ec607/lectures/master/08-regression/08-regression.html#high_dimensional_fes_and_(multiway)_clustering

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


# run this function - critical for workflow for plotting results:
tidy = function(model, ...) {
  data = cbind.data.frame(term = names(coef(model)),
                          coef(summary(model, robust= T)),
                          confint(model))
  names(data) = c("term","estimate","std.error", "statistic","p.value","conf.low","conf.high")
  return(data)
}


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


# sr_non.rare_non.nat.Freq, sr_non.rare_nat.Freq , 

# lagged and change variable removed from above but included here:
coversummaries.robustness = unique(cover[, .(site_code, year,  site_name,  plot,  year_trt , trt, totplotcover.yr.live, LegumePercentcover.yr, cover_nat_dom, cover_nat_sub,
                                  sr_nat_sub, sr_non.nat_sub, cover_tot_non.rare, sr_INT, sr_NAT, sr_domspp, sr_rarespp, sr_subordspp, sr_non_rare_spp, 
                                  sr_non.nat_rare,  sr_nat_rare, sr_non.rare_non.nat, sr_non.rare_nat, sr_nat_dom, sr_non.nat_dom, relabund_sr_domspp,
                                  sr_non_rare_spp.RelA, sr_non_rare_spp.Freq, sr_non.rare_nat.Freq, sr_rare_non.nat.Freq, 
                                  #rare_spp.DI2,
                                  sr_domspp2, sr_rarespp2 , sr_subordspp2, 
                                  relabund_sr_rarespp, relabund_sr_subordspp, sr_non.rare_non.nat.Freq,  sr_nat_dom.Freq, sr_non.nat_dom.Freq,
                                  sr_nat_sub.Freq, sr_non.nat_sub.Freq, sr_non.rare_nat.RelA, sr_non.rare_non.nat.RelA, sr_rare_non.nat.RelA, sr_rare_nat.RelA,
                                  # sr_non.rare_nat.Freq,  sr_rare_non.nat.Freq, sr_rare_nat.Freq, 
                                  non_rare_spp2, non_rare_spp.DI2,
                                  sr_non.rare_nat.Freq2, sr_non.rare_non.nat.Freq2, sr_rare_non.nat.Freq2, sr_rare_nat.Freq2,
                                  sr_non.rare_nat.RelA2, sr_non.rare_non.nat.RelA2, sr_rare_non.nat.RelA2, sr_rare_nat.RelA2, # non_rare_spp2,
                                  sr_non.rare_nat2, sr_non.rare_non.nat2, sr_non.nat_rare2, sr_nat_rare2, sr_Nfixer, 
                                  # sr_rarespp2, 
                                  #rare_spp.DI2,  
                                  sr_non.Nfixer, change_sr_Nfixer, lagged_sr_Nfixer, N_fixer_cover.yr ,change_N_fixer_cover, lagged_N_fixer_cover.yr,
                                  NonNative_cover.yr , Native_cover.yr , sr_annual, sr_peren, sr_null.lspan, sr_biennial ,
                                  sr_indeter, change_sr_annual, change_sr_peren, AnnualPercentcover.yr, PerenPercentcover.yr , sr_graminoid, sr_forbs, sr_woody,
                                  sr_legume, sr_bryophyte , sr_cactus,
                                  change_sr_non.live, GrassPercentcover.yr ,ForbPercentcover.yr, WoodyPercentcover.yr, lagged_LegumePercentCover.yr,
                                  freq_sr_domspp, freq_sr_rarespp, freq_sr_subordspp
                                )])



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


####################################################################################################################################
#### Models for the role of rare and non-native species  ##########################################################################################################
####################################################################################################################################

# This section is organized as follows, and accompanies section S7 in the SM. 
# Analysis for reproducing Figure 5 in the main text. 
# SM analyses for Section 7, where we do seperate analyses for:
#....... *** to fill in ****


######################################################################################################################
###  Figure 5 - Main text. Rare vs Non-Rare and Native vs Invasive  #######################################################################
####################################################################################################################
###########
### Figure 5 - Main text. Grouped based on the Dominance Indicator (DI) and cutoffs of:  breaks=c(0.0,0.2,0.8,1.0),
##########
#. We first present the analyses shown in the main text figure 5, which include 4 groups of species: 
# 1) rare, native: sr_nat_rare
#2) rare non-native: sr_non.nat_rare
#3) non-rare, native: sr_non.rare_nat
#4), non-rare, non-native: 
# the analyses presented in the main text Figure 5 use classify rare versus non-rare groups based 
# on the Dominance Indicator (DI) and cutoffs of:  breaks=c(0.0,0.2,0.8,1.0), where non-rare are
#both subordinate and dominant species (e.g., species with DI greater than the .2 cutoff). 

# run model of main design with the four groups of species richness: 
Mod4A.1 <- felm(log(live_mass) ~ ihs(sr_non.rare_nat) + ihs(sr_non.rare_non.nat)  + ihs(sr_non.nat_rare) +  ihs(sr_nat_rare) | newplotid + site.by.yeardummy | 0 | newplotid, data = mech.data)
summary(Mod4A.1, robust = TRUE)

## Hypothesis tests (F-tests)
# not rare: native vs non-native 
linearHypothesis(Mod4A.1, hypothesis.matrix = "ihs(sr_non.rare_nat) = ihs(sr_non.rare_non.nat)", 
                 test = "F", vcov = Mod4A.1$fevcov,  singular.ok = T)

# Native rare vs non-rare
linearHypothesis(Mod4A.1, hypothesis.matrix = "ihs(sr_nat_rare) = ihs(sr_non.rare_nat)", 
                 test = "F", vcov = Mod4A.1$fevcov,  singular.ok = T)

#non-native rare vs non-rare
linearHypothesis(Mod4A.1, hypothesis.matrix = "ihs(sr_non.nat_rare) = ihs(sr_non.rare_non.nat)", 
                 test = "F", vcov = Mod4A.1$fevcov,  singular.ok = T)

#non-native vs native rare 
linearHypothesis(Mod4A.1, hypothesis.matrix = "ihs(sr_non.nat_rare) = ihs(sr_nat_rare)", 
                 test = "F", vcov = Mod4A.1$fevcov,  singular.ok = T)


###################################################################################################################################
### Plot Figure 5 ######################################################################################################################
#####################################################################################################################################
### Plot Results - Plotting coefficient estimates from felm objects:
# https://raw.githack.com/uo-ec607/lectures/master/08-regression/08-regression.html#high_dimensional_fes_and_(multiway)_clustering

coefs_Mod4A.1 <- tidy(Mod4A.1, conf.int = T, robust = T)

# try to put all models on one line but group them
panelFE.Fig4C.data <-  bind_rows(
  coefs_Mod4A.1 %>% mutate(reg = "Richness Model"),
) 

panelFE.Fig4C.data$term = factor(panelFE.Fig4C.data$term,
                                 levels=c( "ihs(sr_nat_rare)", 
                                           "ihs(sr_non.rare_nat)",
                                           "ihs(sr_non.rare_non.nat)",
                                           "ihs(sr_non.nat_rare)"))

Fig4C <-  # ggplot(panelFE.Fig4C.data, aes(x=term, y=estimate, ymin=conf.low, ymax=conf.high, colour = term)) +
  ggplot(panelFE.Fig4C.data, aes(x=term, y=estimate, ymin=estimate - (1.96*std.error), ymax= estimate + (1.96*std.error), colour = term)) +
  geom_pointrange(size = 1.5, position = position_dodge(width = 0.5)) +
  #  geom_pointrange(aes(col = model), position = position_dodge(width = 0.5)) +
  scale_colour_discrete(name="term") +
  theme_classic() +
  labs(Title = "Marginal effect of richness on live mass") +
  geom_hline(yintercept = 0, col = "black") +
  # geom_hline(yintercept = .2, col = "grey", linetype = "dotdash") +
  ylim(-.7, .8)  +
  scale_y_continuous(name =  "Estimate for log(species richness) effect size", limits=c(-.7, .7), breaks = c(-.8, -.6, -.4, -.2, 0, .2, .4, .6, .8))
Fig4C

#Alternative y-axis label:
Fig4C <- Fig4C + labs(
  caption = "", x = "Type of Species", y = "Estimate for log(species richness) effect size")
Fig4C

#ggplot2 colors to pick from: http://sape.inf.usi.ch/quick-reference/ggplot2/colour
palette <- c("darkslateblue", "green4", "grey69", "maroon4" )

p <- Fig4C  + theme(legend.position = c(0.74, 0.77)) + scale_colour_discrete(name="term") + scale_color_manual(values=palette[c(1,2, 3,4)])   +  labs(
  title = "Effect size of Log Species Richness on Log Productivity",
  caption = "", x = "Type of Species", y = "Estimate for log(species richness) effect size") + 
  theme(legend.title=element_text(size=18), legend.text=element_text(size=18)) + 
  theme(axis.title.y= element_text(size=16)) + theme(axis.title.x= element_text(size=18))
p
# adjusting the legend 
pp <- p + theme(legend.text = element_text(size=14)) +
  theme(legend.title = element_text( size=18,  face="bold")) +
  theme(legend.title = element_blank()) +
  theme(legend.background = element_rect(# fill="lightblue", 
    size=0.5, linetype="solid",
    colour ="black"))
pp 
# adjust the title and the text size:  # https://www.datanovia.com/en/blog/ggplot-title-subtitle-and-caption/
ppp <- pp + theme(axis.text=element_text(size=22),
                  axis.title=element_text(size=20,face="bold")) +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5)) 
ppp

# print final figure:
ppp <- ppp + theme(
  legend.position="none",  #to remove the legend 
  axis.title.x = element_text(size = 20),
  axis.text.x = element_text(size = 16),
  axis.text.y = element_text(size = 16),
  axis.title.y = element_text(size = 16)) + scale_x_discrete(labels = c('ihs(Rare native)','ihs(Non-rare native)', 'ihs(Non-rare Non-native', 'ihs(Rare Non-native)')) 
ppp

#alt simplified variable labels on x-axis
Fig4C <- ppp + theme(
  legend.position="none",  #to remove the legend 
  axis.title.x = element_text(size = 20),
  axis.text.x = element_text(size = 16),
  axis.text.y = element_text(size = 16),
  axis.title.y = element_text(size = 16)) + scale_x_discrete(labels = c('Rare & Native','Non-rare & Native', 'Non-rare & Non-native', 'Rare & Non-native')) 
Fig4C

Fig4C <- Fig4C + labs(title="C. Rare Native vs. Non-rare Native vs. Non-rare non-native vs Rare non-native") +  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0)) 
#print final
Fig4C

#simplify title
Fig4C <- Fig4C + labs(title="C.") +  theme(plot.title = element_text(size = 25, face = "bold", hjust = 0.5)) 
#print final
Fig4C
#title to the left
Fig4C <- Fig4C + labs(title="C.") +  theme(plot.title = element_text(size = 25, face = "bold", hjust = -0.1)) 
#print final
Fig4C




######################################################################################################################
### 1. Native vs Non-Native Richness (in SM only) ################################################################
##################################################################################################################

##### Models with productivity measured as live mass as the Y variable ##########
#Log-log and fixed effects/dummies only.
Mod1A <- felm(log(live_mass) ~ ihs(sr_INT) + ihs(sr_NAT)  | newplotid + site.by.yeardummy | 0 | newplotid, data = mech.data) 
summary(Mod1A, robust = TRUE)

screenreg(list(Mod1A),     # object with results 
          custom.model.names= "Native vs Non-Native",
          omit.coef=c("(site_code)|(newplotid)")) 

linearHypothesis(Mod1A, hypothesis.matrix = "ihs(sr_INT) = ihs(sr_NAT)",
                 test = "F", vcov. =  Mod1A$vcv, 
                 singular.ok = T) # ignore fact that some variables are getting dropped implicitly
               #   white.adjust = "hc1") # heteroskedasticity robust F-Test 
#read about the corrected SEs in this here - https://www.rdocumentation.org/packages/car/versions/3.0-2/topics/hccm
# https://www.econometrics-with-r.org/7-3-joint-hypothesis-testing-using-the-f-statistic.html

#Controlling for evenness: #Log-log and fixed effects/dummies only.
Mod1A.1 <- felm(log(live_mass) ~ ihs(sr_INT) + ihs(sr_NAT) + ihs(even)  | newplotid + site.by.yeardummy | 0 | newplotid, data = mech.data, exactDOF='rM')
summary(Mod1A.1, robust = TRUE, cluster = TRUE) #note the cluster = TRUE isnt a real command. 

linearHypothesis(Mod1A.1, hypothesis.matrix = "ihs(sr_INT) = ihs(sr_NAT)",
                 test = "F", vcov. =  Mod1A.1$vcv, 
                 singular.ok = T) # ignore fact that some variables are getting dropped implicitly
 # white.adjust = "hc1") #

#Log-log and fixed effects/dummies only.
# seperate out sr native into caterogies? 
Mod1A.2 <- felm(log(live_mass) ~ ihs(sr_INT) + ihs(sr_nat_dom) +  ihs(sr_nat_sub) + ihs(sr_nat_rare) | newplotid + site.by.yeardummy | 0 | newplotid, data = mech.data, exactDOF='rM')
summary(Mod1A.2, robust = TRUE, cluster = TRUE)

linearHypothesis(Mod1A.2, hypothesis.matrix = "ihs(sr_INT) = ihs(sr_nat_dom)",
                 test = "F", vcov. =  Mod1A.2$vcv, 
                 singular.ok = T) # ignore fact that some variables are getting dropped implicitly

## Print Results as a table for the SI

#rich as outcome var
ModRich.1A <- felm(log(rich) ~ ihs(sr_INT) + ihs(sr_NAT)  | newplotid + site.by.yeardummy | 0 | newplotid, data = mech.data) 
summary(ModRich.1A, robust = TRUE)

screenreg(list(ModRich.1A),     # object with results 
          custom.model.names= "Native and Non-Native species effects on Richness",
          omit.coef=c("(site_code)|(newplotid)")) 

########################################################################
####### Plot Figure of these results - now in SM only ####################
#########################################################################
coefs_Mod1A <- tidy(Mod1A, conf.int = T, robust = T)
coefs_Mod1A 

Fig4A <- bind_rows(
  coefs_Mod1A %>% mutate(reg = "Native vs Non-Native"),
) %>%
#  ggplot(aes(x=term, y=estimate, ymin=conf.low, ymax=conf.high, colour = term)) +
  ggplot(aes(x=term, y=estimate, ymin=estimate - (1.96*std.error), ymax= estimate + (1.96*std.error), colour = term)) +
  geom_pointrange(size = 1.5, position = position_dodge(width = .25)) +
  scale_colour_discrete(name="term") +
  theme_classic() +
  # labs(Title = "Marginal effect of richness on live mass") +
  geom_hline(yintercept = 0, col = "black") +
  ylim(-.4, .2) +
  labs(
    title = "A. Effect size of Native vs Invasive Species Richness on Productivity",
    caption = ""  ) 
Fig4A 

#Alternative y-axis label:
Fig4A + labs( title =  "A. Effect size of Native vs Invasive Species Richness on Productivity",
  caption = "", x = "Variable", y = "Estimate for log(species richness) effect size")

invcols <- c("gray69", "saddlebrown")
inv <- Fig4A + scale_colour_discrete(name="term") + scale_color_manual(values=invcols[c(1,2)])   +  theme(legend.title=element_text(size=14), legend.text=element_text(size=12)) + theme(axis.title.y= element_text(size=18)) + theme(axis.title.x= element_text(size=18))

# adjust the title and the text size: 
# https://www.datanovia.com/en/blog/ggplot-title-subtitle-and-caption/
inv <- inv + theme(axis.text=element_text(size=22),
            axis.title=element_text(size=22,face="bold")) +
  theme(plot.title = element_text(size = 25, face = "bold", hjust = 0.5)) 
# print final figure:
inv <- inv + theme(
  legend.position="none",   #to remove the legend 
  axis.title.x = element_text(size = 20),
  axis.text.x = element_text(size = 16),
  axis.text.y = element_text(size = 16),
  axis.title.y = element_text(size = 16),
  plot.title = element_text(hjust = -0.1))  
 
#simplified plot labels
inv <-  inv +  labs(caption = "", title =  "A.", x = "Variable", y = "Estimate for log(species richness) effect size") + scale_x_discrete(labels = c('Non-native','Native')) 
inv

Fig4A <- inv + labs( title =  "A. Effect size of Native vs Invasive Species Richness on Productivity",
              caption = "", x = "Type of Species", y = "Estimate for log(species richness) effect size")
Fig4A  
# # adjusting the legend 
# inv <- inv + theme(legend.text = element_text(size=14)) +
#   theme(legend.title = element_text( size=16,  face="bold")) +
#   theme(legend.title = element_blank()) +
#   theme(legend.background = element_rect(# fill="lightblue", 
#     size=0.5, linetype="solid",
#     colour ="black"))
# inv

##################################################################################
## Model 1 Robustness Analyses for SI #############################################
###################################################################################

#control for total cover and evenness
Mod1A.2 <- felm(log(live_mass) ~ ihs(sr_INT) + ihs(sr_NAT) + ihs(even) + total_cover  | newplotid + site.by.yeardummy, data = mech.data, exactDOF='rM')
summary(Mod1A.2, robust = TRUE, cluster = TRUE)

#control for total cover only
Mod1A.3 <- felm(log(live_mass) ~ ihs(sr_INT) + ihs(sr_NAT) + total_cover  | newplotid + site.by.yeardummy, data = mech.data, exactDOF='rM')
summary(Mod1A.3, robust = TRUE, cluster = TRUE)

#B. Log-Level
Mod1B <- felm(log(live_mass) ~ sr_INT + sr_NAT | newplotid + site.by.yeardummy, data = mech.data, exactDOF='rM')
summary(Mod1B, robust = TRUE, cluster = TRUE)

Mod1Bb <- felm(log(live_mass) ~ sr_INT + sr_NAT + ihs(even) | newplotid + site.by.yeardummy, data = mech.data, exactDOF='rM')
summary(Mod1Bb, robust = TRUE, cluster = TRUE)

#C. Levels #B. Log-Level
Mod1C <- felm(live_mass ~ sr_INT + sr_NAT | newplotid + site.by.yeardummy, data = mech.data, exactDOF='rM')
summary(Mod1C, robust = TRUE, cluster = TRUE)

Mod1Cb <- felm(live_mass ~ sr_INT + sr_NAT + ihs(even) | newplotid + site.by.yeardummy, data = mech.data, exactDOF='rM')
summary(Mod1Cb, robust = TRUE, cluster = TRUE)

#print results

screenreg(list(Mod1A, Mod1A.1, Mod1B, Mod1Bb, Mod1C, Mod1Cb),     # object with results 
          #custom.model.names= "Native vs Non-Native",
          omit.coef=c("(site_code)|(newplotid)")) 

#output models into a single table
# screenreg(list(Mod1A, Mod1B ,
#                clus.res.PlotFEs.SiteYear.Log.groundPAR, clust.res.PlotFEs.SiteYear.Levels.groundPAR,
#                clus.res.PlotFEs.SiteYear.LevelLog.groundPAR),       # object with results from clx
#           custom.model.names=c("Log-Log'", "Levels'", "Level-Log'","Log-Log''", "Levels''", "Level-Log''" ),
#           omit.coef=c("(site_code)|(newplotid)"))  # object from estimation (unclustered) for BIC

# linearHypothesis(PlotFEs.SiteYear.Log.Log.Orgin, hypothesis.matrix = "ihs(sr_INT) = ihs(sr_NAT)",
#                  test = "F", 
#                  vcov. = clust.res.PlotFEs.SiteYear.Log.Log.Orgin,
#                  singular.ok = T, # ignore fact that some variables are getting dropped implicitly
#                  white.adjust = "hc1") # heteroskedasticity robust F-Test 
# #read about the corrected SEs in this here - https://www.rdocumentation.org/packages/car/versions/3.0-2/topics/hccm
# # https://www.econometrics-with-r.org/7-3-joint-hypothesis-testing-using-the-f-statistic.html

### Plot Results
#plotting coefficient estimates from felm objects:
# https://raw.githack.com/uo-ec607/lectures/master/08-regression/08-regression.html#high_dimensional_fes_and_(multiway)_clustering
coefs_Mod1A <- tidy(Mod1A, conf.int = T, robust = T)
coefs_Mod1A.1 <- tidy(Mod1A.1, conf.int = T, robust = T)
coefs_Mod1A.2 <- tidy(Mod1A.2, conf.int = T, robust = T)
coefs_Mod1A.3 <- tidy(Mod1A.3, conf.int = T, robust = T)

nvnn <- bind_rows(
   coefs_Mod1A %>% mutate(reg = "Native vs Non-Native"),
   coefs_Mod1A.1 %>% mutate(reg = "Native vs Non-Native Incld Evenness"),
  coefs_Mod1A.2 %>% mutate(reg = "Native vs Non-Native Incld Evenness & Total Cover"),
  coefs_Mod1A.3 %>% mutate(reg = "Native vs Non-Native Incld Total Cover")
) %>%
  #ggplot(aes(x=term, y=estimate, ymin=conf.low, ymax=conf.high, colour = term)) +
  ggplot(aes(x=term, y=estimate, ymin=estimate - (1.96*std.error), ymax= estimate + (1.96*std.error), colour = term)) +
  geom_pointrange(aes(col = reg), position = position_dodge(width = 0.5)) +
  scale_colour_discrete(name="Model") +
  theme_classic() +
 # labs(Title = "Marginal effect of richness on live mass") +
  geom_hline(yintercept = 0, col = "black") +
 # geom_hline(yintercept = .2, col = "grey", linetype = "dotdash") +
  ylim(-.4, .2) +
  labs(
    title = "Effect size of Native vs Invasive Species Richness on Productivitiy",
    caption = ""
  ) 

nvnn + labs(
  # title = "Effect size of Log Species Richness on Log Productivitiy",
  caption = "", x = "Variable", y = "Coefficient Estimate")

######################################################################################################################
### 2. Rare vs Subordinate vs. Dominant (in SM only) ########################################################################
##################################################################################################################
# Do analyses based on DI then based on relative abundance and frequency groups seperately (versus DI in A.)

###########
### A. Grouped based on the Dominance Indicator (DI) and cutoffs of 1
##########

##### Models with productivity as the Y variable ##########
#log-log
Mod2A.1 <- felm(log(live_mass) ~ ihs(sr_domspp) + ihs(sr_rarespp) + ihs(sr_subordspp) | newplotid + site.by.yeardummy | 0 | newplotid, data = mech.data)
summary(Mod2A.1, robust = TRUE)

#hypothesis tests
linearHypothesis(Mod2A.1, hypothesis.matrix = "ihs(sr_rarespp) = ihs(sr_domspp)", 
                 test = "F", vcov = Mod2A.1$fevcov,  singular.ok = T)
linearHypothesis(Mod2A.1, hypothesis.matrix = "ihs(sr_rarespp) = ihs(sr_subordspp)", 
                 test = "F", vcov = Mod2A.1$fevcov,  singular.ok = T)
linearHypothesis(Mod2A.1, hypothesis.matrix = "ihs(sr_subordspp) = ihs(sr_domspp)", 
                 test = "F", vcov = Mod2A.1$fevcov,  singular.ok = T)

################################################################################################
### Plot Figure 4 B ########################################################################
################################################################################################
coefs_Mod2A.1<- tidy(Mod2A.1 , conf.int = T, robust = T)

# try to put all models on one line but group them
panelFE.Fig4b.data <-  bind_rows(
  coefs_Mod2A.1 %>% mutate(reg = "Richness Model"),
) 

panelFE.Fig4b.data$term = factor(panelFE.Fig4b.data$term,
                               levels=c( "ihs(sr_rarespp)", 
                                        "ihs(sr_domspp)" ,
                                         "ihs(sr_subordspp)" ))

Fig4B <-  #ggplot(panelFE.Fig4b.data, aes(x=term, y=estimate, ymin=conf.low, ymax=conf.high, colour = term)) +
  ggplot(panelFE.Fig4b.data, aes(x=term, y=estimate, ymin=estimate - (1.96*std.error), ymax= estimate + (1.96*std.error), colour = term)) +
  geom_pointrange(size = 1.5, position = position_dodge(width = 0.5)) +
  #  geom_pointrange(aes(col = model), position = position_dodge(width = 0.5)) +
  scale_colour_discrete(name="Model") +
  theme_classic() +
  labs(Title = "Marginal effect of richness on live mass") +
  geom_hline(yintercept = 0, col = "black") +
  # geom_hline(yintercept = .2, col = "grey", linetype = "dotdash") +
  ylim(-.7, .8)  +
  scale_y_continuous(name =  "Estimate for log(species richness) effect size", limits=c(-.8, .8), breaks = c(-.8, -.6, -.4, -.2, 0, .2, .4, .6, .8))
labs(
  title = "Effect size of Log Species Richness on Log Productivitiy",
  caption = "", 
) 
Fig4B 

#Alternative y-axis label:
Fig4B < + labs(
  title = "Effect size of Log Species Richness on Log Productivitiy",
  caption = "", x = "Variable", y = "Estimate for log(species richness) effect size")

#panelFE.main.2 +  theme(legend.title=element_text(size=14), legend.text=element_text(size=14)) + theme(axis.title.y= element_text(size=18)) + theme(axis.title.x= element_text(size=18))

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "gray1")
#panelFE.main.2  + scale_color_manual(values=cbPalette[c(7,3,4,8)])   +  theme(legend.title=element_text(size=14), legend.text=element_text(size=12)) + theme(axis.title.y= element_text(size=18)) + theme(axis.title.x= element_text(size=18))

p <- Fig4B  + theme(legend.position = c(0.74, 0.77)) + scale_colour_discrete(name="Model") + scale_color_manual(values=cbPalette[c(7,9,4,8)])   +  labs(
 # title = "Effect size of Log Species Richness on Log Productivity",
  caption = "", x = "Type of Species", y = "Estimate for log(species richness) effect size") + 
  theme(legend.title=element_text(size=18), legend.text=element_text(size=18)) + 
  theme(axis.title.y= element_text(size=16)) + theme(axis.title.x= element_text(size=18))
p
# adjusting the legend 
pp <- p + theme(legend.text = element_text(size=14)) +
  theme(legend.title = element_text( size=18,  face="bold")) +
  theme(legend.title = element_blank()) +
  theme(legend.background = element_rect(# fill="lightblue", 
    size=0.5, linetype="solid",
    colour ="black"))
pp 
# adjust the title and the text size:  # https://www.datanovia.com/en/blog/ggplot-title-subtitle-and-caption/
ppp <- pp + theme(axis.text=element_text(size=22),
                  axis.title=element_text(size=20,face="bold")) +
              theme(plot.title = element_text(size = 25, face = "bold", hjust = 0.5)) 

# print final figure:
ppp + theme(
  legend.position="none",  #to remove the legend 
  axis.title.x = element_text(size = 20),
  axis.text.x = element_text(size = 16),
  axis.text.y = element_text(size = 16),
  axis.title.y = element_text(size = 16)) + scale_x_discrete(labels = c('ihs(rare species)','ihs(dominant species)', 'ihs(subordinate species)')) 

#alt simplified variable labels on x-axis
Fig4b <- ppp + theme(
  legend.position="none",  #to remove the legend 
  axis.title.x = element_text(size = 20),
  axis.text.x = element_text(size = 16),
  axis.text.y = element_text(size = 16),
  axis.title.y = element_text(size = 16)) +
  scale_x_discrete(labels = c('Rare','Dominant', 'Subordinate')) 
Fig4b 
#+  theme(plot.title = element_text(hjust = -0.45, vjust=2.12))


# do model with cut-off 2
Mod2A.DI2 <-felm(log(live_mass) ~ ihs( sr_domspp2) + ihs(sr_rarespp2) + ihs(sr_subordspp2) | newplotid + site.by.yeardummy | 0 | newplotid, data = mech.data)
summary(Mod2A.DI2, robust = TRUE)

#hypothesis tests
linearHypothesis(Mod2A.DI2, hypothesis.matrix = "ihs(sr_rarespp2) = ihs(sr_domspp2)", 
                 test = "F", vcov = Mod2A.DI2$fevcov,  singular.ok = T)
linearHypothesis(Mod2A.DI2, hypothesis.matrix = "ihs(sr_rarespp2) = ihs(sr_subordspp2)", 
                 test = "F", vcov = Mod2A.DI2$fevcov,  singular.ok = T)
linearHypothesis(Mod2A.DI2, hypothesis.matrix = "ihs(sr_subordspp2) = ihs(sr_domspp2)", 
                 test = "F", vcov = Mod2A.DI2$fevcov,  singular.ok = T)

# print results for cut off 1 and 2:
screenreg(list( Mod2A.1, Mod2A.DI2),      # object with results 
          custom.model.names= c("Cutoff 1", "Cutoff 2")) #name the models

#log-level
Mod2B.2 <- felm(log(live_mass) ~ sr_domspp + sr_rarespp + sr_subordspp | newplotid + site.by.yeardummy | 0 | newplotid, data = mech.data, exactDOF='rM')
summary(Mod2B.2, robust = TRUE, cluster = TRUE)

#levels
Mod2B.3 <- felm(live_mass ~ sr_domspp + sr_rarespp + sr_subordspp | newplotid + site.by.yeardummy | 0 | newplotid, data = mech.data, exactDOF='rM')
summary(Mod2B.3, robust = TRUE, cluster = TRUE)

screenreg(list(Mod2B.1, Mod2B.2, Mod2B.3),      # object with results 
          #custom.model.names= "Native vs Non-Native"
          ) 

##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### Models with richness as the Y variable - how do these vars affect overall richness? ##########
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
#log-log
RichMod2A.1 <- felm(log(rich) ~ ihs(sr_domspp) + ihs(sr_rarespp) + ihs(sr_subordspp) | newplotid + site.by.yeardummy | 0 | newplotid, data = mech.data, exactDOF='rM')
summary(RichMod2A.1, robust = TRUE, cluster = TRUE)

#log-log for cut off 2
RichMod2A.DI2 <-  felm(log(rich) ~ ihs(sr_domspp2) + ihs(sr_rarespp2) + ihs(sr_subordspp2) | newplotid + site.by.yeardummy | 0 | newplotid, data = mech.data, exactDOF='rM')
summary(RichMod2A.DI2, robust = TRUE, cluster = TRUE)

# print results for cut off 1 and 2:
screenreg(list(RichMod2A.DI2, RichMod2A.1 ),      # object with results 
          custom.model.names= c("Richness as Response: Cutoff 2", "Richness as Response: Cutoff 1"))

#log-levels
RichMod2A.2 <- felm(log(rich) ~ sr_domspp + sr_rarespp + sr_subordspp | newplotid + site.by.yeardummy  | 0 | newplotid, data = mech.data, exactDOF='rM')
summary(RichMod2A.2, robust = TRUE, cluster = TRUE)

#levels
RichMod2A.3 <- felm(rich ~ sr_domspp + sr_rarespp + sr_subordspp | newplotid + site.by.yeardummy  | 0 | newplotid, data = mech.data, exactDOF='rM')
summary(RichMod2A.3, robust = TRUE, cluster = TRUE)


###########
### AA. Grouped based on the Dominance Indicator (DI) and cutoffs of 2
##########
## do for Cut off 2 # sr_domspp2 sr_rarespp2 sr_subordspp2

##### Models with productivity as the Y variable ##########
#log-log
Mod2AA.1 <- felm(log(live_mass) ~ ihs(sr_domspp2) + ihs(sr_rarespp2) + ihs(sr_subordspp2) | newplotid + site.by.yeardummy | 0 | newplotid, data = mech.data, exactDOF='rM')
summary(Mod2AA.1, robust = TRUE)

#hypothesis tests
linearHypothesis(Mod2AA.1, hypothesis.matrix = "ihs(sr_rarespp) = ihs(sr_domspp)", 
                 test = "F", vcov = Mod2AA.1$fevcov,  singular.ok = T)
linearHypothesis(Mod2AA.1, hypothesis.matrix = "ihs(sr_rarespp) = ihs(sr_subordspp)", 
                 test = "F", vcov = Mod2A.1$fevcov,  singular.ok = T)
linearHypothesis(Mod2AA.1, hypothesis.matrix = "ihs(sr_subordspp) = ihs(sr_domspp)", 
                 test = "F", vcov = Mod2A.1$fevcov,  singular.ok = T)

##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### Models with richness as the Y variable - how do these vars affect overall richness? ##########
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
#log-log
RichMod2AA.1 <- felm(log(rich) ~ ihs(sr_domspp2) + ihs(sr_rarespp2) + ihs(sr_subordspp2) | newplotid + site.by.yeardummy  | 0 | newplotid, data = mech.data, exactDOF='rM')
summary(RichMod2AA.1, robust = TRUE, cluster = TRUE)

########################################################################################
### B. Grouped based on Relative Abundance in year 0 and cutoffs 1 of:  #################
#######################################################################################
## run models with groups defined based on relative abundance
# results are consistent with the results using the DI metrics for groupings of species

Mod2B.1  = felm(log(live_mass) ~ ihs(relabund_sr_domspp) + ihs(relabund_sr_rarespp) + ihs(relabund_sr_subordspp) | newplotid + site.by.yeardummy  | 0 | newplotid, data = mech.data, exactDOF='rM')
summary(Mod2B.1, robust = TRUE, cluster = TRUE)

#how does rare spp richness after SR?
RichMod2B.1 = felm(log(rich) ~ ihs(relabund_sr_domspp) + ihs(relabund_sr_rarespp) + ihs(relabund_sr_subordspp) | newplotid + site.by.yeardummy  | 0 | newplotid, data = mech.data, exactDOF='rM')
summary(RichMod2B.1, robust = TRUE, cluster = TRUE)

# do linear hypothesis tests
linearHypothesis(Mod2B.1 , hypothesis.matrix = "ihs(relabund_sr_domspp) = ihs(relabund_sr_rarespp)",
                 test = "F", 
                 vcov. =  Mod2B.1$vcv,
                 singular.ok = T) # ignore fact that some variables are getting dropped implicitly

linearHypothesis(Mod2B.1 , hypothesis.matrix = "ihs(relabund_sr_subordspp) = ihs(relabund_sr_rarespp)",
                 test = "F", 
                 vcov. =  Mod2B.1$vcv,
                 singular.ok = T) # ignore fact that some variables are getting dropped implicitly

linearHypothesis(Mod2B.1 , hypothesis.matrix = "ihs(relabund_sr_subordspp) = ihs(relabund_sr_domspp)",
                 test = "F", 
                 vcov. =  Mod2B.1$vcv,
                 singular.ok = T) # ignore fact that some variables are getting dropped implicitly

##############################################################################################################
### C. Grouped based on Relative Frequency in year 0 and cutoff 1 of: 0, .2, .8, 1
################################################################################################################
## run models with groups defined based on frequency
# results are consistent with the results using the DI metrics for groupings of species
Mod2C.1  = felm(log(live_mass) ~ ihs(freq_sr_domspp) + ihs(freq_sr_rarespp) + ihs(freq_sr_subordspp) | newplotid + site.by.yeardummy  | 0 | newplotid, data = mech.data, exactDOF='rM')
summary(Mod2C.1 , robust = TRUE, cluster = TRUE)

# Now linear hypotheses tests to assess whether the species' group effect on live biomass is equal
# Do linear hypothesis test - are the effects of rare and dominant species the same?
#linearHypothesis function in car package

linearHypothesis(Mod2C.1 , hypothesis.matrix = "ihs(freq_sr_domspp) = ihs(freq_sr_rarespp)",
                 test = "F", 
                 vcov. =  Mod2C.1$vcv,
                 singular.ok = T) # ignore fact that some variables are getting dropped implicitly

linearHypothesis(Mod2C.1, hypothesis.matrix = "ihs(freq_sr_subordspp) = ihs(freq_sr_rarespp)",
                 test = "F", 
                 vcov. = Mod2C.1$vcv,
                 singular.ok = T) # ignore fact that some variables are getting dropped implicitly

linearHypothesis(Mod2C.1, hypothesis.matrix = "ihs(freq_sr_domspp) = ihs(freq_sr_subordspp)",
                 test = "F", 
                 vcov. = Mod2C.1$vcv,
                 singular.ok = T) # ignore fact that some variables are getting dropped implicitly

## now assess how rare spp affects SR (for mechanisms tests)
#how does rare spp richness after SR?
RichMod2C.1  = felm(log(rich) ~ ihs(freq_sr_domspp) + ihs(freq_sr_rarespp) + ihs(freq_sr_subordspp) | newplotid + site.by.yeardummy,  data = mech.data)
summary(RichMod2C.1 , robust = TRUE, cluster = TRUE)

## putting in other controls
# control for N-fixer cover 
Mod2C.2 = felm(log(live_mass) ~ ihs(freq_sr_domspp) + ihs(freq_sr_rarespp) + ihs(freq_sr_subordspp) + ihs(N_fixer_cover.yr) + ihs(sr_Nfixer) | newplotid + site.by.yeardummy,  data = mech.data)
summary(Mod2C.2 , robust = TRUE, cluster = TRUE)

# are SR rare and SR N-fixing highly correlated? No. r = 0.187
cor(mech.data$sr_Nfixer, mech.data$freq_sr_rarespp)
# 0.1878382
plot(mech.data$sr_Nfixer, mech.data$freq_sr_rarespp, xlab = "SR of N-fixing species", ylab = "SR of rare species")

cor(mech.data$N_fixer_cover.yr, mech.data$freq_sr_rarespp)
plot(mech.data$N_fixer_cover.yr, mech.data$freq_sr_rarespp)

#what about with lagged N-fixers? 
plot(mech.data$lagged_sr_Nfixer, mech.data$freq_sr_rarespp)
cor(mech.data$lagged_sr_Nfixer, mech.data$freq_sr_rarespp)

######################################################################################################################
### 3. Rare vs Non-Rare (Aggregated) - in SM  #############################################################
####################################################################################################

########################################################################################
### A. Grouped based on the Dominance Indicator (DI) and cutoff 1 of:  0, .2, .8, 1
#######################################################################################

##### Models with productivity as the Y variable ##########
#log-log
Mod3A.1 <- felm(log(live_mass) ~ ihs(sr_non_rare_spp) + ihs(sr_rarespp)  | newplotid + site.by.yeardummy | 0 | newplotid, data = mech.data)
summary(Mod3A.1, robust = TRUE)

#hypothesis test 
linearHypothesis(Mod3A.1, hypothesis.matrix = "ihs(sr_rarespp) = ihs(sr_non_rare_spp)", 
                 test = "F", vcov = Mod3A.1$fevcov,  singular.ok = T)

#log-level
Mod3A.2 <- felm(log(live_mass) ~ ihs(sr_non_rare_spp) + ihs(sr_rarespp)  | newplotid + site.by.yeardummy|  0 | newplotid, data = mech.data)
summary(Mod3A.1, robust = TRUE)

##### Models with richness as the Y variable- how do these vars affect overall richness? ##########
#log-log
RichMod3A.1 <- felm(log(rich) ~ ihs(sr_non_rare_spp) + ihs(sr_rarespp)  | newplotid + site.by.yeardummy | 0 | newplotid, data = mech.data)
summary(RichMod3A.1, robust = TRUE)

#print results from richness model
screenreg(list(RichMod3A.1),      # object with results 
          custom.model.names= c("Richness model as outcome, Cut off 1"))

#print results 
screenreg(list(Mod3A.1, Mod3AA.1 ),      # object with results 
          custom.model.names= c("Cut off 1", "Cut off 2"))

###################################################################################################################################
#### Plot Fig. 4D ##################################################################################################################
####################################################################################################################################
coefs_Mod3A.1<- tidy(Mod3A.1 , conf.int = T, robust = T)

# try to put all models on one line but group them
panelFE.Fig4d.data <-  bind_rows(
  coefs_Mod3A.1 %>% mutate(reg = "Richness Model"),
) 

panelFE.Fig4d.data$term = factor(panelFE.Fig4d.data$term,
                                 levels=c( "ihs(sr_rarespp)", 
                                           "ihs(sr_non_rare_spp)" ))

Fig4d <- # ggplot(panelFE.Fig4d.data, aes(x=term, y=estimate, ymin=conf.low, ymax=conf.high, colour = term)) +
  ggplot(panelFE.Fig4d.data , aes(x=term, y=estimate, ymin=estimate - (1.96*std.error), ymax= estimate + (1.96*std.error), colour = term)) +
  geom_pointrange(size = 1.5, position = position_dodge(width = 0.5)) +
  #  geom_pointrange(aes(col = model), position = position_dodge(width = 0.5)) +
  scale_colour_discrete(name="Model") +
  theme_classic() +
  labs(Title = "Marginal effect of richness on live mass") +
  geom_hline(yintercept = 0, col = "black") +
  # geom_hline(yintercept = .2, col = "grey", linetype = "dotdash") +
  ylim(-.7, .8)  +
  scale_y_continuous(name =  "Estimate for log(species richness) effect size", limits=c(-.8, .8), breaks = c(-.8, -.6, -.4, -.2, 0, .2, .4, .6, .8))
Fig4d

#Alternative y-axis label:
Fig4d <- Fig4d + labs(
  title = "Effect size of Log Species Richness on Log Productivitiy",
  caption = "", x = "Variable", y = "Estimate for log(species richness) effect size")

#panelFE.main.2 +  theme(legend.title=element_text(size=14), legend.text=element_text(size=14)) + theme(axis.title.y= element_text(size=18)) + theme(axis.title.x= element_text(size=18))

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "gray1")
#panelFE.main.2  + scale_color_manual(values=cbPalette[c(7,3,4,8)])   +  theme(legend.title=element_text(size=14), legend.text=element_text(size=12)) + theme(axis.title.y= element_text(size=18)) + theme(axis.title.x= element_text(size=18))

#ggplot2 colors to pick from: http://sape.inf.usi.ch/quick-reference/ggplot2/colour
palette <- c("darkslateblue", "green4")
  
p <- Fig4d  + theme(legend.position = c(0.74, 0.77)) + scale_colour_discrete(name="term") + scale_color_manual(values=palette[c(1,2)])   +  labs(
  title = "Effect size of Log Species Richness on Log Productivity",
  caption = "", x = "Type of Species", y = "Estimate for log(species richness) effect size") + 
  theme(legend.title=element_text(size=18), legend.text=element_text(size=18)) + 
  theme(axis.title.y= element_text(size=16)) + theme(axis.title.x= element_text(size=18))
p
# adjusting the legend 
pp <- p + theme(legend.text = element_text(size=14)) +
  theme(legend.title = element_text( size=18,  face="bold")) +
  theme(legend.title = element_blank()) +
  theme(legend.background = element_rect(# fill="lightblue", 
    size=0.5, linetype="solid",
    colour ="black"))
pp 
# adjust the title and the text size:  # https://www.datanovia.com/en/blog/ggplot-title-subtitle-and-caption/
ppp <- pp + theme(axis.text=element_text(size=22),
                  axis.title=element_text(size=20,face="bold")) +
  theme(plot.title = element_text(size = 25, face = "bold", hjust = 0.5)) 

#simplified variable labels on x-axis
Fig4d <- ppp + theme(
  legend.position="none",  #to remove the legend 
  axis.title.x = element_text(size = 20),
  axis.text.x = element_text(size = 16),
  axis.text.y = element_text(size = 16),
  axis.title.y = element_text(size = 16)) + scale_x_discrete(labels = c('Rare species','Non-rare species')) 

Fig4d <-Fig4d + labs(title="B. Rare vs Non-Rare Richness")
Fig4d 
#+  theme(plot.title = element_text(hjust = -0.45, vjust=2.12))

########################################################################################
### AA. Grouped based on the Dominance Indicator (DI) and cutoffs 2
#######################################################################################

##### Models with productivity as the Y variable ##########
#log-log
Mod3AA.1 <- felm(log(live_mass) ~ ihs(non_rare_spp.DI2) + ihs(sr_rarespp2)  | newplotid + site.by.yeardummy | 0 | newplotid, data = mech.data)
summary(Mod3AA.1, robust = TRUE)

#hypothesis test 
linearHypothesis(Mod3AA.1, hypothesis.matrix = "ihs(sr_rarespp2) = ihs(non_rare_spp.DI2)", 
                 test = "F", vcov = Mod3AA.1$fevcov,  singular.ok = T)

#print results 
screenreg(list(Mod3A.1, Mod3AA.1 ),      # object with results 
          custom.model.names= c("Cut off 1", "Cut off 2"))

###########
### B. Grouped based on Relative Abundance in year 0 and cutoffs of:
##########
## run models with groups defined based on relative abundance
# results are consistent with the results using the DI metrics for groupings of species
Mod3B.1  = felm(log(live_mass) ~ ihs(sr_non_rare_spp.RelA) + ihs(relabund_sr_rarespp)  |newplotid + site.by.yeardummy | 0 | newplotid, data = mech.data)
summary(Mod3B.1, robust = TRUE)

#how does rare spp richness after SR?
RichMod3B.1 = felm(log(rich) ~ ihs(sr_non_rare_spp.RelA) + ihs(relabund_sr_rarespp)  |newplotid + site.by.yeardummy | 0 | newplotid, data = mech.data, exactDOF='rM')
summary(RichMod3B.1, robust = TRUE, cluster = TRUE)

# do linear hypothesis tests
linearHypothesis(Mod3B.1 , hypothesis.matrix = "ihs(sr_non_rare_spp.RelA) = ihs(relabund_sr_rarespp)",
                 test = "F", 
                 vcov. =  Mod3B.1$vcv,
                 singular.ok = T) # ignore fact that some variables are getting dropped implicitly

###########
### C. Grouped based on Relative Frequency in year 0 and cutoff 1 of:
##########
## Run models with groups defined based on frequency
# non_rare_spp.Freq
# results are consistent with the results using the DI metrics for groupings of species
Mod3C.1  = felm(log(live_mass) ~ ihs(sr_non_rare_spp.Freq) + ihs(freq_sr_rarespp) | newplotid + site.by.yeardummy | 0 | newplotid, data = mech.data, exactDOF='rM')
summary(Mod3C.1 , robust = TRUE)

# Now linear hypotheses tests to assess whether the species' group effect on live biomass is equal
# Do linear hypothesis test - are the effects of rare and dominant species the same?
linearHypothesis(Mod3C.1 , hypothesis.matrix = "ihs(sr_non_rare_spp.Freq) = ihs(freq_sr_rarespp)",
                 test = "F", 
                 vcov. =  Mod3C.1$vcv,
                 singular.ok = T) # ignore fact that some variables are getting dropped implicitly

######################################################################################################################
### 4. Rare vs Non-Rare and Native vs Invasive  #######################################################################
####################################################################################################################

###########
### A. Grouped based on the Dominance Indicator (DI) and cutoffs of:  breaks=c(0.0,0.2,0.8,1.0),
##########
# these are the variable names/groups:
# sr_non.rare_nat  + sr_non.rare_non.nat + sr_non.nat_rare +  sr_nat_rare
Mod4A.1 <- felm(log(live_mass) ~ ihs(sr_non.rare_nat) + ihs(sr_non.rare_non.nat)  + ihs(sr_non.nat_rare) +  ihs(sr_nat_rare) | newplotid + site.by.yeardummy | 0 | newplotid, data = mech.data)
summary(Mod4A.1, robust = TRUE)

## Hypothesis tests
#Note: fevcov returns a square matrix with the bias corrected covariances. An attribute 'bias' contains
# the biases. The bias corrections have been subtracted from the bias estimates. I.e. vc = vc’ - b,
# where vc’ is the biased variance and b is the bias.

# not rare: native vs non-native 
linearHypothesis(Mod4A.1, hypothesis.matrix = "ihs(sr_non.rare_nat) = ihs(sr_non.rare_non.nat)", 
                 test = "F", vcov = Mod4A.1$fevcov,  singular.ok = T)

# Native rare vs non-rare
linearHypothesis(Mod4A.1, hypothesis.matrix = "ihs(sr_nat_rare) = ihs(sr_non.rare_nat)", 
                 test = "F", vcov = Mod4A.1$fevcov,  singular.ok = T)

#non-native rare vs non-rare
linearHypothesis(Mod4A.1, hypothesis.matrix = "ihs(sr_non.nat_rare) = ihs(sr_non.rare_non.nat)", 
                 test = "F", vcov = Mod4A.1$fevcov,  singular.ok = T)

#non-native vs native rare 
linearHypothesis(Mod4A.1, hypothesis.matrix = "ihs(sr_non.nat_rare) = ihs(sr_nat_rare)", 
                 test = "F", vcov = Mod4A.1$fevcov,  singular.ok = T)

#hypothesis tests using lefm objective and a wald test:
# https://www.rdocumentation.org/packages/lfe/versions/2.8-3/topics/waldtest
#waldtest(object, R, r, type = c("default", "iid", "robust", "cluster"),
#lhs = NULL, df1, df2)
  # hypothesis.matrix = "ihs(sr_INT) = ihs(sr_NAT)"
  # R = hypothesis.matrix 
  # waldtest(Mod4A.1, ihs(sr_non.rare_nat) = ihs(sr_non.rare_non.nat),  type = c("robust", "cluster"))
  #    , lhs = NULL, df1, df2)

#### Richness as Outcome Variable 
ModRich4A.1 <- felm(log(rich) ~ ihs(sr_non.rare_nat) + ihs(sr_non.rare_non.nat)  + ihs(sr_non.nat_rare) +  ihs(sr_nat_rare) | newplotid + site.by.yeardummy | 0 | newplotid, data = mech.data)
summary(ModRich4A.1, robust = TRUE)


###################################################################################################################################
### Plot Figure 4C ######################################################################################################################
#####################################################################################################################################
### Plot Results - Plotting coefficient estimates from felm objects:
# https://raw.githack.com/uo-ec607/lectures/master/08-regression/08-regression.html#high_dimensional_fes_and_(multiway)_clustering

coefs_Mod4A.1 <- tidy(Mod4A.1, conf.int = T, robust = T)

# try to put all models on one line but group them
panelFE.Fig4C.data <-  bind_rows(
  coefs_Mod4A.1 %>% mutate(reg = "Richness Model"),
) 

panelFE.Fig4C.data$term = factor(panelFE.Fig4C.data$term,
                                 levels=c( "ihs(sr_nat_rare)", 
                                           "ihs(sr_non.rare_nat)",
                                           "ihs(sr_non.rare_non.nat)",
                                           "ihs(sr_non.nat_rare)"))

Fig4C <-  # ggplot(panelFE.Fig4C.data, aes(x=term, y=estimate, ymin=conf.low, ymax=conf.high, colour = term)) +
  ggplot(panelFE.Fig4C.data, aes(x=term, y=estimate, ymin=estimate - (1.96*std.error), ymax= estimate + (1.96*std.error), colour = term)) +
  geom_pointrange(size = 1.5, position = position_dodge(width = 0.5)) +
  #  geom_pointrange(aes(col = model), position = position_dodge(width = 0.5)) +
  scale_colour_discrete(name="term") +
  theme_classic() +
  labs(Title = "Marginal effect of richness on live mass") +
  geom_hline(yintercept = 0, col = "black") +
  # geom_hline(yintercept = .2, col = "grey", linetype = "dotdash") +
  ylim(-.7, .8)  +
  scale_y_continuous(name =  "Estimate for log(species richness) effect size", limits=c(-.7, .7), breaks = c(-.8, -.6, -.4, -.2, 0, .2, .4, .6, .8))
Fig4C

#Alternative y-axis label:
Fig4C <- Fig4C + labs(
  caption = "", x = "Type of Species", y = "Estimate for log(species richness) effect size")
Fig4C

#ggplot2 colors to pick from: http://sape.inf.usi.ch/quick-reference/ggplot2/colour
palette <- c("darkslateblue", "green4", "grey69", "maroon4" )

p <- Fig4C  + theme(legend.position = c(0.74, 0.77)) + scale_colour_discrete(name="term") + scale_color_manual(values=palette[c(1,2, 3,4)])   +  labs(
  title = "Effect size of Log Species Richness on Log Productivity",
  caption = "", x = "Type of Species", y = "Estimate for log(species richness) effect size") + 
  theme(legend.title=element_text(size=18), legend.text=element_text(size=18)) + 
  theme(axis.title.y= element_text(size=16)) + theme(axis.title.x= element_text(size=18))
p
# adjusting the legend 
pp <- p + theme(legend.text = element_text(size=14)) +
  theme(legend.title = element_text( size=18,  face="bold")) +
  theme(legend.title = element_blank()) +
  theme(legend.background = element_rect(# fill="lightblue", 
    size=0.5, linetype="solid",
    colour ="black"))
pp 
# adjust the title and the text size:  # https://www.datanovia.com/en/blog/ggplot-title-subtitle-and-caption/
ppp <- pp + theme(axis.text=element_text(size=22),
                  axis.title=element_text(size=20,face="bold")) +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5)) 
ppp

# print final figure:
ppp <- ppp + theme(
  legend.position="none",  #to remove the legend 
  axis.title.x = element_text(size = 20),
  axis.text.x = element_text(size = 16),
  axis.text.y = element_text(size = 16),
  axis.title.y = element_text(size = 16)) + scale_x_discrete(labels = c('ihs(Rare native)','ihs(Non-rare native)', 'ihs(Non-rare Non-native', 'ihs(Rare Non-native)')) 
ppp

#alt simplified variable labels on x-axis
Fig4C <- ppp + theme(
  legend.position="none",  #to remove the legend 
  axis.title.x = element_text(size = 20),
  axis.text.x = element_text(size = 16),
  axis.text.y = element_text(size = 16),
  axis.title.y = element_text(size = 16)) + scale_x_discrete(labels = c('Rare & Native','Non-rare & Native', 'Non-rare & Non-native', 'Rare & Non-native')) 
Fig4C

Fig4C <- Fig4C + labs(title="C. Rare Native vs. Non-rare Native vs. Non-rare non-native vs Rare non-native") +  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0)) 
#print final
Fig4C

#simplify title
Fig4C <- Fig4C + labs(title="C.") +  theme(plot.title = element_text(size = 25, face = "bold", hjust = 0.5)) 
#print final
Fig4C
#title to the left
Fig4C <- Fig4C + labs(title="C.") +  theme(plot.title = element_text(size = 25, face = "bold", hjust = -0.1)) 
#print final
Fig4C


#+  theme(plot.title = element_text(hjust = -0.45, vjust=2.12))

################################################################################################################################################################################################################################
### Put Fig 4 plots in same plot #####################################################################################################
#####################################################################################################
# first make the titles A, B, C, D and left align them

#A invasive vs native
#simplified plot labels
Fig4A <- inv + labs( title =  "A. Native vs Invasive Species Richness",
                     caption = "", x = "Type of Species", y = "Estimate for log(species richness) effect size") 
                    + theme(plot.title = element_text(size = 25, face = "bold", hjust = 0))

Fig4A  <- Fig4A  + labs (title =  "A.") + theme(plot.title = element_text(size = 25, face = "bold", hjust = 0))

#rare vs non rare - B
Fig4d <- Fig4d  + labs(title="B.") + theme(plot.title = element_text(size = 25, face = "bold", hjust = 0))

#  C. dom, sub, rare 
Fig4b <- Fig4b + labs(title="C.") + theme(plot.title = element_text(size = 25, face = "bold", hjust = 0))

#all groups -D
Fig4C <- Fig4C + labs(title="D.") +  theme(plot.title = element_text(size = 25, face = "bold", hjust = 0)) 

##print final
plot_grid(Fig4A, Fig4d, Fig4b, Fig4C)


#########################################
# plot only three of the plots ###########
############################################
Fig4d <- Fig4d  + labs(title="B.") + theme(plot.title = element_text(size = 25, face = "bold", hjust = 0))
Fig4C <- Fig4C + labs(title="C.") +  theme(plot.title = element_text(size = 25, face = "bold", hjust = 0)) 

#print just the Figure 4 C
Fig4C <- Fig4C + labs(title="") +   theme(plot.title = element_text(size = 25, face = "bold", hjust = 0)) 
Fig4C

#print final
plot_grid(Fig4A, Fig4d,  Fig4C)

#Messing around with plot grid to make the plots look better together:
Fig4Apg  <- Fig4A  + labs (title =  "") 
Fig4dpg <- Fig4d  + labs(title="")
Fig4Cpg <- Fig4C + labs(title="")


### FINAL:


## put the three on the same and make
common.ylim = ylim(-0.4, 0.4)
kLabelSize = 25 #Add labels with plotgrid 
common.ylab = ylab("Estimated effect of species richness")  #Estimated coefficient of species richness
plot_grid(
  plot_grid(Fig4Apg + common.ylim + common.ylab, 
            Fig4dpg + common.ylim + common.ylab, 
            labels=c("A","B"), label_size=kLabelSize),
  plot_grid(NULL, Fig4Cpg + common.ylim + common.ylab, NULL, labels = c("", "C", ""), label_size=kLabelSize, rel_widths = c(0.2, 0.6, 0.2), nrow=1), 
nrow = 2)


#tutorial - and way to add joint plot title https://wilkelab.org/cowplot/articles/plot_grid.html

##########################################################################################
#### Robustness analyses for fig 4C analyses ###################################################
##########################################################################################

#need to run some robustness checks with cover variables too
#****cover_tot_dom and cover_tot_sub didn't work in data processing
Mod4A.1 <- felm(log(live_mass) ~ ihs(sr_non.rare_nat) +  cover_tot_INT + cover_tot_NAT + ihs(sr_non.rare_non.nat)  + ihs(sr_non.nat_rare) +  ihs(sr_nat_rare) | newplotid + site.by.yeardummy, data = mech.data, exactDOF='rM')
summary(Mod4A.1, robust = TRUE, cluster = TRUE)

##### Models with richness as the Y variable- how do these vars affect overall richness? ##########
RichMod4A.1 <- felm(log(rich) ~ ihs(sr_non.rare_nat) +  cover_tot_INT + cover_tot_NAT + ihs(sr_non.rare_non.nat)  + ihs(sr_non.nat_rare) 
                    +  ihs(sr_nat_rare) | newplotid + site.by.yeardummy, data = mech.data, exactDOF='rM')
summary(RichMod4A.1, robust = TRUE, cluster = TRUE)

RichMod4A.1 <- felm(log(rich) ~ ihs(sr_non.rare_nat) + ihs(sr_non.rare_non.nat)  + ihs(sr_non.nat_rare) 
                    +  ihs(sr_nat_rare) | newplotid + site.by.yeardummy, data = mech.data, exactDOF='rM')
summary(RichMod4A.1, robust = TRUE, cluster = TRUE)

########## breakng non-rare native into SR of native dom and native subordinate spp.
Mod4A.2 <- felm(log(live_mass) ~ ihs(sr_nat_dom) + ihs(sr_nat_sub) +  ihs(sr_non.rare_non.nat)  + ihs(sr_non.nat_rare) 
                +  ihs(sr_nat_rare) | newplotid + site.by.yeardummy  | 0 | newplotid, data = mech.data, exactDOF='rM')
summary(Mod4A.2, robust = TRUE, cluster = TRUE)


###########
### B. Grouped based on Relative Abundance in year 0 and cutoffs of:  breaks=c(0.0,0.2,0.8,1.0),
##########
# sr_rare_non.nat.RelA sr_rare_nat.RelA sr_rare_non.nat.Freq sr_rare_nat.Freq

Mod4B.1 <- felm(log(live_mass) ~ ihs(sr_non.rare_nat.RelA) + ihs(sr_non.rare_non.nat.RelA)  + ihs(sr_rare_non.nat.RelA) 
                    +  ihs(sr_rare_nat.RelA) | newplotid + site.by.yeardummy | 0 | newplotid, data = mech.data)
summary(Mod4B.1, robust = TRUE)

linearHypothesis(Mod4B.1, hypothesis.matrix = "ihs(sr_non.rare_nat.RelA) = ihs(sr_non.rare_non.nat.RelA)", 
                 test = "F", vcov = Mod4B.1$fevcov,  singular.ok = T)

linearHypothesis(Mod4B.1, hypothesis.matrix = "ihs(sr_non.rare_nat.RelA) = ihs(sr_rare_nat.RelA)", 
                 test = "F", vcov = Mod4B.1$fevcov,  singular.ok = T)

## Richness as the response
RichMod4B.1 <- felm(log(rich) ~ ihs(sr_non.rare_nat.RelA) + ihs(sr_non.rare_non.nat.RelA)  + ihs(sr_rare_non.nat.RelA) 
                      +  ihs(sr_rare_nat.RelA) | newplotid + site.by.yeardummy | 0 | newplotid, data = mech.data, exactDOF='rM')
summary(RichMod4B.1, robust = TRUE, cluster = TRUE)

###########
### C. Grouped based on Relative Frequency in year 0 and cutoffs of:  breaks=c(0.0,0.2,0.8,1.0)
##########
Mod4C.1 <- felm(log(live_mass) ~ ihs(sr_non.rare_nat.Freq) + ihs(sr_non.rare_non.nat.Freq)  + ihs(sr_rare_non.nat.Freq) +  ihs(sr_rare_nat.Freq) | newplotid + site.by.yeardummy | 0 | newplotid, data = mech.data)
summary(Mod4C.1, robust = TRUE)

#make table of results
screenreg(Mod4C.1, custom.model.names= "Rel. Freq Cut off 1")

# hypothesis tests

#native vs non-native, non-rare
linearHypothesis(Mod4C.1, hypothesis.matrix = "ihs(sr_non.rare_nat.Freq) = ihs(sr_non.rare_non.nat.Freq)", 
                 test = "F", vcov = Mod4C.1$fevcov,  singular.ok = T)

#non rare native vs rare native 
linearHypothesis(Mod4C.1, hypothesis.matrix = "ihs(sr_non.rare_nat.Freq) = ihs(sr_rare_nat.Freq)", 
                 test = "F", vcov = Mod4C.1$fevcov,  singular.ok = T)

# rare native vs rare non-native
linearHypothesis(Mod4C.1, hypothesis.matrix = " ihs(sr_rare_non.nat.Freq)  = ihs(sr_rare_nat.Freq)", 
                 test = "F", vcov = Mod4C.1$fevcov,  singular.ok = T)

# non-native rare vs non-rare  
linearHypothesis(Mod4C.1, hypothesis.matrix = " ihs(sr_rare_non.nat.Freq)  = ihs(sr_non.rare_non.nat.Freq)", 
                 test = "F", vcov = Mod4C.1$fevcov,  singular.ok = T)

#richness as response
RichMod4C.1 <- felm(log(rich) ~ ihs(sr_non.rare_nat.Freq) + ihs(sr_non.rare_non.nat.Freq)  + ihs(sr_rare_non.nat.Freq) 
                    +  ihs(sr_rare_nat.Freq) | newplotid + site.by.yeardummy | 0 | newplotid, data = mech.data, exactDOF='rM')
summary(RichMod4C.1, robust = TRUE, cluster = TRUE)

### Plot Results Using the Different ways of defining rarity 

coefs_Mod4A.1 <- tidy(Mod4A.1, conf.int = T, robust = T)
coefs_Mod4B.1 <- tidy(Mod4B.1, conf.int = T, robust = T)
coefs_Mod4C.1 <- tidy(Mod4C.1, conf.int = T, robust = T)

panelMech.all <-  bind_rows(
  coefs_Mod4A.1 %>% mutate(reg = "Model Using Dominance Indicator"),
  coefs_Mod4B.1 %>% mutate(reg = "Model Using Relative Abundance"),
  coefs_Mod4C.1 %>% mutate(reg = "Model Using Frequency")
) %>%
  ggplot(aes(x=term, y=estimate, ymin=conf.low, ymax=conf.high, colour = term)) +
  geom_pointrange() + theme_classic() +
  labs(Title = "Marginal effect of richness on live mass") +
  geom_hline(yintercept = 0, col = "black") +
  geom_hline(yintercept = .2, col = "grey", linetype = "dotdash") +
  scale_colour_discrete(name="Grouping Definition") +
  ylim(-.7, .5) +
  labs(
    title = "Effect size of Log Species Richness By Group on Log Productivitiy",
    caption = "" ) + facet_wrap(~reg) 
#theme(axis.title.x = element_blank())

panelMech.all + labs(
  # title = "Effect size of Log Group Species Richness on Log Productivitiy",
  caption = "", x = "Variable", y = "Coefficient Estimate") 

#####################################################################################################################
#### 5.Run Models from 4 with different grouping cutoffs of breaks=c(0.0,0.4,0.8,1.0).
# Comparing Rare vs Non-rare, Native or Invasive ###################################################################################
#####################################################################################################################

#####################
### A. Grouped based on the Dominance Indicator (DI) and cutoffs of: breaks=c(0.0,0.4,0.8,1.0),
####################
#*for some reason this is now showing up so that the first two variables are 0..
Mod5A.1 <- felm(log(live_mass) ~ ihs(sr_non.rare_nat2) +   ihs(sr_non.rare_non.nat2)  + ihs(sr_non.nat_rare2)  +  ihs(sr_nat_rare2) | newplotid + site.by.yeardummy | 0 | newplotid, data = mech.data)
summary(Mod5A.1, robust = TRUE)

#make table of results
screenreg(Mod5A.1, custom.model.names= "Dominance Indicator: Cut off 2")
#make table of results
screenreg(Mod4A.1, custom.model.names= "Dominance Indicator: cut off 1")
#put both in one table
screenreg(list(Mod4A.1, Mod5A.1), custom.model.names= c("Dominance Indicator: cut off 1","Dominance Indicator: Cut off 2"))

#hypothesis tests
#non-rare native vs non-native 
linearHypothesis(Mod5A.1, hypothesis.matrix = "ihs(sr_non.rare_nat2) = ihs(sr_non.rare_non.nat2)", 
                 test = "F", vcov = Mod5A.1$fevcov,  singular.ok = T)

#native rare vs not rare 
linearHypothesis(Mod5A.1, hypothesis.matrix = "ihs(sr_nat_rare2) = ihs(sr_non.rare_nat2)", 
                 test = "F", vcov = Mod5A.1$fevcov,  singular.ok = T)

#non-native: rare vs non rare 
linearHypothesis(Mod5A.1, hypothesis.matrix = "ihs(sr_non.nat_rare2) = ihs(sr_non.rare_non.nat2)", 
                 test = "F", vcov = Mod5A.1$fevcov,  singular.ok = T)

#rare native vs rare non native
linearHypothesis(Mod5A.1, hypothesis.matrix = "ihs(sr_non.nat_rare2) = ihs(sr_nat_rare2)", 
                 test = "F", vcov = Mod5A.1$fevcov,  singular.ok = T)

###########
### B. Grouped based on Relative Abundance in year 0 and cutoffs of:  breaks=c(0.0,0.4,0.8,1.0),
##########
Mod5B.1 <- felm(log(live_mass) ~ ihs(sr_non.rare_nat.RelA2) + ihs(sr_non.rare_non.nat.RelA2)  + ihs(sr_rare_non.nat.RelA2) 
                +  ihs(sr_rare_nat.RelA2) | newplotid + site.by.yeardummy | 0 | newplotid, data = mech.data)
summary(Mod5B.1, robust = TRUE)

linearHypothesis(Mod5B.1, hypothesis.matrix = "ihs(sr_non.rare_nat.RelA2) = ihs(sr_non.rare_non.nat.RelA2)", 
                 test = "F", vcov = Mod5B.1$fevcov,  singular.ok = T)

linearHypothesis(Mod5B.1, hypothesis.matrix = "ihs(sr_non.rare_nat.RelA2) = ihs(sr_rare_nat.RelA2)", 
                 test = "F", vcov = Mod5B.1$fevcov,  singular.ok = T)

# also control for Nfixers
Mod5B.2 <- felm(log(live_mass) ~ ihs(sr_Nfixer) + ihs(sr_non.rare_nat.RelA2) + ihs(sr_non.rare_non.nat.RelA2)  + ihs(sr_rare_non.nat.RelA2) 
                +  ihs(sr_rare_nat.RelA2) | newplotid + site.by.yeardummy, data = mech.data, exactDOF='rM')
summary(Mod5B.2, robust = TRUE, cluster = TRUE)

#Richness as the response
RichMod5B.1 <- felm(log(rich) ~ ihs(sr_non.rare_nat.RelA2) + ihs(sr_non.rare_non.nat.RelA2)  + ihs(sr_rare_non.nat.RelA2) 
                    +  ihs(sr_rare_nat.RelA2) | newplotid + site.by.yeardummy, data = mech.data, exactDOF='rM')
summary(RichMod5B.1, robust = TRUE, cluster = TRUE)

###########
### C. Grouped based on Relative Frequency in year 0 and cutoffs of:  breaks=c(0.0,0.4,0.8,1.0)
##########
Mod5C.1 <- felm(log(live_mass) ~  ihs(sr_non.rare_nat.Freq2) + ihs(sr_non.rare_non.nat.Freq2)  + ihs(sr_rare_non.nat.Freq2) 
                +  ihs(sr_rare_nat.Freq2) | newplotid + site.by.yeardummy | 0 | newplotid, data = mech.data)
summary(Mod5C.1, robust = TRUE)

# hypothesis tests
linearHypothesis(Mod5C.1, hypothesis.matrix = "ihs(sr_non.rare_nat.Freq2) = ihs(sr_non.rare_non.nat.Freq2)", 
                 test = "F", vcov = Mod5C.1$fevcov,  singular.ok = T)

linearHypothesis(Mod5C.1, hypothesis.matrix = "ihs(sr_non.rare_nat.Freq2) = ihs(sr_rare_nat.Freq2)", 
                 test = "F", vcov = Mod5C.1$fevcov,  singular.ok = T)

# also control for Nfixers
Mod5C.2 <- felm(log(live_mass) ~ ihs(sr_Nfixer) + ihs(sr_non.rare_nat.Freq2) + ihs(sr_non.rare_non.nat.Freq2)  + ihs(sr_rare_non.nat.Freq2) 
                +  ihs(sr_rare_nat.Freq2) | newplotid + site.by.yeardummy, data = mech.data, exactDOF='rM')
summary(Mod5C.2, robust = TRUE, cluster = TRUE)

#richness as response
RichMod5C.1 <- felm(log(rich) ~ ihs(sr_non.rare_nat.Freq2) + ihs(sr_non.rare_non.nat.Freq2)  + ihs(sr_rare_non.nat.Freq2) 
                    +  ihs(sr_rare_nat.Freq2) | newplotid + site.by.yeardummy, data = mech.data, exactDOF='rM')
summary(RichMod5C.1, robust = TRUE, cluster = TRUE)

### Plot Results Using the Different ways of defining rarity 
coefs_Mod5A.1 <- tidy(Mod5A.1, conf.int = T) 
#%>% as_tibble(rownames = "ihs(sr_non.rare_nat2)") -> Non_Rare.Native
coefs_Mod5B.1 <- tidy(Mod5B.1, conf.int = T)
coefs_Mod5C.1 <- tidy(Mod5C.1, conf.int = T)

panelMech.all_cutoff2 <-  bind_rows(
  coefs_Mod5A.1 %>% mutate(reg = "Model Using Dominance Indicator"),
 coefs_Mod5B.1 %>% mutate(reg = "Model Using Relative Abundance"),
 coefs_Mod5C.1 %>% mutate(reg = "Model Using Frequency")
) %>%
  ggplot(aes(x=term, y=estimate, ymin=conf.low, ymax=conf.high, colour = term)) +
  geom_pointrange() + theme_classic() +
  labs(Title = "Marginal effect of richness on live mass") +
  geom_hline(yintercept = 0, col = "black") +
  geom_hline(yintercept = .2, col = "grey", linetype = "dotdash") +
  scale_colour_discrete(name="Grouping Definition") +
  ylim(-.7, .5) +
  labs(
    title = "Effect size of Log Species Richness By Group on Log Productivitiy",
    caption = "" ) + facet_wrap(~reg) 
#theme(axis.title.x = element_blank())

panelMech.all_cutoff2 + labs(
  # title = "Effect size of Log Group Species Richness on Log Productivitiy",
  caption = "", x = "Variable", y = "Coefficient Estimate") 

# save plots my_plot <- stops_facet_plot 
ggsave("Fig4C", plot = panelMech.all_cutoff2, width=NA, height=NA)

#make table of results
screenreg(list(Mod5A.1, Mod5B.1, Mod5C.1),     # object with results 
          custom.model.names= c("Dominance Indicator",  "Relative abundance", "Frequency"))
          # omit.coef=c("(site_code)|(newplotid)")) 

# Print Results Comparing The Cut-Offs 
#DI
screenreg(list(Mod4A.1, Mod5A.1),     # object with results 
          custom.model.names= c("Cut-Off 1", "Cut-Off 2"))

# Print Results Comparing The Cut-Offs -- for relative abundance. 
screenreg(list(Mod4B.1, Mod5B.1),     # object with results 
          custom.model.names= c("Cut-Off 1", "Cut-Off 2"))

# Print Results Comparing The Cut-Offs 
screenreg(list(Mod4C.1, Mod5C.1),     # object with results 
          custom.model.names= c("Cut-Off 1", "Cut-Off 2"))

############################################################################################################################
#### Plot Correlations between all of the SR groupings ###################################################################
#############################################################################################################################
richnessvars.fig4c <- c("sr_non.rare_nat", "sr_non.rare_non.nat", "sr_non.nat_rare", "sr_nat_rare")
richnessvars.to.plot <- mech.data[,richnessvars.fig4c, with=F] #with = F will then use the data.frame conventions, using "" for col names
pairs(richnessvars.to.plot)

######################################################################################################################
#### Visualize Correlations between SR vars ##########################################################################################################
######################################################################################################################
library(corrplot)
sr.metrics <- mech.data[, .(sr_non.nat_rare, sr_nat_rare, non_rare_spp,
                       sr_non.rare_non.nat, sr_non.rare_nat, sr_nat_dom,
                     sr_non.nat_dom, sr_nat_sub, sr_non.nat_sub, cover_tot_non.rare )]
 cor(sr.metrics)
corrplot(cor(sr.metrics), method = "square", tl.cex = .5)

