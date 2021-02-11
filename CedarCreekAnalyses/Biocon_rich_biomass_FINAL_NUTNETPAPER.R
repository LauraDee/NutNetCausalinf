#Realized richness analyses of BioCon
# Laura update Dec 16 2020

rm(list=ls())
#Close graphics windows
graphics.off()

#biocon data downloaded in Nov? 2020
library(tidyr) 
library(here)
library(data.table)
library(lfe)  # to run econometric fixed effect models 
require(ggplot2)
library(texreg)

# When log(0) need to use inverse hyperbolic sine transformation (REF NEEDED)
#https://en.wikipedia.org/wiki/Inverse_hyperbolic_functions#Inverse_hyperbolic_sine
ihs = function(x) {
  return(log(x + sqrt(x^2+1)))
} # v.s. log(x+1) <- is defined for a negative x. 
## it is generally frowned upon, when taking the log, to add small numbers (e.g., log(richness + 0.01))

### LOAD DATA ##########
setwd("~/Documents/GitHub/Causality_BioCON/data/")

#processed by Kaitlin on Dec 16 2020; KK used the biomass files for richness and biomass to be consistent with BigBio.
# The TreatmentSR column = treatment level of spp
# The ObservedSR = observed species based on what was planted in the plot
# Live.mass = biomass
findat <-  fread( "biocon_plantedbiomass_output.csv")

################################################################################
## Process data to prep for models ### #########################################
###############################################################################
findat$Plot <- as.factor(findat$Plot)
findat$Year <- as.factor(findat$Year)

#whats the sample size
table(findat$Year)

#check for any NAs or 0s for Biomass or SR
table(findat$ObservedSR)  
table(findat$TreatmentSR)
summary(findat$live.mass)

hist(findat$ObservedSR)
hist(findat$TreatmentSR)

#Check that we have created a VARIABLE FOR  OF THE SPECIES PLANTED,
# WHICH AND HOW MANY COME UP AND ARE MEASURED IN PLANTED BIOMASS
table(findat[TreatmentSR=="1", ObservedSR])

############################################################################################################
## Plot Data ###############################################################################################
############################################################################################################
## plot the raw realized richness over time, at each richness level - color the points by planted richness level
plot(findat$Year, findat$ObservedSR)
plot(findat$TreatmentSR, findat$ObservedSR)

realrichovertime <- ggplot(findat, aes(x = Year, y = ObservedSR)) + facet_wrap(~TreatmentSR) +
  # geom_smooth(method="lm", se=T) +
  theme_bw() +
  geom_point()+  labs(y = "Realized Planted Species Richness over time") + 
  theme(axis.title.y= element_text(size=14)) + theme(axis.title.x= element_text(size=14)) +
  theme(axis.text.y = element_text(size = 14)) 
realrichovertime 

## plot realized richness over time, for each treatment - color the points by planted richness level
realrichovertime.T <- ggplot(findat, aes(x = Year, y = ObservedSR, colour= TreatmentSR))   + # + facet_wrap(~Trt) +
 #  geom_smooth(method="lm", se=T) +
  theme_bw() +
  geom_point() +  labs(y = "Realized Planted Species Richness over time") + 
  theme(axis.title.y= element_text(size=14)) + theme(axis.title.x= element_text(size=14)) +
  theme(axis.text.y = element_text(size = 14)) 
realrichovertime.T

###########################################################################################################################
### Create copies of the data tables and subset to particular planted diversity levels ###################################
###########################################################################################################################
findat9 = findat[TreatmentSR == "9",]
findat4 = findat[TreatmentSR == "4",]
findat1 = findat[TreatmentSR == "1",]
table(findat9$Year) 
## Forest says to run the analyses on only the 16 species plot since they are more representative, 
# but the sample size is insufficient. Insufficient power for precise estimate
findat16 = findat[TreatmentSR == "16",]
table(findat16$Year)
summary(findat16$ObservedSR)

############################################################################################################################
### Do Fixed Effects Analyses using lfe package ##########################################################################################
##########################################################################################

#ihs transformed because some of the live.mass is 0. Though how is this possible?
mod1 <- felm(ihs(live.mass) ~ ihs(ObservedSR) | Plot  + Year | 0 | Plot,  data = findat)
summary(mod1)

screenreg(mod1,     # object with results 
          custom.model.names= "Main Model for BioCon: Effect on Planted Biomass" ,
          override.se=summary(mod1)$coef[,2],
          override.pval=summary(mod1)$coef[,4])

#do the model in levels (un-transformed variables)
mod2 <- felm(live.mass ~ ObservedSR | Plot  + Year | 0 | Plot,  data = findat)
summary(mod2)

##########################################################################################
## Estimate the treatment effect of the experiment, using the lfe package 
 # and controlling for year effects
####################################################################################

findat$SR.level <- as.factor(findat$TreatmentSR)

#as factor
mod1.treat <- felm(ihs(live.mass) ~ SR.level | Year | 0 | 0,  data = findat)
summary(mod1.treat)

# as continuous (though it should be as a factor)


#### Run models on subset of the treatment levels: starting with the most diverse plots
mod.16 <- felm(ihs(live.mass) ~ ihs(ObservedSR) | Plot  + Year | 0 | Plot,  data = findat16)
summary(mod.16)

mod.1 <- felm(ihs(live.mass) ~ ihs(ObservedSR) | Plot  + Year | 0 | Plot,  data = findat1)
summary(mod.1)

mod.4 <- felm(ihs(live.mass) ~ ihs(ObservedSR) | Plot  + Year | 0 | Plot,  data = findat4)
summary(mod.4)

mod.9 <- felm(ihs(live.mass) ~ ihs(ObservedSR) | Plot  + Year | 0 | Plot,  data = findat9)
summary(mod.9)



#older - 
################################################################################################################
### Process Data to mimic Econometric Fixed Effect Approach: First DIfference and Substract Group Means  #######
################################################################################################################
# Rationale: 
# Baseline levels differ by plot 
# Variation through time after removing important differences.

# changerich and changelive_mass are first differences in rich # and live_mass at the plot level. 

# My NutNet anaysis has site, plot, year, so in that model, the Preferred models use plot FE and
# site x year effects, so take the changes and demean by plots within same site and year. In contrast here we only have plots and years.


### demean by plot - what we are doing with the plot fixed effects (i think we will want to do this for year too!)
findat[,dm.log.RR  := l.RR -mean(l.RR, na.rm=T), by=.(Plot)]
findat[, dm.log.totbio := l.biomass - mean(l.biomass, na.rm = T), by = .(Plot)]

##***** Should we also be de-meaning by Year seperately - cant put in the above line bc it will make it a plot by year mean so need to do seperately ******
##*

realrich.biomass.amb <- ggplot(findat[Trt == "Namb.Camb",], aes(x = dm.log.RR, y = dm.log.totbio, colour= TreatmentSR)) + facet_wrap(~Trt) +
  geom_smooth(method="lm", se=T) +
  theme_bw() +
  geom_point() +
  labs(y = "De-meaned Tot  Tot Biomass per Year", x = "de-meaned (plot mean) realized richness") + 
  theme(axis.title.y= element_text(size=14)) + theme(axis.title.x= element_text(size=14)) +
  theme(axis.text.y = element_text(size = 14)) 
realrich.biomass.amb

# we should only be looking within treatment or de-meaning from treartment too But:
realrich.biomass <- ggplot(findat, aes(x = dm.log.RR, y = dm.log.totbio, colour= TreatmentSR)) + facet_wrap(~Trt) +
  geom_smooth(method="lm", se=T) +
  theme_bw() +
  geom_point() +
  labs(y = "De-meaned Tot Biomass per Year") + 
  theme(axis.title.y= element_text(size=14)) + theme(axis.title.x= element_text(size=14)) +
  theme(axis.text.y = element_text(size = 14)) 
realrich.biomass



### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### Do a first difference -- i.e., the change. ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

#as for changes you can use the shift function in data.table; but things need to be sorted first
# if you want changes, rows need to be sorted by time. Sort your data table by whatever the time column is. then use the shift function.
#if you already computed laggedrich, you don't need the by part or the order part. 
# they are only needed for constructing the lags   #dat[, changerich := rich-laggedrich,]
findat[order(Year), changerealrich := ObservedSR - shift(ObservedSR), by =.(Plot)]
# do the same for biomass (untransformed)
findat[order(Year), changetotbio := TotalBiomass -shift(TotalBiomass), by =.(Plot)]

# most year-to-year changes in realized richness in the most diversity plots
plot(findat$TreatmentSR, findat$changerealrich)

## Make a lag of the differences ##
findat[order(Year), lagged_changetotbio := shift(changetotbio ), by =.(Plot)]
findat[order(Year), lagged_changerealrich:= shift(changerealrich), by =.(Plot)]

## other data processing examples from the NutNet for reference

# We'll also do (double) demeaning of logged values to mimic 
# what's done in log-log model
# comb[,dm.changerich:=changerich-mean(changerich, na.rm=T), by=.(site,year)]
# comb[,dm.changelive_mass:=changelive_mass-mean(changelive_mass, na.rm=T), by=.(site,year)]
# 
# comb[,`:=`(log.rich=log(rich), log.live_mass=log(live_mass))]
# comb[order(year), change.log.rich := log(rich)-shift(log(rich)), by =.(plot, site_code)]
# comb[order(year), change.log.live_mass := log(live_mass)-shift(log(live_mass)), by =.(plot, site_code)]
# comb[,dm.change.log.rich:=change.log.rich-mean(change.log.rich, na.rm=T), by=.(site,year)]
# comb[,dm.change.log.live_mass:=change.log.live_mass-mean(change.log.live_mass, na.rm=T), by=.(site,year)]

######## PLOT ########
plot(findat$changetotbio ~ findat$changerealrich)

## plot realized richness over time, for each treatment - color the points by planted richness level
realrich.biomass.change <- ggplot(findat, aes(x = changerealrich, y = changetotbio, colour= TreatmentSR)) + facet_wrap(~Trt) +
  geom_smooth(method="lm", se=T) +
  theme_bw() +
  geom_point() +
  labs(y = "Change in Tot Biomass per YEar") + 
  theme(axis.title.y= element_text(size=14)) + theme(axis.title.x= element_text(size=14)) +
  theme(axis.text.y = element_text(size = 14)) 
realrich.biomass.change 

realrich.biomass.change2 <- ggplot(findat, aes(x = changerealrich, y = changetotbio, colour= Trt)) + facet_wrap(~TreatmentSR) +
  geom_smooth(method="lm", se=T) +
  theme_bw() +
  geom_point() +
  labs(y = "Change in Tot Biomass per YEar") + 
  theme(axis.title.y= element_text(size=14)) + theme(axis.title.x= element_text(size=14)) +
  theme(axis.text.y = element_text(size = 14)) 
realrich.biomass.change2 


realrich.biomass.change.amb <- ggplot(findat[Trt == "Namb.Camb",], aes(x = changerealrich, y = changetotbio, colour= Trt)) + facet_wrap(~TreatmentSR) +
  geom_smooth(method="lm", se=T) +
  theme_bw() +
  geom_point() +
  labs(y = "Change in Tot Biomass per YEar") + 
  theme(axis.title.y= element_text(size=14)) + theme(axis.title.x= element_text(size=14)) +
  theme(axis.text.y = element_text(size = 14)) 
realrich.biomass.change.amb


## change models 
mod2 <- felm(ihs(changetotbio) ~ ihs(changerealrich) + TreatmentSR | 0 | 0 | Plot, data = findat)
summary(mod2)


## Older

# options - get rid of species with 0 biomass or add 0.01 to values of 0. 
# I'm currently adding 0.01 so we can take the log
# findat$Biomass0.01 <- findat$Aboveground.Biomass..g.m.2. + 0.01
nat.non <- lm(log(Biomass0.01) ~ Native_status + Year, data = findat)
summary(nat.non)
nat.non1 <- lm(log(Biomass0.01) ~ poanat + Year, data = tbnat)
summary(nat.non1)
