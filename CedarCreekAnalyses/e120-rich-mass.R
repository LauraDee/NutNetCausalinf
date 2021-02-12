# Two way fixed effect panel models applied to Big Bio Data
# Laura Dee Dec 3 2020, 
# modifying ode by Eric Seabloom sent Nov 19 2020 that processed the data, made plots, and applied random effects models

#Clear all existing data
rm(list=ls())
#Close graphics windows
graphics.off()

library(ggplot2)
library(lme4)
library(lmerTest)
library(data.table)
library(lfe)

# When log(0) need to use inverse hyperbolic sine transformation (Bellemare & Wichman 2020 Oxford Bulletin of Economics and Statistics)
#https://en.wikipedia.org/wiki/Inverse_hyperbolic_functions#Inverse_hyperbolic_sine
ihs = function(x) {
  return(log(x + sqrt(x^2+1)))
} # v.s. log(x+1) <- is defined for a negative x. 


#######    READ IN DATA    #######
setwd("~/Documents/GitHub/Causality_BioCON/bigbio/")

#d<-read.csv("e120mega96-19.csv", strip.white=T)
#from eric: 
  # d$biomass <- d$TotalBiomassPlantedSpNoOak.g.m2.
  # d$biomass.lg <- log(d$biomass)
  # d$s.obs <- d$SrPlantedSp.Biomass.
  # d$d.obs <- d$ShanWPlantedSpNoOak.Biomass.

# read in file that Kaitlin Kimmel processed using the R code e120-biomass-data-2020-12-04.R 
# KK: e120-plantedbiomass-data-output.csv is the one you want to use because it only includes planted species.
d <-fread("e120-plantedbiomass-data-output.csv", strip.white=T)   

# numsp is the treatment and
# rich is the observed richness of planted species and
# mass.live is the biomass of the observed planted species

#check to make sure rich is the observed richness of planted species for THAT PLOT vs the whole experiment
table(d[numsp=="1", rich])
table(d$rich)

#whats the sample size
table(d$year)

#planted richness variable SpNum 
class(d$SpNum)
d$SpNum <- as.factor(d$SpNum)

# Delete missing data
d <- d[!is.na(d$biomass),]
summary(d)

#also make a variable for SpNum as numeric 
# d$SpNum_numeric <- as.integer(d$SpNum) # this command seemed to only make it 11 spp. so i did the following line:
d$SpNum_numeric <- as.integer(d$NumSp)

#Eric - Random Effect Models
# code by Eric Seabloom sent Nov 19 2020 
# 1. Effects of observed richness on biomass in all plots
summary(lme1<-lmer(biomass.lg~s.obs + (1|Plot) +(1|Year), data=d))

# 2. Effects of observed richness on biomass in plots with 8 or more species
summary(lme1<-lmer(biomass.lg~s.obs + (1|Plot) +(1|Year), data=d[d$SpNum >= 8,]))    ## *** ERIC NOT SURE HOW YOU RAN THIS BECAUSE SpNum is a factor??? **

# 3. Effects of observed richness on biomass in 16 specie splots
summary(lme1<-lmer(biomass.lg~s.obs + (1|Plot) +(1|Year), data=d[d$SpNum >= 16,]))  ## *** ERIC NOT SURE HOW YOU RAN THIS BECAUSE SpNum is a factor??? **

#########################################################################################################################
## Laura adding these models but with two-way fixed effects design #####################################################
########################################################################################################################
plot(d$NumSp, d$s.obs, xlab = "planted species", ylab = "realized richness", main = "BigBio species richness")
#1.***** do the same as above but log richness -- have a econometric fixed and plot effects  ***** This is the one to focus on***
mod_main <- felm(log(biomass) ~ log(s.obs) | Year + Plot | 0 | Plot, data =d)
summary(mod_main)

#other analysis, testing robustness to modeling decisions about which fixed effects, data, or clusters to include:
# do with just the year effect
mod2 <- felm(log(biomass) ~ log(s.obs) | Year  | 0 | Plot, data =d)
summary(mod2)
#do with year effect and plot effect but don't use cluster robust SEs 
mod3 <- felm(log(biomass) ~ log(s.obs) | Year + Plot  | 0 | 0, data =d)
summary(mod3)
# same as above but year fixed effect only 
mod4 <- felm(log(biomass) ~ log(s.obs) | Year   | 0 | 0, data =d)
summary(mod4)

#2. in species with 8 or more 
mod1_over8 <- felm(log(biomass) ~ log(s.obs) | Year + Plot | 0 | Plot, data =d[d$SpNum_numeric >= 8,])
summary(mod1_over8)

#3. in species with 16 or more 
mod1_over16 <- felm(log(biomass) ~ log(s.obs) | Year + Plot | 0 | Plot, data =d[d$SpNum_numeric >= 16,])
summary(mod1_over16)

### do within SpNum - within treatment, with the the treatment as the fixed effect:
withinmod <- felm(log(biomass) ~ log(s.obs) |  SpNum + Year | 0 | Plot, data =d)
summary(withinmod)

###############################################################################
## Compare to analyzing as an experimental treatment of planted richness #####
##############################################################################

#### 1. FIRST RUN ANALYSES WITH SR AS NUMERIC, no fixed effects of SE clustering
# SR as a integer, both transformed on on a log-log scale - no fixed effects or clustering of the standard error
treat.effect <- felm(log(biomass + .01) ~ log(SpNum_numeric) | 0 | 0 | 0, data =d)
summary(treat.effect)
# instead, transformed using an ihs (inverse hyperbolic sign approximation of a log)  -- no fixed effects or clustering of the error
treat.effect2 <- felm(ihs(biomass + .01) ~ ihs(SpNum_numeric) | 0 | 0 | 0, data =d)
summary(treat.effect2)
#no transformation:
treat.effect3 <- felm(biomass  ~ SpNum_numeric | 0 | 0 | 0, data =d)
summary(treat.effect3)

#### 2. Second, RUN ANALYSES WITH SR AS NUMERIC, no fixed effects, but cluster SE on the plot
treat.effect.cluster <- felm(log(biomass + .01) ~ log(SpNum_numeric) | 0 | 0 | Plot, data =d)
summary(treat.effect.cluster)
#ihs transformaton:
treat.effect.cluster2 <- felm(ihs(biomass + .01) ~ ihs(SpNum_numeric) | 0 | 0 | Plot, data =d)
summary(treat.effect.cluster2)
#untransformed
treat.effect.cluster3 <- felm(biomass  ~ SpNum_numeric | 0 | 0 | Plot, data =d)
summary(treat.effect.cluster3)

##### 3. Add year as a fixed effect ##### It's unchanged:
treat.effect.cluster.yr <- felm(log(biomass + .01) ~ log(SpNum_numeric) | Year | 0 | Plot, data =d)
summary(treat.effect.cluster.yr)
#ihs transformaton:
treat.effect.cluster2.yr <- felm(ihs(biomass + .01) ~ ihs(SpNum_numeric) | Year | 0 | Plot, data =d)
summary(treat.effect.cluster2.yr)
#untransformed
treat.effect.cluster3.yr <- felm(biomass  ~ SpNum_numeric | Year | 0 | Plot, data =d)
summary(treat.effect.cluster3.yr)

###########################################################################
### RUN THE ABOVE WITH LM for comparison -- ##############################
########################################################################
treat.effectlm <- lm(log(biomass + 0.01) ~ log(SpNum_numeric), data =d)
summary(treat.effectlm)
treat.effectlm1 <- lm(ihs(biomass) ~ ihs(SpNum_numeric), data = d)
summary(treat.effectlm1)
treat.effectlm2 <- lm(biomass ~ SpNum_numeric, data =d)
summary(treat.effectlm2)

#add year as fixed effect in lm:
treat.effectlm.yr <- lm(log(biomass + 0.01) ~ log(SpNum_numeric) + Year, data =d)
summary(treat.effectlm.yr)
treat.effectlm1.yr <- lm(ihs(biomass) ~ ihs(SpNum_numeric) +Year, data = d)
summary(treat.effectlm1.yr)
treat.effectlm2.yr <- lm(biomass ~ SpNum_numeric +Year, data =d)
summary(treat.effectlm2.yr)

####################################################################################
### RUN models estimating the treatment effect with Planted Richness as a factor ######
#######################################################################################

#*** Not sure why there is a planted 17 ??

#no fixed effects or clustering
treat.effect.factor <- felm(log(biomass +0.01) ~ SpNum | 0 | 0 | 0, data =d)
summary(treat.effect.factor)
#ihs transformed
treat.effect.factor <- felm(ihs(biomass) ~ SpNum | 0 | 0 | 0, data =d)
summary(treat.effect.factor)
#untransformed
treat.effect.factor <- felm(biomass ~ SpNum | 0 | 0 | 0, data =d)
summary(treat.effect.factor)

#add clustering of SE on the plot level as the cluster
#as a factor, log-levels 
treat.effect <- felm(log(biomass +0.01) ~ SpNum |  0 | 0 | Plot, data =d)
summary(treat.effect)
#as a factor, ihs transformed biomass 
treat.effect <- felm(ihs(biomass) ~ SpNum |  0 | 0 | Plot, data =d)
summary(treat.effect)
#as a factor, levels-levels (untransformed ) 
treat.effect.untransformed <- felm(biomass ~ SpNum |  0 | 0 | Plot, data =d)
summary(treat.effect.untransformed)

## Add year as a fixed effect - and clustered SEs
  #as a factor, log-levels with a year fixed effect
treat.effect <- felm(log(biomass +0.01) ~ SpNum |  Year | 0 | Plot, data =d)
summary(treat.effect)
  #as a factor, ihs transformed biomass with a year fixed effect
treat.effect <- felm(ihs(biomass) ~ SpNum |  Year | 0 | Plot, data =d)
summary(treat.effect)
  #as a factor, levels-levels (untransformed ) with a year fixed effect
treat.effect.untransformed <- felm(biomass ~ SpNum |  Year | 0 | Plot, data =d)
summary(treat.effect.untransformed)

#######################################################################################################################################
###### Run the same models as above with treatment of planted richness levels (factor) as a regular LM 
##########################################################################################################
treat.effectlm <- lm(log(biomass + 0.01) ~ SpNum, data =d)
summary(treat.effectlm)
treat.effectlm1 <- lm(ihs(biomass) ~ SpNum, data = d)
summary(treat.effectlm1)
treat.effectlm2 <- lm(biomass ~ SpNum, data =d)
summary(treat.effectlm2)

#add year as fixed effect in lm:
treat.effectlm.yr <- lm(log(biomass + 0.01) ~ SpNum + Year, data =d)
summary(treat.effectlm.yr)
treat.effectlm1.yr <- lm(ihs(biomass) ~ SpNum + Year, data = d)
summary(treat.effectlm1.yr)
treat.effectlm2.yr <- lm(biomass ~ SpNum + Year, data =d)
summary(treat.effectlm2.yr)



# Eric's plots of data
pdf("e120-rich-mass.pdf", height=6, width = 9)
ggplot(aes(x=Year, y=s.obs, pch=as.factor(Plot), col=as.factor(SpNum)), data=d) +
scale_color_brewer(palette="Spectral", name="Planted Richness") +
geom_smooth(se=FALSE)+ labs(y="Species Richness in Biomass Sample")

ggplot(aes(x=Year, y=d.obs, pch=as.factor(Plot), col=as.factor(SpNum)), data=d) + 
scale_color_brewer(palette="Spectral", name="Planted Richness") +
geom_smooth(se=FALSE) + labs(y="Shannon Diversity in Biomass Sample")

ggplot(aes(x=Year, y=biomass, pch=as.factor(Plot), col=as.factor(SpNum)), data=d) + 
scale_color_brewer(palette="Spectral", name="Planted Richness") +
geom_smooth(se=FALSE) + labs(y="Planted Species Biomass")


dev.off()
system("open e120-rich-mass.pdf")

