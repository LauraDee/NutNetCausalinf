########################################################################
# Two way fixed effect panel models applied to Big Bio Data #############
# Laura Dee Dec 3 2020    ; cleaned Feb 12, 2021            ###############
# Supplemental Materials: Section 8                         ################
########################################################################
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

# read in file that Kaitlin Kimmel processed using the R code e120-biomass-data-2020-12-04.R 
# KK: e120-plantedbiomass-data-output.csv is the one you want to use because it only includes planted species.
d <-fread("e120-plantedbiomass-data-output.csv", strip.white=T)   
# numsp is the treatment -- Planted number of species: 1, 2, 4, 8, or 16 
# rich is the observed richness ofthe  planted species
# mass.live is the biomass of the observed planted species

## to run for the models
d$plot <- as.factor(d$plot)
d$year <- as.factor(d$year)

#### process data and prep for models ### 
#check to make sure rich is the observed richness of planted species for THAT plot vs the whole experiment
table(d[numsp=="1", rich])
 table(d$rich)
 summary(d$rich) # is this right that there is always one species in the monoculture that grows (vs some with 0 in biocon?)

## **to run for the models **
d$plot <- as.factor(d$plot)
d$year <- as.factor(d$year)

#planted richness variable - create a factor version of numsp  
class(d$numsp)
d$Treat.numsp <- as.factor(d$numsp)

# Delete missing data
d <- d[!is.na(d$mass.live),]
summary(d)

plot(d$numsp, d$rich, xlab = "planted species", ylab = "realized richness", main = "BigBio species richness")


####################################################################################################################
## Implement the two-way fixed effects design #####################################################################
#####################################################################################################################
mod_main <- felm(ihs(mass.live ) ~ ihs(rich) | year + plot | 0 | plot, data =d)
summary(mod_main, robust = T)

#print results for SM Table, as in SM Section S8.
screenreg(mod_main,     # object with results 
          custom.model.names= "Main Model for BigBio: Effect on Planted Biomass" ,
          override.se=summary(mod_main,)$coef[,2],
          override.pval=summary(mod_main)$coef[,4])

### Implement on just plots with 16 species planted, because these plots are 
# found to be more representative of nautral communities in some aspects 
   # (e.g., see Jochum, M., et al. (2020). The results of biodiversity–ecosystem functioning experiments are realistic. Nat. Ecol. Evol., 4, 1485–1494.)
d16 = d[Treat.numsp == "16",]
mod16 <- felm(ihs(mass.live) ~ ihs(rich) | year + plot | 0 | plot, data = d16) 
summary(mod16, robust = T)  # we find the results are still consistent with our conjecture.

screenreg(mod16,  # object with results 
          custom.model.names= " 16 species only plots",
          override.se=summary(mod16,)$coef[,2],
          override.pval=summary(mod16)$coef[,4])

##****SHOULD THE REST OF THIS BE CUT ****? ####
##*
###############################################################################################################
## Comparisons to other types of models ######################################################################
###############################################################################################################

#Eric - Random Effect Models -  code by Eric Seabloom sent Nov 19 2020 
#to run Eric's models create mass.live.lg
d$mass.live.lg <- log(d$mass.live + .001)

# 1. Effects of observed richness on biomass in all plots
summary(lme1<-lmer(mass.live.lg~rich + (1|plot) +(1|year), data=d))

# 2. Effects of observed richness on biomass in plots with 8 or more species
summary(lme1<-lmer(mass.live.lg~rich + (1|plot) +(1|year), data=d[d$numsp >= 8,]))    ## *** ERIC NOT SURE HOW YOU RAN THIS BECAUSE numsp is a factor??? **

# 3. Effects of observed richness on biomass in 16 specie splots
summary(lme1<-lmer(mass.live.lg~rich + (1|plot) +(1|year), data=d[d$numsp >= 16,]))  ## *** ERIC NOT SURE HOW YOU RAN THIS BECAUSE numsp is a factor??? **



#4. check for lower spp plots too:
#in only plots with 8 species 
mod8 <- felm(log(mass.live) ~ log(rich) | year + plot | 0 | plot, data =d[d$numsp == 8,])
summary(mod8)

mod4 <- felm(log(mass.live) ~ log(rich) | year + plot | 0 | plot, data =d[d$numsp == 4,])
summary(mod4)

mod2 <- felm(ihs(mass.live) ~ log(rich) | year + plot | 0 | plot, data =d[d$numsp == 2,])
summary(mod2)


### do within numsp - within treatment, with the the treatment as the fixed effect:
withinmod <- felm(ihs(mass.live) ~ log(rich) |  numsp + year | 0 | plot, data =d)
summary(withinmod)

###############################################################################
## Compare to analyzing as an experimental treatment of planted richness #####
##############################################################################

#### 1. FIRST RUN ANALYSES WITH SR AS NUMERIC, no fixed effects of SE clustering
# SR as a integer, both transformed on on a log-log scale - no fixed effects or clustering of the standard error
treat.effect <- felm(log(mass.live + .001) ~ log(numsp) | 0 | 0 | 0, data =d)
summary(treat.effect)
# instead, transformed using an ihs (inverse hyperbolic sign approximation of a log)  -- no fixed effects or clustering of the error
treat.effect2 <- felm(ihs(mass.live) ~ ihs(numsp) | 0 | 0 | 0, data =d)
summary(treat.effect2)

#no transformation:
treat.effect3 <- felm(mass.live  ~ numsp | 0 | 0 | 0, data =d)
summary(treat.effect3)

#### 2. Second, RUN ANALYSES WITH SR AS NUMERIC, no fixed effects, but cluster SE on the plot
treat.effect.cluster <- felm(log(mass.live + .01) ~ log(numsp) | 0 | 0 | plot, data =d)
summary(treat.effect.cluster)
#ihs transformaton:
treat.effect.cluster2 <- felm(ihs(mass.live + .01) ~ ihs(numsp) | 0 | 0 | plot, data =d)
summary(treat.effect.cluster2)
#untransformed
treat.effect.cluster3 <- felm(mass.live  ~ numsp | 0 | 0 | plot, data =d)
summary(treat.effect.cluster3)

##### 3. Add year as a fixed effect ##### It's unchanged:
treat.effect.cluster.yr <- felm(log(mass.live + .01) ~ log(numsp) | year | 0 | plot, data =d)
summary(treat.effect.cluster.yr)
#ihs transformaton:
treat.effect.cluster2.yr <- felm(ihs(mass.live + .01) ~ ihs(numsp) | year | 0 | plot, data =d)
summary(treat.effect.cluster2.yr)
#untransformed
treat.effect.cluster3.yr <- felm(mass.live  ~ numsp | year | 0 | plot, data =d)
summary(treat.effect.cluster3.yr)

###########################################################################
### RUN THE ABOVE WITH LM for comparison -- ##############################
##########################################################################
treat.effectlm <- lm(log(mass.live + 0.01) ~ log(numsp), data =d)
summary(treat.effectlm)
treat.effectlm1 <- lm(ihs(mass.live) ~ ihs(numsp), data = d)
summary(treat.effectlm1)
treat.effectlm2 <- lm(mass.live ~ numsp, data =d)
summary(treat.effectlm2)

#add year as fixed effect in lm:
treat.effectlm.yr <- lm(log(mass.live + 0.01) ~ log(numsp) + year, data =d)
summary(treat.effectlm.yr)
treat.effectlm1.yr <- lm(ihs(mass.live) ~ ihs(numsp) +year, data = d)
summary(treat.effectlm1.yr)
treat.effectlm2.yr <- lm(mass.live ~ numsp +year, data =d)
summary(treat.effectlm2.yr)

####################################################################################
### RUN models estimating the treatment effect with Planted Richness as a factor ######
#######################################################################################

#no fixed effects or clustering
treat.effect.factor <- felm(log(mass.live + 0.001) ~ Treat.numsp | 0 | 0 | 0, data =d)
summary(treat.effect.factor)
#ihs transformed
treat.effect.factor <- felm(ihs(mass.live) ~ Treat.numsp | 0 | 0 | 0, data =d)
summary(treat.effect.factor)
#untransformed
treat.effect.factor <- felm(mass.live ~ Treat.numsp | 0 | 0 | 0, data =d)
summary(treat.effect.factor)

#add clustering of SE on the plot level as the cluster
#as a factor, log-levels 
treat.effect <- felm(log(mass.live + 0.001) ~ Treat.numsp |  0 | 0 | plot, data =d)
summary(treat.effect)
#as a factor, ihs transformed mass.live 
treat.effect <- felm(ihs(mass.live) ~ Treat.numsp |  0 | 0 | plot, data =d)
summary(treat.effect)
#as a factor, levels-levels (untransformed ) 
treat.effect.untransformed <- felm(mass.live ~ Treat.numsp |  0 | 0 | plot, data =d)
summary(treat.effect.untransformed)

## Add year as a fixed effect - and clustered SEs
  #as a factor, log-levels with a year fixed effect
treat.effect <- felm(log(mass.live +0.01) ~ Treat.numsp |  year | 0 | plot, data =d)
summary(treat.effect)
  #as a factor, ihs transformed biomass (mass.live) with a year fixed effect
treat.effect <- felm(ihs(mass.live) ~ Treat.numsp |  year | 0 | plot, data =d)
summary(treat.effect)
  #as a factor, levels-levels (untransformed ) with a year fixed effect
treat.effect.untransformed <- felm(mass.live ~ Treat.numsp |  year | 0 | plot, data =d)
summary(treat.effect.untransformed)

#######################################################################################################################################
###### Run the same models as above with treatment of planted richness levels (factor) as a regular LM 
##########################################################################################################
treat.effectlm <- lm(log(mass.live + 0.01) ~ Treat.numsp, data =d)
summary(treat.effectlm)
treat.effectlm1 <- lm(ihs(mass.live) ~ Treat.numsp, data = d)
summary(treat.effectlm1)
treat.effectlm2 <- lm(mass.live ~ Treat.numsp, data =d)
summary(treat.effectlm2)

#add year as fixed effect in lm:
treat.effectlm.yr <- lm(log(mass.live + 0.01) ~ Treat.numsp + year, data =d)
summary(treat.effectlm.yr)
treat.effectlm1.yr <- lm(ihs(mass.live) ~ Treat.numsp + year, data = d)
summary(treat.effectlm1.yr)
treat.effectlm2.yr <- lm(mass.live ~ Treat.numsp+ year, data =d)
summary(treat.effectlm2.yr)



