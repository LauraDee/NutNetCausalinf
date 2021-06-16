#Realized richness analyses of BioCon
# Laura update Dec 16 2020; checked by Kaitlin Kimmel and Tim O.

rm(list=ls())
#Close graphics windows
graphics.off()

#biocon data downloaded in Nov? 2020
library(tidyr) 
library(here)
library(data.table)
library(fixest)  # to run econometric fixed effect models 
require(ggplot2)
library(texreg)
library(here)

# When log(0) need to use inverse hyperbolic sine transformation (REF NEEDED)
#https://en.wikipedia.org/wiki/Inverse_hyperbolic_functions#Inverse_hyperbolic_sine
ihs = function(x) {
  return(log(x + sqrt(x^2+1)))
} # v.s. log(x+1) <- is defined for a negative x. 
## it is generally frowned upon, when taking the log, to add small numbers (e.g., log(richness + 0.01))

### LOAD DATA ##########

#processed by Kaitlin on Dec 16 2020; KK used the biomass files for richness and biomass to be consistent with BigBio.
#setwd("~/Documents/GitHub/NutNetCausalinf/code/r/SMSection8_CedarCreekAnalyses/data/")
findat <-  fread(here::here("code", "r", "SMSection8_CedarCreekAnalyses","data","biocon_plantedbiomass_output.csv"))
# The TreatmentSR column = treatment level of spp
# The ObservedSR = observed species based on what was planted in the plot
# Live.mass = biomass

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

#plot richness 
hist(findat$ObservedSR)
hist(findat$TreatmentSR)

#Check that we have created a VARIABLE FOR  OF THE SPECIES PLANTED,
# WHICH AND HOW MANY COME UP AND ARE MEASURED IN PLANTED BIOMASS
table(findat[TreatmentSR=="1", ObservedSR])  # this should be 1 and 0 

############################################################################################################
## Plot Data ###############################################################################################
############################################################################################################
## plot the raw realized richness over time, at each richness level - color the points by planted richness level
plot(findat$Year, findat$ObservedSR)
plot(findat$TreatmentSR, findat$ObservedSR)
hist(findat$live.mass, main = "Live biomass of planted species")

# species richness in the plots can go down but not up due to weeding of species colonizating the plot
realrichovertime <- ggplot(findat, aes(x = Year, y = ObservedSR)) + facet_wrap(~TreatmentSR) +
  theme_bw() +
  geom_point()+  labs(y = "Realized Planted Species Richness over time") + 
  theme(axis.title.y= element_text(size=14)) + theme(axis.title.x= element_text(size=14)) +
  theme(axis.text.y = element_text(size = 14)) 
realrichovertime 

## plot realized richness over time, for each treatment - color the points by planted richness level
realrichovertime.T <- ggplot(findat, aes(x = Year, y = ObservedSR, colour= TreatmentSR))   + # + facet_wrap(~Trt) +
  theme_bw() +
  geom_point() +  labs(y = "Realized Planted Species Richness over time") + 
  theme(axis.title.y= element_text(size=14)) + theme(axis.title.x= element_text(size=14)) +
  theme(axis.text.y = element_text(size = 14)) 
realrichovertime.T

############################################################################################################################
### Do two-way fixed fffects Analyses using lfe package - see Section SM S8 ######################################################
##############################################################################################################################

#ihs transformed because some of the live.mass is 0. Though how is this possible?
# also we do this because of support for this relation (the effect of SR on productivity) to follow a power relationship,
# or a saturating relationship, see: 
#  Reich,  et al. (2012). Impacts of biodiversity loss escalate through time as redundancy fades. Science, 336, 589–92.
#  Cardinale, B.J., et al. (2011). The functional role of producer diversity in ecosystems. Am. J. Bot., 98, 572–592.
biocon.mod  <- feols(ihs(live.mass) ~ ihs(ObservedSR)  | Year + Plot,  findat, cluster = "Plot") 
# 
# mod1 <- felm(ihs(live.mass) ~ ihs(ObservedSR) | Plot  + Year | 0 | Plot,  data = findat)
# summary(mod1)


########################################### 
#### Print Results for Table S16 ########
########################################### 

esttex(biocon.mod, 
       coefstat = "se",  replace = TRUE,
       file = "Table_S16_R_se.tex")
