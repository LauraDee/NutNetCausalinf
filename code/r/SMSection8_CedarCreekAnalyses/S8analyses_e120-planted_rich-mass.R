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
library(fixest)  # v 0.8.2
library(data.table)
library(lfe)
library(here)
library(texreg)

# When log(0) need to use inverse hyperbolic sine transformation (Bellemare & Wichman 2020 Oxford Bulletin of Economics and Statistics)
#https://en.wikipedia.org/wiki/Inverse_hyperbolic_functions#Inverse_hyperbolic_sine
ihs = function(x) {
  return(log(x + sqrt(x^2+1)))
} # v.s. log(x+1) <- is defined for a negative x. 

#######    READ IN DATA    #######
# read in file that Kaitlin Kimmel processed using the R code e120-biomass-data-2020-12-04.R 
 # d <-fread(here("CedarCreekAnalyses", "data", "e120-plantedbiomass-data-output.csv"), strip.white=T)  

setwd("~/Documents/GitHub/NutNetCausalinf/code/r/SMSection8_CedarCreekAnalyses/data/")
d <- fread("e120-plantedbiomass-data-output.csv")

## meta-data and variables used: ##
  # numsp is the treatment -- Planted number of species: 1, 2, 4, 8, or 16 
  # rich is the observed richness ofthe  planted species
  # mass.live is the biomass of the observed planted species

## check the data:
#check to make sure rich is the observed richness of planted species for THAT plot vs the whole experiment
table(d[numsp=="1", rich])
table(d$rich)
summary(d$rich) 

#planted richness variable - create a factor version of numsp  
class(d$numsp)
d$Treat.numsp <- as.factor(d$numsp)

# Delete missing data
d <- d[!is.na(d$mass.live),]
summary(d)

# plot the data comparing the number of species planted in a treatment versus that are found
plot(d$numsp, d$rich, xlab = "planted species", ylab = "realized richness", main = "BigBio species richness")

##########################################
## #### process data and prep for models ###   
##########################################
# to run for the models to create a dummery variable for each plot and for each year 
d$plot <- as.factor(d$plot)
d$year <- as.factor(d$year)

####################################################################################################################
## Implement the two-way fixed effects design #####################################################################
#####################################################################################################################
mod_main <- feols(ihs(mass.live) ~ ihs(rich)  | year + plot,  d, cluster = "plot") 


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

mod16<- feols(ihs(mass.live) ~ ihs(rich)  | year + plot,  d16, cluster = "plot") 

#mod16 <- felm(ihs(mass.live) ~ ihs(rich) | year + plot | 0 | plot, data = d16) 
summary(mod16, robust = T)  # we find the results are still consistent with our conjecture.

screenreg(mod16,  # object with results 
          custom.model.names= " 16 species only plots",
          override.se=summary(mod16,)$coef[,2],
          override.pval=summary(mod16)$coef[,4])

