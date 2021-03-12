###############################################################
## Create SM Data figures for Dee et al ###################
###############################################################
#Close graphics and clear local memory
#graphics.off()
rm(list=ls())

#load libraries
require(ggplot2)
library(plyr)
library(data.table)

setwd("~/Dropbox/IV in ecology/NutNet")
##Load processed Data, processed from version 'comb-by-plot-clim-soil-diversity-28-Apr-2017.csv'
comb <- fread("NutNetControlPlotDataToUseApril2018.csv",na.strings='NA')

####################################################################################################
## Filter cover data to only control plots used  ###############################
##################################################################################################
comb = comb[is.control==T & has.5.yrs.data==T,]

##########################################################################################
### Process comb data to prep to use it in the models & merge #############################
############################################################################################
## Filter data to records with non-NA live_mass and non-NA richness.
## Models can generally handle NAs by implicitly dropping, but clustered SE
# # function currently does not.
 comb = comb[!is.na(live_mass) & !is.na(rich)]

# # make year a character, to be a dummy variable: 
comb$year <- as.character(comb$year)
cover$year <- as.character(cover$year)
# same with plot
comb$plot <- as.character(comb$plot)

## Figure S1
FigS1 <- ggplot(data = comb, aes(x = changerich)) + geom_histogram()+ facet_wrap(~site_code) + theme_bw() +
  geom_vline(xintercept=c(0,0), color = "blue", linetype="dashed") +
  labs(x = "Plot-level change in species richness year to year") +  theme_bw() +
  theme(axis.title.y= element_text(size=14)) + theme(axis.title.x= element_text(size=12)) +
  theme(axis.text.y = element_text(size = 14)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(size=14)) 
FigS1





