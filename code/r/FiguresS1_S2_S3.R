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
merged.data <- fread("merged_NutNet_Oct2019.csv", na.strings = 'NA')

####################################################################################################
## Filter  data to only control plots used  ###############################
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
merged.data$year <- as.character(merged.data$year)
# same with plot
comb$plot <- as.character(comb$plot)

###############################################################
## Figure S1  Plot-level (a &b) and site-level (c & d) biodiversity (as species richness) and productivity (as above-ground live mass) from between 2007-2017. 
###############################################################
par(mfrow=c(2,2))
plot(live_mass ~ rich, data = comb, main = "a)", xlab = "species richness", ylab = "live mass")
plot(log(live_mass) ~ log(rich), data = comb, main = "b)", xlab = "log(species richness)", ylab = "log(live mass)")
## also plot site-level data 
plot(comb$site_richness, comb$ave_site_live_mass, main = "c)", xlab = "species richness", ylab = "average site live mass")
plot(log(comb$site_richness), log(comb$ave_site_live_mass), main = "d)", ylab = "log(average site live mass)", xlab = "log(site richness)") 
mtext("Plot-level", side=3, outer=TRUE, line=-3)
mtext("Site-level", side=3, outer=TRUE, line=-17)

###############################################################
# Figure S2 ###################################################
###############################################################
FigS2 <- ggplot(data = comb, aes(x = changerich)) + geom_histogram()+ facet_wrap(~site_code) + theme_bw() +
  geom_vline(xintercept=c(0,0), color = "blue", linetype="dashed") +
  labs(x = "Plot-level change in species richness year to year") +  theme_bw() +
  theme(axis.title.y= element_text(size=14)) + theme(axis.title.x= element_text(size=12)) +
  theme(axis.text.y = element_text(size = 14)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(size=14)) 
FigS2

############################################################### 
## Figure S3. Change in species evenness per plot (overall data)
###############################################################
FigS3 <- ggplot(data = comb, aes(x = changeEvenness)) + geom_histogram() +
  geom_vline(xintercept=c(0,0), color = "blue", linetype="dashed") +
  labs(x = "Plot-level change in evenness") +  theme_bw() +
  theme(axis.title.y= element_text(size=14)) + theme(axis.title.x= element_text(size=12)) +
  theme(axis.text.y = element_text(size = 14)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(size=14)) 
FigS3






